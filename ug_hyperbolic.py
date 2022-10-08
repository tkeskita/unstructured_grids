# ##### BEGIN GPL LICENSE BLOCK #####
#
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation; either version 3
#  of the License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software Foundation,
#  Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
#
# ##### END GPL LICENSE BLOCK #####

# <pep8 compliant>

# Operations related to hyperbolic extrusion of new cells

import bpy
from . import ug
from . import ug_op
from .ug_extrude import *
import logging
l = logging.getLogger(__name__)

LARGE = 100.0  # Cut-off for large values for weight function
ug_props = bpy.context.scene.ug_props

def extrude_cells_hyperbolic(niter, bm, bmt, speeds, new_ugfaces, \
                             initial_face_areas, is_last_layer):
    '''Extrude new cells from current face selection. Initial faces
    argument provides optional list of initial UGFaces whose normal direction
    is reversed at the end.
    '''

    import bmesh
    bm.verts.ensure_lookup_table()
    bm.faces.ensure_lookup_table()

    # Selected faces are base faces for extrusion. New cell index number is
    # index number of mesh face in faces list. Actual UGCell index
    # number is cell number plus initial number of prior ugcells.
    base_faces = [f for f in bm.faces if f.select]
    if fulldebug: l.debug("Face count at beginning: %d" % len(base_faces))

    ugci0 = len(ug.ugcells) # Number of UGCells before extruding
    ugfi0 = len(ug.ugfaces) # Number of UGFaces before extruding

    base_verts, base_vis_of_fis, base_fis_of_vis, ibasevertmap = \
        get_verts_and_relations(base_faces)

    # Populate initial verts and faces to trajectory bmesh
    if ug_props.extrusion_create_trajectory_object and len(bmt.verts) == 0:
        for v in base_verts:
            bmt.verts.new(v.co)
        bmt.verts.ensure_lookup_table()
        bmt.verts.index_update()
        for vis in base_vis_of_fis:
            verts = [bmt.verts[vi] for vi in vis]
            bmt.faces.new(verts)

    base_edges, edge2sideface_index, fis2edges = get_edges_and_face_map(base_faces)
    new_ugfaces = create_UG_verts_faces_and_cells(base_verts, base_edges, base_faces, new_ugfaces)

    # Generate vertex classification information
    if ug_props.extrusion_method == "hyperbolic":
        df = 1e-3  # Initial step size for hyperbolic extrusion
        is_corners, is_boundaries, bnvis_of_vi, anvis_of_vi = \
            classify_verts(base_verts, base_edges, base_faces, \
                           base_fis_of_vis, ibasevertmap)
        if fulldebug:
            l.debug("is_corners %s" % str(is_corners))
            l.debug("is_boundaries %s" % str(is_boundaries))
            l.debug("bnvis_of_vi %s" % str(bnvis_of_vi))
            l.debug("anvis_of_vi %s" % str(anvis_of_vi))

        neighbour_vis_of_vi, fils_of_neighbour_vis = \
            get_face_vis_of_vi(base_verts, base_fis_of_vis, base_vis_of_fis)


    # Initialization of state variables
    if len(speeds) == 0:
        # Initial speeds from vertex normal speeds
        speeds = get_vertex_normal_speeds(base_verts, base_faces, \
                                          base_fis_of_vis)

    bm, top_verts, vert_map = cast_vertices(\
        bm, base_verts, speeds, df)
    top_faces = create_mesh_faces(bm, edge2sideface_index, vert_map, \
                                  base_faces, fis2edges, ugci0, ugfi0)
    correct_face_normals(bm, base_faces, ugci0)

    if ug_props.extrusion_method == "hyperbolic":
        fi_areas = calculate_face_areas(base_faces)
        area_coeffs = calculate_mean_vertex_area_change( \
            fi_areas, initial_face_areas, base_fis_of_vis)
        intvis, intvpairs0, intvpairs1 = \
            calculate_intvpairs(anvis_of_vi, base_fis_of_vis, is_boundaries)

        # Save vertex normal speeds for cone limitation
        vnspeeds0 = get_vertex_normal_speeds(top_verts, top_faces, \
                                             base_fis_of_vis)

        # Save top vertex locations for cone limitation
        coords0 = [v.co for v in top_verts]

        # Carry out iterations until one whole layer is done
        layer_frac = 0.0 # Fraction of simulated layer
        while layer_frac < 1.0 - EPS:
            if fulldebug: l.debug("===== Extrusion iteration %d =====" % niter)

            bm, speeds, layer_frac, df, target_cos = evolve_iteration( \
                bm, top_verts, speeds, is_boundaries, is_corners, anvis_of_vi, \
                intvpairs0, intvpairs1, top_faces, base_fis_of_vis, \
                coords0, vnspeeds0, area_coeffs, layer_frac, \
                neighbour_vis_of_vi)

            if ug_props.extrusion_create_trajectory_object:
                top_vert_cos = [v.co for v in top_verts]
                bmt = update_trajectory_mesh(bmt, top_vert_cos) # Current coordinates
                #bmt = update_trajectory_mesh(bmt, target_cos) # Target coordinates

            if print_iterations:
                max_vel = max([x.length for x in speeds])
                l.debug("iter %d: frac %f, df %f, max_vel %f" % (niter, layer_frac, df, max_vel))

            niter += 1

        # Add faces to trajectory object at last layer
        if is_last_layer and ug_props.extrusion_create_trajectory_object:
            nv0 = len(bmt.verts) - len(top_verts)
            bmt = add_faces_to_trajectory_mesh(bmt, nv0, base_vis_of_fis)

    add_base_face_to_cells(base_faces, ugci0)
    thickness_update()

    # End of extrude_cells()
    return niter, bm, bmt, len(base_faces), speeds, new_ugfaces


def evolve_iteration(bm, top_verts, speeds, is_boundaries, is_corners, anvis_of_vi, \
                     intvpairs0, intvpairs1, top_faces, base_fis_of_vis, \
                     coords0, vnspeeds0, area_coeffs, layer_frac, \
                     neighbour_vis_of_vi):
    '''Evolve hyperbolic extrusion by one iteration'''

    min_velocity = ug_props.extrusion_thickness # Absolute minimum velocity

    # Calculate vertex normal speeds
    vnspeeds = get_vertex_normal_speeds(top_verts, top_faces, \
                                        base_fis_of_vis)

    # Current coordinates
    current_cos = [v.co for v in top_verts]

    # Calculate convexity_sums = sum of convex edge convexity values per vertex
    convexity_sums = calculate_convexity_sums(bm, top_verts, top_faces)

    ### First round: Calculate targets, accelerations and step sizes ###

    dfs = dict() # step sizes for each vertex
    accs = dict() # accelerations for each vertex
    target_cos = [] # target coordinates

    for vi, v in enumerate(top_verts):
        # Skip corners and boundaries
        if is_corners[vi] or is_boundaries[vi]:
            target_cos.append(v.co)
            continue

        if fulldebug:
            l.debug("1. Starting on internal vi %d index %d" % (vi, v.index))

        # Generate bmesh and bvhtree of projected faces surrounding this
        # vertex, with vertices located at vertex normal projected coordinates
        bmp, bt = get_bmp(vi, top_faces, top_verts, base_fis_of_vis, \
                          anvis_of_vi, current_cos)

        # Calculate projected geometric mean target coordinates
        # TODO: Restore option for using anvis or neighbour_vis?
        # pvgm_target_co = get_pvgm_target_co(vi, top_verts, anvis_of_vi[vi], \
        pvgm_target_co = get_pvgm_target_co(vi, top_verts, neighbour_vis_of_vi[vi], \
                                            speeds, bmp, bt, vnspeeds, \
                                            convexity_sums[vi])

        # Calculate pvs = project neighbour vertices to direction normal plane
        # Calculate pvspeeds = project neighbour speeds to direction normal plane
        pvs, pvspeeds = project_pvs_pvspeeds(vi, top_verts, anvis_of_vi[vi], vnspeeds)

        # Calculate convex_target_cos = target coordinates for all vpairs
        # is_aboves = booleans for original target being located above faces
        # vplengths = lengths between vpairs
        convex_target_cos, is_aboves, vplengths = \
            get_convex_target_cos(v, vi, intvpairs0[vi], intvpairs1[vi], \
                                  current_cos, bmp, bt, vnspeeds[vi], \
                                  convexity_sums[vi])

        # Calculate layers until vpairs intersect
        vimap = {} # Map from vi index to anvis index
        for i, x in enumerate(anvis_of_vi[vi]):
            vimap[x] = i
        cut_layers = get_cut_layers(intvpairs0[vi], intvpairs1[vi], vimap, pvs, pvspeeds)

        # Calculate average neighbour convexity sum
        anc = get_ave_neigh_convexity(convexity_sums, anvis_of_vi[vi])

        # Calculate convex target weights
        convex_weights = get_convex_weights(is_aboves, cut_layers, vplengths)

        # Calculate convex target
        convex_target_co = get_convex_target(v, convex_target_cos, convex_weights, pvgm_target_co)

        # Calculate target coordinates
        target_co = get_target_co(v, convexity_sums[vi], anc, v.co, \
                                  pvgm_target_co, convex_target_co)
        target_cos.append(target_co)

        # Calculate acceleration
        acc = get_acceleration(v, speeds[vi], target_co, min_velocity, convexity_sums[vi])
        accs[vi] = acc

        # Calculate step size
        min_len = get_min_len_from_faces(v, top_verts, anvis_of_vi[vi])
        df = get_step_size(layer_frac, min_len, speeds[vi].length, acc, min_velocity)
        dfs[vi] = df

        if fulldebug:
            l.debug("pvgm_target_co %s" % str(pvgm_target_co))
            l.debug("convex_target_co %s" % str(convex_target_co))
            l.debug("-- is_aboves %s" % str(is_aboves))
            l.debug("-- convex_target_cos %s" % str(convex_target_cos))
            l.debug("-- cut_layers %s" % str(cut_layers))
            l.debug("-- convex_weights %s" % str(convex_weights))
            l.debug("convexity_sum %f" % convexity_sums[vi] \
                    + ", anc %f" % anc)

        bmp.free()
        del bt

    ### Second round: Make the step ###

    # Use smallest step size
    if len(dfs) == 0:
        df = 0.01  # 1.0
    else:
        df = min(dfs.values())
    if fulldebug:
        l.debug("df %f layer_frac %f" % (df, layer_frac))
    layer_frac += df

    new_speeds = [] # new speed vectors, to be calculated
    unscaled_new_speeds = [] # new speed vectors (without area change scaling)

    for vi, v in enumerate(top_verts):
        # Speed of corner vertices are not changed
        if is_corners[vi]:
            new_speeds.append(speeds[vi])
            unscaled_new_speeds.append(speeds[vi])
            continue

        # Boundary vertices use plain vertex normal speed
        if is_boundaries[vi]:
            # TODO: Restore boundary neighbour interpolation method?
            new_speeds.append(vnspeeds[vi])
            unscaled_new_speeds.append(vnspeeds[vi])
            continue

        if fulldebug:
            l.debug("2. internal vi %d index %d" % (vi, v.index))

        # Limit target
        target_speed, unscaled_target_speed = limit_target( \
            v, speeds[vi], target_cos[vi], convexity_sums[vi], min_velocity, \
            coords0[vi], vnspeeds0[vi], area_coeffs[vi], df, accs[vi])

        if fulldebug:
            l.debug("FINAL target_co %s " % str(v.co + target_speed) \
                    + "velocity %f" % target_speed.length)

        new_speeds.append(target_speed)
        unscaled_new_speeds.append(unscaled_target_speed)


    # Evolve vertex positions and return new speeds
    for v, s in zip(top_verts, new_speeds):
        v.co += s * df

    # TODO: Not needed? Remove
    for f in top_faces:
        f.normal_update()

    # End of evolve_iteration()
    return bm, unscaled_new_speeds, layer_frac, df, target_cos


def control_acceleration(ndv, ndl):
    '''Calculate acceleration from normalized velocity difference (ndv) and
    normalized length difference (ndl) from source to target
    '''
    slope = 1.0 # Slope for control line # TODO: Parametrize?
    dacc = 0.1 # Acceleration increase factor # TODO: Parametrize?
    ndv0 = slope * ndl
    acc = (ndv0 - ndv) * dacc
    return acc


def get_acceleration(oldv, oldspeed, target_co, tgvel, convexity_sum):
    '''Calculate acceleration for this vertex using current location,
    target and speed data
    '''
    thickness = ug_props.extrusion_thickness
    max_acc = 20.0 # Maximum acceleration
    
    vel = oldspeed.length # Current velocity
    vdir = oldspeed.normalized()
    u = target_co - oldv.co # Vector to target coordinates
    
    # Project u to vdir normal plane
    cos_alpha_u = u @ vdir
    q = cos_alpha_u * vdir # u component aligned with vdir
    target_len = q.length # Length to projected target

    # Deduce is vertex ahead or behind of target coordinates:
    vert_is_ahead = (cos_alpha_u >= 0.0)

    # Normalized velocity difference
    ndv = (vel - tgvel) / tgvel

    # Normalized length difference to target
    if cos_alpha_u >= 0.0:
        ndl = target_len / thickness
    else:
        ndl = -target_len / thickness

    # Calculate acceleration
    acc = control_acceleration(ndv, ndl)
    acc = max(-max_acc, min(max_acc, acc))

    return acc


def get_min_len_from_faces(v, top_verts, anvis_of_vi):
    '''Calculate minimum edge length from vertex to it's neighbour
    vertices
    '''
    from sys import float_info
    min_len = float_info.max
    for nvi in anvis_of_vi:
        vec = top_verts[nvi].co - v.co
        elen = vec.length
        min_len = min(min_len, elen)
    return min_len


def get_step_size(layer_frac, min_len, vel, acc, min_velocity):
    '''Calculate iteration step size (unit 1/layer) using Courant
    number. layer_frac is the current fraction value,
    top_faces are BMesh faces, and max_velocity is the maximum
    velocity in units m/layer.
    '''

    # from .ug_checks import get_edge_stats_from_bmesh_faces

    # Option 1. Maximum normalized step size
    df1 = 0.2 # TODO: Parametrize?

    # Option 2. Courant limitation for distance
    courant = ug_props.extrusion_courant_number
    df2 = courant * min_len / vel

    # Option 3. Courant limitation for acceleration
    df3 = courant * vel / abs(acc) # TODO: CHECKME

    df = min(df1, df2, df3) # Choose smallest
    #df = max(0.05, df) # But clamp at a minimum step size # TODO: Parametrize?
    # Do not overshoot layer thickness
    if (layer_frac + df) > 1.0:
        df = 1.0 - layer_frac + EPS
        
    if df < 0.0:
        raise Exception("df %f < 0.0" % df)

    return df


# TODO: Remove?
def get_max_vec_length(vectors):
    '''Calculate maximum length of vectors'''

    max_len = 0.0
    for vec in vectors:
        if vec.length > max_len:
            max_len = vec.length
    return max_len


def get_estimated_cos(verts, speeds, df):
    '''Calculate locations of vertices from speeds and step size
    '''
    estimated_cos = []
    for v, s in zip(verts, speeds):
        estimated_cos.append(v.co + df * s)
    return estimated_cos


def get_bmp(vi, top_faces, top_verts, fis_of_vis, anvis_of_vi, \
            new_cos):
    '''Generate bmesh and bvhtree with projected estimated surrounding
    faces, with verts at new coordinates
    '''
    import bmesh
    bmp = bmesh.new()
    vertmap = {} # map from old vertex to projected vertex

    # Add base vertex
    v = bmp.verts.new(new_cos[vi])
    vertmap[top_verts[vi]] = v

    # Add vertices
    for oldvi in anvis_of_vi[vi]:
        v = bmp.verts.new(new_cos[oldvi])
        vertmap[top_verts[oldvi]] = v

    # Add triangle faces generated from neighbour verts to bmesh
    oldfaces = [top_faces[x] for x in fis_of_vis[vi]]
    for oldf in oldfaces:
        verts = []
        for oldv in oldf.verts:
            if oldv in vertmap:
                verts.append(vertmap[oldv])
        f = bmp.faces.new(verts)
        f.normal_update()

    # Generate bvhtree
    from mathutils.bvhtree import BVHTree
    bt = BVHTree.FromBMesh(bmp)
    return bmp, bt


def is_above_planes(v, speed0, co, bmp, bt):
    '''Return True if co is located above normal plane (vector speed0)
    at vertex v location, and above (=on the normal direction side)
    of all faces, all of which connect at v. Otherwise return False.
    '''

    # Return False if co is located below vector normal plane
    u = co - v.co
    m = u.normalized()
    vdir = speed0.normalized()
    cos_alpha = m @ vdir
    if cos_alpha < 0.0:
        return False

    # Cast rays from co towards face centers and check the
    # normal direction. Rays are cast until counter results in
    # difference of two. Default is not above.

    hits_from_above = 0 # Counter
    for f in bmp.faces:
        center = f.calc_center_median()
        ray_dir = center - co
        ray_dir.normalize()

        hit_co, hit_nor, hit_index, hit_length = \
            bt.ray_cast(co, ray_dir)
        if not hit_co:
            continue

        cos_angle = hit_nor @ ray_dir
        if cos_angle < 0.0:
            hits_from_above += 1
        else:
            hits_from_above -= 1

        if fulldebug:
            l.debug("co %s hit_co %s" % (str(co), str(hit_co)))
            l.debug("hit_nor %s ray_dir %s" % (str(hit_nor), str(ray_dir)))
            l.debug("hit_nor @ ray_dir cos_angle %f" % cos_angle)

        if abs(hits_from_above) > 1:
            break

    if fulldebug: l.debug("hits_from_above %d" % hits_from_above)
    if hits_from_above > 1:
        return True
    else:
        return False


def project_co_above_plane(baseco, co, vec):
    '''Project coordinate co above vector normal plane located
    at baseco.
    '''
    from mathutils import Vector
    newco = Vector(co)
    u = co - baseco
    m = u.normalized()
    vdir = vec.normalized()
    cos_alpha = m @ vdir
    # Project u to normal plane if needed
    if cos_alpha < 0.0:
        q = (u @ vdir) * vdir # u component aligned with vdir
        p = u - q # u component normal to vdir
        newco = baseco + p
    return newco


def project_co_to_planes(baseco, co, bt, speed0):
    '''Project coordinates co to faces in bvhtree bt and above v speed
    normal plane. First return value is True if projection was
    successful (False otherwise). Second return value is hit
    coordinates.
    '''
    from mathutils import Vector
    speed = Vector(speed0)
    
    # Use ray casting to get collision point with projected mesh
    hit_co, hit_nor, hit_index, hit_length = bt.ray_cast(co, speed)
    if not hit_co:
        speed.negate()
        hit_co, hit_nor, hit_index, hit_length = bt.ray_cast(co, speed)
    if not hit_co:
        hit_co = Vector(co)

    # Project hit co above speed normal plane at baseco (+ speed)
    u = hit_co - baseco
    m = u.normalized()
    vdir = speed0.normalized()
    cos_alpha = m @ vdir

    # Project u above vdir normal plane
    if cos_alpha < 0.0:
        q = (u @ vdir) * vdir # u component aligned with vdir
        hit_co += speed0 - q

    return hit_co


def get_pvgm_target_co(vi, top_verts, anvis, speeds, bmp, bt, vnspeeds, convexity_sum):
    '''Calculate geometric mean from vertex neighbours' locations and
    speeds. Target coordinate is projected up from planes.
    '''
    cv = top_verts[vi]
    from mathutils import Vector
    co = Vector((0, 0, 0))
    for i in anvis:
        v = top_verts[i]
        s = speeds[i]
        co += v.co + EPS * s
    co /= float(len(anvis))

    # If target is below planes, project it above
    if not is_above_planes(cv, vnspeeds[vi], co, bmp, bt):
        co = project_co_to_planes(cv.co, co, bt, EPS * vnspeeds[vi])

    # Extra projection for convex verts
    cfac = ug_props.extrusion_convex_speed_factor
    co += cfac * convexity_sum * (co - cv.co)

    return co


def project_pvs_pvspeeds(vi, verts, anvis, speeds):
    '''Project neighbour vertex locations and speeds to vertex vi speed
    normal plane
    '''
    from mathutils import Vector
    cv = verts[vi] # center vertex
    cvdir = Vector(speeds[vi]) # copy of center vertex normal vector
    cvdir.normalize()

    nvs = [verts[x] for x in anvis] # neighbour vertices
    nvspeeds = [speeds[x] for x in anvis] # neighbour speeds

    pvs = [] # projected neighbour coordintes (center vertex at origin)
    pvspeeds = [] # projected neighbour speeds

    i = 0
    for v, speed in zip(nvs, nvspeeds):
        # Project coordinates
        u = v.co - cv.co # vector from center vertex to neighbour vertex
        q = (u @ cvdir) * cvdir # u component aligned with cvdir
        p = u - q # u component normal to cvdir
        pvs.append(p)

        # Project speed
        u = speed
        q = (u @ cvdir) * cvdir
        p = u - q
        pvspeeds.append(p)

    return pvs, pvspeeds


def get_convex_target_cos(v, vi, vpairs0, vpairs1, estimated_cos, \
                          bmp, bt, vnspeed, convexity_sum):
    '''First return value is target coordinates (middle coordinates
    which are projected above planes) for all vpairs.
    Second return value is booleans for middle coordinates
    being originally above faces.
    Third return value is distances between vpairs.
    '''
    from mathutils import Vector
    mid_cos = [] # Middle coordinates
    is_aboves = [] # Booleans for above planes
    vplengths = [] # Distances between vpairs
    
    for v0, v1 in zip(vpairs0, vpairs1):
        # Check if co is above planes
        co = Vector((estimated_cos[v0] + estimated_cos[v1]) / 2.0)
        is_above = is_above_planes(v, vnspeed, co, bmp, bt)
        is_aboves.append(is_above)

        # If needed, project co above planes and vertex normal plane
        if not is_above:
            co = project_co_to_planes(v.co, co, bt, EPS * vnspeed)

        # Extra projection for convex verts
        cfac = ug_props.extrusion_convex_speed_factor
        co += cfac * convexity_sum * (co - v.co)

        mid_cos.append(co)
        vec = estimated_cos[v1] - estimated_cos[v0]
        vplengths.append(vec.length)
    return mid_cos, is_aboves, vplengths


def get_cut_layers(vpairs0, vpairs1, vimap, pvs, pvspeeds):
    '''Calculate how many layer extrusions it would take for
    intersection to happen for each neighbour vertex pair if
    neighbours continued with their current pspeeds
    '''
    cut_layers = [] # estimated layers until intersection

    # Process each vertex pair
    for i, j in zip(vpairs0, vpairs1):
        co0 = pvs[vimap[i]]
        co1 = pvs[vimap[j]]
        s0 = pvspeeds[vimap[i]]
        s1 = pvspeeds[vimap[j]]

        # Calculate layers until intersection
        vec0 = co1 - co0 # vector between verts
        vdir = vec0.normalized()
        q0 = (s0 @ vdir) * vdir # speed component aligned with vdir
        q1 = (s1 @ vdir) * vdir # speed component aligned with vdir

        vec1 = (co1 + q1) - (co0 + q0) # vector between moved verts

        # If distance is increasing (not going to intersect),
        # then set cut_layer to LARGE
        SAFETY = 0.999 # avoid singularity
        if vec1.length > SAFETY * vec0.length:
            cut_layers.append(LARGE)

        # Otherwise calculate layers until intersection
        else:
            cut_step = vec0.length / (vec0.length - vec1.length)
            cut_step = min(LARGE, cut_step)
            cut_layers.append(cut_step)
    return cut_layers


def get_ave_neigh_convexity(convexity_sums, anvis):
    '''Calculate average convexity among neighbours'''
    cs = [convexity_sums[x] for x in anvis]
    return sum(cs) / float(len(cs))


def get_convex_weights(is_aboves, cut_layers, vplengths):
    '''Calculate weights for target coordinates. Weights are used to
    calculate a target coordinate and speed.
    '''

    def weight(x):
        '''Weight function'''
        z = ug_props.extrusion_weight_smoothing_coefficient
        wmax = 1.0 / z # Normalisation factor
        w = 1.0 / (x + z) / wmax
        return w * w

    weights = []

    # Weights for target coordinates
    if len(is_aboves) > 0:
        max_vplen = max(vplengths)
        for ia, cs, vl in zip(is_aboves, cut_layers, vplengths):
            # Weight option 1: Closeness of vpair
            vplen_fac = (1.0 + 10.0 * (vl / max_vplen)) # TODO: Parametrize
            wo1 = weight(vplen_fac)
            # Weight option 2: Value of cut layers
            wo2 = weight(cs)
            w = max(wo1, wo2)
            weights.append(w)
            if fulldebug:
                l.debug("vfrac %f " % (vl / max_vplen) \
                        + "closeness %f " % wo1 \
                        + "cut_layers %f " % wo2 \
                        + "w %f" % w)

    return weights


def get_convex_target(v, convex_target_cos, convex_weights, pvgm_target_co):
    '''Calculate convex target from target coordinates and weights'''

    # Default to geometric mean if there is no other data
    if not convex_target_cos:
        return pvgm_target_co

    from mathutils import Vector
    convex_target_co = Vector((0, 0, 0))

    # Target is sum of weighted coordinates
    for co, w in zip(convex_target_cos, convex_weights):
        convex_target_co += co * w
    convex_target_co /= sum(convex_weights)
    return convex_target_co


# TODO: Remove
def get_min_neighbour_vel(speeds, anvis_of_vi):
    '''Calculate minimum neighbour velocity'''
    vels = [x.length for x in get_items_from_list(speeds, anvis_of_vi)]
    return min(vels)


def steering_f(alpha):
    '''Calulate steering factor from angle between current direction and
    target direction. Zero results in no direction change,
    one in full steering.
    '''
    alphamin = -0.7
    alphamax = 0.2
    slope = 1.0 / (alphamax - alphamin)
    val = slope * (alpha - alphamin)
    val = min(1.0, max(0.0, val))
    return val


def limit_co_by_angle_deviation(co, baseco, speed, min_cos_alpha):
    '''Limit coordinates co by angle deviation by moving
    coordinate inside a conical volume. Cone starting point is
    baseco, and cone axis normal vector (vdir) is calculated from
    speed vector. min_cos_alpha defines the angle of the cone.
    '''
    from math import sqrt

    # Vector from vertex to new coordinates
    u = co - baseco
    m = u.normalized()
    vdir = speed.normalized()
    cos_alpha = m @ vdir

    # Target for angle deviation is limited to a mix between
    # co0 (small distance directly ahead) and
    # co1 (point somewhere above speed normal plane).
    EPSC = 1e-4 # Small fraction for jump distance
    jump = 2.0 * EPSC * speed # Small jump vector
    co0 = baseco + jump
    co1 = project_co_above_plane(baseco, co, speed) + jump
    sf = steering_f(cos_alpha)
    target_co = (1-sf) * co0 + sf * co1

    # New target vector
    u = target_co - baseco

    # Project u to vdir normal plane
    q = (u @ vdir) * vdir # u component aligned with vdir
    p = u - q # u component normal to vdir

    # Decrease angle by moving co on vdir normal plane
    if cos_alpha < min_cos_alpha:
        plen = sqrt(1.0/(min_cos_alpha * min_cos_alpha) - 1.0)
        p.normalize()
        p = plen * q.length * p
        u = p + q

    # New coordinates
    return baseco + u


def get_target_co(v, convexity_sum, anc, vn_target, pvgm_target, convex_target):
    '''Calculate an unconstrained target coordinate from several options'''

    def mix_f(x, cut_off):
        '''Linear mixing function for [0, cut_off]'''
        slope = 1.0 / cut_off
        y = slope * x
        y = max(0.0, min(1.0, y)) # Limit to [0, 1]
        return y

    # target1: Mix vertex normal target with geometric mean target by

    # Average neighbour convexity sum
    par1 = ug_props.extrusion_cut_off_anc
    fac1 = mix_f(anc, par1)

    # Convexity factor: Make convex vertices always favor
    # geometric mean over vertex normal target
    cf = mix_f(convexity_sum, 1e-3) # TODO: Need to tune value?

    # Minimum geometric target mix fraction
    mgf = ug_props.extrusion_minimum_geometric_frac

    fac1 = max(mgf, fac1, cf)
    target1 = fac1 * pvgm_target + (1 - fac1) * vn_target

    # target2: Mix target1 with convex target by convexity_sum
    par2 = ug_props.extrusion_cut_off_convexity
    fac2 = mix_f(convexity_sum, par2)
    target2 = fac2 * convex_target + (1 - fac2) * target1

    if fulldebug:
        l.debug("fac1 %f, fac2 %f" % (fac1, fac2))
        l.debug("mgf %f, cf %f" % (mgf, cf))
        l.debug("target2 %s" % str(target2))

    return target2


def limit_target(oldv, oldspeed, target_co, \
                 convexity_sum, min_velocity, coord0, \
                 vnspeed0, area_coeff, df, acc):
    '''Constrain and adjust velocity and direction of speed'''
    from math import sqrt
    from mathutils import Vector

    # Update velocity
    vel = oldspeed.length
    vel += acc * df

    ### Velocity limitations ###
    
    # Maximum allowed velocity limit
    #vel = min(vel, max_vel)

    # Minimum allowed velocity limit
    vel = max(vel, min_velocity)

    ### Geometric limitations ###

    # Limit the change of direction angle compared to previous direction
    min_cos_alpha = ug_props.extrusion_deviation_angle_min
    limited_co = limit_co_by_angle_deviation \
        (target_co, oldv.co, oldspeed, min_cos_alpha)

    # Limit coordinates within accepted region cone
    cone_cos_alpha = ug_props.extrusion_cone_angle
    limited_co = limit_co_by_angle_deviation \
        (limited_co, coord0, vnspeed0, cone_cos_alpha)

    # limited_co = target_co # TEST
    vec = limited_co - oldv.co
    if vec.length < EPS:
        vec = Vector(oldspeed)
    vec.normalize()

    # Speed-up factor from surrounding top face area change
    #speedup = max(1.0, sqrt(area_coeff))
    speedup = 1.0 # TODO: Disabled area scaling for now
    unscaled_target_speed = vel * vec
    target_speed = unscaled_target_speed * speedup

    return target_speed, unscaled_target_speed


def calculate_face_areas(faces):
    '''Calculate areas of faces'''
    areas = []
    for f in faces:
        areas.append(f.calc_area())
    return areas


def calculate_mean_vertex_area_change(fi_areas, initial_face_areas, \
                                      base_fis_of_vis):
    '''Calculate mean relative change of areas of faces (compared to
    initial face areas) around each vertex
    '''
    mean_changes = []
    for vi in range(len(base_fis_of_vis)):
        fis = base_fis_of_vis[vi]
        mean_change = 0.0
        for fi in fis:
            mean_change += initial_face_areas[fi] / fi_areas[fi]

        mean_change /= float(len(fis))
        mean_changes.append(mean_change)
    return mean_changes
