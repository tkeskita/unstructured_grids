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

# Operations related to shell extrusion of new cells

import bpy
from . import ug
from . import ug_op
from .ug_extrude import *
import logging
l = logging.getLogger(__name__)
ic_ob_name = 'Interactive Correction'  # Name for interactive correction object


def extrude_cells_shell(n, niter, bm, bmt, speeds, new_ugfaces, \
                        is_last_layer):
    '''Extrude new cells from current face selection using shell extrusion
    method
    '''

    import bmesh
    info_text = ""
    ug_props = bpy.context.scene.ug_props

    # Selected faces are base faces for extrusion. New cell index number is
    # index number of mesh face in faces list. Actual UGCell index
    # number is cell number plus initial number of prior ugcells.
    base_faces = [f for f in bm.faces if f.select]
    if fulldebug: l.debug("Face count at beginning: %d" % len(base_faces))

    # Get topology information from mesh
    base_verts, base_vis_of_fis, base_fis_of_vis, ibasevertmap = \
        get_verts_and_relations(base_faces)
    base_edges, edge2sideface_index, fis2edges = get_edges_and_face_map(base_faces)
    neighbour_vis_of_vi, fils_of_neighbour_vis = \
        get_face_vis_of_vi(base_verts, base_fis_of_vis, base_vis_of_fis)
    is_corners, is_boundaries, bnvis_of_vi, anvis_of_vi = \
        classify_verts(base_verts, base_edges, base_faces, \
                       base_fis_of_vis, ibasevertmap)

    ugci0 = len(ug.ugcells) # Number of UGCells before extruding
    ugfi0 = len(ug.ugfaces) # Number of UGFaces before extruding
    new_ugfaces = create_UG_verts_faces_and_cells(base_verts, base_edges, base_faces, new_ugfaces)

    if ug_props.interactive_correction_mode and len(speeds) == 0:
        speeds = get_speeds_from_interactive_correction(base_verts)
        delete_object(ic_ob_name)
    else:
        if len(speeds) == 0:
            speeds = get_shell_speeds( \
                bm, base_verts, base_faces, base_fis_of_vis, neighbour_vis_of_vi,
                is_corners, is_boundaries
            )
            if ug_props.shell_ensure_thickness:
                speeds = adjust_speeds(bm, base_verts, speeds)

    df = calculate_layer_thickness(n)
    bm, top_verts, vert_map = cast_vertices(bm, base_verts, speeds, df=df)
    top_faces = create_mesh_faces(bm, edge2sideface_index, vert_map, \
                                  base_faces, fis2edges, ugci0, ugfi0)
    correct_face_normals(bm, base_faces, ugci0)
    add_base_face_to_cells(base_faces, ugci0)

    # Trajectory bmesh
    if ug_props.extrusion_create_trajectory_object and len(bmt.verts) == 0:
        for v, speed in zip(base_verts, speeds):
            v0 = bmt.verts.new(v.co)
            v1 = bmt.verts.new(v.co + speed * ug_props.extrusion_thickness)
            bmt.edges.new([v0, v1])
        bmt.verts.ensure_lookup_table()
        bmt.verts.index_update()

    return niter + 1, bm, bmt, len(base_faces), speeds, new_ugfaces, info_text


def get_shell_speeds(bm, base_verts, base_faces, base_fis_of_vis, \
                     neighbour_vis_of_vi, is_corners, is_boundaries):
    '''Calculate vertex speeds (direction + velocity) to be applied in
    shell extrusion. Vertex speeds are obtained iteratively, starting
    from vertex normal direction. The direction is iteratively updated
    from neighbor vertices. Vertex convexity is used as a weight when
    direction is determined.
    '''

    speeds = get_vertex_normal_speeds(base_verts, base_faces, base_fis_of_vis)
    convexity_sums = calculate_convexity_sums(bm, base_verts, base_faces)
    ug_props = bpy.context.scene.ug_props

    # Weight calculation
    minimum_weight = 0.1
    max_weight = 2.0
    weights = []  # Weight factors for verts
    for c in convexity_sums:
        if c > minimum_weight:
            weights.append(min(c, max_weight))
        else:
            weights.append(minimum_weight)
    for i in range(len(base_verts)):
        if is_corners[i] or is_boundaries[i]:
            weights[i] = max_weight

    niter = 12  # Number of direction propagation iterations
    nw_factor = 0.2  # New weight increase factor

    for iter in range(niter):
        l.debug("Direction propagation iteration %d" % iter)
        new_speeds = list()
        new_weights = list()

        # Calculate new speed
        for i in range(len(base_verts)):
            # Corners and boundaries use initial speed
            if is_corners[i] or is_boundaries[i]:
                new_speeds.append(speeds[i])
                new_weights.append(weights[i])
                continue

            # New speed calculation
            own_w = weights[i]
            own_weight_factor = (1.0 + own_w)**6.0 - 1.0 + own_w
            new_speed = own_weight_factor * speeds[i]
            max_nw = 0.0  # Maximum neighbour weight
            # Neighbour contribution
            for nvi in neighbour_vis_of_vi[i]:
                w = weights[nvi]
                new_speed += w * speeds[nvi]
                max_nw = max(max_nw, w)
            # Final scaling
            new_speed.normalize()
            new_speeds.append(new_speed)

            # Calculate new weight based on maximum neigbour weight
            nw = nw_factor * max_nw + (1.0 - nw_factor) * own_w
            nw = max(own_w, nw)
            new_weights.append(nw)

        speeds = list(new_speeds)
        weights = list(new_weights)

    return speeds


def adjust_speeds(bm, base_verts, speeds):
    '''Scale up speed vectors to ensure layer thickness
    '''

    ug_props = bpy.context.scene.ug_props
    thickness = ug_props.extrusion_thickness
    scaled_speeds = []
    from mathutils.bvhtree import BVHTree
    from .ug_extrude import EPS
    bt = BVHTree.FromBMesh(bm)

    # For each top point, cast rays towards inverse neighbor face
    # normal direction and find the smallest ray length, which is used
    # to scale the speed.
    for v, speed in zip(base_verts, speeds):
        min_len = 1.0e+38  # Practical single precision float max
        co = v.co + speed * thickness
        for ray_dir in [-1.0 * f.normal for f in v.link_faces]:
            hit_co, hit_nor, hit_index, hit_length = \
                bt.ray_cast(co, ray_dir)
            if not hit_co:
                continue
            ray_len = (co - hit_co).length
            if ray_len < min_len:
                min_len = ray_len
        if min_len < EPS:
            # min_len = EPS
            raise Exception("Internal error: min_len < EPS")
        if min_len < thickness:
            scale = thickness/min_len
            scale = min(scale, 2.0)  # Should this be an option?
            scaled_speeds.append(scale * speed)
        else:
            scaled_speeds.append(speed)

    return scaled_speeds


def interactively_correct_speeds(bm, base_verts, base_faces, speeds):
    '''Create a temporary Interactive Correction object to allow user to
    manually modify top face vertex locations before extruding cells.
    '''

    import bmesh
    ug_props = bpy.context.scene.ug_props
    info_text = ""

    # Create top vertices with edges, and top faces
    vertmap = {}
    top_verts = []
    for i, bv in enumerate(base_verts):
        v = bm.verts.new(bv.co + speeds[i])
        top_verts.append(v)
        vertmap[bv] = v
        bm.edges.new([bv, v])
    bm.verts.ensure_lookup_table()
    bm.verts.index_update()
    for i, bf in enumerate(base_faces):
        verts = [vertmap[bv] for bv in bf.verts]
        bm.faces.new(verts)
    bm.faces.ensure_lookup_table()
    bm.faces.index_update()

    # TODO: This intersection search is missing side faces
    # TODO: Code duplication
    intersecting_verts = []
    if ug_props.check_for_intersections:
        intersecting_verts = check_for_intersections(bm, top_verts)
        if intersecting_verts:
            bm.select_mode = {'VERT'}
            for v in bm.verts:
                v.select = False
            for v in intersecting_verts:
                v.select = True
            info_text = "WARNING: Highlighted %d intersecting vertices." % len(intersecting_verts)
        else:
            info_text = "No intersections detected."
        bm.select_flush_mode()

    # Select top verts only
    if not intersecting_verts:
        bm.select_mode = {'VERT'}
        for v in bm.verts:
            v.select = False
        for v in top_verts:
            v.select = True
        bm.select_flush_mode()


    # # FIXME: Highlight only extrusion side edges
    # bm.select_mode = {'VERT'}
    # for v in bm.verts:
    #     v.select = False
    # bm.select_flush_mode()
    # bm.edges.ensure_lookup_table()
    # bm.edges.index_update()
    # bm.select_mode = {'EDGE'}
    # bm.edges.ensure_lookup_table()
    # for i in range(len(base_verts)):
    #     bm.edges[i].select = True
    # bm.select_flush_mode()

    # Delete possibly existing old Interactive Correction object and
    # hide other objects
    bpy.ops.object.mode_set(mode='OBJECT')
    bpy.ops.object.select_all(action='DESELECT')
    if ic_ob_name in bpy.data.objects:
        bpy.data.objects[ic_ob_name].select_set(True)
        mesh = bpy.data.objects[ic_ob_name].data
        bpy.ops.object.delete()
        bpy.data.meshes.remove(mesh)
    for ob in bpy.data.objects:
        ob.hide_set(True)

    # Create new Interactive correction object and make it active
    mesh_data = bpy.data.meshes.new(ic_ob_name)
    bm.to_mesh(mesh_data)
    bm.free()
    ob = bpy.data.objects.new(ic_ob_name, mesh_data)
    bpy.context.scene.collection.objects.link(ob)
    bpy.context.view_layer.objects.active = bpy.data.objects[ic_ob_name]
    ob.select_set(True)
    ob.hide_set(False)
    bpy.ops.object.mode_set(mode='EDIT')
    bpy.ops.mesh.select_mode(type='VERT')
    return info_text


def prepare_interactive_correction(source_ob):
    '''Top level routine for creating Interactive Correction object'''

    import bmesh
    ug_props = bpy.context.scene.ug_props
    bm = create_bmesh_from_selection(source_ob)

    base_faces = [f for f in bm.faces if f.select]
    if fulldebug: l.debug("Face count at beginning: %d" % len(base_faces))

    # Get topology information from mesh
    base_verts, base_vis_of_fis, base_fis_of_vis, ibasevertmap = \
        get_verts_and_relations(base_faces)
    base_edges, edge2sideface_index, fis2edges = get_edges_and_face_map(base_faces)
    neighbour_vis_of_vi, fils_of_neighbour_vis = \
        get_face_vis_of_vi(base_verts, base_fis_of_vis, base_vis_of_fis)
    is_corners, is_boundaries, bnvis_of_vi, anvis_of_vi = \
        classify_verts(base_verts, base_edges, base_faces, \
                       base_fis_of_vis, ibasevertmap)

    speeds = get_shell_speeds( \
        bm, base_verts, base_faces, base_fis_of_vis, neighbour_vis_of_vi,
        is_corners, is_boundaries
    )
    if ug_props.shell_ensure_thickness:
        speeds = adjust_speeds(bm, base_verts, speeds)

    # Scale by extrusion thickness to show end result correctly
    speeds = [x * ug_props.extrusion_thickness for x in speeds]

    info_text = interactively_correct_speeds(bm, base_verts, base_faces, speeds)
    ug_props.interactive_correction = False
    info_text += " Interactive Correction Object created for editing. "
    return info_text


def get_speeds_from_interactive_correction(base_verts):
    '''Update top vertex positions from Interactive Correction object'''

    ob = bpy.data.objects[ic_ob_name]
    speeds = []
    nv = len(base_verts)
    for i, v in enumerate(base_verts):
        top_vert = ob.data.vertices[i - nv]
        speeds.append(top_vert.co - v.co)
    # Divide by extrusion thickness to scale back to normalized speeds
    ug_props = bpy.context.scene.ug_props
    speeds = [x / ug_props.extrusion_thickness for x in speeds]
    return speeds


def delete_object(ob_name):
    '''Delete an object and it's mesh'''

    bpy.ops.object.select_all(action='DESELECT')
    bpy.data.objects[ob_name].select_set(True)
    mesh = bpy.data.objects[ob_name].data
    bpy.ops.object.delete()
    bpy.data.meshes.remove(mesh)
