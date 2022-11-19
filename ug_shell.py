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

ug_props = bpy.context.scene.ug_props

def extrude_cells_shell(niter, bm, bmt, speeds, new_ugfaces, \
                             is_last_layer):
    '''Extrude new cells from current face selection using shell extrusion
    method
    '''

    import bmesh
    bm.verts.ensure_lookup_table()
    bm.faces.ensure_lookup_table()
    info_text = ""

    # Selected faces are base faces for extrusion. New cell index number is
    # index number of mesh face in faces list. Actual UGCell index
    # number is cell number plus initial number of prior ugcells.
    base_faces = [f for f in bm.faces if f.select]
    if fulldebug: l.debug("Face count at beginning: %d" % len(base_faces))

    ugci0 = len(ug.ugcells) # Number of UGCells before extruding
    ugfi0 = len(ug.ugfaces) # Number of UGFaces before extruding

    base_verts, base_vis_of_fis, base_fis_of_vis, ibasevertmap = \
        get_verts_and_relations(base_faces)
    base_edges, edge2sideface_index, fis2edges = get_edges_and_face_map(base_faces)
    new_ugfaces = create_UG_verts_faces_and_cells(base_verts, base_edges, base_faces, new_ugfaces)
    neighbour_vis_of_vi, fils_of_neighbour_vis = \
        get_face_vis_of_vi(base_verts, base_fis_of_vis, base_vis_of_fis)
    is_corners, is_boundaries, bnvis_of_vi, anvis_of_vi = \
        classify_verts(base_verts, base_edges, base_faces, \
                       base_fis_of_vis, ibasevertmap)

    if len(speeds) == 0:
        speeds = get_shell_speeds( \
            bm, base_verts, base_faces, base_fis_of_vis, neighbour_vis_of_vi,
            is_corners, is_boundaries
        )
    else:
        speeds = scale_speeds(speeds)
    if ug_props.shell_ensure_thickness:
        speeds = adjust_speeds(bm, base_verts, speeds)

    bm, top_verts, vert_map = cast_vertices(bm, base_verts, speeds, df=1.0)
    top_faces = create_mesh_faces(bm, edge2sideface_index, vert_map, \
                                  base_faces, fis2edges, ugci0, ugfi0)
    correct_face_normals(bm, base_faces, ugci0)
    add_base_face_to_cells(base_faces, ugci0)
    thickness_update()

    if ug_props.check_for_intersections:
        intersecting_verts = check_for_intersections(bm, top_verts)
        if intersecting_verts:
            bm.select_mode = {'VERT'}
            for f in bm.verts:
                f.select = False
            for v in intersecting_verts:
                v.select = True
            bm.select_flush_mode()
            info_text = "WARNING: Highlighted %d intersecting vertices." % len(intersecting_verts)

    # Trajectory bmesh
    if ug_props.extrusion_create_trajectory_object and len(bmt.verts) == 0:
        for v, speed in zip(base_verts, speeds):
            v0 = bmt.verts.new(v.co)
            v1 = bmt.verts.new(v.co + speed)
            bmt.edges.new([v0, v1])
        bmt.verts.ensure_lookup_table()
        bmt.verts.index_update()

    return niter, bm, bmt, len(base_faces), speeds, new_ugfaces, info_text


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
    ext_len = ug_props.extrusion_thickness

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
            new_speed *= ext_len  # Set velocity to minimum velocity
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

    thickness = ug_props.extrusion_thickness
    scaled_speeds = []
    from mathutils.bvhtree import BVHTree
    bt = BVHTree.FromBMesh(bm)

    # For each top point, cast rays towards inverse neighbor face
    # normal direction and find the smallest ray length, which is used
    # to scale the speed.
    for v, speed in zip(base_verts, speeds):
        min_len = 1.0e+38  # Practical single precision float max
        co = v.co + speed
        for ray_dir in [-1.0 * f.normal for f in v.link_faces]:
            hit_co, hit_nor, hit_index, hit_length = \
                bt.ray_cast(co, ray_dir)
            if not hit_co:
                continue
            ray_len = (co - hit_co).length
            if ray_len < min_len:
                min_len = ray_len

        if min_len < thickness:
            scale = thickness/min_len
            scale = min(scale, 2.0)  # Should this be an option?
            scaled_speeds.append(scale * speed)
        else:
            scaled_speeds.append(speed)

    return scaled_speeds
