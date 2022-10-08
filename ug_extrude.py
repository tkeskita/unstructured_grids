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

# Operations related to extrusion of new cells

import bpy
from . import ug
from . import ug_op
import logging
l = logging.getLogger(__name__)

fulldebug = False # Set to True if you wanna see walls of logging debug
print_iterations = True # Set to True to print iteration stats for hyperbolic extrusion

# Blender coordinates are stored in single precision floats,
# so single precision must be used for zero tolerance.
# https://en.wikipedia.org/wiki/Machine_epsilon
EPS_FLOAT = 1.2e-7 # Single precision float machine epsilon
EPS = 2.0 * EPS_FLOAT # Precision with a safety margin

# Overall description of all extrusion methods
#
# Extrusion starts from selected mesh faces (base faces). Each face
# produces one new cell (in matching order and numbering). Extrusion
# is started by casting one new vertex (top vertex) from each vertex
# of base faces (base vertex) towards a normal direction calculated
# from surrounding boundary faces. Each edge of base faces produces
# one side face by connecting old and new vertices.  Finally top faces
# are added by creating faces connecting new vertices, to close each
# new cell.

# Extrusion methods
#
# Fixed Extrusion Method - uses constant vertex normal direction and
#    constant extrusion length.
# Hyperbolic Extrusion Method - an experimental and iterative method
#    which adjusts extrusion direction and extrusion length of individual
#    vertices according to surroundings.
#    Unit of extrusion "velocity" is meters per cell layer.
#    Iteration step size unit is fraction of cell layer.
#    One way to think about it: If one cell layer is extruded in one
#    second, then velocity unit would be meter per second.

# Terminology
#
# iteration = internal iteration step in hyperbolic extrusion of a cell layer
# velocity = vertex extrusion speed, which is directly proportional to
#            extrusion thickness (scalar value, in unit of m/layer)
# speed = velocity * direction normal vector (vector)

# Ideas to optimize extrusion speed further:
# - don't create internal faces into mesh
# - don't create top faces into mesh until last layer
# TODO: Autosplit twisted faces after some threshold?


class UG_OT_ExtrudeCells(bpy.types.Operator):
    '''Extrude new cells from current face selection'''
    bl_idname = "unstructured_grids.extrude_cells"
    bl_label = "Extrude Cells (UG)"

    @classmethod
    def poll(cls, context):
        return context.mode in {'OBJECT', 'EDIT_MESH'}

    def execute(self, context):
        mode = context.active_object.mode

        # Initialize from selected faces if needed
        is_ok, text, initial_ugfaces = initialize_extrusion()
        if not is_ok:
            self.report({'ERROR'}, "Initialization failed: " + text)
            return {'FINISHED'}

        # Initialize stuff
        ug_props = bpy.context.scene.ug_props
        n = 0 # new cell count
        niter = 0 # iteration counter
        speeds = [] # Extrusion speed vectors, can be updated per extruded layer
        new_ugfaces = list(initial_ugfaces) # List of created ugfaces

        import bmesh
        import time
        ob = ug.get_ug_object()
        bpy.ops.object.mode_set(mode='EDIT')
        bm = bmesh.from_edit_mesh(ob.data)  # Actual mesh
        bmt = bmesh.new()  # Trajectory mesh

        # Save initial face areas (used for scaling extrusion length)
        initial_face_areas = \
            [f.calc_area() for f in bm.faces if f.select]

        # Extrude layers
        for i in range(ug_props.extrusion_layers):
            is_last_layer = (i == (ug_props.extrusion_layers - 1))

            t0 = time.time()
            if ug_props.extrusion_method == "fixed":
                bm, nf, speeds, new_ugfaces = extrude_cells_fixed( \
                    bm, speeds, new_ugfaces)
            elif ug_props.extrusion_method == "hyperbolic":
                from . import ug_hyperbolic
                niter, bm, bmt, nf, speeds, new_ugfaces = \
                    ug_hyperbolic.extrude_cells_hyperbolic(\
                        niter, bm, bmt, speeds, new_ugfaces, \
                        initial_face_areas, is_last_layer)
            t1 = time.time()

            l.debug("Extruded layer %d, " % (i + 1) \
                    + "cells: %d " % n \
                    + "iters: %d " % niter \
                    + "(%d cells/s)" % int(nf/(t1-t0)))
            n += nf

        flip_initial_faces(bm, initial_ugfaces)
        bm.normal_update()
        bmesh.update_edit_mesh(mesh=ob.data)
        bm.free()
        recreate_trajectory_object(bmt)
        bmt.free()
        ug_op.set_faces_boundary_to_default(new_ugfaces)
        ug.update_ug_all_from_blender()
        # Return to original mode
        bpy.ops.object.mode_set(mode=mode)

        self.report({'INFO'}, "Extruded %d new cells" % n)
        return {'FINISHED'}


def initialize_extrusion():
    '''Initialize UG data for extrusion. For a new unstructured grid,
    create UG object and UGFaces from faces of active object.
    Return values are boolean for successful initialization and
    list of initial faces (if initializing from faces when no cells exist).
    '''

    initial_ugfaces = [] # List of new UGFaces

    # Do nothing if there is already an UG state
    if ug.exists_ug_state():
        return True, "UG state exists, OK to continue", initial_ugfaces

    source_ob = bpy.context.active_object
    if source_ob.name == ug.obname:
        return False, "Source object name can't be " + ug.obname, initial_ugfaces

    # Mode switch is needed to make sure mesh is saved to original object
    bpy.ops.object.mode_set(mode='OBJECT')
    #bpy.ops.object.mode_set(mode='EDIT')

    import bmesh
    # Couldn't use bmesh.from_edit_mesh() here because original mesh
    # was modified after changes when bailing out.
    bm = bmesh.new()

    # Create selected faces to bmesh
    new_verts = {}
    new_faces = []
    for p in source_ob.data.polygons:
        if not p.select:
            continue
        verts = []
        for vi in p.vertices:
            # New vert exists
            if vi in new_verts:
                verts.append(new_verts[vi])
                continue
            # Create new vert
            v = bm.verts.new(source_ob.data.vertices[vi].co)
            verts.append(v)
            new_verts[vi] = v
        new_faces.append(verts)
    bm.verts.ensure_lookup_table()
    bm.verts.index_update()

    # Create faces
    for nf in new_faces:
        f = bm.faces.new(nf)
        f.normal_update()
        f.select_set(True)
    bm.faces.ensure_lookup_table()
    bm.faces.index_update()

    # Bail out if there are no faces
    if len(bm.faces) == 0:
        return False, "No faces selected", initial_ugfaces

    # Check mesh for hanging verts
    vertlist = check_hanging_face_verts(bm)
    if vertlist:
        for f in bm.faces:
            f.select_set(False)
        for v in vertlist:
            v.select_set(True)
        bm.to_mesh(source_ob.data)
        bm.free()
        return False, "Found %d hanging vert(s)" % len(vertlist), initial_ugfaces

    ob = ug.initialize_ug_object()

    # Generate ugverts
    for i in range(len(bm.verts)):
        ug.UGVertex()
    l.debug("Initial vertex count: %d" % len(ug.ugverts))

    # Generate ugfaces
    for i in range(len(bm.faces)):
        verts_ind = [v.index for v in bm.faces[i].verts]
        ugf = ug.UGFace(verts_ind)
        ugf.add_mesh_face(i)
        for vi in verts_ind:
            ug.ugverts[vi].add_face(ugf)
        initial_ugfaces.append(ugf)
    l.debug("Initial Face count: %d" % len(ug.ugfaces))

    bm.to_mesh(ob.data)
    bm.free()

    # Hide everything else than UG object
    bpy.ops.object.mode_set(mode = 'OBJECT')
    ug.hide_other_objects()
    bpy.ops.object.mode_set(mode = 'EDIT')

    return True, "initialization done", initial_ugfaces


def recreate_trajectory_object(bm):
    '''Replace trajectory object mesh with argument mesh'''

    # Do nothing in fixed extrusion mode
    ug_props = bpy.context.scene.ug_props
    if ug_props.extrusion_method == "fixed":
        return None

    obname = "Extrusion Trajectory"
    bpy.ops.object.mode_set(mode='OBJECT')

    # Delete old object
    if obname in bpy.data.objects:
        for ob in bpy.data.objects:
            ob.select_set(False)
        bpy.data.objects[obname].select_set(True)
        mesh = bpy.data.objects[obname].data
        bpy.ops.object.delete()
        bpy.data.meshes.remove(mesh)

    # Create new object and add mesh
    if len(bm.verts) > 0:
        mesh_data = bpy.data.meshes.new(obname)
        bm.to_mesh(mesh_data)
        ob = bpy.data.objects.new(obname, mesh_data)
        bpy.context.scene.collection.objects.link(ob)
        ob.select_set(False)


def check_hanging_face_verts(bm):
    '''Check selected faces for hanging verts.
    Returns list of hanging vert coordinate vectors.
    '''

    vertlist = []
    faces = [f for f in bm.faces if f.select]
    for f in faces:
        for v in f.verts:
            if len(v.link_edges) == 2 and len(v.link_faces) == 2:
                if v not in vertlist:
                    vertlist.append(v)
    return vertlist


def extrude_cells_fixed(bm, speeds, new_ugfaces):
    '''Extrude new cells from current face selection using a simple
    extrusion method, where extrusion direction vector (speeds)
    point towards initial vertex normal direction.
    '''

    bm.verts.ensure_lookup_table()
    bm.faces.ensure_lookup_table()
    ug_props = bpy.context.scene.ug_props

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

    if len(speeds) == 0:
        # Initial speeds from vertex normal speeds
        speeds = get_vertex_normal_speeds(base_verts, base_faces, \
                                          base_fis_of_vis)

    bm, top_verts, vert_map = cast_vertices(bm, base_verts, speeds, df=1.0)
    top_faces = create_mesh_faces(bm, edge2sideface_index, vert_map, \
                                  base_faces, fis2edges, ugci0, ugfi0)
    correct_face_normals(bm, base_faces, ugci0)
    add_base_face_to_cells(base_faces, ugci0)
    thickness_update()

    return bm, len(base_faces), speeds, new_ugfaces


def get_verts_and_relations(faces):
    '''Generate face and vertex topology relationships.
    First return value is list of bmesh verts that are part of argument
    bmesh faces.
    Second return value is list of face vertex lists (to map from face to
    it's vertex indices).
    Third return value is list of vertex face lists (to map
    from vertex to it's face indices).
    Fourth return value maps face vertices to vertex index number.
    '''

    verts = [] # vertex list
    ivertmap = {} # dictionary of vert indices
    vis_of_fis = [] # vertex index list for faces (starting from 0)
    fis_of_vis = [] # face index list for verts (starting from 0)
    vi = 0 # vertex index
    fi = 0 # face index
    for f in faces:
        face_iverts = []
        for v in f.verts:
            if v not in verts:
                verts.append(v)
                ivertmap[v] = vi
                vi += 1
            face_iverts.append(ivertmap[v])
            add_entry(fis_of_vis, ivertmap[v], fi)
        vis_of_fis.append(face_iverts)
        fi += 1
    return verts, vis_of_fis, fis_of_vis, ivertmap


def get_edges_and_face_map(faces):
    '''Generate face and edge topology relationships.
    First return value is list of unique bmesh edges that are part
    of argument bmesh faces. Index number of edges is used to
    generate indices for new faces.
    Second return value is map from edges to side face index.
    Third return value is map from face index to face edges.
    '''

    edges = [] # edge list
    edge2sideface_index = dict()
    fis2edges = [] # face index to edges list of lists
    fi = 0 # side face index
    for f in faces:
        face_edges = []
        for e in f.edges:
            if e not in edges:
                edges.append(e)
                edge2sideface_index[e] = fi
                fi += 1
            face_edges.append(e)
        fis2edges.append(face_edges)
    return edges, edge2sideface_index, fis2edges


def create_UG_verts_faces_and_cells(base_verts, base_edges, base_faces, new_ugfaces):
    '''Create new bare bone UGCells, UGFaces and UGVerts for extrusion of
    a single cell layer
    '''

    # Vertices
    for i in range(len(base_verts)):
        ug.UGVertex()

    # Faces (sides and top)
    for i in range(len(base_edges) + len(base_faces)):
        ugf = ug.UGFace()
        new_ugfaces.append(ugf)

    # Cells
    for i in range(len(base_faces)):
        ug.UGCell()

    return new_ugfaces


def get_face_vis_of_vi(verts, fis_of_vis, vis_of_fis):
    '''Find surrounding face vertices for all vertices.
    First return value is list which contains, for each vertex, a list
    of vertex indices of neighbour faces excluding the base vertex.
    Second return value is list of face indices for the same vertices.
    '''
    neighbour_vis_of_vi = [] # neighbour vertex indices
    fils_of_neighbour_vis = [] # face index lists of neighbour vis
    for vi in range(len(verts)):
        vis = []
        fis = []
        for fi in fis_of_vis[vi]:
            for i in vis_of_fis[fi]:
                if i == vi:
                    continue # Skip self
                if i in vis:
                    # Vertex index is already in vis list, add fi
                    ind = vis.index(i)
                    fis[ind].append(fi)
                else:
                    # New vertex, initalize fils with this index
                    vis.append(i)
                    fis.append([fi])
        neighbour_vis_of_vi.append(vis)
        fils_of_neighbour_vis.append(fis)
    return neighbour_vis_of_vi, fils_of_neighbour_vis


def get_face_vertex_cos_angle(v, f):
    '''Calculate cos(angle) between the two edges connected to vertex v,
    where both edges are part of face f
    '''

    # Find edges
    edges = []
    for e in v.link_edges:
        if e not in f.edges:
            continue
        edges.append(e)

    # Generate edge vectors and calculate cos(angle)
    vec1 = edges[0].other_vert(v).co - v.co
    vec1.normalize()
    vec2 = edges[1].other_vert(v).co - v.co
    vec2.normalize()
    cos_angle = vec1 @ vec2
    return min(1.0, max(-1.0, cos_angle))


def get_items_from_list(alist, ilist):
    '''Return items in alist located at index locations in ilist'''
    res = []
    for i in ilist:
        res.append(alist[i])
    return res


def get_vertex_normal_speeds(verts, faces, fis_of_vis):
    '''Calculate minimum vertex normal speeds by averaging face
    normals surrounding each vertex, weighted by face vertex angle
    (the angle between the two edges of the face connected at vertex).
    '''

    from mathutils import Vector
    import math
    ug_props = bpy.context.scene.ug_props

    # Layer thickness
    ext_len = ug_props.extrusion_thickness

    # Calculate new extrusion directions
    speeds = []
    for vi in range(len(verts)):
        angle_factors = []
        neigh_faces = get_items_from_list(faces, fis_of_vis[vi])

        # Calculate weight coefficients from face vertex angles
        for f in neigh_faces:
            # Convert cos_angle to angle. Scale to [0, 1] for weight factor.
            cos_angle = get_face_vertex_cos_angle(verts[vi], f)
            angle_factor = math.acos(cos_angle) / math.pi
            angle_factors.append(angle_factor)

        # Calculate final speeds from normals and weights
        norvecs = [f.normal for f in neigh_faces]
        norvec = Vector((0, 0, 0))
        for nv, factor in zip(norvecs, angle_factors):
            norvec += factor * nv
        norvec.normalize()

        norvec *= ext_len # Set velocity to minimum velocity
        speeds.append(norvec)

    return speeds


def calculate_convexity_sums(bm, verts, faces):
    '''Calculate sum of convexities for verts.
    Convexity is a measure of face-face angles:
    Convexity value near 0 means extremely sharp concave angle.
    Convexity value 0.25 means 90 degree concave angle.
    Convexity calue 0.5 means flat terrain (neither concave or convex).
    Convexity value 0.75 means 90 degree convex angle.
    Convexity value near 1 means extremely convex hole.
    Convexity sum is the sum of (convexity - 0.5) > 0 values
    calculated for all convex edges surrounding each vertex.
    '''

    # Calculate convexity sums
    convexity_sums = []

    for vi, v in enumerate(verts):
        convexity_sum = 0.0
        if fulldebug: l.debug("Vert %d:" % v.index)
        for e in v.link_edges:
            efaces = [f for f in e.link_faces if f in faces]
            if len(efaces) == 2:
                if fulldebug: l.debug("  edge %s: " % str(e))
                f0 = efaces[0]
                f1 = efaces[1]

                # face-face-edge angle:
                cos_epsilon = face_face_cos_angle(e, f0, f1)
                if fulldebug: l.debug("  cos_epsilon = %f" % cos_epsilon)

                # Angle between face center to face center and first face normal
                vec_f2f = f0.calc_center_median() - f1.calc_center_median()
                vec_f2f.normalize()
                cos_theta = vec_f2f @ f0.normal
                if fulldebug: l.debug("  cos_theta = %f" % cos_theta)

                if cos_theta <= 0.0:
                    convexity = (cos_epsilon + 1.0) / 4.0 + 0.5
                else:
                    convexity = (-cos_epsilon + 1.0) / 4.0
                if fulldebug: l.debug("  convexity = %f" % convexity)

                if convexity > 0.5:
                    convexity_sum += (convexity - 0.5)
        convexity_sums.append(convexity_sum)
    return convexity_sums


def face_face_cos_angle(e, f, f_neighbor):
    '''Calculate cosine of angle between two connected bmesh faces, which
    share a common edge e. Return None if edge does not connect two
    faces.
    '''

    # Make sure edge connects given two faces before continuing
    etest = [etest for etest in f.edges if etest in f_neighbor.edges]
    if e not in etest:
        return None

    # face center coordinates
    f_center = f.calc_center_median()
    f_neighbor_center = f_neighbor.calc_center_median()
    return edge_vec_vec_cos_angle(e, f_center, f_neighbor_center)


def bmesh_edge_center(edge):
    '''Calculates coordinates for edge center using edge vertex
    coordinates
    '''
    return (edge.verts[0].co + edge.verts[1].co) / 2


def edge_vec_vec_cos_angle(e, vec1, vec2):
    """Calculates cosine of angle between two vectors which are
    orthogonal projections for edge e. This is used to calculate
    cos(angle) for two faces which share edge e, and vec1 and vec2
    represent the faces' center coordinate vectors.
    """

    # Calculate cos(epsilon) where epsilon is the angle between the two
    # faces connected by edge e. cos(epsilon) is mathematically angle between
    # vectors e_ortho_f -> f_center and e_ortho_f_neighbor -> f_neighbor_center.
    # e_ortho_f is point along edge e which forms 90 degree angle between
    # edge point 1 -> edge point 2 and e_ortho_f -> f_center.
    # Similarly e_ortho_f_neighbor is point along edge e which forms 90 degree
    # angle between edge point 1 -> edge point 2 and
    # e_ortho_f_neighbor -> f_neighbor_center.
    #
    # Diagram example with two triangular faces: + marks face boundary edges,
    # x = e_ortho_f, y = e_ortho_f_neighbor, c = edge or face center
    #
    #      +++
    #     +   +++++
    #    +  c       +++++            <-- f
    #   +   I            +++++++
    #  1----x------c---y--------2    <-- edge e
    #   +++++          I      ++
    #        +++++     c   +++       <-- f_neighbor
    #             ++++   ++
    #                 +++

    e_center = bmesh_edge_center(e) # edge center
    vec_e = e.verts[1].co - e.verts[0].co # edge vector

    # Calculate orthogonal vector e_ortho_f -> f_center and normalize it
    f_center = vec1 # face center coordinates
    vec_f_center = f_center - e_center
    project_f = vec_f_center.project(vec_e) # project vec_f to vec_e
    e_ortho_f = e_center + project_f # coordinates for x
    vec_ortho_f = f_center - e_ortho_f # orthogonal vector
    vec_ortho_f.normalize() # normalize it

    # Similarly to above, calculate orthogonal vector
    # e_ortho_f_neighbor -> f_neighbor_center and normalize it
    f_neighbor_center = vec2
    vec_f_neighbor_center = f_neighbor_center - e_center
    project_f_neighbor = vec_f_neighbor_center.project(vec_e)
    e_ortho_f_neighbor = e_center + project_f_neighbor
    vec_ortho_f_neighbor = f_neighbor_center - e_ortho_f_neighbor
    vec_ortho_f_neighbor.normalize()

    # Finally calculate cos(angle) between faces
    cos_epsilon = vec_ortho_f @ vec_ortho_f_neighbor
    # Limit -1.0 <= cos_epsilon <= 1.0 for physical correctness
    cos_epsilon = max(-1.0, min(1.0, cos_epsilon))
    return cos_epsilon


def classify_vert(v, ibasevertmap, edges, vfaces):
    '''Classify argument vertex v.
    First return value is True if v is corner vertex.
    Second return value is True if v is boundary vertex.
    Third return value is list of two vertex indices (vis) pointing to
    neighbour boundary vertices (only for boundary vertices).
    Fourth return value is list of all neighbour vertices,
    in face connection order.
    '''
    bn_vis = [] # boundary neighbour vis (vertex indices)
    bn_faces = [] # boundary neighbour faces
    an_vis = [] # all neighbour vis (connected by edges)
    processed_faces = [] # list of traversed faces
    processed_edges = [] # list of traversed edges

    # Find an edge to start from. Must be a boundary edge for
    # boundary verts, for face traversal to work correctly.
    v_edges = [x for x in v.link_edges if x in edges]
    for e in v_edges:
        e_faces = [x for x in e.link_faces if x in vfaces]
        if len(e_faces) == 1:
            break
    if e not in edges:
        raise ValueError("error edge %d" % e.index)

    # Go through all edges in face connection order
    while e not in processed_edges and v in e.verts:
        # Process new edge
        an_vis.append(ibasevertmap[e.other_vert(v)])
        link_faces = [x for x in e.link_faces if x in vfaces]
        if len(link_faces) == 1:
            bn_vis.append(ibasevertmap[e.other_vert(v)])
            bn_faces.append(link_faces[0])
        processed_edges.append(e)

        # Get a connected unprocessed vface
        for f in link_faces:
            if f not in processed_faces:
                break
        processed_faces.append(f)

        # Get the next connected edge
        for e in f.edges:
            if e in processed_edges:
                continue
            if v not in e.verts:
                continue
            break

    # Return classification information
    if len(bn_vis) == 2:
        if bn_faces[0] == bn_faces[1]:
            return True, True, [bn_vis[0], bn_vis[1]], an_vis
        else:
            return False, True, [bn_vis[0], bn_vis[1]], an_vis
    return False, False, [None, None], an_vis


def classify_verts(verts, edges, faces, fis_of_vis, ibasevertmap):
    '''Generate vertex classification information.
    First return value is True for corner vertices.
    Second return value is True for boundary vertices.
    Third return value is list of two neighbour vertex indices for
    boundary vertices (or Nones if vertex is not boundary vertex).
    Fourth return value is list of neighbour vertices connected by
    edges at vertex.
    '''

    # Turn lists into sets for fast testing
    edgeset = set(edges)
    is_corners = []
    is_boundaries = []
    bn_vis = [] # boundary vertex neighbours
    an_vis = [] # all vertex neighbours
    for vi, v in enumerate(verts):
        vfaces = [faces[x] for x in fis_of_vis[vi]] # faces around vertex
        is_corner, is_boundary, bn_verts, an_verts = \
            classify_vert(v, ibasevertmap, edgeset, vfaces)
        is_corners.append(is_corner)
        is_boundaries.append(is_boundary)
        bn_vis.append(bn_verts)
        an_vis.append(an_verts)
    return is_corners, is_boundaries, bn_vis, an_vis


def cast_vertices(bm, base_verts, speeds, df):
    '''Create new top vertices from base vertices in argument bmesh, by
    casting each vertex towards speeds (length elen). Returns updated bmesh,
    top vertices and map from base verts to new top verts.
    '''

    vert_map = {} # Dictionary for mapping original face verts to new verts
    ug_props = bpy.context.scene.ug_props
    top_verts = []

    # Cast new vertices very small distance ahead
    for i in range(len(base_verts)):
        if ug_props.extrusion_method == "fixed":
            newco = base_verts[i].co + speeds[i]
        else:
            newco = base_verts[i].co + df * speeds[i]
        v2 = bm.verts.new(newco)
        vert_map[base_verts[i]] = v2
        top_verts.append(v2)

    bm.verts.ensure_lookup_table()
    bm.verts.index_update()

    if fulldebug: l.debug("Cast %d vertices" % len(speeds))
    return bm, top_verts, vert_map


def create_mesh_faces(bm, edge2sideface_index, vert_map, base_faces, \
                      fis2edges, ugci0, ugfi0):
    '''Create bmesh side and top faces, and link mesh face with
    UGFace. Return top faces.
    '''

    # Face creation help functions

    ugfi = ugfi0 # UGFace index
    fi = len(bm.faces) # mesh face index

    # Create side faces one cell at a time
    processed_edges = []
    for ci in range(len(base_faces)):
        ugci = ugci0 + ci # UGCell index
        for e in fis2edges[ci]:
            if e not in processed_edges:
                # New face (boundary or internal)
                create_side_face_from_edge(bm, e, vert_map, ugci, fi, ugfi)
                processed_edges.append(e)
                fi += 1
                ugfi += 1
            else:
                # Existing face (internal). Add existing UGFace to UGCell
                ugf = ug.ugfaces[ugfi0 + edge2sideface_index[e]]
                c = ug.ugcells[ugci]
                c.add_face_and_verts(ugf)
                ugf.neighbour = c

    # Create top faces
    top_faces = []
    for ci in range(len(base_faces)):
        ftop = create_top_face_from_base_face(bm, base_faces[ci], vert_map, \
                                              ugci0 + ci, fi, ugfi)
        top_faces.append(ftop)
        fi += 1
        ugfi += 1

    return top_faces


def create_face_from_verts(bm, verts, ugci, fi, ugfi):
    '''Creates BMFace from verts and adds it to cell index ugci.
    fi is mesh face index and ugfi UGFace index for this new face.
    '''
    f = bm.faces.new(verts)
    f.normal_update()
    ugf = ug.ugfaces[ugfi]

    # Add vertices of mesh face f to ugface index ugfi
    for v in f.verts:
        ugv = ug.ugverts[v.index]
        ugf.add_verts([ugv])
        ugv.add_face(ugf)

    # Link UGFace and mesh face in UG data
    ugf.add_mesh_face(fi)

    # Add UGFace to UGCell
    c = ug.ugcells[ugci]
    c.add_face_and_verts(ugf)
    ugf.owner = c
    return f


def create_side_face_from_edge(bm, e, vert_map, ugci, fi, ugfi):
    '''Creates a face from edge e'''
    # Generate vertex list for face creation
    e0 = e.verts[0]
    e1 = e.verts[1]
    verts = [e0, vert_map[e0], vert_map[e1], e1]
    create_face_from_verts(bm, verts, ugci, fi, ugfi)


def create_top_face_from_base_face(bm, f, vert_map, ugci, fi, ugfi):
    '''Creates a top face from base face f'''
    verts = []
    for i in f.verts:
        verts.append(vert_map[i])
    ftop = create_face_from_verts(bm, verts, ugci, fi, ugfi)

    # Deselect base face and select top face
    f.select_set(False)
    ftop.select_set(True)
    return ftop


def correct_face_normals(bm, faces, ugci0):
    '''Ensure that face normal directions point out of cells'''

    # Note: For internal faces normal should point from lower index cell
    # to higher index cell. That is achieved by the reverse loop below.

    bm.faces.ensure_lookup_table()
    bm.faces.index_update()

    # Loop through cell faces in reverse and flip
    # face normals that point inside
    for i in range(len(faces)-1, -1, -1):
        for ugf in ug.ugcells[ugci0 + i].ugfaces:
            # Flip face normal if face normal points inside of new cell.
            # Uses face center to base face center as reference vector.
            # Cells are required to be convex for this to work correctly.
            f = bm.faces[ugf.bi]
            facevec = f.calc_center_median()
            refvec = faces[i].calc_center_median() - facevec
            refvec.normalize()
            cos_epsilon = f.normal @ refvec
            if fulldebug:
                l.debug("f.normal:%s, refvec:%s" %(str(f.normal), str(refvec)))
                l.debug("cos_epsilon = %f" % cos_epsilon)
            if (cos_epsilon > 0.0):
                f.normal_flip()
                f.normal_update()
                ugf.invert_face_dir(False) # Flip vertex list only
                if fulldebug: l.debug("Flipped mesh face %d normal" % f.index)


def calculate_intvpairs(anvis_of_vi, fis_of_vis, is_boundaries):
    '''Calculate index list of vpairs indices for verts. vpairs are
    neighbour vertex pairs of a vertex, which don't share
    faces. First return value is list of lists of first vertices,
    and second return value similar list of lists of paired
    vertices.
    '''

    intvis = [] # Internal vertex indices
    intvpairs0 = [] # Internal vertex indices, first of pair
    intvpairs1 = [] # Internal vertex indices, second of pair
    for vi in range(len(anvis_of_vi)):
        # Boundary vertices are not needed
        if is_boundaries[vi]:
            intvpairs0.append([])
            intvpairs1.append([])
            continue

        # Internal vertex found
        intvis.append(vi)

        # Build lists of vpairs for vertex
        ivp0 = []
        ivp1 = []
        anvis = anvis_of_vi[vi]
        for i in range(len(anvis)):
            for j in range(i + 1, len(anvis)):
                fsi = set(fis_of_vis[anvis[i]])
                fsj = set(fis_of_vis[anvis[j]])

                # If any face is common to both sets skip this combination
                if len(fsi.intersection(fsj)) > 0:
                    continue
                # Otherwise, vpair has been found, add 'em
                ivp0.append(anvis[i])
                ivp1.append(anvis[j])
        intvpairs0.append(ivp0)
        intvpairs1.append(ivp1)
    return intvis, intvpairs0, intvpairs1


def update_trajectory_mesh(bmt, cos):
    '''Add coordinates to trajectory mesh'''

    bmt.verts.ensure_lookup_table()
    bmt.verts.index_update()
    nverts = len(cos)
    nv0 = len(bmt.verts) - nverts
    oldverts = [bmt.verts[x] for x in range(nv0, nv0 + nverts)]
    for i,co in enumerate(cos):
        nv = bmt.verts.new(co)
        v0 = oldverts[i]
        bmt.edges.new([v0, nv])
    return bmt


def add_faces_to_trajectory_mesh(bmt, nv0, vis_of_fis):
    '''Add final faces to trajectory mesh'''

    bmt.verts.ensure_lookup_table()
    bmt.verts.index_update()
    for vis in vis_of_fis:
        verts = [bmt.verts[nv0 + vi] for vi in vis]
        bmt.faces.new(verts)
    return bmt


def add_base_face_to_cells(faces, ugci0):
    '''Add base faces to cells as neighbours'''

    for i in range(len(faces)):
        # Add existing UGFace to UGCell
        ugf = ug.facemap[faces[i].index]
        ugci = ugci0 + i
        c = ug.ugcells[ugci]
        c.add_face_and_verts(ugf)
        ugf.neighbour = c


def thickness_update():
    '''Update layer thickness using expression specified by user'''

    ug_props = bpy.context.scene.ug_props
    x = ug_props.extrusion_thickness
    expr = ug_props.extrusion_scale_thickness_expression
    try:
        rval = eval(expr)
        if fulldebug: l.debug("Expression returned %s" % str(rval))
        ug_props.extrusion_thickness = float(rval)
    except:
        l.error("Error in evaluating: %r" % expr)
    # TODO: Add error notification to user if expression fails


def flip_initial_faces(bm, initial_ugfaces):
    '''Flip the normal direction of initial faces in bmesh'''

    bm.faces.ensure_lookup_table()
    for ugf in initial_ugfaces:
        if fulldebug: l.debug("Final flipping face %d" % ugf.bi)
        bm.faces[ugf.bi].normal_flip()
        bm.faces[ugf.bi].normal_update()
        ugf.invert_face_dir()


##### Generic help functions #####


def add_entry(lol, i, val):
    '''Add entry val to list with index i into list of lists lol'''
    if i == len(lol):
        lol.append([])
    if i < len(lol):
        lol[i].append(val)
        return lol
    else:
        raise ValueError("Illegal index %d, list length %d" % (i, len(lol)))


