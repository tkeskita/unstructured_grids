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
LARGE = 100.0 # Cut-off for large values for weight function

# Blender coordinates are stored in single precision floats,
# so single precision must be used for zero tolerance.
# https://en.wikipedia.org/wiki/Machine_epsilon
EPS_FLOAT = 1.2e-7 # Single precision float machine epsilon
EPS = 2.0 * EPS_FLOAT # Precision with a safety margin

# Summary of extrusion method: Extrusion starts from selected mesh
# faces (base faces). Each face produces one new cell (in matching
# order and numbering). Extrusion is started by casting one new vertex
# (top vertex) from each vertex of base faces (base vertex) towards a
# normal direction calculated from surrounding boundary faces. Each edge
# of base faces produces one side face by connecting old and new vertices.
# Finally top faces are added by creating faces connecting new vertices,
# to close each new cell.
#
# There are two extrustion modes:
# 1. Fixed Extrusion Method - uses vertex normal direction and
#    extrusion velocities constant. Simple, mostly uncontrollable and
#    relatively fast method.
# 2. Hyperbolic Extrusion Method - adjusts extrusion direction and
#    velocities of individual vertices according to surrounding topology.
#    Adaptive but slow and experimental method, with lots of controls.
#
# Ideas to optimize extrusion speed further:
# - don't create internal faces into mesh
# - don't create top faces into mesh until last layer

# TODO: Autosplit twisted faces after some threshold?

# Notes for Hyperbolic Extrusion:
# Unit of velocity is meters per cell layer. Iteration step size
# unit is fraction of cell layer. One way to think about it: If one
# cell layer was extruded in one second, then velocity unit would be
# in meters per second.

# Nomenclature: (TODO: Add more terms)
#
# iteration = internal iteration step in hyperbolic extrusion of a cell layer
# velocity = measure of speed (scalar value, in unit of m/layer)
# speed = velocity * direction normal vector (vector)

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
        is_ok, text, initial_faces = initialize_extrusion()
        if not is_ok:
            self.report({'ERROR'}, "Initialization failed: " + text)
            return {'FINISHED'}

        # Layer extrusion, initialize stuff
        ug_props = bpy.context.scene.ug_props
        n = 0 # new cell count
        niter = 0 # iteration counter
        speeds = [] # Extrusion speed vectors, can be updated per layer
        new_ugfaces = [] # List of new ugfaces created in extrusion

        import bmesh
        import time
        ob = ug.get_ug_object()
        bpy.ops.object.mode_set(mode='EDIT')
        bm = bmesh.from_edit_mesh(ob.data)
        bmt = bmesh.new() # Extrusion trajectory mesh

        # Save initial face areas (used for scaling extrusion length)
        initial_face_areas = \
            [f.calc_area() for f in bm.faces if f.select]

        # Extrude layers
        for i in range(ug_props.extrusion_layers):
            t0 = time.time()
            is_last_layer = (i == (ug_props.extrusion_layers - 1))
            if i == 0:
                niter, bm, bmt, nf, speeds, new_ugfaces = \
                    extrude_cells(niter, bm, bmt, initial_faces, speeds, \
                                  new_ugfaces, initial_face_areas, \
                                  is_last_layer)
            else:
                niter, bm, bmt, nf, speeds, new_ugfaces = \
                    extrude_cells(niter, bm, bmt, [], speeds, \
                                  new_ugfaces, initial_face_areas, \
                                  is_last_layer)
            t1 = time.time()
            l.debug("Extruded layer %d, " % (i + 1) \
                    + "cells: %d " % n \
                    + "iters: %d " % niter \
                    + "(%d cells/s)" % int(nf/(t1-t0)))
            n += nf

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


def recreate_trajectory_object(bm):
    '''Replace trajectory object mesh with argument mesh'''

    # Do nothing in fixed extrusion mode
    ug_props = bpy.context.scene.ug_props
    if ug_props.extrusion_uses_fixed_initial_directions:
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


def initialize_extrusion():
    '''Initialize UG data for extrusion. For a new unstructured grid,
    create UG object and UGFaces from faces of active object.
    Return values are boolean for successful initialization and
    list of initial faces (if initializing from faces when no cells exist).
    '''

    initial_faces = [] # List of new UGFaces

    # Do nothing if there is already an UG state
    if ug.exists_ug_state():
        return True, "UG state exists, OK to continue", initial_faces

    source_ob = bpy.context.active_object
    if source_ob.name == ug.obname:
        return False, "Source object name can't be " + ug.obname, initial_faces

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
        return False, "No faces selected", initial_faces

    # Check mesh for hanging verts
    vertlist = check_hanging_face_verts(bm)
    if vertlist:
        for f in bm.faces:
            f.select_set(False)
        for v in vertlist:
            v.select_set(True)
        bm.to_mesh(source_ob.data)
        bm.free()
        return False, "Found %d hanging vert(s)" % len(vertlist), initial_faces

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
        initial_faces.append(ugf)
    l.debug("Initial Face count: %d" % len(ug.ugfaces))

    bm.to_mesh(ob.data)
    bm.free()

    # Hide everything else than UG object
    bpy.ops.object.mode_set(mode = 'OBJECT')
    ug.hide_other_objects()
    bpy.ops.object.mode_set(mode = 'EDIT')

    return True, "initialization done", initial_faces


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


def extrude_cells(niter, bm, bmt, initial_faces, speeds, new_ugfaces, \
                  initial_face_areas, is_last_layer):
    '''Extrude new cells from current face selection. Initial faces
    argument provides optional list of initial UGFaces whose direction
    is reversed at the end.
    '''

    import bmesh
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
    if len(initial_faces) > 0:
        new_ugfaces = list(initial_faces) # List of created UGFaces

    def add_entry(lol, i, val):
        '''Add entry val to list with index i into list of lists lol'''
        if i == len(lol):
            lol.append([])
        if i < len(lol):
            lol[i].append(val)
            return lol
        else:
            raise ValueError("Illegal index %d, list length %d" % (i, len(lol)))

    def get_verts_and_relations(faces):
        '''Generate face and vertex topology relationships.
        First return value is list of bmesh verts that are part of argument
        bmesh faces.
        Second return value is list of face vertex lists (to map from face to
        it's vertex indices).
        Third return value is list of vertex face lists (to map
        from vertex to it's face indices).
        Fourth return value is maps face vertices to their index number.
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

    base_edges, edge2sideface_index, fis2edges = get_edges_and_face_map(base_faces)


    def classify_verts(verts, edges, faces, fis_of_vis, ibasevertmap):
        '''Generate vertex classification information.
        First return value is True for corner vertices.
        Second return value is True for boundary vertices.
        Third return value is list of two neighbour vertex indices for
        boundary vertices (or Nones if vertex is not boundary vertex).
        Fourth return value is list of neighbour vertices connected by
        edges at vertex.
        '''

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

    # Generate vertex classification information
    if not ug_props.extrusion_uses_fixed_initial_directions:
        is_corners, is_boundaries, bnvis_of_vi, anvis_of_vi = \
            classify_verts(base_verts, base_edges, base_faces, \
                           base_fis_of_vis, ibasevertmap)
        if fulldebug:
            l.debug("is_corners %s" % str(is_corners))
            l.debug("is_boundaries %s" % str(is_boundaries))
            l.debug("bnvis_of_vi %s" % str(bnvis_of_vi))
            l.debug("anvis_of_vi %s" % str(anvis_of_vi))


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

    if not ug_props.extrusion_uses_fixed_initial_directions:
        neighbour_vis_of_vi, fils_of_neighbour_vis = \
            get_face_vis_of_vi(base_verts, base_fis_of_vis, base_vis_of_fis)


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

    new_ugfaces = create_UG_verts_faces_and_cells(base_verts, base_edges, base_faces, new_ugfaces)


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

    # Initialization of state variables
    if len(speeds) == 0:
        # Initial speeds from vertex normal speeds
        speeds = get_vertex_normal_speeds(base_verts, base_faces, \
                                          base_fis_of_vis)

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

        # Help geometry functions

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

    def calculate_face_areas(faces):
        '''Calculate areas of faces'''
        areas = []
        for f in faces:
            areas.append(f.calc_area())
        return areas

    fi_areas = calculate_face_areas(base_faces)


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

    area_coeffs = calculate_mean_vertex_area_change( \
        fi_areas, initial_face_areas, base_fis_of_vis)


    # Initial step size for hyperbolic extrusion
    df = 1e-3

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
            if ug_props.extrusion_uses_fixed_initial_directions:
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

    bm, top_verts, vert_map = cast_vertices(\
        bm, base_verts, speeds, df)


    def create_mesh_faces(bm, edge2sideface_index, vert_map, base_faces, \
                          faces2edge, ugci0, ugfi0):
        '''Create bmesh side and top faces, and link mesh face with
        UGFace. Return top faces.
        '''

        # Face creation help functions

        def create_face_from_verts(verts, ugci, fi, ugfi):
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

        def create_side_face_from_edge(e, vert_map, ugci, fi, ugfi):
            '''Creates a face from edge e'''
            # Generate vertex list for face creation
            e0 = e.verts[0]
            e1 = e.verts[1]
            verts = [e0, vert_map[e0], vert_map[e1], e1]
            create_face_from_verts(verts, ugci, fi, ugfi)

        def create_top_face_from_base_face(f, vert_map, ugci, fi, ugfi):
            '''Creates a top face from base face f'''
            verts = []
            for i in f.verts:
                verts.append(vert_map[i])
            ftop = create_face_from_verts(verts, ugci, fi, ugfi)

            # Deselect base face and select top face
            f.select_set(False)
            ftop.select_set(True)
            return ftop


        ugfi = ugfi0 # UGFace index
        fi = len(bm.faces) # mesh face index

        # Create side faces one cell at a time
        processed_edges = []
        for ci in range(len(base_faces)):
            ugci = ugci0 + ci # UGCell index
            for e in fis2edges[ci]:
                if e not in processed_edges:
                    # New face (boundary or internal)
                    create_side_face_from_edge(e, vert_map, ugci, fi, ugfi)
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
            ftop = create_top_face_from_base_face(base_faces[ci], vert_map, \
                                                  ugci0 + ci, fi, ugfi)
            top_faces.append(ftop)
            fi += 1
            ugfi += 1

        return top_faces

    top_faces = create_mesh_faces(bm, edge2sideface_index, vert_map, \
                                  base_faces, fis2edges, ugci0, ugfi0)


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

    correct_face_normals(bm, base_faces, ugci0)


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

    if not ug_props.extrusion_uses_fixed_initial_directions:
        intvis, intvpairs0, intvpairs1 = \
            calculate_intvpairs(anvis_of_vi, base_fis_of_vis, is_boundaries)


    ########################################
    # Formation Flying Extrusion Algorithm #
    ########################################

    def evolve_iteration(bm, top_verts, speeds, is_boundaries, anvis_of_vi, \
                         intvpairs0, intvpairs1, top_faces, fis_of_vis, \
                         coords0, vnspeeds0, area_coeffs, layer_frac, \
                         neighbour_vis_of_vi):
        '''Evolve hyperbolic extrusion by one iteration'''

        # Help functions

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
            ug_props = bpy.context.scene.ug_props
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

            from .ug_checks import get_edge_stats_from_bmesh_faces

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
            ug_props = bpy.context.scene.ug_props

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

            ug_props = bpy.context.scene.ug_props
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
            ug_props = bpy.context.scene.ug_props
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
            ug_props = bpy.context.scene.ug_props

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


        #########################################################
        # Main Iteration Loop of the Formation Flying Algorithm #
        #########################################################

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
            bmp, bt = get_bmp(vi, top_faces, top_verts, fis_of_vis, \
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
            df = 1.0
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

    # Call the extrusion iteration in loop
    ug_props = bpy.context.scene.ug_props
    if not ug_props.extrusion_uses_fixed_initial_directions:

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
                bm, top_verts, speeds, is_boundaries, anvis_of_vi, \
                intvpairs0, intvpairs1, top_faces, base_fis_of_vis, \
                coords0, vnspeeds0, area_coeffs, layer_frac, \
                neighbour_vis_of_vi)

            if ug_props.extrusion_create_trajectory_object:
                top_vert_cos = [v.co for v in top_verts]
                bmt = update_trajectory_mesh(bmt, top_vert_cos) # Current coordinates
                #bmt = update_trajectory_mesh(bmt, target_cos) # Target coordinates

            if print_iterations:
                max_vel = max([x.length for x in speeds])
                l.debug("iter %d: df %f, max_vel %f" % (niter, df, max_vel))

            niter += 1


    def add_faces_to_trajectory_mesh(bmt, nv0, vis_of_fis):
        '''Add final faces to mesh'''
        bmt.verts.ensure_lookup_table()
        bmt.verts.index_update()
        for vis in vis_of_fis:
            verts = [bmt.verts[nv0 + vi] for vi in vis]
            bmt.faces.new(verts)
        return bmt

    # Add faces to trajectory object at last layer
    if not ug_props.extrusion_uses_fixed_initial_directions:
        if is_last_layer and ug_props.extrusion_create_trajectory_object:
            nv0 = len(bmt.verts) - len(top_verts)
            bmt = add_faces_to_trajectory_mesh(bmt, nv0, base_vis_of_fis)


    def add_base_face_to_cells(faces, ugci0):
        '''Add base faces to cells as neighbours'''

        for i in range(len(faces)):
            # Add existing UGFace to UGCell
            ugf = ug.facemap[faces[i].index]
            ugci = ugci0 + i
            c = ug.ugcells[ugci]
            c.add_face_and_verts(ugf)
            ugf.neighbour = c

    add_base_face_to_cells(base_faces, ugci0)


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

    thickness_update()


    # Reverse direction of initial faces (initial extrusion only)
    bm.faces.ensure_lookup_table()
    for ugf in initial_faces:
        if fulldebug: l.debug("Final flipping face %d" % ugf.bi)
        bm.faces[ugf.bi].normal_flip()
        bm.faces[ugf.bi].normal_update()
        ugf.invert_face_dir()

    # End of extrude_cells()
    return niter, bm, bmt, len(base_faces), speeds, new_ugfaces
