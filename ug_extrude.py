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
LARGE = 100.0 # Cut-off for large values for weight function

# Summary of extrusion method: Extrusion starts from selected mesh
# faces (base faces). Each face produces one new cell (in matching
# order and numbering). Extrusion is started by casting one new vertex
# (top vertex) from each vertex of base faces (base vertex) towards a
# normal direction calculated from surrounding boundary faces. Each edge
# of base faces produces one side face by connecting old and new vertices.
# Finally top faces are added by creating faces connecting new vertices,
# to close each new cell.

# Ideas to optimize extrusion speed further:
# - don't create internal faces into mesh
# - don't create top faces into mesh until last layer

class UG_OT_ExtrudeCells(bpy.types.Operator):
    '''Extrude new cells from current face selection'''
    bl_idname = "unstructured_grids.extrude_cells"
    bl_label = "Extrude Cells (UG)"

    @classmethod
    def poll(cls, context):
        return context.mode in {'OBJECT', 'EDIT_MESH'}

    def execute(self, context):
        # Initialize from selected faces if needed
        initialization_ok, initial_faces = initialize_extrusion()
        if not initialization_ok:
            self.report({'ERROR'}, "Initialization failed. Maybe " \
                        + "no faces were selected, or object name is " \
                        + "%r?" % ug.obname)
            return {'FINISHED'}

        # Layer extrusion, initialize stuff
        ug_props = bpy.context.scene.ug_props
        n = 0 # new cell count
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
            t0 = time.clock()
            if i == 0:
                bm, bmt, nf, speeds, new_ugfaces = \
                    extrude_cells(bm, bmt, initial_faces, speeds, \
                                  new_ugfaces, initial_face_areas)
            else:
                bm, bmt, nf, speeds, new_ugfaces = \
                    extrude_cells(bm, bmt, [], speeds, \
                                  new_ugfaces, initial_face_areas)
            t1 = time.clock()
            l.debug("Extruded layer %d, " % (i + 1) \
                    + "cells added: %d " % n \
                    + "(%d cells/s)" % int(nf/(t1-t0)))
            n += nf

        bm.normal_update()
        bmesh.update_edit_mesh(mesh=ob.data)
        bm.free()
        recreate_trajectory_object(bmt)
        bmt.free()
        ug_op.set_faces_boundary_to_default(new_ugfaces)
        ug.update_ug_all_from_blender()

        self.report({'INFO'}, "Extruded %d new cells" % n)
        return {'FINISHED'}


def recreate_trajectory_object(bm):
    '''Replace trajectory object mesh with argument mesh'''

    obname = "Extrusion Trajectory"
    bpy.ops.object.mode_set(mode='OBJECT')

    # Delete old object
    if obname in bpy.data.objects:
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
        bpy.context.view_layer.objects.active = bpy.data.objects[obname]


def initialize_extrusion():
    '''Initialize UG data for extrusion. For a new unstructured grid,
    create UG object and UGFaces from faces of active object.
    Return values are boolean for successful initialization and
    list of initial faces (if initializing from faces when no cells exist).
    '''

    initial_faces = [] # List of new UGFaces

    # Do nothing if there is already an UG state
    if ug.exists_ug_state():
        return True, initial_faces

    source_ob = bpy.context.active_object
    if source_ob.name == ug.obname:
        return False, initial_faces

    # Mode switch is needed to make sure mesh is saved to original object
    bpy.ops.object.mode_set(mode='OBJECT')
    bpy.ops.object.mode_set(mode='EDIT')
    ob = ug.initialize_ug_object()

    import bmesh
    bm = bmesh.from_edit_mesh(source_ob.data)

    # Delete unselected faces
    facelist = []
    for f in bm.faces:
        if f.select == False:
            facelist.append(f)
    bmesh.ops.delete(bm, geom=facelist, context='FACES_ONLY')

    # Delete leftover verts which are not part of any faces
    vertlist = []
    for v in bm.verts:
        if len(v.link_faces) == 0:
            vertlist.append(v)
    bmesh.ops.delete(bm, geom=vertlist, context='VERTS')

    bm.verts.ensure_lookup_table()
    bm.faces.ensure_lookup_table()

    # Bail out if no faces are left
    if len(bm.faces) == 0:
        ug.delete_ug_object()
        return False, initial_faces

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

    return True, initial_faces


def extrude_cells(bm, bmt, initial_faces, speeds, new_ugfaces, initial_face_areas):
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

    # Populate initial verts to trajectory bmesh
    if ug_props.extrusion_create_trajectory_object and len(bmt.verts) == 0:
        for v in base_verts:
            bmt.verts.new(v.co)


    def get_edges_and_face_map(faces):
        '''Generate face and edge topology relationshops.
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

    #else: # TODO: remove
    #    # Dummy default for fixed extrusion
    #    is_corners = [False] * len(base_verts)


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


    def get_vertex_normal_speeds(oldspeeds, verts, faces, fis_of_vis):
        '''Calculate vertex normal speeds by averaging face
        normals surrounding each vertex, weighted by face vertex angle
        (the angle between the two edges of the face connected at vertex).
        '''

        from mathutils import Vector
        ug_props = bpy.context.scene.ug_props

        # Base extension length of a substep
        ext_len = ug_props.extrusion_thickness / float(ug_props.extrusion_substeps)

        # Do nothing if initial direction option is enabled
        if len(oldspeeds) > 0 and ug_props.extrusion_uses_fixed_initial_directions:
            return oldspeeds

        # Calculate new extrusion directions
        speeds = []
        for vi in range(len(verts)):
            angle_factors = []
            neigh_faces = get_items_from_list(faces, fis_of_vis[vi])

            # Calculate weight coefficients from face vertex angles
            for f in neigh_faces:
                # Invert cos angle and make range from 0 to 2:
                angle_factor = (-1.0 * get_face_vertex_cos_angle(verts[vi], f)) + 1.0
                angle_factors.append(angle_factor)

            # Calculate final speeds from normals and weights
            norvecs = [f.normal for f in neigh_faces]
            norvec = Vector((0, 0, 0))
            for nv, factor in zip(norvecs, angle_factors):
                norvec += factor * nv
            norvec.normalize()

            if len(oldspeeds) > 0:
                ext_len = oldspeeds[vi].length

            norvec *= ext_len # Set velocity
            speeds.append(norvec)

        return speeds

    if len(speeds) == 0:
        speeds = get_vertex_normal_speeds(speeds, base_verts, base_faces, \
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


    # TODO: Remove?
    def calculate_face_areas(faces):
        '''Calculate areas of faces'''
        areas = []
        for f in faces:
            areas.append(f.calc_area())
        return areas

    #fi_areas = calculate_face_areas(base_faces)

    # TODO: Remove?
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

    #area_coeffs = calculate_mean_vertex_area_change( \
    #    fi_areas, initial_face_areas, base_fis_of_vis)

    # TODO: Remove?
    def get_convexity_factor(max_convexity):
        '''Calculated convexity factor from maximum convexity.
        Used both for extension length scaling and angle deviation
        maximum length scaling.
        '''

        ug_props = bpy.context.scene.ug_props
        convexity_factor = ug_props.extrusion_convexity_scale_factor

        if max_convexity > 0.5:
            return (max_convexity - 0.5) * 2.0 * convexity_factor + 1.0
        else:
            return 1.0


    # TODO: Remove?
    def calculate_extension_length(area_coeffs, is_corners, max_convexities):
        '''Calculate length how far vertices are extended towards extension
        direction during casting or extension steps
        '''

        ug_props = bpy.context.scene.ug_props
        # Base extension length of a substep
        ext_len = ug_props.extrusion_thickness / float(ug_props.extrusion_substeps)
        # Additional length factor for corners
        corner_factor = ug_props.extrusion_corner_factor
        # Factor used in scaling length according to area
        area_factor = ug_props.extrusion_area_factor
        # Factor used to limit length when growing length
        growth_scale_factor = ug_props.extrusion_growth_scale_factor

        elens = [] # Extension lengths, to be calculated
        for a, ic, mc in zip(area_coeffs, is_corners, max_convexities):
            factor = 1.0 # Cumulative length factor

            # Convexity scaling is done always
            factor *= get_convexity_factor(mc)

            # Other scalings are done only for non-fixed extrusions
            if not ug_props.extrusion_uses_fixed_initial_directions:
                # Area change scaling
                if a > 1.0:
                    # Area change suggests extension of step legth. Scale
                    # area coefficient by growth scaling factor < 1.0 to reduce
                    # zigzagging in concave extrusion fronts.
                    afac = (a - 1.0) * growth_scale_factor * area_factor + 1.0
                else:
                    afac = (a - 1.0) * area_factor + 1.0
                factor *= afac

                # Special scaling for corners
                if ic:
                    factor = corner_factor

            step_length = ext_len * factor
            elens.append(step_length)

        return elens

    # elens = calculate_extension_length(area_coeffs, is_corners, max_convexities)

    def cast_vertices(bm, base_verts, speeds):
        '''Create new top vertices from base vertices in argument bmesh, by
        casting each vertex towards speeds (length elen). Returns updated bmesh,
        top vertices and map from base verts to new top verts.
        '''

        vert_map = {} # Dictionary for mapping original face verts to new verts
        ug_props = bpy.context.scene.ug_props
        top_verts = []

        # Full step for fixed extrusion, small initial step for hyperbolic
        if ug_props.extrusion_uses_fixed_initial_directions:
            speed_scale = 1.0
        else:
            speed_scale = 0.001

        # Cast new vertices very small distance ahead
        for i in range(len(base_verts)):
            newco = base_verts[i].co + speed_scale * speeds[i]
            v2 = bm.verts.new(newco)
            vert_map[base_verts[i]] = v2
            top_verts.append(v2)

        bm.verts.ensure_lookup_table()
        bm.verts.index_update()

        if fulldebug: l.debug("Cast %d vertices" % len(speeds))
        return bm, top_verts, vert_map

    bm, top_verts, vert_map = cast_vertices(\
        bm, base_verts, speeds)


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

    def evolve_substep(bm, top_verts, speeds, is_boundaries, anvis_of_vi, \
                       intvpairs0, intvpairs1, top_faces, fis_of_vis):
        '''Evolve extrusion by one substep'''

        # Help functions

        def get_estimated_cos(verts, speeds):
            '''Estimate next locations of vertices based on current locations and
            speeds
            '''
            estimated_cos = []
            for v, s in zip(verts, speeds):
                estimated_cos.append(v.co + s)
            return estimated_cos

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

        def is_above_planes(vco, faces, co):
            '''Return True if co is located "above" (=on the normal direction
            side) of all faces, all of which connect at
            coordinates vco. Otherwise return False.
            '''
            u = co - vco
            u.normalize()
            test = True
            for f in faces:
                if (u @ f.normal) <= 0.0:
                    test = False
            return test

        def get_target_cos(v, vi, vpairs0, vpairs1, estimated_cos, \
                           faces, pvgm_target_co):
            '''First return value is middle coordinates for all vpairs
            (except pvgm_target_co if middle coordinate is below faces).
            Second return value is booleans mid_cos being above faces.
            Third return value is distances between vpairs.
            '''
            from mathutils import Vector
            mid_cos = [] # Middle coordinates
            is_aboves = [] # Booleans for above planes
            vplengths = [] # Distances between vpairs

            for v0, v1 in zip(vpairs0, vpairs1):
                co = (estimated_cos[v0] + estimated_cos[v1]) / 2.0
                is_above = is_above_planes(v.co, faces, co)
                is_aboves.append(is_above)
                # Use pvgm_target_co for targets that are below planes
                if not is_above:
                    co = Vector(pvgm_target_co)
                mid_cos.append(co)
                vec = estimated_cos[v1] - estimated_cos[v0]
                vplengths.append(vec.length)
            return mid_cos, is_aboves, vplengths

        def get_cut_substeps(vpairs0, vpairs1, vimap, pvs, pvspeeds):
            '''Calculate how many substeps it would take for intersection to
            happen for each neighbour vertex pair if neighbours
            continued with their current pspeeds
            '''
            cut_substeps = [] # estimated substeps until intersection

            # Process each vertex pair
            for i, j in zip(vpairs0, vpairs1):
                co0 = pvs[vimap[i]]
                co1 = pvs[vimap[j]]
                s0 = pvspeeds[vimap[i]]
                s1 = pvspeeds[vimap[j]]

                # Calculate substeps until intersection
                vec0 = co1 - co0 # vector between verts
                vdir = vec0.normalized()
                q0 = (s0 @ vdir) * vdir # speed component aligned with vdir
                q1 = (s1 @ vdir) * vdir # speed component aligned with vdir

                vec1 = co1 + q1 - (co0 + q0) # vector between moved verts

                # If distance is increasing (not going to intersect),
                # then set cut_substep to LARGE
                if vec1.length > vec0.length:
                    cut_substeps.append(LARGE)

                # Otherwise calculate substeps until intersection
                else:
                    cut_step = vec0.length / (vec0.length - vec1.length)
                    cut_step = min(LARGE, cut_step)
                    cut_substeps.append(cut_step)
            return cut_substeps

        def project_co_above_planes(v, co, faces):
            '''Project coordinates co above faces, all of which connect at vertex
            v, by length speed
            '''
            for f in faces:
                u = co - v.co
                unorm = u.normalized()
                cos_angle = unorm @ f.normal
                if cos_angle < 0.0:
                    # u component aligned with f.normal
                    q = (u @ f.normal) * f.normal
                    co -= q # TODO: Does this need extra factor like 1.05?
            return co

        def get_pvgm_target_co(vi, top_verts, anvis, speeds, faces):
            '''Calculate geometric mean from vertex neighbours' locations and
            speeds. Target coordinate is projected up from planes.
            '''
            cv = top_verts[vi]
            from mathutils import Vector
            co = Vector((0, 0, 0))
            for i in anvis:
                v = top_verts[i]
                s = speeds[i]
                co += v.co + s
            co /= float(len(anvis))
            return (project_co_above_planes(cv, co, faces))

        def limit_co_by_angle_deviation(co, v, speed):
            '''Limit coordinates co by angle deviation by moving
            coordinate inside a conical volume. Cone starting point is
            base vertex, and normal vector (vdir) is calculated from
            speed vector.
            '''

            from math import sqrt
            ug_props = bpy.context.scene.ug_props

            # Minimum allowed cosine of angle
            min_cos_alpha = ug_props.extrusion_deviation_angle_min

            # Vector from vertex to new coordinates
            u = co - v.co
            m = u.normalized()
            vdir = speed.normalized()
            cos_alpha = m @ vdir

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
            return v.co + u

        def get_weights(is_aboves, cut_substeps, vplengths, \
                        pvgm_target_co, convexity_sums, anvis, vi):
            '''Calculate weights for target coordinates. Weights are used to
            calculate a target coordinate and speed.
            '''
            ug_props = bpy.context.scene.ug_props

            def weight(x):
                '''Calculate weight from argument'''
                z = ug_props.extrusion_weight_smoothing_coefficient
                wmax = 1.0 / z # Normalisation factor
                w = 1.0 / (x + z) / wmax
                return w * w

            def get_ave_neigh_convexity(convexity_sums, anvis):
                '''Calculate average convexity among neighbours'''
                cs = [convexity_sums[x] for x in anvis]
                return sum(cs) / float(len(cs))

            gmf = ug_props.extrusion_geometric_mean_factor
            weights = []

            # Weights for target coordinates
            if len(is_aboves) > 0:
                max_vplen = max(vplengths)
                for ia, cs, vl in zip(is_aboves, cut_substeps, vplengths):
                    # Weight option 1: Closeness of vpair
                    vplen_fac = (1.0 + 10.0 * (vl / max_vplen)) # TODO: Parametrize
                    wo1 = weight(vplen_fac)
                    # Weight option 2: Value of cut substeps
                    wo2 = weight(cs)
                    # Weight option 3: Minimum convex weight # TODO: Remove?
                    wo3 = weight(gmf)
                    w = max(wo1, wo2, wo3)
                    weights.append(w)
                    if fulldebug:
                        l.debug("vfrac %f " % (vl/max_vplen) \
                                + "closeness %f " % wo1 \
                                + "cut_substeps %f " % wo2 \
                                + "min weight %f " % wo3 \
                                + "w %f" % w)

            # Average neighbour convexity
            anv = get_ave_neigh_convexity(convexity_sums, anvis)

            # Append geometric mean weight to end
            gmw = 4.0 + anv * 80.0 # TODO: Parametrize
            weights.append(weight(gmw))

            # Append vertex normal weight to end
            vnw = 1.0 + anv * 200.0 # TODO: Parametrize
            weights.append(weight(vnw))

            return weights, anv

        def get_speeds(v, target_cos, is_aboves, cut_substeps, pvgm_target_co, \
                       min_velocity, convexity_sum, vnspeed):
            '''Calculate target speeds for target_cos'''

            from mathutils import Vector
            ug_props = bpy.context.scene.ug_props
            # Speed factor for convex vertices
            csf = ug_props.extrusion_convex_speed_factor
            cfac = (1.0 + convexity_sum * csf)

            speeds = []
            # Otherwise ivpairs data exists, get velocity
            if len(is_aboves) > 0:
                for co, ia, cs in zip(target_cos, is_aboves, cut_substeps):
                    # Calculate velocity such that target_co
                    # is reached in fraction of cut_substep, but
                    # always at least with minimum convex velocity.
                    vec = co - v.co
                    vel = max(vec.length, min_velocity) * cfac
                    vel = vec.length * cfac
                    speeds.append(vel * vec.normalized())

            # Append geometric mean coordinate vector to end
            pvgm_vec = pvgm_target_co - v.co
            if pvgm_vec.length < min_velocity:
                pvgm_vec.normalize()
                pvgm_vec *= min_velocity
            speeds.append(pvgm_vec)

            # Append vertex normal speed to end
            speeds.append(cfac * vnspeed)

            return speeds, cfac

        def get_min_neighbour_vel(speeds, anvis_of_vi):
            '''Calculate minimum neighbour velocity'''
            vels = [x.length for x in get_items_from_list(speeds, anvis_of_vi)]
            return min(vels)

        def get_target_speed(oldv, oldspeed, target_speeds, weights, \
                             min_neighbour_vel, min_velocity, anv, convexity_sum):
            '''Calculate new speed vector from target_cos and
            weights
            '''
            ug_props = bpy.context.scene.ug_props
            from mathutils import Vector
            target_speed = Vector((0, 0, 0))

            # Calculate target speed
            for s, w in zip(target_speeds, weights):
                target_speed += s * w
            target_speed /= sum(weights)

            #return target_speed # for debugging

            # Maximum allowed relative velocity change
            mrv = 1.0 + convexity_sum * 0.5 # TODO: Parametrize acceleration factor
            mrv *= ug_props.extrusion_max_relative_velocity

            # Limit velocity change rate
            vel = target_speed.length
            oldvel = oldspeed.length
            vel = min(vel, oldvel * mrv) # Clamp upper limit
            vel = max(vel, oldvel / mrv * 0.5) # Clamp lower limit # TODO: Parametrize braking factor
            vel = max(vel, min_velocity) # Ensure minimum velocity

            # Limit the change of direction angle.
            limited_co = limit_co_by_angle_deviation\
                         (oldv.co + target_speed, oldv, oldspeed)
            vec = limited_co - oldv.co
            vec.normalize()
            target_speed = vel * vec

            # Limit by maximum allowed relative velocity among neighbours
            #if target_speed.length > (mrv * min_neighbour_vel):
            #    target_speed = (mrv * min_neighbour_vel) * vec.normalized()

            return target_speed


        #######################################################
        # Main Substep Loop of the Formation Flying Algorithm #
        #######################################################

        # Calculate vertex normal speeds
        vnspeeds = get_vertex_normal_speeds(speeds, top_verts, top_faces, \
                                            base_fis_of_vis)

        # Calculate estimated_cos = estimate where vertices would be
        # after this substep
        estimated_cos = get_estimated_cos(top_verts, speeds)

        # Calculate convexity_sums = sum of convex edge convexity values per vertex
        convexity_sums = calculate_convexity_sums(bm, top_verts, top_faces)

        new_speeds = [] # new speed vectors, to be calculated
        for vi, v in enumerate(top_verts):
            # Speed of corner vertices are not changed
            if is_corners[vi]:
                new_speeds.append(speeds[vi])
                continue

            # List of faces surrounding this vertex
            vfaces = [top_faces[x] for x in fis_of_vis[vi]]

            # Boundary vertices use plain vertex normal speed
            if is_boundaries[vi]:
                new_speeds.append(vnspeeds[vi])
                continue

            if fulldebug:
                l.debug("Starting on internal vi %d index %d" % (vi, v.index))

            # Calculate projected geometric mean target coordinates
            min_velocity = ug_props.extrusion_thickness / float(ug_props.extrusion_substeps)
            min_speed = min_velocity * speeds[vi].normalized()
            pvgm_target_co = get_pvgm_target_co(vi, top_verts, anvis_of_vi[vi], \
                                                speeds, vfaces)

            # Calculate pvs = project neighbour vertices to direction normal plane
            # Calculate pvspeeds = project neighbour speeds to direction normal plane
            pvs, pvspeeds = project_pvs_pvspeeds(vi, top_verts, anvis_of_vi[vi], speeds)

            # Calculate target_cos = target coordinates for all vpairs
            # is_aboves = booleans for original target being located above faces
            # vplengths = lengths between vpairs
            target_cos, is_aboves, vplengths = \
                get_target_cos(v, vi, intvpairs0[vi], intvpairs1[vi], \
                               estimated_cos, vfaces, pvgm_target_co)

            # Calculate amount of substeps until vpairs intersect
            vimap = {} # Map from vi index to anvis index
            for i, x in enumerate(anvis_of_vi[vi]):
                vimap[x] = i
            cut_substeps = get_cut_substeps(intvpairs0[vi], intvpairs1[vi], vimap, pvs, pvspeeds)

            # Calculate weights for target_cos + geometric mean
            weights, anv = get_weights(is_aboves, cut_substeps, vplengths, \
                                       pvgm_target_co, convexity_sums, \
                                       anvis_of_vi[vi], vi)

            # Calculate target_speeds for target_cos + geometric mean + vertex normal
            target_speeds, cfac = get_speeds(\
                v, target_cos, is_aboves, cut_substeps, pvgm_target_co, \
                min_velocity, convexity_sums[vi], vnspeeds[vi])

            # Calculate minimum neighbour velocity to limit speed
            min_neighbour_vel = get_min_neighbour_vel(speeds, anvis_of_vi[vi])

            # Calculate a new target speed vector for this vertex
            target_speed = get_target_speed(v, speeds[vi], target_speeds, \
                                            weights, min_neighbour_vel, \
                                            min_velocity, pvgm_target_co, convexity_sums[vi])
            if fulldebug:
                l.debug("is_aboves %s" % str(is_aboves))
                l.debug("target_cos %s" % str(target_cos))
                l.debug("cut_substeps %s" % str(cut_substeps))
                l.debug("weights %s" % str(weights))
                l.debug("convexity_sum %f" % convexity_sums[vi] \
                        + ", anv %f" % anv \
                        + ", cfac %f" % cfac)
                l.debug("target_speeds lengths %s" % str([x.length for x in target_speeds]))
                l.debug("target_co %s " % str(v.co + target_speed) \
                        + "velocity %f" % target_speed.length)

            new_speeds.append(target_speed)

        # Evolve vertex positions and return new speeds
        for v, s in zip(top_verts, new_speeds):
            v.co += s
        return new_speeds


    def update_trajectory_mesh(bmt, verts):
        '''Add top verts to trajectory mesh'''
        bmt.verts.ensure_lookup_table()
        bmt.verts.index_update()
        nverts = len(verts)
        nv0 = len(bmt.verts) - nverts
        oldverts = [bmt.verts[x] for x in range(nv0, nv0 + nverts)]
        for i,v in enumerate(top_verts):
            nv = bmt.verts.new(v.co)
            v0 = oldverts[i]
            bmt.edges.new([v0, nv])
        return bmt

    # Call the extrusion substep in loop
    ug_props = bpy.context.scene.ug_props
    if not ug_props.extrusion_uses_fixed_initial_directions:
        # Carry out substeps
        for j in range(ug_props.extrusion_substeps):
            if fulldebug: l.debug("Extrusion substep %d" % j)
            speeds = evolve_substep(bm, top_verts, speeds, \
                                    is_boundaries, anvis_of_vi, \
                                    intvpairs0, intvpairs1, \
                                    top_faces, base_fis_of_vis)
            if ug_props.extrusion_create_trajectory_object:
                bmt = update_trajectory_mesh(bmt, top_verts)


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

    return bm, bmt, len(base_faces), speeds, new_ugfaces
