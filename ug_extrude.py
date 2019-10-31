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
smoothing_factor = 0.5 # TODO: parametrize

# Summary of extrusion method: Extrusion starts from selected mesh
# faces (base faces). Each face produces one new cell (in matching
# order and numbering). Extrusion is started by casting one new vertex
# from each vertex of base faces towards a normal direction calculated
# from surrounding boundary faces. Each edge of base faces produces
# one side face by connecting old and new vertices. Finally top faces
# are added by creating faces connecting new vertices, to close each
# new cell.

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

        # Layer extrusion
        ug_props = bpy.context.scene.ug_props
        n = 0 # new cell count
        vdir = [] # Extrusion directions, can be updated per layer
        coeffs = [] # Extrusion lengths, can be updated per layer
        new_ugfaces = [] # List of new ugfaces created in extrusion
        import bmesh
        import time
        ob = ug.get_ug_object()
        bm = bmesh.from_edit_mesh(ob.data)
        for i in range(ug_props.extrusion_layers):
            t0 = time.clock()
            if i == 0:
                bm, nf, vdir, coeffs, new_ugfaces = \
                    extrude_cells(bm, initial_faces, vdir, coeffs, new_ugfaces)
            else:
                bm, nf, vdir, coeffs, new_ugfaces = \
                    extrude_cells(bm, [], vdir, coeffs, new_ugfaces)
                t1 = time.clock()
                l.debug("Extruding layer %d, " % i \
                + "cells added: %d " % n \
                + "(%d cells/s)" % int(nf/(t1-t0)))
            n += nf
        bm.normal_update()
        bmesh.update_edit_mesh(mesh=ob.data)
        bm.free()
        ug_op.set_faces_boundary_to_default(new_ugfaces)
        ug.update_ug_all_from_blender()

        self.report({'INFO'}, "Extruded %d new cells" % n)
        return {'FINISHED'}


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



def extrude_cells(bm, initial_faces, vdir, coeffs, new_ugfaces):
    '''Extrude new cells from current face selection. Initial faces
    argument provides optional list of initial UGFaces whose direction
    is reversed at the end. vdir and coeffs are directions and length
    coefficients for extrusion in multilayer extrusion, and new_ugfaces
    is a list of newly created ugfaces.
    '''

    import bmesh
    bm.verts.ensure_lookup_table()
    bm.faces.ensure_lookup_table()

    # List selected base faces for extrusion. New cell index number is
    # index number of mesh face in faces list. Actual UGCell index
    # number is cell number plus initial number of prior ugcells.
    faces = [f for f in bm.faces if f.select]
    if fulldebug: l.debug("Face count at beginning: %d" % len(faces))

    ugci0 = len(ug.ugcells) # Number of UGCells before extruding
    ugfi0 = len(ug.ugfaces) # Number of UGFaces before extruding
    if len(initial_faces) > 0:
        new_ugfaces = list(initial_faces) # List of created UGFaces

    def get_verts(faces):
        '''Return list of bmesh verts that are part of argument bmesh
        faces. These verts will be used for casting vertices in
        first part of extrusion.
        '''

        verts = [] # vertex list
        vi = 0
        for f in faces:
            for v in f.verts:
                if v not in verts:
                    verts.append(v)
        return verts

    verts = get_verts(faces)


    def get_edges_and_face_map(faces):
        '''Return list of unique bmesh edges that are part of argument bmesh
        faces. Index number of edges in edges list is used to generate
        indices for new faces.
        '''

        edges = [] # edge list

        # map dictionary from edge to cell index
        # (= same as base face index = same as edge index)
        edge2cell_index = dict()
        faces2edges = [] # face to edges list
        fi = 0 # face index
        for f in faces:
            face_edges = []
            for e in f.edges:
                if e not in edges:
                    edges.append(e)
                    edge2cell_index[e] = fi
                    fi += 1
                face_edges.append(e)
            faces2edges.append(face_edges)
        return edges, edge2cell_index, faces2edges

    edges, edge2cell_index, faces2edges = get_edges_and_face_map(faces)


    def create_UG_verts_faces_and_cells(verts, edges, faces, new_ugfaces):
        '''Create new bare bone UGCells, UGFaces and UGVerts for extrusion of
        a single cell layer
        '''

        # Vertices
        for i in range(len(verts)):
            ug.UGVertex()

        # Faces (sides and top)
        for i in range(len(edges) + len(faces)):
            ugf = ug.UGFace()
            new_ugfaces.append(ugf)

        # Cells
        for i in range(len(faces)):
            ug.UGCell()

        return new_ugfaces

    new_ugfaces = create_UG_verts_faces_and_cells(verts, edges, faces, new_ugfaces)


    def calculate_extrusion_dir_and_coeffs(verts, vdir, coeffs):
        '''Calculate normalized extrusion direction vector and length
        scale coefficients based on angles of surrounding faces.
        Return dictionaries of directions and coeffs for each argument
        mesh vertex.
        '''

        def get_vdir(norvecs):
            '''Calculate vector direction from argument normal vectors'''
            vec = Vector((0, 0, 0))
            for norvec in norvecs:
                vec += norvec
            return (vec / float(len(norvecs)))


        def get_coeffs(v, faces):
            '''Calculate vector length coefficient'''
            #return 1.0

            # Note: Idea here is to push convex vertices further out
            # to fill cavities. Alone, it works only partially.
            # If pushing effect is small, there will be intersections.
            # If pushing effect is large, then extrusion becomes
            # unstable. Therefore smoothing is required additionally.

            coeff = 0.0
            for e in v.link_edges:
                efaces = [f for f in e.link_faces if f in faces]
                if len(efaces) == 2:
                    if fulldebug: l.debug("edge %s: " % str(e))
                    f0 = efaces[0]
                    f1 = efaces[1]

                    # face-face edge angle
                    cos_epsilon = face_face_cos_angle(e, f0, f1)
                    if fulldebug: l.debug("  cos_epsilon = %f" % cos_epsilon)

                    # Angle between face center-to-face center and first face normal
                    vec_f2f = f0.calc_center_median() - f1.calc_center_median()
                    vec_f2f.normalize()
                    cos_theta = vec_f2f @ f0.normal
                    if fulldebug: l.debug("  cos_theta = %f" % cos_theta)

                    if cos_theta < 0.0:
                        coeff = max(coeff, cos_epsilon + 1) # TODO: Improve this somehow?

            # Amplify and shift
            coeff *= 2.0 # TODO: parametrize?
            coeff += 1.0
            return coeff


        from mathutils import Vector
        ug_props = bpy.context.scene.ug_props

        # Do nothing if initial direction and coeffs is to be used
        if len(vdir) > 0 and len(coeffs) > 0:
            if ug_props.extrusion_uses_fixed_initial_directions:
                return vdir, coeffs

        vdir = [] # new extrusion directions to be calculated
        coeffs = [] # new extrusion length coefficients to be calculated

        for v in verts:
            norvecs = []
            faces = []
            for f in v.link_faces:
                fi = f.index
                uf = ug.facemap[fi]
                if uf.deleted:
                    continue
                if uf.neighbour != None:
                    continue
                if ug_props.extrusion_ignores_unselected_face_normals:
                    if f.select == False:
                        continue
                norvecs.append(f.normal)
                faces.append(f)

            # Extrusion direction
            vdir.append(get_vdir(norvecs))

            # Length coefficient
            coeffs.append(get_coeffs(v, faces))
            if fulldebug: l.debug("vert %d coeff %f" % (v.index, coeffs[-1]))

        return vdir, coeffs


    def cast_vertices(bm, faces, verts, vdir, coeffs):
        '''Create new vertices from vertices of faces in argument bmesh, by
        casting each vertex towards initial (vdir) or updated average
        face normal direction. Return updated bmesh, vertex mapping
        dictionary, and initial extrusion direction dictionary.
        '''

        # Dictionary for mapping original face verts to new verts
        vert_map = {}
        save_vdir = dict() # saved direction vector for next layer
        bvi = bm.verts[-1].index # Index of last vertex
        ug_props = bpy.context.scene.ug_props

        # Extrusion length
        extrude_len = ug_props.extrusion_thickness

        # Calculate updated extrusion direction and length
        # coefficients for vertices based on current face normals
        vdir, coeffs = calculate_extrusion_dir_and_coeffs(verts, vdir, coeffs)

        # Cast new vertices
        for i in range(len(verts)):
            newco = verts[i].co + extrude_len * vdir[i] * coeffs[i]
            v2 = bm.verts.new(newco)
            vert_map[verts[i]] = v2

        bm.verts.ensure_lookup_table()
        bm.verts.index_update()

        if fulldebug: l.debug("Cast %d vertices" % len(save_vdir))
        return bm, vert_map, vdir, coeffs

    bm, vert_map, vdir, coeffs = cast_vertices(bm, faces, verts, vdir, coeffs)


    def create_mesh_faces(bm, edges, edge2cell_inex, vert_map, faces, \
                          faces2edge, ugci0, ugfi0):
        '''Create bmesh side faces and top and links mesh face with
        UGFace. Return map from cell index to list of cell faces.
        '''

        # First define face creation help functions

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


        ugfi = ugfi0 # UGFace index
        fi = len(bm.faces) # mesh face index

        # Create side faces one cell at a time
        processed_edges = []
        for ci in range(len(faces)):
            ugci = ugci0 + ci # UGCell index
            for e in faces2edges[ci]:
                if e not in processed_edges:
                    create_side_face_from_edge(e, vert_map, ugci, fi, ugfi)
                    processed_edges.append(e)
                    fi += 1
                    ugfi += 1
                else:
                    # Internal face. Add existing UGFace to UGCell
                    ugf = ug.ugfaces[ugfi0 + edge2cell_index[e]]
                    c = ug.ugcells[ugci]
                    c.add_face_and_verts(ugf)
                    ugf.neighbour = c

        # Create top faces
        for ci in range(len(faces)):
            create_top_face_from_base_face(faces[ci], vert_map, ugci0 + ci, fi, ugfi)
            fi += 1
            ugfi += 1


    create_mesh_faces(bm, edges, edge2cell_index, vert_map, faces, faces2edges, ugci0, ugfi0)


    def correct_face_normals(bm, faces, ugci0):
        '''Ensure face normal directions point out of cells'''

        # Note: For internal faces normal should point from lower
        # index cell to higher index cell.

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

    correct_face_normals(bm, faces, ugci0)


    def smoothen_verts(bm, vdir):
        '''Smoothen boundary vertex locations'''

        def get_edges_of_faces(faces):
            '''Return list of unique edges of argument faces'''
            edges = []
            for f in faces:
                for e in f.edges:
                    if e not in edges:
                        edges.append(e)
            return edges

        def get_faces_of_vert(v, faces):
            '''Get the faces from faces list which are faces of vertex v'''
            vertfaces = []
            for vf in v.link_faces:
                if vf in faces:
                    vertfaces.append(vf)
            return vertfaces

        def get_verts_of_faces(vx, faces):
            '''Get list of vertices of argument faces, excluding vertex vx'''
            tot_area = 0.0 # total area of all faces
            for f in faces:
                tot_area += f.calc_area()

            verts = []
            weights = []
            for f in faces:
                weight = f.calc_area() / tot_area
                nverts = 0
                for v in f.verts:
                    if v == vx:
                        continue
                    if v not in verts:
                        verts.append(v)
                        weights.append(weight)
            return verts, weights

        def new_internal_co_from_verts(v, verts, weights):
            '''Calculate new location for argument vertex v from verts and weights'''
            from mathutils import Vector
            co = Vector((0, 0, 0))
            nverts = len(verts)
            for i in range(nverts):
                co += (verts[i].co - v.co) * weights[i] / nverts * smoothing_factor
            if fulldebug: l.debug("co %s %f" % (str(co), co.length))
            if fulldebug: l.debug("max weight %f" % max(weights))
            return v.co + co

        def new_boundary_co_from_verts(v, v1, v2):
            '''Calculate new vertex v location from neighbour boundary vertices
            v1 and v2
            '''
            vec1 = (v1.co - v.co)
            vec2 = (v2.co - v.co)

            # Do nothing for corner vertices
            cos_alpha = vec1.normalized() @ vec2.normalized()
            if fulldebug: l.debug("cos_alpha %f" % cos_alpha)
            max_cos_angle = 2.5e-1 # maximum angle cosine to detect corners
            if abs(cos_alpha) < max_cos_angle:
                return v.co

            sqlen1 = vec1.length * vec1.length
            sqlen2 = vec2.length * vec2.length
            sqtot = sqlen1 + sqlen2
            co = 0.5 * (vec1 * sqlen1 + vec2 * sqlen2) / sqtot * smoothing_factor
            return v.co + co

        def is_boundary_vertex_among_faces(v, edges, faces):
            '''Return True if argument vertex v is a boundary vertex among
            argument edges and faces. Vertex is boundary vertex if it
            is connected to an edge that is part of only one face.
            Returns also neighbour vertices (or None for no boundary edges)
            '''
            neighbour_verts = []
            for e in v.link_edges:
                if e not in edges:
                    continue
                link_faces = [x for x in e.link_faces if x in faces]
                if len(link_faces) == 1:
                    neighbour_verts.append(e.other_vert(v))
            if len(neighbour_verts) == 2:
                return True, neighbour_verts[0], neighbour_verts[1]
            return False, None, None

        # Only selected faces are processed
        faces = [f for f in bm.faces if f.select]
        edges = get_edges_of_faces(faces)
        verts = [] # vertices
        coords = {} # Dictionary to map vertex to new coordinate

        # Generate new coordinates
        for f in faces:
            for v in f.verts:
                if v in verts:
                    continue
                verts.append(v)

                is_boundary_vertex, v1, v2 = \
                    is_boundary_vertex_among_faces(v, edges, faces)

                # Boundary vertices use boundary edges
                if is_boundary_vertex:
                    if fulldebug: l.debug("Boundary vertex %d" % v.index)
                    co = new_boundary_co_from_verts(v, v1, v2)

                # Internal vertices use faces
                else:
                    if fulldebug: l.debug("Internal vertex %d" % v.index)
                    vertfaces = get_faces_of_vert(v, faces)
                    nbverts, weights = get_verts_of_faces(v, vertfaces)
                    co = new_internal_co_from_verts(v, nbverts, weights)

                if fulldebug: l.debug("Move vertex %d " % v.index \
                                      + "from %s to %s" % (str(v.co), str(co)))
                coords[v] = co

        # Move vertices to new coordinates
        for v in coords:
            v.co = coords[v]

    # Loop vertex smoothing
    for i in range(3): # TODO: parametrize
        smoothen_verts(bm, vdir)

    def add_base_face_to_cells(faces, ugci0):
        '''Add base faces to cells as neighbours'''

        for i in range(len(faces)):
            # Add existing UGFace to UGCell
            ugf = ug.facemap[faces[i].index]
            ugci = ugci0 + i
            c = ug.ugcells[ugci]
            c.add_face_and_verts(ugf)
            ugf.neighbour = c

    add_base_face_to_cells(faces, ugci0)


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

    return bm, len(faces), vdir, coeffs, new_ugfaces


def bmesh_edge_center(edge):
    '''Calculates coordinates for edge center using edge vertex
    coordinates
    '''
    return (edge.verts[0].co + edge.verts[1].co) / 2


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
