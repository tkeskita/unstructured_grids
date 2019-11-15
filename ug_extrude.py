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

        # Layer extrusion, initialize stuff
        ug_props = bpy.context.scene.ug_props
        n = 0 # new cell count
        vdir = [] # Extrusion directions, can be updated per layer
        prev_vlens = [] # Previous extrusion length for vertices
        new_ugfaces = [] # List of new ugfaces created in extrusion
        import bmesh
        import time
        ob = ug.get_ug_object()
        bm = bmesh.from_edit_mesh(ob.data)

        # Save initial face areas (used for scaling extrusion length)
        # TODO: Implement scaling by area changes
        initial_face_areas = \
            [f.calc_area() for f in bm.faces if f.select]

        # Do the layers
        for i in range(ug_props.extrusion_layers):
            t0 = time.clock()
            if i == 0:
                bm, nf, vdir, prev_vlens, new_ugfaces = \
                    extrude_cells(bm, initial_faces, vdir, prev_vlens, \
                                  new_ugfaces, initial_face_areas)
            else:
                bm, nf, vdir, prev_vlens, new_ugfaces = \
                    extrude_cells(bm, [], vdir, prev_vlens, \
                                  new_ugfaces, initial_face_areas)
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


def extrude_cells(bm, initial_faces, vdir, prev_vlens, new_ugfaces, \
                  initial_face_areas):
    '''Extrude new cells from current face selection. Initial faces
    argument provides optional list of initial UGFaces whose direction
    is reversed at the end.
    new_ugfaces is a list of newly created ugfaces.
    '''

    import bmesh
    bm.verts.ensure_lookup_table()
    bm.faces.ensure_lookup_table()
    ug_props = bpy.context.scene.ug_props

    # List selected base faces for extrusion. New cell index number is
    # index number of mesh face in faces list. Actual UGCell index
    # number is cell number plus initial number of prior ugcells.
    base_faces = [f for f in bm.faces if f.select]
    if fulldebug: l.debug("Face count at beginning: %d" % len(base_faces))

    ugci0 = len(ug.ugcells) # Number of UGCells before extruding
    ugfi0 = len(ug.ugfaces) # Number of UGFaces before extruding
    if len(initial_faces) > 0:
        new_ugfaces = list(initial_faces) # List of created UGFaces

    def add_entry(lol, i, val):
        '''Add entry val to list with index i in list of lists lol'''
        if i == len(lol):
            lol.append([])
        if i < len(lol):
            lol[i].append(val)
            return lol
        else:
            raise ValueError("Illegal index %d, list length %d" % (i, len(lol)))

    def get_verts_and_relations(faces):
        '''Return list of bmesh verts that are part of argument bmesh faces,
        list of face vertex lists (to map from face to it's vertex
        indices), and list of vertex face lists (to map from vertex to
        it's face indices)
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


    def get_edges_and_face_map(faces):
        '''Return list of unique bmesh edges that are part of argument bmesh
        faces. Index number of edges in edges list is used to generate
        indices for new faces.
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
        '''Classify verts to corner verts, boundary verts and search boundary
        vertex indices for boundary verts
        '''

        def classify_vert(v, ibasevertmap, edges, faces):
            '''Classify argument vertex v.
            First return value is True if v is corner vertex.
            Second return value is True if v is boundary vertex.
            Third return value is list of two vertex indices (vis) pointing to
            neighbour boundary vertices (only for boundary vertices)
            '''
            neighbour_vis = []
            neighbour_faces = []
            for e in v.link_edges:
                if e not in edges:
                    continue
                link_faces = [x for x in e.link_faces if x in faces]
                if len(link_faces) == 1:
                    neighbour_vis.append(ibasevertmap[e.other_vert(v)])
                    neighbour_faces.append(link_faces[0])

            if len(neighbour_vis) == 2:
                if neighbour_faces[0] == neighbour_faces[1]:
                    return True, True, [neighbour_vis[0], neighbour_vis[1]]
                else:
                    return False, True, [neighbour_vis[0], neighbour_vis[1]]
            return False, False, [None, None]

        # Turn lists into sets for fast search
        edgeset = set(edges)
        faceset = set(faces)

        is_corners = []
        is_boundaries = []
        boundary_vert_neighbours = []
        for v in verts:
            is_corner, is_boundary, neigh_verts = \
                classify_vert(v, ibasevertmap, edgeset, faceset)
            is_corners.append(is_corner)
            is_boundaries.append(is_boundary)
            boundary_vert_neighbours.append(neigh_verts)
        return is_corners, is_boundaries, boundary_vert_neighbours

    if not ug_props.extrusion_uses_fixed_initial_directions:
        is_corners, is_boundaries, boundary_vert_neighbours = \
            classify_verts(base_verts, base_edges, base_faces, base_fis_of_vis, ibasevertmap)


    def get_vis_of_vi(verts, fis_of_vis, vis_of_fis):
        '''First return value is list which contains, for each vertex, a list
        of other vertex indices of neighbour vertices.
        Second return value list of face indices for the same vertices.
        '''
        neighbour_vis_of_vi = []
        fis_of_neighbour_vis = []
        for vi in range(len(verts)):
            vis = []
            fis = []
            for fi in fis_of_vis[vi]:
                for i in vis_of_fis[fi]:
                    if i in vis:
                        continue
                    if i != vi:
                        vis.append(i)
                        fis.append(fi)
            neighbour_vis_of_vi.append(vis)
            fis_of_neighbour_vis.append(fis)
        return neighbour_vis_of_vi, fis_of_neighbour_vis

    if not ug_props.extrusion_uses_fixed_initial_directions:
        neighbour_vis_of_vi, fis_of_neighbour_vis = \
            get_vis_of_vi(base_verts, base_fis_of_vis, base_vis_of_fis)


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


    def get_mean_dir(vecs):
        '''Calculate average vector direction from argument vectors'''
        from mathutils import Vector
        resvec = Vector((0, 0, 0))
        for vec in vecs:
            resvec += vec
        return resvec.normalized()


    def get_items_from_list(alist, ilist):
        '''Return items in alist located at index locations in ilist'''
        res = []
        for i in ilist:
            res.append(alist[i])
        return res


    def calculate_initial_extrusion_dir(vdir, verts, faces, fis_of_vis):
        '''Calculate initial extrusion direction vectors by averaging face
        normals surrounding each vertex
        '''

        ug_props = bpy.context.scene.ug_props

        # Do nothing if initial direction option is enabled
        if len(vdir) > 0 and ug_props.extrusion_uses_fixed_initial_directions:
            return vdir

        vdir = [] # new extrusion directions to be calculated
        for vi in range(len(verts)):
            neigh_faces = get_items_from_list(faces, fis_of_vis[vi])
            norvecs = [f.normal for f in neigh_faces]
            vdir.append(get_mean_dir(norvecs))

        return vdir

    vdir = calculate_initial_extrusion_dir(vdir, base_verts, base_faces, \
                                           base_fis_of_vis)


    def cast_vertices(bm, base_verts, vdir):
        '''Create new vertices from base vertices in argument bmesh, by
        casting each vertex towards vdir
        '''

        vert_map = {} # Dictionary for mapping original face verts to new verts
        ug_props = bpy.context.scene.ug_props
        top_verts = []
        vlens = [] # Extrusion lengths

        # Extrusion length
        extrude_len = ug_props.extrusion_thickness / float(ug_props.extrusion_substeps)

        # Cast new vertices
        for i in range(len(base_verts)):
            newco = base_verts[i].co + extrude_len * vdir[i]

            v2 = bm.verts.new(newco)
            vert_map[base_verts[i]] = v2
            top_verts.append(v2)
            vlens.append(extrude_len)

        bm.verts.ensure_lookup_table()
        bm.verts.index_update()

        if fulldebug: l.debug("Cast %d vertices" % len(vdir))
        return bm, top_verts, vert_map, vlens

    bm, top_verts, vert_map, vlens = cast_vertices(bm, base_verts, vdir)

    if len(prev_vlens) == 0:
        prev_vlens = vlens


    def create_mesh_faces(bm, edge2sideface_index, vert_map, base_faces, \
                          faces2edge, ugci0, ugfi0):
        '''Create bmesh side faces and top and links mesh face with
        UGFace. Return map from top BMFace to original face area used
        in extrusion length scaling by face area.
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

    top_faces = create_mesh_faces(bm, edge2sideface_index, vert_map, base_faces, \
                                  fis2edges, ugci0, ugfi0)


    def correct_face_normals(bm, faces, ugci0):
        '''Ensure face normal directions point out of cells'''

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


    ###############################
    ### Extending and smoothing ###
    ###############################

    def calculate_max_convexities(bm, verts, faces):
        '''Calculate minimum convexity limitation factors for verts'''

        ug_props = bpy.context.scene.ug_props
        # minimum max_convexity to allow smoothing
        a = ug_props.extrusion_convexity_min
        # maximum max_convexity to clamp smoothing
        b = ug_props.extrusion_convexity_max
        # maximum clamp value
        c = ug_props.extrusion_convexity_clamp

        slope = c / (b - a)
        root = -a * slope

        max_convexities = []
        for vi, v in enumerate(verts):
            max_convexity = 0.0
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

                    # Calculate convexity: value 0 means extremely
                    # sharp concave point, value 0.5 means flat
                    # terrain (neither concave or convex) and value 1
                    # means extremely convex deep hole.
                    if cos_theta <= 0.0:
                        convexity = (cos_epsilon + 1.0) / 4.0 + 0.5
                    else:
                        convexity = (-cos_epsilon + 1.0) / 4.0
                    if fulldebug: l.debug("  convexity = %f" % convexity)
                    max_convexity = max(max_convexity, convexity)

            # Calculate limitation coefficient based on minimum convexity
            coeff = root + slope * max_convexity
            coeff = min(c, max(0.0, coeff))
            max_convexities.append(coeff)
            if fulldebug: l.debug("  max_convexity = %f" % max_convexity)
            if fulldebug: l.debug("  max_convexity_coeff %f" % coeff)
        return max_convexities

    def propagate_max_convexities(verts, max_convexities, neighbour_vis_of_vi):
        '''Propagate maximum convexity value from neighbouring vertices'''

        max_vals = []
        radius = 0.07 # allowed radius for propagation
        for i, v in enumerate(verts):
            nvals = get_items_from_list(max_convexities, neighbour_vis_of_vi[i])
            nverts = get_items_from_list(verts, neighbour_vis_of_vi[i])
            for nvi, nv in enumerate(nverts):
                vec = nv.co - v.co
                if vec.length < radius:
                    scale = 1.0
                else:
                    scale = 0.0
                nvals[nvi] *= scale
            nvals.append(max_convexities[i]) # own value is possible, too
            max_vals.append(max(nvals))
        return max_vals

    def calculate_face_areas(faces):
        '''Calculate areas of faces'''
        areas = []
        for f in faces:
            areas.append(f.calc_area())
        return areas

    def calculate_mean_vertex_area_change(fi_areas, initial_face_areas, \
                                          base_fis_of_vis):
        '''Calculate mean relative change of areas of faces around each vertex'''
        mean_changes = []
        for vi in range(len(base_fis_of_vis)):
            fis = base_fis_of_vis[vi]
            mean_change = 0.0
            for fi in fis:
                mean_change += initial_face_areas[fi] / fi_areas[fi]

            mean_change /= float(len(fis))
            mean_changes.append(mean_change)
        return mean_changes

    def extend_verts(bm, top_verts, base_verts, fi_areas, initial_face_areas, \
                     base_fis_of_vis, prev_vlens):
        '''Move top vertices outwards to extend cell in a substep'''

        ug_props = bpy.context.scene.ug_props
        nsteps = float(ug_props.extrusion_substeps)
        ext_len = ug_props.extrusion_thickness
        corner_factor = ug_props.extrusion_corner_factor
        area_factor = ug_props.extrusion_area_factor
        growth_damping_factor = ug_props.extrusion_growth_damping_factor

        # Calculate mean area change per vertex
        area_coeffs = calculate_mean_vertex_area_change( \
            fi_areas, initial_face_areas, base_fis_of_vis)

        for tv, bv, a, prev_len, ic in \
            zip(top_verts, base_verts, area_coeffs, prev_vlens, is_corners):

            vdir = tv.co - bv.co
            vdir.normalize()

            # Extension length scaling
            if a > 1.0:
                # Area change suggests extension of step legth. Dampen
                # area coefficient by damping factor to reduce
                # zigzagging in concave extrusion fronts.
                a2 = (a - 1.0) * growth_damping_factor + 1.0
                step_len = ext_len * a2 / nsteps
            else:
                step_len = ext_len * a / nsteps

            # Area scaling (non-corners and corners)
            if not ic:
                tv.co = tv.co + vdir * step_len * area_factor
            else:
                tv.co = tv.co + vdir * step_len * area_factor * corner_factor


    def smoothen_verts(bm, vdir, top_verts, base_verts, top_faces, \
                       base_faces, neighbour_vis_of_vi, fis_of_neighbour_vis, \
                       fi_areas, boundary_vert_neighbours, base_fis_of_vis, \
                       convexities):
        '''Smoothen top vertex locations'''

        def new_internal_co_from_vert(v, vi, verts, neighbour_vis_of_vi, fis_of_neighbour_vis, fi_areas):
            '''Calculate new location for argument vertex v from verts and weights'''
            from mathutils import Vector
            ug_props = bpy.context.scene.ug_props
            smoothing_factor = ug_props.extrusion_smoothing_factor
            neighbour_verts = [verts[i] for i in neighbour_vis_of_vi[vi]]
            areas = [fi_areas[i] for i in fis_of_neighbour_vis[vi]]
            tot_area = sum(areas)

            co = Vector((0, 0, 0))
            for nv, area in zip(neighbour_verts, areas):
                co += (nv.co - v.co) * area
            co = co / tot_area * smoothing_factor

            if fulldebug: l.debug("co %s %f" % (str(co), co.length))
            if fulldebug: l.debug("tot_area %f" % tot_area)
            return v.co + co

        def new_boundary_co_from_verts(v, vi, verts, boundary_vert_neighbours):
            '''Calculate new vertex v location from neighbour boundary vertices
            v1 and v2
            '''
            ug_props = bpy.context.scene.ug_props
            smoothing_factor = ug_props.extrusion_smoothing_factor
            v1 = verts[boundary_vert_neighbours[vi][0]]
            v2 = verts[boundary_vert_neighbours[vi][1]]
            vec1 = (v1.co - v.co)
            vec2 = (v2.co - v.co)
            sqlen1 = vec1.length * vec1.length
            sqlen2 = vec2.length * vec2.length
            sqtot = sqlen1 + sqlen2
            co = 0.5 * (vec1 * sqlen1 + vec2 * sqlen2) / sqtot * smoothing_factor
            return v.co + co

        def limit_co_by_angle_deviation(co, v, vdir):
            '''Limit (smoothened) coordinates co by the angle it is deviating from
            vdir. Additionally limit the length
            '''

            from math import sqrt
            ug_props = bpy.context.scene.ug_props

            # Minimum allowed cosine of angle between vdir and u (vector from
            # base vertex to co)
            min_cos_alpha = ug_props.extrusion_deviation_angle_min
            # Minimum length coefficient for u
            min_len_coeff = ug_props.extrusion_deviation_length_min
            # Maximum length coefficient for u
            max_len_coeff = ug_props.extrusion_deviation_length_max

            extrude_len = ug_props.extrusion_thickness \
                / float(ug_props.extrusion_substeps + 1) # TODO: +1 enough, or need more?
            # cvec = vdir * extrude_len

            u = co - v.co # Vector from base vertex to new coordinates
            q = (u @ vdir) * vdir # u component aligned with vdir
            p = u - q # u component normal to vdir
            m = u.normalized() # Normalized u
            cos_alpha = m @ vdir
            if cos_alpha < 0.0:
                raise ValueError("cos_alpha negative: %f" % cos_alpha)

            # Decrease angle by moving co on vdir normal plane
            if cos_alpha < min_cos_alpha:
                plen = sqrt(1.0/(min_cos_alpha * min_cos_alpha) - 1.0)
                p.normalize()
                p = plen * q.length * p
                u = p + q

            # Increase length if needed
            if u.length < min_len_coeff * extrude_len:
                u.normalize()
                u = min_len_coeff * extrude_len * u

            # Decrease length if needed
            if u.length > max_len_coeff * extrude_len:
                u.normalize()
                u = max_len_coeff * extrude_len * u

            return v.co + u

        def limit_co_by_plane(co, v, vdir):
            '''Limit (smoothened) coordinates co by a plane whose normal is vdir
            and plane location is a small distance from vertex v towards vdir.
            If new coordinates co are on wrong side of the plane,
            return new coordinates projected to the plane.
            '''
            ug_props = bpy.context.scene.ug_props
            extrude_len = ug_props.extrusion_thickness \
                / float(ug_props.extrusion_substeps + 1) # TODO: +1 enough, or need more?
            cvec = vdir * extrude_len

            # Center coordinates of limitation plane
            projection_fraction = 0.2 # TODO: Parametrize or remove
            center = v.co + projection_fraction * cvec

            testvec = co - center

            # Return current coordinates if plane is not cut
            cos_alpha = testvec.normalized() @ vdir
            if cos_alpha > 0.0:
                return co

            # Otherwise, return projection to plane
            scale = testvec @ vdir
            newco = co - scale * vdir
            return newco

        def limit_by_faces(co, v, faces):
            '''Limit coordinate co by faces around vertex v'''
            for f in faces:
                co = limit_co_by_plane(co, v, f.normal)
            return co


        new_coords = [] # new coordinate list

        for vi, v in enumerate(top_verts):
            # Corners are not smoothened
            if is_corners[vi]:
                new_coords.append(v.co)
                continue

            # Smoothing of other boundary vertices
            if is_boundaries[vi]:
                newco = new_boundary_co_from_verts(v, vi, top_verts, boundary_vert_neighbours)

            # Smoothing of internal vertices
            else:
                newco = new_internal_co_from_vert(v, vi, top_verts, neighbour_vis_of_vi, fis_of_neighbour_vis, fi_areas)

            if fulldebug:
                l.debug("Propose move vertex %d " % v.index \
                        + "from %s to %s" % (str(v.co), str(newco)))

            # Limitations
            oldv = base_verts[vi]
            old_faces = [base_faces[i] for i in base_fis_of_vis[vi]]
            # Limit by plane calculated from extrusion vdir
            newco = limit_co_by_plane(newco, oldv, vdir[vi])
            # Limit by angle deviation
            if ug_props.extrusion_uses_angle_deviation:
                newco = limit_co_by_angle_deviation(newco, oldv, vdir[vi])
            # Limit by factor calculated from convexity
            if ug_props.extrusion_uses_convexity:
                convexity_coeff = convexities[vi]
                newco = v.co + (newco - v.co) * convexity_coeff
            # Limit by surrounding faces
            newco = limit_by_faces(newco, oldv, old_faces)
            # Limit again by vdir plane
            newco = limit_co_by_plane(newco, oldv, vdir[vi])

            new_coords.append(newco)

        # Move vertices to new coordinates after all new positions
        # have been calculated
        for v, c in zip(top_verts, new_coords):
            v.co = c


    # Main vertex extension + smoothing loops

    ug_props = bpy.context.scene.ug_props
    if not ug_props.extrusion_uses_fixed_initial_directions and \
        ug_props.extrusion_smoothing_iterations > 0:

        convexities = []
        if ug_props.extrusion_uses_convexity:
            convexities = calculate_max_convexities(bm, base_verts, base_faces)
            for i in range(ug_props.extrusion_convexity_propagations):
                convexities = propagate_max_convexities(top_verts, \
                                                        convexities, \
                                                        neighbour_vis_of_vi)

        # First run smoothing rounds for first substep
        for i in range(ug_props.extrusion_smoothing_iterations):
            fi_areas = calculate_face_areas(top_faces)
            smoothen_verts(bm, vdir, top_verts, base_verts, top_faces, \
                           base_faces, neighbour_vis_of_vi, \
                           fis_of_neighbour_vis, fi_areas, \
                           boundary_vert_neighbours, base_fis_of_vis, \
                           convexities)

        # Carry out the rest of substeps (extend + smoothings)
        for j in range(ug_props.extrusion_substeps - 1):
            fi_areas = calculate_face_areas(top_faces)
            extend_verts(bm, top_verts, base_verts, fi_areas, \
                         initial_face_areas, base_fis_of_vis, prev_vlens)

            for i in range(ug_props.extrusion_smoothing_iterations):
                fi_areas = calculate_face_areas(top_faces)
                smoothen_verts(bm, vdir, top_verts, base_verts, top_faces, \
                               base_faces, neighbour_vis_of_vi, \
                               fis_of_neighbour_vis, fi_areas, \
                               boundary_vert_neighbours, base_fis_of_vis, \
                               convexities)



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

    # Update extrusion lengths for next round
    prev_vlens=[]
    for vold, v in zip(base_verts, top_verts):
        vec = v.co - vold.co
        prev_vlens.append(vec.length)

    return bm, len(base_faces), vdir, prev_vlens, new_ugfaces


# TODO: Remove obsolete help functions?

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
