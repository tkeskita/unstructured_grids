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
import logging
l = logging.getLogger(__name__)
fulldebug = True # Set to True if you wanna see walls of logging debug

class UG_OT_ExtrudeCells(bpy.types.Operator):
    '''Extrude new cells from current face selection'''
    bl_idname = "unstructured_grids.extrude_cells"
    bl_label = "Extrude Cells (UG)"

    @classmethod
    def poll(cls, context):
        return context.mode == 'EDIT_MESH' and ug.exists_ug_state()

    def execute(self, context):
        n = extrude_cells()
        if not n:
            self.report({'ERROR'}, "No object %r" % ug.obname)
        self.report({'INFO'}, "Extruded %d new cells" % n)
        return {'FINISHED'}


def extrude_cells():
    '''Extrude new cells from current face selection'''

    import bmesh
    ob = ug.get_ug_object()
    bm = bmesh.from_edit_mesh(ob.data)
    bm.verts.ensure_lookup_table()
    bm.faces.ensure_lookup_table()

    # Get selected faces
    faces = [f for f in bm.faces if f.select]
    l.debug("Face count: %d" % len(faces))

    # Extrusion length, TODO
    extrude_len = 0.01

    # Create new vertices (cast extrude length towards vertex normal direction)
    # Create list of boundary faces and cells to be created
    # Create new faces and hide/delete old faces
    # Create new cells

    def cast_vertices(bm, faces):
        '''Create new vertices from vertices of faces in argument bmesh, by
        casting each vertex towards vertex normal direction. Return
        updated bmesh and vertex mapping dictionary.
        '''

        processed_verts = [] # List of processed vertices
        vert_map = {} # Dictionary mapping original face vertices to new vertices
        bverts = [] # List of boundary vertices
        for f in faces:
            for v in f.verts:
                if v in processed_verts:
                    continue
                processed_verts.append(v)

                newco = v.co + extrude_len * v.normal
                v2 = bm.verts.new(newco)
                vert_map[v] = v2

        return bm, vert_map

    bm, vert_map = cast_vertices(bm, faces)


    def find_boundary_edges(bm, faces):
        '''Find boundary edges within faces'''

        edges = [] # Boundary edge list
        processed_edges = [] # List of already processed edges
        for f in faces:
            for e in f.edges:
                if e in edges:
                    continue
                if e.select == False:
                    continue
                if e in processed_edges:
                    continue
                processed_edges.append(e)
                edge_faces = [f2 for f2 in faces if e in f2.edges]
                if len(edge_faces) < 2:
                    edges.append(e)
                    if fulldebug: l.debug("Edge %d is boundary" % e.index)
                else:
                    if fulldebug: l.debug("Edge %d is internal" % e.index)

        return edges

    b_edges = find_boundary_edges(bm, faces)


    def create_faces(bm, faces, vert_map, b_edges):
        '''Create faces to boundary sides and extrudate top'''

        # Create side faces at boundary edges
        for f in faces:
            for e in f.edges:
                if e not in b_edges:
                    continue
                # Create new side face
                e0 = e.verts[0]
                e1 = e.verts[1]
                side_verts = [e0, vert_map[e0], vert_map[e1], e1]
                f2 = bm.faces.new(side_verts)
                f2.normal_update()
                # Flip face normal if face normal points "inside"
                # Use edge center to face center as reference vector
                edgevec = 0.5*(e0.co+e1.co)
                refvec = f.calc_center_median() - edgevec
                refvec.normalize()
                cos_epsilon = f2.normal @ refvec
                if fulldebug:
                    l.debug("f2.normal:%s, refvec:%s" %(str(f2.normal), str(refvec)))
                    l.debug("cos_epsilon of %d is %f" % (f.index, cos_epsilon))
                if (cos_epsilon > 0.0):
                    f2.normal_flip()
                    f2.normal_update()
                    if fulldebug:
                        l.debug("Flipped boundary face normal")

        # Create face at top of extrusion
        for f in faces:
            topverts = []
            for i in f.verts:
                topverts.append(vert_map[i])
            ftop = bm.faces.new(topverts)
            ftop.normal_update()
            # Flip top face normal if it is opposite to original face normal
            # TODO: Need to implement or not?

            # Select new faces and deselect original faces
            f.select_set(False)
            ftop.select_set(True)

        return bm

    bm = create_faces(bm, faces, vert_map, b_edges)

    bmesh.update_edit_mesh(mesh=ob.data)
    bm.free()
    bpy.ops.object.mode_set(mode = 'OBJECT')
    bpy.ops.object.mode_set(mode = 'EDIT')
    return len(faces)
