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

# Operations related to Unstructured Grids

import bpy
from . import ug
import logging
l = logging.getLogger(__name__)
fulldebug = False # Set to True if you wanna see walls of logging debug
print_interval = 100000 # Debug print progress interval

class UG_OT_SelectCellsInclusive(bpy.types.Operator):
    '''Extend current vertex selection to include whole cells'''
    bl_idname = "unstructured_grids.select_cells_inclusive"
    bl_label = "Select UG Cells (Inclusive)"

    @classmethod
    def poll(cls, context):
        return context.mode == 'EDIT_MESH' and ug.exists_ug_state()

    def execute(self, context):
        n = select_cells_inclusive()
        if not n:
            self.report({'ERROR'}, "No cells found. Maybe selection was " \
                        + "empty, cells don't exist there, or object name " \
                        + "is not %r?" % ug.obname)
        self.report({'INFO'}, "Selected vertices of %d cells" % n)
        return {'FINISHED'}


def select_cells_inclusive():
    '''Extend vertex selection to include all cells that currently
    selected vertices are part of
    '''

    ob = ug.get_ug_object()
    mode = ob.mode # Save original mode
    # Return to object mode to update selection
    if mode == 'EDIT':
        bpy.ops.object.mode_set(mode='OBJECT')

    # First get cells that are part of selected vertices
    verts = [v.index for v in ob.data.vertices if v.select]
    l.debug("Initially selected vertex count: %d" % len(verts))
    clist = get_ugcells_from_vertices_inclusive(verts)

    # Select cell vertices
    n = select_vertices_from_ugcells(ob, clist)
    l.debug("Finally selected vertex count: %d" % n)
    # Return to original mode
    bpy.ops.object.mode_set(mode=mode)

    return len(clist)


def get_ugcells_from_vertices_inclusive(vilist):
    '''Return list of UGCells that are part of argument vertex index list'''

    clist = []
    i = 0
    viset = set(vilist) # Convert to set for fast search speed

    for c in ug.ugcells:
        if c.deleted:
            continue
        for ugv in c.ugverts:
            if ugv.bi in viset:
                if not clist or clist[-1] != c:
                    clist.append(c)
            if i % print_interval == 0:
                l.debug("... processed vertex count: %d" % i)
            i += 1
    return clist


def select_vertices_from_ugcells(ob, clist):
    '''Select mesh object vertices that are part of cells in argument UGCell
    list. Return number of selected vertices.
    '''

    n = 0
    for c in clist:
        for ugf in c.ugfaces:
            for ugv in ugf.ugverts:
                if ob.data.vertices[ugv.bi].select == False:
                    ob.data.vertices[ugv.bi].select = True
                    n += 1
    return n


class UG_OT_SelectCellsExclusive(bpy.types.Operator):
    '''Reduce current vertex selection to include only whole cells'''

    bl_idname = "unstructured_grids.select_cells_exclusive"
    bl_label = "Select UG Cells (Exclusive)"

    @classmethod
    def poll(cls, context):
        return context.mode == 'EDIT_MESH' and ug.exists_ug_state()

    def execute(self, context):
        n = len(select_cells_exclusive())
        self.report({'INFO'}, "Selected vertices of %d cells" % n)
        return {'FINISHED'}


def select_cells_exclusive():
    '''Reduce vertex selection to include only whole cells.
    Return list of cells.
    '''

    ob = ug.get_ug_object()
    mode = ob.mode # Save original mode
    # Return to object mode to update selection
    if mode == 'EDIT':
        bpy.ops.object.mode_set(mode='OBJECT')

    # First get cells that are part of selected vertices
    verts = [v.index for v in ob.data.vertices if v.select]
    l.debug("Initially selected vertex count: %d" % len(verts))
    clist = get_ugcells_from_vertices_exclusive(verts)

    # Deselect all vertices, edges and faces
    bpy.ops.object.mode_set(mode="EDIT")
    bpy.ops.mesh.select_all(action='DESELECT')
    bpy.ops.object.mode_set(mode="OBJECT")

    # Select only whole cell vertices
    n = select_vertices_from_ugcells(ob, clist)
    l.debug("Finally selected vertex count: %d" % n)
    # Return to original mode
    bpy.ops.object.mode_set(mode=mode)

    return clist

def get_ugcells_from_vertices_exclusive(vilist):
    '''Return list of UGCells that are completely defined by vertices in
    vertex index list vilist
    '''

    clist = []
    i = 0
    viset = set(vilist) # Convert to set for fast search speed

    for c in ug.ugcells:
        if c.deleted:
            continue
        test = True
        for ugv in c.ugverts:
            if ugv.bi not in viset:
                test = False
        if test:
            clist.append(c)
        if i % print_interval == 0:
            l.debug("... processed cell count: %d" % i)
        i += 1
    return clist


def get_ugfaces_from_vertices_exclusive(vilist):
    '''Return list of UGFaces that are completely defined by vertices in
    vertex index list vilist
    '''

    ugflist = []
    i = 0
    viset = set(vilist) # Convert to set for fast search speed

    for ugf in ug.ugfaces:
        if ugf.deleted:
            continue
        test = True
        for ugv in ugf.ugverts:
            if ugv.bi not in viset:
                test = False
        if test:
            ugflist.append(ugf)
        if i % print_interval == 0:
            l.debug("... processed face count: %d" % i)
        i += 1
    return ugflist


def select_vertices_from_ugfaces(ob, ugflist):
    '''Select mesh object vertices that are part of faces in argument UGFace
    list. Return number of selected vertices.
    '''

    n = 0
    for ugf in ugflist:
        for ugv in ugf.ugverts:
            if ob.data.vertices[ugv.bi].select == False:
                ob.data.vertices[ugv.bi].select = True
                n += 1
    return n


class UG_OT_ResetView(bpy.types.Operator):
    '''Reset view in Edit Mode. Show boundary faces and hide deleted faces and verts'''

    bl_idname = "unstructured_grids.reset_view"
    bl_label = "Reset View (UG)"

    @classmethod
    def poll(cls, context):
        return context.mode in {'OBJECT','EDIT_MESH'} and ug.exists_ug_state()

    def execute(self, context):
        reset_view()
        self.report({'INFO'}, "View has been reset")
        return {'FINISHED'}


def reset_view():
    '''Reset view in Edit Mode. Show boundary faces and hide deleted and
    internal faces and deleted verts.
    '''

    # Note: Hidden elements show up correctly only in Edit Mode.

    import bmesh
    ob = ug.get_ug_object()
    bpy.ops.object.mode_set(mode = 'OBJECT')
    bpy.ops.object.mode_set(mode = 'EDIT')

    bm = bmesh.from_edit_mesh(ob.data)
    bm.faces.ensure_lookup_table()
    bm.verts.ensure_lookup_table()

    nf = 0 # number of hidden faces
    nv = 0 # number of hidden vertices

    # Hiding order seems to matter. First vertices, then
    # faces. Otherwise faces are not hidden correctly.

    # TODO: Reset view still does not work correctly when first
    # extrude new cells (internal faces formed in mesh), then delete
    # some cells. Deletion works if mesh is exported and imported
    # after extrusion prior to cell deletions.

    for ugv in ug.ugverts:
        if ugv.deleted:
            bm.verts[ugv.bi].hide_set(True)
            nv += 1
        else:
            bm.verts[ugv.bi].hide_set(False)

    for ugf in ug.ugfaces:
        if ugf.bi == -1:
            continue
        if ugf.deleted or ugf.neighbour != None:
            bm.faces[ugf.bi].hide_set(True)
            nf += 1
        else:
            bm.faces[ugf.bi].hide_set(False)

    bmesh.update_edit_mesh(mesh=ob.data)
    bm.free()
    bpy.ops.object.mode_set(mode = 'OBJECT')
    bpy.ops.object.mode_set(mode = 'EDIT')
    l.debug("Hidden faces: %d, hidden verts: %d" % (nf, nv))


###########################
##### CELL OPERATIONS #####
###########################


class UG_OT_DeleteCells(bpy.types.Operator):
    '''Delete cells covered by current vertex selection'''

    bl_idname = "unstructured_grids.delete_cells"
    bl_label = "Delete Cells (UG)"

    @classmethod
    def poll(cls, context):
        return context.mode == 'EDIT_MESH' and ug.exists_ug_state()

    def execute(self, context):
        n = delete_cells_from_vertex_selection()
        self.report({'INFO'}, "Deleted %d cells" % n)
        return {'FINISHED'}


def delete_cells_from_vertex_selection():
    '''Delete cells in current vertex selection'''

    clist = select_cells_exclusive()

    bfacelist = [] # new boundary face list
    delfacelist = [] # deletion face list

    # Delete each cell and collect list of new boundary faces and
    # faces to be deleted
    for c in clist:
        bfacelist, delfacelist = delete_cell(c, bfacelist, delfacelist)
        if fulldebug: l.debug("Deleted cell %d" % c.ii)

    l.debug("New boundary faces: %d" % len(bfacelist))
    l.debug("Old boundary faces for deletion: %d" % len(delfacelist))

    # Add faces to geometry for new boundary faces
    add_faces_i2b(bfacelist)

    # Hide old boundary faces and deleted vertices
    reset_view()

    return len(clist)


def delete_cell(c, bfacelist, delfacelist):
    '''Delete argument UGCell c. Append list of boundary faces for
    deletion (bfacelist) and list of internal faces for transformation
    into boundary faces (delfacelist). Return both lists updated.
    '''

    ugflist = [] # UGFaces to be deleted

    # Handle cell faces
    for ugf in c.ugfaces:
        if fulldebug: l.debug("Checking face %d" % ugf.bi)
        if ugf.owner == c and ugf.neighbour != None:
            if fulldebug: l.debug("internal face, c is owner")
            ugf.invert_face_dir()
            ugf.neighbour = None
            bfacelist.append(ugf)

        elif ugf.neighbour == c:
            if fulldebug: l.debug("internal face, c is neighbour")
            ugf.neighbour = None
            bfacelist.append(ugf)

        elif ugf.neighbour == None:
            if fulldebug: l.debug("boundary face")
            # Mark face. Can't use ugf.delete() right away, maybe
            # because it affects c.ugfaces loop? Delete faces after loop.
            ugf.deleted = True
            ugf.owner = None

            # Remove face from new boundary face list if needed
            if ugf in bfacelist:
                bfacelist.remove(ugf)
            delfacelist.append(ugf)

    # Delete faces
    for ugf in ugflist:
        ugf.delete()

    # Delete cell
    c.delete()

    return bfacelist, delfacelist


def add_faces_i2b(ugflist):
    '''Add faces to geometry for internal UGFaces turned into boundary faces'''

    import bmesh

    # Add new face to geometry
    ob = ug.get_ug_object()
    bpy.ops.object.mode_set(mode = 'OBJECT')
    bpy.ops.object.mode_set(mode = 'EDIT')

    bm = bmesh.from_edit_mesh(ob.data)
    bm.faces.ensure_lookup_table()
    bm.edges.ensure_lookup_table()
    bm.verts.ensure_lookup_table()

    fi = len(bm.faces) # face index

    for ugf in ugflist:
        if ugf.bi > -1:
            continue
        bvlist = []
        for ugv in ugf.ugverts:
            bvlist.append(bm.verts[ugv.bi])
        bf = bm.faces.new(bvlist)
        if not bf:
            l.error("Could not create face " + str(bvlist))
        ugf.bi = fi
        ug.facemap[fi] = ugf
        fi += 1

    bm.verts.index_update()
    bm.edges.index_update()
    bm.faces.index_update()
    bmesh.update_edit_mesh(mesh=ob.data)
    bm.free()

    # Refresh modified mesh
    bpy.ops.object.mode_set(mode = 'OBJECT')
    bpy.ops.object.mode_set(mode = 'EDIT')

    set_faces_boundary_to_default(ugflist)
    l.debug("Converted faces: %d" % len(ugflist))


def set_faces_boundary_to_default(ugflist):
    '''Set material (boundary patch) of argument UGFaces to default'''

    matname = 'default'
    ob = ug.get_ug_object()

    patch = None
    for p in ug.ugboundaries:
        if p.patchname == matname:
            patch = p
            patch.deleted = False
            break
    if not patch:
        patch = ug.UGBoundary(matname)

    if 'default' in bpy.data.materials:
        mat = bpy.data.materials['default']
    else:
        mat = bpy.data.materials.new(name=matname)
        bpy.ops.object.material_slot_add()
    ob.active_material = mat

    # Find slot index for default material:
    for mati in range(len(ob.material_slots)):
        if ob.material_slots[mati].name == matname:
            break

    bpy.ops.object.mode_set(mode = 'OBJECT')
    for ugf in ugflist:
        if ugf.is_boundary_face():
            ob.data.polygons[ugf.bi].material_index = mati
    bpy.ops.object.mode_set(mode = 'EDIT')
