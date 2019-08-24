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
        return context.mode in {'OBJECT','EDIT_MESH'} and ug.exists_ug_state()

    def execute(self, context):
        n = select_cells_inclusive()
        if not n:
            self.report({'ERROR'}, "No object %r" % ug.obname)
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
        for v in c.ugverts:
            if v.bi in viset:
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
        for f in c.ugfaces:
            for v in f.ugverts:
                if ob.data.vertices[v.bi].select == False:
                    ob.data.vertices[v.bi].select = True
                    n += 1
    return n


class UG_OT_SelectCellsExclusive(bpy.types.Operator):
    '''Reduce current vertex selection to include only whole cells'''

    bl_idname = "unstructured_grids.select_cells_exclusive"
    bl_label = "Select UG Cells (Exclusive)"

    @classmethod
    def poll(cls, context):
        return context.mode in {'OBJECT','EDIT_MESH'} and ug.exists_ug_state()

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
        for v in c.ugverts:
            if v.bi not in viset:
                test = False
        if test:
            clist.append(c)
        if i % print_interval == 0:
            l.debug("... processed cell count: %d" % i)
        i += 1
    return clist


def select_vertices_from_ugfaces(ob, flist):
    '''Select mesh object vertices that are part of faces in argument UGFace
    list. Return number of selected vertices.
    '''

    n = 0
    for f in flist:
        for v in f.ugverts:
            if ob.data.vertices[v.bi].select == False:
                ob.data.vertices[v.bi].select = True
                n += 1
    return n


###########################
##### CELL OPERATIONS #####
###########################


class UG_OT_DeleteCells(bpy.types.Operator):
    '''Delete whole cells in current vertex selection'''

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
    for c in clist:
        delete_cell(c)
        if fulldebug: l.debug("Deleted cell %d" % c.ii)
    return len(clist)


def delete_cell(c):
    '''Delete argument UGCell'''

    # Process faces: Delete boundary faces and transform internal
    # faces to boundary faces
    for f in c.ugfaces:
        # internal face, c is owner
        if f.owner == c and f.neighbour != None:
            f.invert_face_dir()
            convert_face_i2b(f)

        # internal face, c is neighbour
        elif f.neighbour == c:
            convert_face_i2b(f)

        # c is boundary face
        elif f.neighbour == None:
            delete_boundary_face(f)

    # Mark cell as deleted
    c.deleted = True

    # Mark unused vertices as deleted
    delete_vertices_of_deleted_cell(c)


def convert_face_i2b(f):
    '''Convert internal face into boundary face'''

    import bmesh

    # Debug safety check
    if fulldebug:
        if f.owner == None or f.neighbour == None:
            l.error("Not internal face")
        if f.owner.deleted:
            l.error("Cell %d is deleted" % f.owner.ii)

    # Mark as boundary face
    f.neighbour = None

    # Add new face to geometry
    ob = ug.get_ug_object()
    bpy.ops.object.mode_set(mode = 'OBJECT')
    bpy.ops.object.mode_set(mode = 'EDIT')
    bm = bmesh.from_edit_mesh(ob.data)
    bm.faces.ensure_lookup_table()
    bm.edges.ensure_lookup_table()
    bm.verts.ensure_lookup_table()
    bvlist = []
    for v in f.ugverts:
        bvlist.append(bm.verts[v.bi])
    bf = bm.faces.new(bvlist)
    bm.faces.index_update()
    bmesh.update_edit_mesh(mesh=ob.data)
    bm.free()
    # Refresh modified mesh
    bpy.ops.object.mode_set(mode = 'OBJECT')
    bpy.ops.object.mode_set(mode = 'EDIT')

    # Set face index
    f.bi = ob.data.polygons[-1].index
    if fulldebug: l.debug("New face bi %d" % f.bi)
    set_face_boundary_to_default(f)


def delete_boundary_face(f):
    '''Delete argument boundary UGFace'''

    # Debug safety check
    if fulldebug:
        if f.neighbour != None:
            l.error("Not boundary face")
        if f.owner.deleted:
            l.error("Cell %d is deleted" % f.owner.ii)

    # Hide face in mesh. Hiding seems to work only in Object Mode.
    # And it seems to only show up correctly in Edit Mode.
    ob = ug.get_ug_object()
    bpy.ops.object.mode_set(mode = 'OBJECT')
    ob.data.polygons[f.bi].hide = True
    bpy.ops.object.mode_set(mode = 'EDIT')

    # Mark face as deleted
    f.deleted = True


def delete_vertices_of_deleted_cell(c):
    '''Delete unused UGVertices after deletion of argument UGCell'''

    for v in c.ugverts:
        cellFound = False
        for vc in v.ugcells:
            if vc.deleted:
                continue
            if vc == c:
                continue
            cellFound = True
            break
        if not cellFound:
            v.deleted = True


def set_face_boundary_to_default(f):
    '''Set material (boundary patch) of argument UGFace to default'''

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
        mat = bpy.data.materials.new(name=matname)
        bpy.ops.object.material_slot_add()
        ob.active_material = mat

    # Find slot index for default material:
    for mati in range(len(ob.material_slots)):
        if ob.material_slots[mati].name == matname:
            break

    bpy.ops.object.mode_set(mode = 'OBJECT')
    ob.data.polygons[f.bi].material_index = mati
    bpy.ops.object.mode_set(mode = 'EDIT')
