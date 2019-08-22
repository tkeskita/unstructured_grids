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
        n = select_cells_exclusive()
        self.report({'INFO'}, "Selected vertices of %d cells" % n)
        return {'FINISHED'}


def select_cells_exclusive():
    '''Reduce vertex selection to include only whole cells'''

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

    return len(clist)

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

