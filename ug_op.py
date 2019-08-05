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

class UG_OT_Select_Cells_Inclusive(bpy.types.Operator):
    '''Operator to extend vertex selection to include all cells that currently
    selected vertices are part of
    '''
    bl_idname = "unstructured_grids.select_cells_inclusive"
    bl_label = "UG Select Cells (Inclusive)"

    @classmethod
    def poll(cls, context):
        ob = ug.get_ug_object()
        return (ob and ob.type == 'MESH' and \
                context.mode in {'OBJECT','EDIT_MESH'})

    def execute(self, context):
        n = select_cells_inclusive()
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

    verts = [v for v in ob.data.vertices if v.select]
    l.debug("Initially selected vertex count: %d" % len(verts))

    # Find cells whose vertices are to be selected 
    clist = []
    for v in verts:
        ugvert = ug.ugverts[v.index]
        iis = ''
        for c in ugvert.ugcells:
            iis += str(c.ii) + ' '
            if c not in clist:
                clist.append(c)
        if fulldebug: l.debug("Vert %d " % v.index + "is part of cell(s): " + iis)

    # Select all cell vertices
    for c in clist:
        for f in c.ugfaces:
            for v in f.ugverts:
                if ob.data.vertices[v.bi].select == False:
                    ob.data.vertices[v.bi].select=True
    
    # Return to original mode
    bpy.ops.object.mode_set(mode=mode)

    return len(clist)
