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

# Input/Output routines for OpenFOAM PolyMesh unstructured grids

import bpy
import logging
l = logging.getLogger(__name__)


class UG_OT_PolyMeshToUG(bpy.types.Operator):
    '''Generate UG data and mesh object from OpenFOAM PolyMesh text blocks'''
    bl_idname = "unstructured_grids.polymesh_to_ug"
    bl_label = "Generate UG from polyMesh texts"

    def execute(self, context):
        polymesh_to_ugdata(self)
        return {'FINISHED'}


def polymesh_to_ugdata(self):
    '''Convert OpenFOAM polyMesh data from Blender text blocks
    into UG data structures and Blender mesh
    '''

    ob = initialize_ug_object()
    verts = polymesh_get_verts('points')
    edges = [] # TODO
    faces = [] # TODO
    # Create vertices and faces into mesh object
    ob.data.from_pydata(verts, edges, faces)


def initialize_ug_object():
    '''Creates and returns an initialized and empty UG mesh object'''
    
    name = "Unstructured Grid"   
    if name in bpy.data.objects:
        l.debug("Delete existing object " + name)
        bpy.ops.object.select_all(action='DESELECT')
        bpy.data.objects[name].select_set(True)
        bpy.ops.object.delete()

    l.debug("Create and activate new mesh object " + name)
    mesh_data = bpy.data.meshes.new(name)
    ob = bpy.data.objects.new(name, mesh_data)
    bpy.context.scene.collection.objects.link(ob)
    bpy.context.view_layer.objects.active = bpy.data.objects[name]
    ob.select_set(True)
    return ob


def polymesh_get_verts(text_name):
    '''Creates list of vertex triplets from polymesh points text block'''

    import re
    verts = [] # list of x, y, z point coordinate triplets

    text_points = bpy.data.texts[text_name]
    for tl in text_points.lines:
        line = tl.body
        regex = re.search(r'^\(([dDeE\d\.\-]+)\s+([dDeE\d\.\-]+)\s+([dDeE\d\.\-]+)\)', line, re.M)
        if regex:
            x = float(regex.group(1))
            y = float(regex.group(2))
            z = float(regex.group(3))
            verts.append(tuple([x, y, z]))

    l.debug("Number of triplets read from %s: %d" % (text_name, len(verts)))
    return verts

