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

# Input/Output routines for OpenFOAM PolyMesh unstructured grids.
# More information about PolyMesh:
# https://cfd.direct/openfoam/user-guide/mesh-description/

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
    [edges, faces] = polymesh_get_faces('owner', 'neighbour', 'faces')
    # TODO: boundary
    # Create vertices and faces into mesh object
    ob.data.from_pydata(verts, edges, faces)
    ob.data.validate()


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


def polymesh_get_verts(name):
    '''Creates list of vertex triplets from PolyMesh points text block'''

    import re
    verts = [] # list of x, y, z point coordinate triplets

    if not name in bpy.data.texts:
        l.error("text not found: %s" % name)
        return verts

    text_name = bpy.data.texts[name]
    for tl in text_name.lines:
        line = tl.body
        regex = re.search(r'^\(([dDeE\d\.\-]+)\s+([dDeE\d\.\-]+)\s+([dDeE\d\.\-]+)\)', line, re.M)
        if regex:
            x = float(regex.group(1))
            y = float(regex.group(2))
            z = float(regex.group(3))
            verts.append(tuple([x, y, z]))

    l.debug("Number of triplets read from %s: %d" % (name, len(verts)))
    return verts


def polymesh_get_faces(owner_name, neighbour_name, faces_name):
    '''Creates edge and face list from PolyMesh owner, neighbour and
    faces text blocks
    '''

    edges = [] # List of edge vertex index pairs, to be generated
    faces = [] # List of face vertex index lists, to be generated

    # Read in owner and neighbour lists
    owner = polymesh_get_intlist(owner_name)
    neighbour = polymesh_get_intlist(neighbour_name)
    face_verts = polymesh_get_list_intlist(faces_name)

    # Create faces at boundary and only edges for internal faces
    i = 0
    for verts in face_verts:
        if i < len(neighbour):
            # Internal face, add edges
            for j in range(len(face_verts[i])):
                edges.append(tuple([face_verts[i][j-1], face_verts[i][j]]))
        else:
            # Boundary face, add faces
            faces.append(tuple(face_verts[i]))
        i += 1

    l.debug("Number of edge index pairs generated from %s: %d" % (faces_name, len(edges)))
    l.debug("Number of boundary face index lists generated from %s: %d" % (faces_name, len(faces)))

    return edges, faces


def polymesh_get_intlist(name):
    '''Creates integer list from argument PolyMesh integer text block'''

    import re
    iList = [] # list of integers to be generated
    inside = False # boolean for marking integer list in text

    if not name in bpy.data.texts:
        l.error("text not found: %s" % name)
        return iList

    text_name = bpy.data.texts[name]
    for tl in text_name.lines:
        line = tl.body

        # Opening of integer list by single parenthesis
        regex = re.search(r'^\(', line, re.M)
        if regex:
            inside = True

        # Closing of integer list by single parenthesis
        regex2 = re.search(r'^\)', line, re.M)
        if regex2:
            inside = False

        # Integer, at start of line
        regex3 = re.search(r'^(\d+)', line, re.M)
        if inside and regex3:
            iList.append(int(regex3.group(1)))

    l.debug("Number of integers read from %s: %d" % (name, len(iList)))
    return iList


def polymesh_get_list_intlist(name):
    '''Creates list of integer lists from argument PolyMesh
    text block
    '''

    # TODO: Get rid of code duplication

    import re
    iList = [] # list of integers lists to be generated
    inside = False # boolean for marking integer list in text

    if not name in bpy.data.texts:
        l.error("text not found: %s" % name)
        return iList

    text_name = bpy.data.texts[name]
    for tl in text_name.lines:
        line = tl.body

        # Opening of integer list by single parenthesis
        regex = re.search(r'^\(', line, re.M)
        if regex:
            inside = True

        # Closing of integer list by single parenthesis
        regex2 = re.search(r'^\)', line, re.M)
        if regex2:
            inside = False

        # List of integer list within parenthesis
        regex3 = re.search(r'^\d+\(([\d\s]+)\)', line, re.M)
        if inside and regex3:
            vals = regex3.group(1).split()
            valList = []
            for val in vals:
                valList.append(int(val))
            iList.append(valList)

    l.debug("Number of integer lists read from %s: %d" % (name, len(iList)))
    return iList
