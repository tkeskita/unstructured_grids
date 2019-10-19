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

# Input/Output routines for VTK Unstructured Grid XML (VTU) format.

# More information about VTK: https://vtk.org/
# More information about VTK file formats:
# https://lorensen.github.io/VTKExamples/site/VTKFileFormats/

# Note: VTK Unstructured Grid's don't seem to support patches
# (named boundaries) or zones.

import bpy
from bpy_extras.io_utils import (
    ImportHelper,
    ExportHelper,
    orientation_helper,
    axis_conversion,
)
from . import ug
from . import ug_op
import logging
l = logging.getLogger(__name__)

# Global variables
print_interval = 100000 # Debug print progress interval
fulldebug = True # Set to True if you wanna see walls of logging debug


##################
##### IMPORT #####
##################

class UG_OT_ImportVtu(bpy.types.Operator, ImportHelper):
    '''Import VTK Unstructured Grid (.vtu) Files into Blender as UG Data'''
    bl_idname = "unstructured_grids.import_vtu"
    bl_label = "Import VTK Unstructured Grid (.vtu) (UG)"
    filename_ext = ".vtu"

    @classmethod
    def poll(cls, context):
        return context.mode in {'OBJECT','EDIT_MESH'}

    def execute(self, context):
        text = read_vtu_files(self)
        res = validate_vtu(text)
        if res:
             self.report({'ERROR'}, "Validation failed: " + res)
             return {'FINISHED'}
        n = vtu_to_ugdata(text)
        ug.update_ug_all_from_blender(self)
        self.report({'INFO'}, "Imported %d cells" % n)
        return {'FINISHED'}


def read_vtu_files(self):
    '''Read VTU file contents into string variable'''

    import os
    ug_props = bpy.context.scene.ug_props
    filepath = self.filepath

    if not (os.path.isfile(filepath)):
        self.report({'ERROR'}, "Could not find %r" \
                    % filepath)
        return None

    with open(filepath, 'r') as infile:
        text = infile.read()

    return text


def validate_vtu(text):
    '''Validate VTU data is suitable for parsing'''

    import re
    rec1 = re.compile(r'VTKFile\ type\=\"(.*?)\"', re.M)
    rec2 = re.compile(r'\ format\=\"(.*?)\"', re.M)

    for line in text.splitlines():
        regex1 = rec1.search(line)
        if regex1:
            value = str(regex1.group(1))
            if value != "UnstructuredGrid":
                return "VTK file type is not 'UnstructuredGrid' but " + value
        regex2 = rec2.search(line)
        if regex2:
            value = str(regex2.group(1))
            if value != "ascii":
                return "VTK file format is not 'ascii' but " + value


def vtu_to_ugdata(text):
    '''Convert VTU data from text into UG data structures and Blender
    mesh
    '''

    import bmesh
    ug_props = bpy.context.scene.ug_props
    ug.hide_other_objects()
    ob = ug.initialize_ug_object()
    bm = bmesh.new()

    # Get lists of data from text
    points = get_data_array_block("Points", "float", text)
    connectivities = get_data_array_block("connectivity", "int", text)
    offsets = get_data_array_block("offsets", "int", text)
    celltypes = get_data_array_block("types", "int", text)
    cellfaces = get_data_array_block("faces", "int", text)
    cellfaceoffsets = get_data_array_block("faceoffsets", "int", text)

    l.debug("VTU data contains %d coordinates" % len(points) \
            + ", %d connectivities" % len(connectivities) \
            + ", %d offsets" % len(offsets) \
            + ", %d celltypes" % len(celltypes) \
            + ", %d cellfaces" % len(cellfaces) \
            + " and %d cellfaceoffsets" % len(cellfaceoffsets))

    create_points(bm, points)
    vtu_datalists_to_ugdata(connectivities, offsets, celltypes, \
                            cellfaces, cellfaceoffsets)
    bpy.ops.object.mode_set(mode='OBJECT')
    bm.to_mesh(ob.data)
    return len(celltypes)

def get_data_array_block(name, vartype, text):
    '''Return list of items (type vartype) in DataArray name from text'''

    import re
    rec1 = re.compile(r'\<DataArray\ .*Name\=\"(.*?)\"', re.M)
    rec2 = re.compile(r'\<\/DataArray\>', re.M)
    inside = False
    datablock = ""

    for line in text.splitlines():
        # Detect start of wanted section
        regex1 = rec1.search(line)
        if regex1:
            if str(regex1.group(1)) == name:
                inside = True
                continue
            else:
                inside = False

        # Detect end of wanted section
        regex2 = rec2.search(line)
        if regex2:
            inside = False

        # Add line if it is inside wanted section
        if inside:
            datablock += line

    return get_list_from_text(datablock, vartype)


def get_list_from_text(text, vartype):
    '''Creates list from argument text block by taking anything separated
    by spaces and then converting it to type vartype
    '''

    valuelist = []
    for value in text.split( ):
        value = value.strip()
        if value:
            command = vartype + "(" + value + ")"
            valuelist.append(eval(command))

    return valuelist


def create_points(bm, points):
    '''Create UGVerts and bmesh vertices from points coordinate list'''

    for i in range(0, len(points), 3):
        x = points[i]
        y = points[i+1]
        z = points[i+2]
        ug.UGVertex()
        bm.verts.new([x, y, z])
        if i % print_interval == 0:
            l.debug("... processed vertex count: %d" % i)

    bm.verts.ensure_lookup_table()
    bm.verts.index_update()
    l.debug("Number of vertices: %d" % (i/3 + 1))


def vtu_datalists_to_ugdata(connectivities, offsets, celltypes, cellfaces, cellfaceoffsets):
    '''Generate UGFaces and UGCells from VTU datalists'''

    conn_index = 0 # index for connectivities
    facemap = dict() # dictionary to map from vertex list string to UGFace index

    # Loop through all cells
    for ci in range(len(celltypes)):
        vtk_cell_type = celltypes[ci]
        conn_end = offsets[ci]
        vilist = connectivities[conn_index:conn_end] # Vertex index list
        if fulldebug: l.debug("Cell %d vertices: " % ci + str(vilist))
        # TODO: Add polyhedron stuff here

        facemap = vtu_add_cell(vtk_cell_type, vilist, facemap)
        conn_index += conn_end - conn_index # Increment connectivities index


def vtu_add_cell(vtk_cell_type, vilist, facemap):
    '''Add new cell of argument type number and vertex index list. New
    faces are added to face map
    '''

    def add_cell_faces(c, fis, vilist, facemap):
        '''Create and/or add faces to cell c using face index list (fis) and
        vertex indices in vilist
        '''
        for faceverts in fis:
            real_vilist = [vilist[v] for v in faceverts] # Actual vertex indices
            string = get_vert_string(real_vilist) # Get facemap key

            # Create new UGFace or use existing. Map to owner/neighbour.
            # If ugf is new, add to facemap.
            if string in facemap:
                ugf = facemap[string]
                ugf.neighbour = c
            else:
                ugf = ug.UGFace(real_vilist)
                ugf.owner = c
                facemap[string] = ugf
            # Add to cell
            c.add_face_and_verts(ugf)
        return facemap


    # Create new cell, and depending on VTK cell type, add faces
    c = ug.UGCell()

    if vtk_cell_type == 10: # VTK_TETRA
        fis = [[0,2,1], [0,1,3], [1,2,3], [0,3,2]]
        facemap = add_cell_faces(c, fis, vilist, facemap)
        if fulldebug: l.debug("Created cell %d: tetra" % c.ii)

    elif vtk_cell_type == 12: # VTK_HEXAHEDRON
        fis = [[0,4,5,1], [0,3,2,1], [5,6,2,1], [4,7,6,5], [0,3,7,4], [7,3,2,6]]
        facemap = add_cell_faces(c, fis, vilist, facemap)
        if fulldebug: l.debug("Created cell %d: hexa" % c.ii)

    elif vtk_cell_type == 13: # VTK_WEDGE (=prism)
        fis = [[0,1,2], [0,3,4,1], [1,4,5,2], [2,5,3,0], [3,5,4]]
        facemap = add_cell_faces(c, fis, vilist, facemap)
        if fulldebug: l.debug("Created cell %d: prism" % c.ii)

    elif vtk_cell_type == 14: # VTK_PYRAMID
        fis = [[0,3,2,1], [0,4,3], [3,4,2], [2,4,1], [1,4,0]]
        facemap = add_cell_faces(c, fis, vilist, facemap)
        if fulldebug: l.debug("Created cell %d: pyramid" % c.ii)

    elif vtk_cell_type == 15: # VTK_PENTAGONAL_PRISM
        fis = [[0,1,2,3,4], [0,5,6,1], [1,6,7,2], [2,7,8,3], [3,8,9,4], [4,9,5,0], [9,8,7,6,5]]
        facemap = add_cell_faces(c, fis, vilist, facemap)
        if fulldebug: l.debug("Created cell %d: pentaprism" % c.ii)

    elif vtk_cell_type == 16: # VTK_HEXAGONAL_PRISM
        fis = [[0,1,2,3,4,5], [0,6,7,1], [1,7,8,2], [2,8,9,3], [3,9,10,4], [4,10,11,5], [5,11,6,0], [11,10,9,8,7,6]]
        facemap = add_cell_faces(c, fis, vilist, facemap)
        if fulldebug: l.debug("Created cell %d: hexaprism" % c.ii)

    # TODO: elif vtk_cell_type == 42: # VTK_POLYHEDRON

    else:
        raise ValueError("Unsupported VTK cell type %d" % vtk_cell_type)

    return facemap

def get_vert_string(vilist):
    '''Return ordered vertex index list converted to a string'''

    vilist.sort()
    string = ""
    for i in vilist:
        string += str(i) + ","
    return string

##################
##### EXPORT #####
##################


class UG_OT_ExportVtu(bpy.types.Operator, ExportHelper):
    '''Export UG Data as Vtk Vtu Files'''
    bl_idname = "unstructured_grids.export_vtu"
    bl_label = "Export VTK Unstructured Grid (.vtu) (UG)"
    filename_ext = ".vtu"

    @classmethod
    def poll(cls, context):
        return context.mode in {'OBJECT','EDIT_MESH'} and ug.exists_ug_state()

    def execute(self, context):
        ug.update_ug_all_from_blender(self)
        # TODO: write_vtu(self)
        return {'FINISHED'}
