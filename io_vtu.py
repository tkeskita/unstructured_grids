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
        self.report({'INFO'}, "TODO: Import %d cells" % n)
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

    #ug_props = bpy.context.scene.ug_props
    #ug.hide_other_objects()
    #ob = ug.initialize_ug_object()

    # Get lists of data from text
    points = get_data_array_block("Points", "float", text)
    connectivity = get_data_array_block("connectivity", "int", text)
    offsets = get_data_array_block("offsets", "int", text)
    celltypes = get_data_array_block("types", "int", text)
    cellfaces = get_data_array_block("faces", "int", text)
    cellfaceoffsets = get_data_array_block("faceoffsets", "int", text)

    l.debug("VTU data contains %d coordinates" % len(points) \
            + ", %d connectivities" % len(connectivity) \
            + ", %d offsets" % len(offsets) \
            + ", %d celltypes" % len(celltypes) \
            + ", %d cellfaces" % len(cellfaces) \
            + " and %d cellfaceoffsets" % len(cellfaceoffsets))

    #bpy.ops.object.mode_set(mode='OBJECT')
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
