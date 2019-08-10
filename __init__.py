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

bl_info = {
    "name": "Unstructured Grids for Blender",
    "author": "Tuomo Keskitalo",
    "version": (0, 1, 0),
    "blender": (2, 80, 0),
    "location": "TBA",
    "description": "Import, Editing and Export of Unstructured Grids (3D Volume Meshes)",
    "warning": "WIP",
    "wiki_url": "https://github.com/tkeskita/unstructured_grids",
    "support": 'COMMUNITY',
    "category": "Mesh",
}


if "bpy" in locals():
    import importlib
    importlib.reload(ug)
    importlib.reload(io_polymesh)
    importlib.reload(ug_op)
else:
    import math
    import bpy
    from . import(
        ug,
        io_polymesh,
        ug_op,
        )
    from bpy.app.handlers import persistent

# Set up logging of messages using Python logging
# Logging is nicely explained in:
# https://code.blender.org/2016/05/logging-from-python-code-in-blender/
# To see debug messages, configure logging in file
# $HOME/.config/blender/{version}/scripts/startup/setup_logging.py
# add there something like:
# import logging
# logging.basicConfig(format='%(funcName)s: %(message)s', level=logging.DEBUG)
import logging
l = logging.getLogger(__name__)

# Common settings as property group
class UGProperties(bpy.types.PropertyGroup):
    export_path: bpy.props.StringProperty(
        name="Export Path",
        description="Path to Export Unstructured Grid",
        default="//",
        maxlen=1024,
        subtype="DIR_PATH",
    )
    text_boundary: bpy.props.StringProperty(
        name="PolyMesh Boundary File Contents",
        description="PolyMesh Boundary File Contents",
        default="",
        maxlen=0,
    )
    text_faces: bpy.props.StringProperty(
        name="PolyMesh Faces File Contents",
        description="PolyMesh Faces File Contents",
        default="",
        maxlen=0,
    )
    text_neighbour: bpy.props.StringProperty(
        name="PolyMesh Neighbour File Contents",
        description="PolyMesh Neighbour File Contents",
        default="",
        maxlen=0,
    )
    text_owner: bpy.props.StringProperty(
        name="PolyMesh Owner File Contents",
        description="PolyMesh Owner File Contents",
        default="",
        maxlen=0,
    )
    text_points: bpy.props.StringProperty(
        name="PolyMesh Points File Contents",
        description="PolyMesh Points File Contents",
        default="",
        maxlen=0,
    )
    generate_internal_edges: bpy.props.BoolProperty(
        name="Generate Edges for Internal Faces",
        description="Boolean for Generating Internal Face Edges",
        default=False,
    )

def menu_import(self, context):
    self.layout.operator(io_polymesh.UG_OT_ImportPolyMesh.bl_idname, \
                         text="OpenFOAM PolyMesh (UG)"
    )

def menu_export(self, context):
    self.layout.operator(io_polymesh.UG_OT_ExportPolyMesh.bl_idname, \
                         text="OpenFOAM PolyMesh (UG)"
    )

@persistent
def load_handler(dummy):
    '''Updates UG data from string variables after loading Blend file'''

    ug_props = bpy.context.scene.ug_props
    if len(ug_props.text_points) == 0:
        return None
    l.debug("Executing load_post handler")
    bpy.ops.unstructured_grids.polymesh_to_ug()

@persistent
def save_handler(dummy):
    '''Updates string variables from UG data before saving Blend file'''

    ug_props = bpy.context.scene.ug_props
    if len(ug.ugcells) == 0:
        return None
    l.debug("Executing pre_save handler")
    bpy.ops.unstructured_grids.ug_to_polymesh()



classes = (
    UGProperties,
    ug.UG_OT_UpdateBoundariesFromFaceMaterials,
    io_polymesh.UG_OT_ImportPolyMesh,
    io_polymesh.UG_OT_ExportPolyMesh,
    io_polymesh.UG_OT_PolyMeshToUG,
    io_polymesh.UG_OT_UGToPolyMesh,
    ug_op.UG_OT_Select_Cells_Inclusive,
    ug_op.UG_OT_Select_Cells_Exclusive,
)

def register():
    for cls in classes:
        bpy.utils.register_class(cls)

    bpy.types.Scene.ug_props = \
        bpy.props.PointerProperty(type = UGProperties)

    bpy.types.TOPBAR_MT_file_import.append(menu_import)
    bpy.types.TOPBAR_MT_file_export.append(menu_export)

    bpy.app.handlers.load_post.append(load_handler)
    bpy.app.handlers.save_pre.append(save_handler)

def unregister():
    for cls in classes:
        bpy.utils.unregister_class(cls)

    bpy.types.TOPBAR_MT_file_import.remove(menu_import)
    bpy.types.TOPBAR_MT_file_export.remove(menu_export)

    bpy.app.handlers.load_post.remove(load_handler)
    bpy.app.handlers.save_pre.remove(save_handler)

if __name__ == "__main__":
    register()
