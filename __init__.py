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
    "version": (0, 2, 0),
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
    importlib.reload(ug_op_extrude)
else:
    import math
    import bpy
    from . import(
        ug,
        io_polymesh,
        ug_op,
        ug_op_extrude,
        )
    from bpy.app.handlers import persistent
    from sys import float_info

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
    text_cellZones: bpy.props.StringProperty(
        name="PolyMesh cellZones File Contents",
        description="PolyMesh cellZones File Contents",
        default="",
        maxlen=0,
    )
    text_faceZones: bpy.props.StringProperty(
        name="PolyMesh faceZones File Contents",
        description="PolyMesh faceZones File Contents",
        default="",
        maxlen=0,
    )
    generate_internal_edges: bpy.props.BoolProperty(
        name="Generate Edges for Internal Faces",
        description="Boolean for Generating Internal Face Edges",
        default=False,
    )
    extrusion_thickness: bpy.props.FloatProperty(
        name="Extrusion Thickness",
        description="Extrusion Thickness (Cell Side Length Perpendicular to Surface)",
        default=0.05,
        precision=4,
        min=float_info.min, max=float_info.max
    )
    extrusion_ignores_unselected_face_normals: bpy.props.BoolProperty(
        name="Ignore Unselected Face Normals",
        description="Ignore Unselected Neighbour Face Normals in Calculation of Extrusion Direction",
        default=True,
    )
    extrusion_uses_fixed_initial_directions: bpy.props.BoolProperty(
        name="Use Fixed Extrusion Directions",
        description="Use Initial Normal Directions for All Layers in Extrusion",
        default=True,
    )
    extrusion_layers: bpy.props.IntProperty(
        name="Extrusion Layers",
        description="Number of Layers to Extrude",
        default=1,
        min=1, max=10000000
    )
    extrusion_scale_thickness_expression: bpy.props.StringProperty(
        name="Layer Thickness (x) Scaling Expression",
        description="Python Expression to Scale Layer Thickness After Layer Addition",
        default="x*1.0",
        maxlen=0,
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
    bpy.ops.unstructured_grids.update_all_from_blender()


class VIEW3D_PT_UG_GUI:
    '''UG Tool Bar, common for object and edit modes'''
    bl_label = "Unstructured Grid"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'

    @classmethod
    def poll(cls, context):
        obj = context.active_object
        return obj and obj.type == 'MESH'

    def draw(self, context):
        layout = self.layout
        ug_props = context.scene.ug_props

        row = layout.row()
        row.label(text=ug.ug_print_stats())

        row = layout.row()
        row.operator("unstructured_grids.import_openfoam_polymesh", text="Import PolyMesh")
        row = layout.row()
        row.operator("unstructured_grids.update_all_from_blender", text="Update to Storage")
        row = layout.row()
        row.operator("unstructured_grids.polymesh_to_ug", text="Restore from Storage")
        row = layout.row()
        row.operator("unstructured_grids.reset_view", text="Reset View")
        row = layout.row()
        row.operator("unstructured_grids.export_openfoam_polymesh", text="Export PolyMesh")

        row = layout.row()
        row.label(text="Select Cells:")
        col = layout.column()
        rowsub = col.row(align=True)
        rowsub.operator("unstructured_grids.select_cells_exclusive", text="Exclusive")
        rowsub.operator("unstructured_grids.select_cells_inclusive", text="Inclusive")

        row = layout.row()
        row.operator("unstructured_grids.delete_cells", text="Delete Cells")

        row = layout.row()
        row.label(text="Extrusion Settings:")
        col = layout.column()
        rowsub = col.row(align=True)
        rowsub.prop(ug_props, "extrusion_layers", text="Layers")
        rowsub.prop(ug_props, "extrusion_ignores_unselected_face_normals",
                    icon='NORMALS_FACE', text="")
        rowsub.prop(ug_props, "extrusion_uses_fixed_initial_directions",
                    icon='NORMALS_VERTEX_FACE', text="")

        row = layout.row()
        row.prop(ug_props, "extrusion_thickness", text="Thickness")

        row = layout.row()
        row.label(text="Expression for Scaling Thickness:")
        row = layout.row()
        row.prop(ug_props, "extrusion_scale_thickness_expression", text="")

        row = layout.row()
        row.operator("unstructured_grids.extrude_cells", text="Extrude Cells")

        row = layout.row()
        row.label(text="Debug:")
        row = layout.row()
        row.operator("unstructured_grids.print_info_of_selected_cells", text="Print Cell Info")
        row = layout.row()
        row.operator("unstructured_grids.print_info_of_selected_faces", text="Print Face Info")

        # Object Mode warning
        if context.mode == 'OBJECT':
            box = layout.box()
            col = box.column(align=True)
            row = col.row(align=True)
            row.label(text="Note: Deletions are not shown")
            row = col.row(align=True)
            row.label(text="correctly in Object Mode")


class VIEW3D_PT_UG_GUI_Object(bpy.types.Panel, VIEW3D_PT_UG_GUI):
    '''UG Panel in Object Mode'''
    bl_category = "UG"
    bl_idname = "VIEW3D_PT_ug_object_mode"
    bl_context = "objectmode"


class VIEW3D_PT_UG_GUI_Edit(bpy.types.Panel, VIEW3D_PT_UG_GUI):
    '''UG Panel in Edit Mode'''
    bl_category = "UG"
    bl_idname = "VIEW3D_PT_ug_edit_mode"
    bl_context = "mesh_edit"


classes = (
    VIEW3D_PT_UG_GUI_Object,
    VIEW3D_PT_UG_GUI_Edit,
    UGProperties,
    ug.UG_OT_UpdateBoundariesFromFaceMaterials,
    ug.UG_OT_UpdateZonesFromVertexGroups,
    ug.UG_OT_UpdateUGAllFromBlender,
    ug.UG_OT_PrintSelectedCellsInfo,
    ug.UG_OT_PrintSelectedFacesInfo,
    ug.UG_OT_PrintSelectedVertexIndices,
    io_polymesh.UG_OT_ImportPolyMesh,
    io_polymesh.UG_OT_ExportPolyMesh,
    io_polymesh.UG_OT_PolyMeshToUG,
    ug_op.UG_OT_SelectCellsInclusive,
    ug_op.UG_OT_SelectCellsExclusive,
    ug_op.UG_OT_ResetView,
    ug_op.UG_OT_DeleteCells,
    ug_op_extrude.UG_OT_ExtrudeCells,
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
