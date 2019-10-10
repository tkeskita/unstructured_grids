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

# Operations related to cell and face zones

import bpy
from . import ug
from . import ug_op
import logging
l = logging.getLogger(__name__)
fulldebug = False # Set to True if you wanna see walls of logging debug


def exist_face_zones():
    '''Return True if any face zones exist in UG data'''

    for z in ug.ugzones:
        if z.deleted:
            continue
        if z.zonetype == 'face':
            return True
    return False


class UG_OT_EditFaceZoneOrientations(bpy.types.Operator):
    '''Start Editing Face Orientations of a Face Zone'''
    bl_idname = "unstructured_grids.facezone_edit_face_orientations"
    bl_label = "Edit Face Orientations of a Face Zone (UG)"

    @classmethod
    def poll(cls, context):
        return context.mode in {'OBJECT','EDIT_MESH'} and ug.exists_ug_state()

    def execute(self, context):
        zone = get_selected_face_zone(self)
        if zone == None:
            return {'FINISHED'}
        edit_face_zone(zone)
        self.report({'INFO'}, "Orientation object created. " \
                    + "Click 'Finish Editing' to save changes.")
        return {'FINISHED'}


def get_selected_face_zone(self):
    '''Return UGZone corresponding to face zone selection UI option'''

    ug_props = bpy.context.scene.ug_props
    i = 0
    for z in ug.ugzones:
        if z.deleted:
            continue
        if z.zonetype == 'face':
            i += 1
            if i == ug_props.facezone_selection:
                return z

    self.report({'ERROR'}, "No face zone order number " \
                + "%d" % ug_props.facezone_selection)
    return None


def edit_face_zone(zone):
    '''Create object for face zone orientation editing'''

    # Create mesh from faces in face zone
    import bmesh
    bm = bmesh.new()
    ob = ug.get_ug_object()
    bpy.ops.object.mode_set(mode='OBJECT')

    # First generate index dictionary for new face vertices (vmap)
    # and list of vertices (vlist)
    vmap = dict()
    vlist = []
    vi = 0
    for uf in zone.ugfaces:
        verts = [] # vertex index list
        for ugv in uf.ugverts:
            verts.append(ugv.bi)
        for v in verts:
            if v not in vlist:
                vmap[v] = vi
                vlist.append(v)
                vi += 1

    # Generate BMVerts
    for v in vlist:
        bm.verts.new(ob.data.vertices[v].co)
    bm.verts.ensure_lookup_table()
    bm.verts.index_update()

    # Generate BMFaces
    for i in range(len(zone.ugfaces)):
        uf = zone.ugfaces[i]
        verts = [] # vertex index list
        for ugv in uf.ugverts:
            verts.append(ugv.bi)
        bmverts = [] # bmesh (new) vertex list
        for v in verts:
            bmverts.append(bm.verts[vmap[v]])

        if fulldebug:
            text = "UGFace %s verts: " % str(uf)
            for v in verts:
                text += "%s " % str(v)
            l.debug(text)

        f = bm.faces.new(bmverts)
        # Take flipMap into account
        if zone.flipMap[i] == 1:
            f.normal_flip()
        f.normal_update()

    # Initialize Editing object
    ed_ob_name = 'Face Zone Orientation'
    bpy.ops.object.select_all(action='DESELECT')

    if ed_ob_name in bpy.data.objects:
        bpy.data.objects[ed_ob_name].select_set(True)
        mesh = bpy.data.objects[ed_ob_name].data
        bpy.ops.object.delete()
        bpy.data.meshes.remove(mesh)

    mesh_data = bpy.data.meshes.new(ed_ob_name)
    ed_ob = bpy.data.objects.new(ed_ob_name, mesh_data)
    bpy.context.scene.collection.objects.link(ed_ob)
    bpy.context.view_layer.objects.active = bpy.data.objects[ed_ob_name]
    ed_ob.select_set(True)
    ob.hide_set(True)
    bm.to_mesh(mesh_data)
    bpy.ops.object.mode_set(mode='EDIT')
    # Set "Face Orientation" overlay on
    bpy.context.area.spaces[0].overlay.show_face_orientation = True
    bpy.ops.mesh.select_mode(type="FACE") # Face selection mode
    bm.free()

class UG_OT_FinishFaceZoneOrientations(bpy.types.Operator):
    '''Finish Editing Face Orientations of a Face Zone'''
    bl_idname = "unstructured_grids.facezone_finish_face_orientations"
    bl_label = "Finish Editing Face Orientations of a Face Zone (UG)"

    @classmethod
    def poll(cls, context):
        return context.mode in {'OBJECT','EDIT_MESH'} and ug.exists_ug_state()

    def execute(self, context):
        zone = get_selected_face_zone(self)
        if zone == None:
            return {'FINISHED'}

        self.report({'INFO'}, "TODO")
        return {'FINISHED'}
