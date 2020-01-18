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

# Operations related to UG data integrity checking and statistics

import bpy
from . import ug
from . import ug_op
import logging
l = logging.getLogger(__name__)


class UG_OT_CheckCells(bpy.types.Operator):
    '''Check selected cells for integrity, statistics etc'''

    bl_idname = "unstructured_grids.check_cells"
    bl_label = "Check Selected Cells (UG)"

    @classmethod
    def poll(cls, context):
        return context.mode == 'EDIT_MESH' and ug.exists_ug_state()


    def execute(self, context):
        nerr = 0 # number of errors
        ctext = "" # cell check result text
        tot_vol = 0.0 # Total volume of selected cells

        clist = ug_op.select_cells_exclusive()
        for c in clist:
            rval, rtext, vol, area = check_cell_integrity(c)
            nerr += rval
            ctext += rtext
            tot_vol += vol

        text = "Unstructured Grid - result of Check Selected Cells:\n\n"
        text += "Total cell count=%d" % len(clist)
        text += ", issue count=%d" % nerr
        text += ", total volume=%.5g" % tot_vol
        text += "\n\n"
        text += ctext
        set_text_to_text_block(text)
        bpy.ops.object.mode_set(mode='EDIT')

        self.report({'INFO'}, "Issue count: %d /" % nerr \
                    + " %d cells." % len(clist) \
                    + " See 'UG' Text Block.")
        return {'FINISHED'}


def check_cell_integrity(c):
    '''Check integrity of argument cell and print out issues and stats.
    First return value is 0 if no issues were found and 1 if something 
    is seriously wrong. Second return value is the text string of the 
    analysis that was made. Rest are cell statistics.
    '''

    import bmesh

    if c.deleted:
        return 1, "ERROR: Cell %d has been deleted" % c.ii

    ob = ug.get_ug_object()
    bpy.ops.object.mode_set(mode='OBJECT')
    bm = bmesh.new()    
    bm = add_cell_faces_to_bmesh(c, bm, ob)

    test1, text1, vol, area = check_cell_closedness_area_volume(c, bm)

    # TODO: Add intersection and other checks?

    # Concatenate tests and texts and return
    bm.free()
    errors_found = test1
    text = text1
    l.debug(text.strip("\n"))

    if errors_found:
        return 1, text, vol, area
    else:
        return 0, text, vol, area


def add_cell_faces_to_bmesh(c, bm, ob):
    '''Add argument cell c faces to argument bmesh. 
    Faces are flipped so that normals point out of cell.
    '''

    # Create bmverts
    vertmap = dict() # dictionary to map from UGVert to new bmvert
    for ugv in c.ugverts:
        # Note: Assumes that there is no vertex in this location already
        v = bm.verts.new(ob.data.vertices[ugv.bi].co)
        vertmap[ugv] = v

    # Create bmfaces
    for ugf in c.ugfaces:
        bmverts = [vertmap[ugv] for ugv in ugf.ugverts]
        f = bm.faces.new(bmverts)
        if ugf.neighbour == c:
            f.normal_flip()
        f.normal_update()

    bm.faces.ensure_lookup_table()
    return bm


def check_cell_closedness_area_volume(c, bm):
    '''Calculate area and volume of cell c in UG object ob.
    First return value is True if errors were found, False otherwise.
    Second return value is text of the test results.
    Third return value is volume and fourth is total area of faces.
    '''

    from mathutils import Vector

    face_areas = [] # Face areas [m2]
    face_awns = [] # Face area weighted normal vectors
    for f in bm.faces:
        area = f.calc_area()
        face_areas.append(area)
        face_awns.append(area * f.normal)

    vol = bm.calc_volume()
    area = sum(face_areas)

    text = "Cell %d:" % c.ii
    text += " volume=%.5f" % vol
    text += " face count=%d" % len(c.ugfaces)
    text += " area=%.5f" % area

    # Calculated cell closedness
    csum = Vector((0, 0, 0))
    for v in face_awns:
        csum = csum + v
    cl = abs(csum.x) + abs(csum.y) + abs(csum.z)

    TOL = 1e-6 # tolerance for closedness
    if cl < TOL:
        text += " closedness=OK\n"
        return False, text, vol, area
    else:
        text += " closedness=ERROR (%.3f)\n" % cl
        return True, text, vol, area


def set_text_to_text_block(text):
    '''Set argument text into Blender Text Block, to be available in
    Blender Text Editor.
    '''

    name = "UG" # name of text

    # Get or create text block, set text
    if name not in bpy.data.texts.keys():
        t = bpy.data.texts.new(name)
    else:
        t = bpy.data.texts[name]
    t.from_string(text)

    # Set data into Text Editor
    areas = bpy.context.screen.areas
    for area in areas:
        if area.type != 'TEXT_EDITOR':
            continue
        for space in area.spaces:
            if space.type != 'TEXT_EDITOR':
                continue
            space.text = t
            space.top = 0


class UG_OT_PrintSelectedCellsInfo(bpy.types.Operator):
    '''Print information about selected cells'''
    bl_idname = "unstructured_grids.print_info_of_selected_cells"
    bl_label = "Print Selected UG Cell Info"

    @classmethod
    def poll(cls, context):
        return context.mode == 'EDIT_MESH' and ug.exists_ug_state()

    def execute(self, context):
        clist = ug_op.select_cells_exclusive()
        text = "Info for selected %d cells\n\n" % len(clist)
        for c in clist:
            text += ug_print_cell_info(c)
        set_text_to_text_block(text)
        l.debug(text.strip("\n"))
        self.report({'INFO'}, "%d cell infos printed." % len(clist) \
                    + " See 'UG' Text Block.")
        return {'FINISHED'}


def ug_print_cell_info(c):
    '''Print information about argument cell'''

    text = "Cell %d " % c.ii
    if c.deleted:
        text += "(DELETED) "
    text += "contains %d UGFaces " % len(c.ugfaces)
    text += "and %d UGVerts\n" % len(c.ugverts)

    text += "  mesh face indices: "
    for f in c.ugfaces:
        text += "%d " % f.bi
    text += "\n"

    text += "  face ownership: "
    for f in c.ugfaces:
        if f.owner == c:
            text += "O"
        if f.neighbour == c:
            text += "N"
        text += " "
    text += "\n"

    text += "  mesh vertex indices: "
    for v in c.ugverts:
        text += "%d " % v.bi
    text += "\n"
    return text


class UG_OT_PrintSelectedFacesInfo(bpy.types.Operator):
    '''Print information about selected faces'''
    bl_idname = "unstructured_grids.print_info_of_selected_faces"
    bl_label = "Print Selected UG Face Info (via Python Logging)"

    @classmethod
    def poll(cls, context):
        return context.mode == 'EDIT_MESH' and ug.exists_ug_state()

    def execute(self, context):
        ob = ug.get_ug_object()
        mode = ob.mode # Save original mode
        # Visit object mode to update selection
        bpy.ops.object.mode_set(mode='OBJECT')
        bpy.ops.object.mode_set(mode='EDIT')

        vilist = [v.index for v in ob.data.vertices if v.select]
        ugflist = ug_op.get_ugfaces_from_vertices_exclusive(vilist)
        text = "Info for selected %d faces\n\n" % len(ugflist)
        for ugf in ugflist:
            text += ug_print_face_info(ugf)

        set_text_to_text_block(text)
        l.debug(text.strip("\n"))
        self.report({'INFO'}, "%d face infos printed." % len(ugflist) \
                    + " See 'UG' Text Block.")
        return {'FINISHED'}


def ug_print_face_info(ugf):
    '''Print information about argument UGFace'''

    if ugf.bi == -1:
        text = "Internal face "
        if ugf.ei == -1:
            text += "(no indices)"
        else:
            text += "(export index %d)" % ugf.ei
    else:
        text = "Boundary face "
        if ugf.ei == -1:
            text += "(mesh index %d)" % ugf.bi
        else:
            text += "(mesh index %d, export index %d)" % (ugf.bi, ugf.ei)

    if ugf.deleted:
        text += " (DELETED)"
    text += " contains %d UGVerts: " % len(ugf.ugverts)

    for v in ugf.ugverts:
        text += "%d " % v.bi

    text += "owner: "
    if ugf.owner != None:
        text += "%d " % ugf.owner.ii
    else:
        text += "None "

    text += "neighbour: "
    if ugf.neighbour != None:
        text += "%d" % ugf.neighbour.ii
    else:
        text += "None"
    text += "\n"

    if ugf.bi in ug.facemap:
        if ug.facemap[ugf.bi] != ugf:
            text += "  ERROR: wrong facemap[%d] point to %d\n" % (ugf.bi, ug.facemap[ugf.bi].bi)

    return text


class UG_OT_PrintSelectedVertexIndices(bpy.types.Operator):
    '''Debug print indices of selected vertices'''
    bl_idname = "unstructured_grids.print_selected_vertex_indices"
    bl_label = "Print Selected Vertex Indices"

    @classmethod
    def poll(cls, context):
        return context.mode in {'EDIT_MESH'}

    def execute(self, context):
        n = print_selected_vertex_indices()
        self.report({'INFO'}, "%d vertex infos printed." % n \
                    + " See 'UG' Text Block.")
        return {'FINISHED'}


def print_selected_vertex_indices():
    '''Debug print indices of selected vertices'''

    ob = bpy.context.active_object
    bpy.ops.object.mode_set(mode = 'OBJECT')
    bpy.ops.object.mode_set(mode = 'EDIT')
    verts = [v for v in ob.data.vertices if v.select]

    text = "Selected vertices: "
    for v in verts:
        text += "%d " % v.index

    set_text_to_text_block(text)
    l.debug(text.strip("\n"))
    return len(verts)


class UG_OT_PrintEdgeStatsText(bpy.types.Operator):
    '''Generate Edge Statistics of Mesh Faces (Minimum / Average / Maximum Length)'''
    bl_idname = "unstructured_grids.update_edge_stats_text"
    bl_label = "Update Edge Statistics Text"

    @classmethod
    def poll(cls, context):
        return context.mode in {'OBJECT','EDIT_MESH'}

    def execute(self, context):
        ug_props = bpy.context.scene.ug_props
        min_len, mean_len, max_len = get_edge_stats()

        if min_len == 0.0:
            text = "Error: No faces selected."
            set_text_to_text_block(text)
            self.report({'ERROR'}, text)
            return {'FINISHED'}

        text = "%.6e / %.6e / %.6e" % (min_len, mean_len, max_len)
        text = "Edge Stats of Selected Mesh Faces: " + text
        set_text_to_text_block(text)
        l.debug(text.strip("\n"))
        self.report({'INFO'}, "Edge stats ready. See 'UG' Text Block.")
        return {'FINISHED'}


def get_edge_stats():
    '''Return minimum, average and maximum lengths of
    selected edges in currently selected object mesh faces
    '''

    import bmesh
    if not bpy.context.active_object:
        return 0.0, 0.0, 0.0

    ob = bpy.context.active_object
    mode = ob.mode

    if mode == 'OBJECT':
        bpy.ops.object.mode_set(mode='EDIT')

    bm = bmesh.from_edit_mesh(ob.data)
    if not bm.faces:
        return 0.0, 0.0, 0.0

    faces = [f for f in bm.faces if f.select]
    min_len, mean_len, max_len = get_edge_stats_from_bmesh_faces(faces)
    # Return to original mode
    bpy.ops.object.mode_set(mode=mode)

    return min_len, mean_len, max_len


def get_edge_stats_from_bmesh_faces(faces):
    '''Return minimum, average and maximum lengths of
    selected edges in argument bmesh faces
    '''

    from sys import float_info
    min_len = float_info.max # minimum length
    max_len = 0.0 # maximum length
    mean_len = 0.0 # average length

    processed_edges = set()
    for f in faces:
        for e in f.edges:
            if e in processed_edges:
                continue
            processed_edges.add(e)

            # Edge length
            co0 = e.verts[0].co
            co1 = e.verts[1].co
            vec = co1 - co0
            elen = vec.length

            if elen < min_len:
                min_len = elen
            if elen > max_len:
                max_len = elen
            mean_len += elen

    if len(processed_edges) == 0:
        return 0.0, 0.0, 0.0

    mean_len /= float(len(processed_edges))
    return min_len, mean_len, max_len
