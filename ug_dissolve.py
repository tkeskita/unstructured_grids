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

# Operations related to dissolving or reducing topology

import bpy
from . import ug
from . import ug_op
import logging
l = logging.getLogger(__name__)
fulldebug = False # Set to True if you wanna see walls of logging debug


class UG_OT_DissolveEdges(bpy.types.Operator):
    '''Dissolve Edges. Merge Selected Vertices Connected by Edges'''
    bl_idname = "unstructured_grids.dissolve_edges"
    bl_label = "Dissolve Edges (UG)"

    @classmethod
    def poll(cls, context):
        return context.mode == 'EDIT_MESH' and ug.exists_ug_state()

    def execute(self, context):
        n = dissolve_selected_edges()
        self.report({'INFO'}, "Dissolved %d edges" % n)
        return {'FINISHED'}


def dissolve_selected_edges():
    '''Dissolve Edges. Merge Selected Vertices Connected by Edges'''

    import bmesh
    bpy.ops.object.mode_set(mode='OBJECT')
    bpy.ops.object.mode_set(mode='EDIT')
    ob = ug.get_ug_object()
    bm = bmesh.from_edit_mesh(ob.data)
    verts = [v for v in bm.verts if v.select]

    processed_verts = [] # list of processed vertices
    n = 0
    for v in verts:
        # Check all edges for selected neighbour vertex
        for e in v.link_edges:
            ov = e.other_vert(v)
            if ov in verts and not ug.ugverts[ov.index].deleted:
                dissolve_vertex_pair(bm, v, ov)
                n += 1
    bm.normal_update()
    bmesh.update_edit_mesh(mesh=ob.data)
    bm.free()
    return n


def dissolve_vertex_pair(bm, v1, v2):
    '''Dissolve argument bmesh vertices in argument bmesh and UG data'''

    # Create new vertex at middle point
    bm.verts.new((v1.co + v2.co) / 2)
    bm.verts.ensure_lookup_table()
    bm.verts.index_update()

    ugv = ug.UGVertex()
    ugv1 = ug.ugverts[v1.index]
    ugv2 = ug.ugverts[v2.index]
    if fulldebug: l.debug("ugv1=%d ugv2=%d" % (ugv1.bi, ugv2.bi))

    def get_ugfaces_ugcells_from_vertices(vilist):
        '''Return lists of non-deleted UGFaces and UGCells that are part of
        argument list of vertex indices
        '''

        ugfaces = []
        for vi in vilist:
            for f in ug.ugverts[vi].ugfaces:
                if f.deleted:
                    continue
                if f in ugfaces:
                    continue
                ugfaces.append(f)

        ugcells = []
        for vi in vilist:
            for c in ug.ugverts[vi].ugcells:
                if c.deleted:
                    continue
                if c in ugcells:
                    continue
                ugcells.append(c)

        return ugfaces, ugcells

    ugfaces, ugcells = get_ugfaces_ugcells_from_vertices([v1.index, v2.index])

    if fulldebug:
        text = "Dissolve affects ugfaces: "
        for ugf in ugfaces:
            text += "%d " % ugf.bi
        text += "  and ugcells: "
        for c in ugcells:
            text += "%d " % c.ii
        l.debug(text)

    def replace_ugvertex_in_ugface(bm, ugf, old_ugv, ugv=None):
        '''Replace old UGVertex old_ugv with new UGVertex ugv for UGFace
        ugf. If no new UGvertex ugv is given, remove old UGVertex. New
        mesh face is created for boundary faces, and old face is
        hidden.
        '''

        i = ugf.ugverts.index(old_ugv) # Index of old vertex

        # Replace with new vertex
        if ugv:
            # Update face
            ugf.ugverts[i] = ugv
            # Update vertices
            old_ugv.remove_face(ugf)
            ugv.add_face(ugf)
            # Update cells
            if ugf.owner:
                ugf.owner.remove_vert(old_ugv)
                ugf.owner.add_vert(ugv)
            if ugf.neighbour:
                ugf.neighbour.remove_vert(old_ugv)
                ugf.neighbour.add_vert(ugv)

        # Remove old vertex
        else:
            # Update face
            ugf.ugverts.pop(i)
            # Update vertices
            old_ugv.remove_face(ugf)
            # Update cells
            if ugf.owner:
                ugf.owner.remove_vert(old_ugv)
            if ugf.neighbour:
                ugf.neighbour.remove_vert(old_ugv)

        # Create new boundary face if needed
        if ugf.is_boundary_face():
            old_fi = ugf.bi
            bmverts = [bm.verts[ugv.bi] for ugv in ugf.ugverts]
            fi = len(bm.faces)
            f = bm.faces.new(bmverts)
            f.normal_update()
            ugf.add_mesh_face(fi)
            bm.faces.ensure_lookup_table() # TODO: Circumvent somehow?
            bm.faces[old_fi].hide_set(True) # Hide old face

    def replace_face_vertices(bm, ugf, ugv, ugv1, ugv2):
        '''Replace old UGVerts ugv1 and/or ugv2 with new vertex ugv
        corresponding to BMVert v in UGFace ugf
        '''

        # Both vertices are part of this face, handle both
        if ugv1 in ugf.ugverts and ugv2 in ugf.ugverts:

            # If number of vertices of a face is 3, then delete face
            if len(ugf.ugverts) == 3:
                if fulldebug:
                    l.debug("Delete face (%d ugverts)" % len(ugf.ugverts))
                bm.faces[ugf.bi].hide_set(True) # Hide old face
                ugf.delete()
                return True

            # Otherwise replace first vertex with new and remove second
            if fulldebug: 
                l.debug("Face %d vertex" % ugf.bi \
                        + " %d->%d, remove %d" % (ugv1.bi, ugv.bi, ugv2.bi))
            replace_ugvertex_in_ugface(bm, ugf, ugv1, ugv)
            replace_ugvertex_in_ugface(bm, ugf, ugv2, None)

        # Only one vertex on this face, vertex gets replaced by new vertex
        else:
            old_ugv = None
            if ugv1 in ugf.ugverts:
                old_ugv = ugv1
            elif ugv2 in ugf.ugverts:
                old_ugv = ugv2
            if fulldebug: 
                l.debug("Face %d vertex" % ugf.bi \
                        + " %d->%d" % (old_ugv.bi, ugv.bi))
            replace_ugvertex_in_ugface(bm, ugf, old_ugv, ugv)

    for ugf in ugfaces:
        replace_face_vertices(bm, ugf, ugv, ugv1, ugv2)
