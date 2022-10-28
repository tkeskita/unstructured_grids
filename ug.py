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

# Core classes and routines for Unstructured Grids

import bpy
from mathutils import Vector
from . import ug_op
from . import io_polymesh
import logging
l = logging.getLogger(__name__)

obname = "Unstructured Grid" # Name for the Blender object
ugcells = [] # global list of all UGCells
ugfaces = [] # global list of all UGFaces
ugverts = [] # global list of all UGVerts
ugboundaries = [] # global list of all UGBoundaries
ugzones = [] # global list of all UGFaceZones and UGCellZones
facemap = dict() # global dictionary to get from mesh face index to an UGFace
fulldebug = False # Set to True if you wanna see walls of logging debug

# Note about facemap: Every UGVertex has a corresponding mesh vertex,
# and their indices match. This makes it trivial to map between UGVert
# and mesh vertex by the index number. However, same is NOT true for
# faces. Internal UGFaces usually don't have (but may in some cases
# have) a corresponding mesh face with same index number, since only
# boundary faces are created to the mesh. That makes mesh handling and
# mesh import faster, but it created a need to map from mesh face
# index to UGFace. Searching through all ugfaces for right index
# number takes too long for large meshes. Facemap returns UGFace much
# faster.

class UGCell:
    '''Class for Unstructured Grid cell data objects'''

    deleted = False # boolean for marking deleted cells
    ugfaces = [] # list ugfaces that make up this cell
    ugverts = [] # list of ugverts that make up this cell
    ei = -1 # Cell index (used at export only)
    ii = -1 # Internal index (used for debugging)

    # Cell operations:

    def __init__(self):
        '''Initialize a new UGCell'''
        global ugcells
        self.ugfaces = []
        self.ugverts = []
        self.ii = len(ugcells)
        ugcells.append(self)

    def add_face_and_verts(self, ugf):
        '''Add argument UGFace ugf and it's UGVerts to this cell'''
        if ugf in self.ugfaces:
            return None
        self.ugfaces.append(ugf)
        for ugv in ugf.ugverts:
            self.add_vert(ugv)

    def remove_face_and_verts(self, ugf):
        '''Remove argument UGFace ugf and it's UGVerts from this cell'''
        if ugf in self.ugfaces:
            self.ugfaces.remove(ugf)
            for ugv in ugf.ugverts:
                self.remove_vert(ugv)

    def add_vert(self, ugv):
        '''Add argument UGVertex ugv to this cell'''
        if ugv in self.ugverts:
            return None
        self.ugverts.append(ugv)

    def remove_vert(self, ugv):
        '''Remove argument UGVertex ugv from this cell. Remove only if no
        faces are using UGVertex
        '''
        if ugv in self.ugverts:
            ugflist = [ugf for ugf in self.ugfaces if ugv in ugf.ugverts]
            if len(ugflist) == 0:
                self.ugverts.remove(ugv)
        if len(self.ugverts) == 0:
            self.deleted = True

    def delete(self):
        '''Delete this cell'''
        for ugf in self.ugfaces:
            ugf.remove_cell(self)
        for ugv in self.ugverts:
            ugv.remove_cell(self)
        self.ugfaces = []
        self.ugverts = []
        self.deleted = True


class UGFace:
    '''Class for Unstructured Grid face data objects'''

    deleted = False # boolean for marking deleted faces
    owner = None # UGCell owning this face
    neighbour = None # UGCell neighbouring this face
    bi = -1 # Blender face index corresponding to this UGFace
            # Note: mapping from mesh face to ugface does NOT work
            # like ugfaces[face.index]. Instead, use facemap[face.index]
    ei = -1 # Face index (used at export only)
    mati = -1 # Blender material slot index number
    ugverts = [] # list of ordered ugverts used by this face

    def __init__(self, verts=[]):
        '''Initialize UGFace with optional vertex index list'''
        global ugfaces
        self.ugverts = []
        for v in verts:
            self.ugverts.append(ugverts[v])
        ugfaces.append(self)

    def add_verts(self, ugverts):
        '''Add argument UGVerts to this UGFace'''
        for ugv in ugverts:
            if ugv not in self.ugverts:
                self.ugverts.append(ugv)

    def remove_verts(self, ugverts):
        '''Add argument UGVerts to this UGFace'''
        for ugv in ugverts:
            if ugv in self.ugverts:
                self.ugverts.remove(ugv)

    # Note: no add_owner/neighbour, just access those directly.

    def remove_cell(self, c):
        '''Remove argument UGCell c from this UGFace'''
        if self.neighbour == c:
            self.neighbour = None
        elif self.owner == c:
            self.owner = None
            self.deleted = True

    def add_mesh_face(self, fi):
        '''Add mesh face index fi to this UGFace. Also update facemap'''
        global facemap
        self.bi = fi
        facemap[fi] = self

    # Note: No remove_mesh_face, that is not supported (no mesh
    # elements are really deleted)

    def delete(self):
        '''Delete this UGFace'''
        # Remove face from owner and neighbour cells
        if self.owner:
            self.owner.remove_face_and_verts(self)
        if self.neighbour:
            self.neighbour.remove_face_and_verts(self)
        # Remove face from UGVerts
        for ugv in self.ugverts:
            ugv.remove_face(self)
        self.ugverts = []
        self.owner = None
        self.neighbour = None
        self.deleted = True

    def invert_face_dir(self, switch=True):
        '''Invert face direction. Swap face owner and neighbour cells
        (done only if switch argument is true), and reverse face ugverts
        to mirror PolyMesh face normal vector.
        '''
        if switch:
            c = self.owner
            self.owner = self.neighbour
            self.neighbour = c
        self.ugverts.reverse() # Nicer than list(reversed(self.ugverts))

    def is_boundary_face(self):
        '''Return True if this UGFace is a boundary face, otherwise False'''
        if self.neighbour == None:
            return True
        else:
            return False


class UGVertex:
    '''Class for Unstructured Grid vertex data object'''

    deleted = False # boolean for marking deleted vertex
    bi = -1 # Blender vertex index corresponding to this UGVertex
            # Note: mapping between mesh vertices and ugverts should
            # always work in both ways, like:
            # ug.ugverts[vert.index] or: bm.verts[ugvert.bi]
    ei = -1 # Vertex index (used at export only)
    ugfaces = [] # List of ugfaces that this vertex is part of
    ugcells = [] # List of ugcells that this vertex is part of

    def __init__(self):
        '''Initialize a new UGVertex'''
        global ugverts
        self.bi = len(ugverts)
        self.ugfaces = []
        self.ugcells = []
        ugverts.append(self)

    def add_face(self, ugf):
        '''Add argument UGFace ugf to this UGVertex'''
        if ugf in self.ugfaces:
            return None
        self.ugfaces.append(ugf)
        self.deleted = False

    def remove_face(self, ugf):
        '''Remove argument UGFace ugf from this UGVertex'''
        if ugf in self.ugfaces:
            self.ugfaces.remove(ugf)
        if len(self.ugfaces) == 0:
            self.deleted = True

    def add_cell(self, c):
        '''Add argument UGCell c to this UGVertex'''
        if c in self.ugcells:
            return None
        self.ugcells.append(c)
        self.deleted = False

    def remove_cell(self, c):
        '''Remove argument UGCell c from this UGVertex'''
        if c in self.ugcells:
            self.ugcells.remove(c)
        if len(self.ugcells) == 0:
            self.deleted = True



class UGBoundary:
    '''Class for Unstructured Grid patch (boundary face) objects'''

    deleted = False # boolean for patches which contain zero boundary faces
    patchname = 'default' # patch name
    typename = 'default' # patch type name
    inGroups = '' # name for patch group this patch is part of
    nFaces = 0 # number of faces in patch
    startFace = -1 # index to first face for this patch (used at export only)
    mati = -1 # Blender object material slot index number of this boundary patch
    # TODO: Remove mati? This is redunant since index in ugboundaries matches slot index
    ugfaces = [] # list of UGFaces that are part of this boundary patch (updated after changes)

    def __init__(self, patchname):
        '''Initialize UGBoundary with name patchname'''
        global ugboundaries
        self.ugfaces = []
        self.patchname = patchname
        self.mati = len(ugboundaries)
        ugboundaries.append(self)


class UGZone:
    '''Class for Unstructured Grid face and cell zone objects.
    Face zones use ugfaces list and cell zones ugcells list.
    '''

    deleted = False # boolean for patches which contain zero faces
    zonetype = 'cell' # zone type name: cell or face
    zonename = 'default' # zone name
    ugfaces = [] # list of UGFaces that are part of face zone
    ugcells = [] # list of UGCells that are part of cell zone
    flipMap = [] # storage for face flipMap

    def __init__(self, zonetype, zonename):
        '''Initialize new zone with name zone type and zone name'''
        global ugzones
        self.ugfaces = []
        self.ugcells = []
        self.flipMap = []
        self.zonetype = zonetype
        self.zonename = zonename
        ugzones.append(self)


##########################
##### HELP FUNCTIONS #####
##########################

def hide_other_objects():
    '''Hides and deselects all other objects except UG object from view'''

    for ob in bpy.data.objects:
        ob.select_set(False)
        ob.hide_set(True)
    ob = get_ug_object()
    if ob:
        ob.hide_set(False)
        ob.select_set(True)


def delete_ug_object():
    '''Delete UG mesh object'''

    bpy.ops.object.mode_set(mode='OBJECT')
    l.debug("Delete existing object " + obname)
    bpy.ops.object.select_all(action='DESELECT')
    bpy.data.objects[obname].select_set(True)
    mesh = bpy.data.objects[obname].data
    bpy.ops.object.delete()
    bpy.data.meshes.remove(mesh)


def remove_ug_object_and_data():
    '''Remove UG mesh object and clear UG data'''

    global ugverts, ugfaces, ugcells, ugboundaries, ugzones

    # Zero UG data
    ugverts = []
    ugfaces = []
    ugcells = []
    ugboundaries = []
    ugzones = []

    # Delete UG object
    ob = get_ug_object()
    if ob:
        delete_ug_object()


def initialize_ug_object():
    '''Creates and returns an initialized and empty UG mesh object and
    initializes UG data
    '''

    remove_ug_object_and_data()
    l.debug("Create and activate new mesh object " + obname)
    mesh_data = bpy.data.meshes.new(obname)
    ob = bpy.data.objects.new(obname, mesh_data)
    bpy.context.scene.collection.objects.link(ob)
    bpy.context.view_layer.objects.active = bpy.data.objects[obname]
    ob.select_set(True)

    return ob


def get_ug_object():
    '''Returns UG Blender object or None if it does not exist'''

    if obname in bpy.data.objects:
        return bpy.data.objects[obname]
    return None


def exists_ug_state():
    '''Return True if UG Data (up-to-date or not) exists in memory,
    False otherwise.
    '''

    global ugcells
    test = False
    for c in ugcells:
        if c.deleted == False:
            test = True
            break
    ug_props = bpy.context.scene.ug_props
    if (len(ug_props.text_points) < 1):
        test = False
    if not get_ug_object():
        test = False
    return test


def ug_print_stats():
    '''Print short text with statistics of UG Data'''

    global ugcells

    if not exists_ug_state():
        return "No Cells. Import or Extrude to begin."
    clist = [c for c in ugcells if c.deleted == False]
    text = "Cells: " + str(len(clist))
    return text


class UG_OT_UpdateBoundariesFromFaceMaterials(bpy.types.Operator):
    '''Run after any changes to material slot or face material assignments.
    Materials are used to define boundary patch name and face assignments.
    '''
    bl_idname = "unstructured_grids.update_ugboundaries"
    bl_label = "Update UG Boundaries From Materials"

    @classmethod
    def poll(cls, context):
        return context.mode in {'OBJECT','EDIT_MESH'} and exists_ug_state()

    def execute(self, context):
        update_ugboundaries()
        return {'FINISHED'}


def update_ugboundaries():
    '''Updates Boundary definitions and properties according to Materials
    for all boundary faces.
    '''

    global ugboundaries

    ob = get_ug_object()
    # Initialize ugboundaries
    for b in ugboundaries:
        b.deleted = True
        b.ugfaces = []
        b.nFaces = 0
        b.startFace = -1

    # Process all boundary faces
    for i in range(len(ugfaces)):
        ugf = ugfaces[i]
        if ugf.deleted:
            continue
        if ugf.neighbour:
            continue
        if ugf.bi == -1:
            continue

        mati = ob.data.polygons[ugf.bi].material_index

        # Create default slot and material if nothing exists
        if len(ob.material_slots) == 0:
            if "default" in bpy.data.materials:
                mat = bpy.data.materials["default"]
            else:
                mat = bpy.data.materials.new(name="default")
            bpy.ops.object.material_slot_add()
            ob.active_material = mat

        mat = ob.material_slots[mati]

        while (mati >= len(ugboundaries)):
            UGBoundary("new")

        patch = ugboundaries[mati]
        patch.patchname = mat.name

        patch.deleted = False
        patch.ugfaces.append(ugf)

    # Update number of faces in boundaries. Note: startFace is updated at export.
    for b in ugboundaries:
        b.nFaces = len(b.ugfaces)
        l.debug("Faces on patch %s: %d" %(b.patchname, b.nFaces))

    return None


class UG_OT_UpdateZonesFromVertexGroups(bpy.types.Operator):
    '''Update Zone Data from Vertex Groups After Zone Changes'''
    bl_idname = "unstructured_grids.update_ugzones"
    bl_label = "Update UG Zones From Vertex Groups"

    @classmethod
    def poll(cls, context):
        return context.mode in {'OBJECT','EDIT_MESH'} and exists_ug_state()

    def execute(self, context):
        n = update_ugzones()
        self.report({'INFO'}, "Updated %d zone(s) from vertex groups" % n)
        return {'FINISHED'}


def update_ugzones():
    '''Update face and cell zones from vertex group assignments'''

    global ugzones

    ob = get_ug_object()
    # Initialize zones
    for z in ugzones:
        z.deleted = True
        z.ugfaces = []
        z.ugcells = []

    # Process all vertex groups:
    i = 0
    n = 0
    for vg in ob.vertex_groups:
        vgname = vg.name

        # Determine zonetype
        if vgname.startswith('cell'):
            zonetype = 'cell'
        elif vgname.startswith('face'):
            zonetype = 'face'
        else:
            continue

        zonename = vgname.replace(zonetype + 'Zone_', '')
        l.debug("Updating %s zone %r" % (zonetype, zonename))

        # Get existing zone or create new
        if len(ugzones) <= i:
            z = UGZone(zonetype, zonename)
        else:
            z = ugzones[i]
        z.deleted = False

        # Get vertices belonging to this vertex group
        verts = [v.index for v in ob.data.vertices if \
                 vg.index in [vg.group for vg in v.groups]]

        if zonetype == 'cell':
            l.debug("Search cells among %d verts" % len(verts))
            z.ugcells = ug_op.get_ugcells_from_vertices_exclusive(verts)
            l.debug("Cell zone: " + z.zonename + " cell count: %d" % len(z.ugcells))
            n += 1
        elif zonetype == 'face':
            l.debug("Search cells among %d verts" % len(verts))
            z.ugfaces = ug_op.get_ugfaces_from_vertices_exclusive(verts)
            l.debug("Face zone: " + z.zonename + " face count: %d" % len(z.ugfaces))
            n += 1
        i += 1

    return n


class UG_OT_UpdateUGAllFromBlender(bpy.types.Operator):
    '''Update All Changes Made in Blender into UG Data'''
    bl_idname = "unstructured_grids.update_all_from_blender"
    bl_label = "Update UG Data and Storage From Blender"

    @classmethod
    def poll(cls, context):
        return context.mode in {'OBJECT','EDIT_MESH'} and exists_ug_state()

    def execute(self, context):
        update_ug_all_from_blender(self)
        self.report({'INFO'}, "Updated UG Data and Storage")
        return {'FINISHED'}


def update_ug_all_from_blender(self=None):
    '''Update all changes made in Blender into UG data model.
    Returns True if Unstructured Grid is OK, False otherwise.
    '''

    global ugverts, ugfaces, ugcells, ugboundaries, ugzones

    ug_props = bpy.context.scene.ug_props
    found_cells = False
    for c in ugcells:
        if c.deleted == False:
            found_cells = True

    # Zero all if there are no cells
    if not found_cells:
        ugcells = []
        ugfaces = []
        ugverts = []
        ugboundaries = []
        ugzones = []
        ug_props.text_boundary = ''
        ug_props.text_faces = ''
        ug_props.text_neighbour = ''
        ug_props.text_owner = ''
        ug_props.text_points = ''
        ug_props.text_cellZones = ''
        ug_props.text_faceZones = ''
        if self:
            self.report({'ERROR'}, "No cells detected in %r" % obname)
        return False

    # Force refresh if in Edit Mode
    ob = get_ug_object()
    if ob.mode == 'EDIT':
        bpy.ops.object.mode_set(mode = 'OBJECT')
        bpy.ops.object.mode_set(mode = 'EDIT')

    # Force update boundaries
    update_ugboundaries()
    # Force update zones
    update_ugzones()
    # Update PolyMesh text strings
    io_polymesh.ugdata_to_polymesh()

    return True


######################
##### SCRAP YARD #####
######################


def get_next_undeleted_cell(clist, ind):
    '''Finds next undeleted UGCell in clist starting from index ind+1.
    Returns UGCell and it's index in clist.
    '''
    # TODO: Delete function if not needed

    for i in range(ind + 1, len(clist)):
        c = clist[i]
        if not c.deleted:
            return c, i
    return None, None


def order_ugcells_by_BFS():
    '''Breadth-first search (BFS) algorithm to sort cells. See
    https://en.wikipedia.org/wiki/Breadth-first_search.
    Return list of ugcells ordered so that each cell is neighbour to
    at least one of the previous cells. Ordered cell list is used
    for exporting UG data to ensure that cell editing result is
    a single closed region and to order faces and cells.

    Note: This function uses UGFace.ei to track processed faces.
    '''

    global ugverts, ugfaces, ugcells, ugboundaries

    from collections import deque

    # Initialize
    print_interval = 1000 # Debug print progress interval
    i = 0 # Number of processed cells
    d = deque() # deque of cells to be processed
    clist = [] # ordered list of cells, to be generated here

    # Initialize face ei to -1 (0 will mark processed face)
    for f in ugfaces:
        f.ei = -1

    # Start from first undeleted cell
    for c in ugcells:
        if not c.deleted:
            d.append(c)
            break

    # Process cells in deque
    while (len(d) > 0):
        # Get last unprocessed cell and add to list
        c = d.popleft()
        clist.append(c)
        if fulldebug: l.debug("Processing cell %d", c.ii)

        # Go through all cell faces:
        for f in c.ugfaces:
            # Skip boundary faces
            if not f.neighbour:
                continue
            # Skip processed faces
            if f.ei == 0:
                continue
            # Make c owner of face
            if f.neighbour == c:
                f.invert_face_dir()
            # Add cell arcross unprocessed face to deque unless already there
            if not f.neighbour in d:
                d.append(f.neighbour)
            # Mark face as processed
            f.ei = 0

        if i % print_interval == 0:
            l.debug("Ordered cells: %d" % i)
        i += 1

    l.debug("Total ordered cells: %d" % (i-1))
    return clist


def order_ugcells_by_internal_face_search():
    '''An ad-hoc internal face search algorithm to sort cells.
    Return list of ugcells ordered so that each cell is neighbour to
    at least one of the previous cells. Ordered cell list is used
    for exporting UG data to ensure that cell editing result is
    a single closed region and to order faces and cells.

    Note: This function uses UGCell.ei to store number of unprocessed
    neighbour cells, to track which cells to include next in the list.
    '''

    global ugverts, ugfaces, ugcells, ugboundaries

    def init_ncells():
        '''Initialize ugcells ei property with number of neighbour cells
        and ugfaces ei property with zero.
        '''

        n = 0 # Number of cells
        # Initialize cells
        for c in ugcells:
            if c.deleted:
                c.ei = 0
                continue
            n += 1
            # Initialize with all faces, then decrease number of boundary faces
            c.ei = len(c.ugfaces)
            for f in c.ugfaces:
                if not f.neighbour:
                    c.ei -= 1

        return n

    def get_next_neighbour_cell(c, clist):
        '''Return first neighbour cell next to argument cell c
        which is not in argument cell list.
        '''

        if fulldebug: l.debug("Searching through faces of " + str(c.ii) + " whose c.ei=%d" % c.ei)
        for f in c.ugfaces:
            if fulldebug: l.debug("- checking face " + str(f))
            if not f.neighbour:
                if fulldebug: l.debug("- boundary face")
                continue
            if f.owner == c:
                nc = f.neighbour
            else:
                nc = f.owner
            # Find new path to an unprocessed cell
            if nc in clist:
                if fulldebug: l.debug("- Already in list " + str(nc.ii))
                continue
            if nc.ei > 0:
                # If neighbour is unprocessed, decrease unprocessed cell counts and return
                if fulldebug: l.debug("- Found new cell " + str(nc.ii))
                return nc
            else:
                if fulldebug: l.debug("- nc.ei=%" % nc.ei)

        # No neighbour found
        return None

    def get_next_unfinished_cell(clist):
        '''Return first cell in argument clist which contains neighbours
        not in clist according to ei variable.
        '''

        for c in clist:
            if c.ei > 0:
                return c
        return None

    def decrease_cell_counters(c, clist):
        '''Decreases cell ei index for every cell in clist connected to cell c'''

        for f in c.ugfaces:
            if not f.neighbour:
                continue
            if f.owner == c:
                nc = f.neighbour
            else:
                nc = f.owner
            if nc in clist and nc.ei > 0:
                c.ei -= 1
                nc.ei -= 1
                if fulldebug: l.debug("-- Reduced neighbor count for " + str(nc.ii))


    # Intialize
    print_interval = 1000 # Debug print progress interval
    i = 0
    n = init_ncells()
    clist = [] # list of ordered ugcells, to be populated
    for c in ugcells:
        if c.ei > 0:
            clist.append(c)
            break

    # Process all cells
    for i in range(n-1):
        if fulldebug: l.debug("i=%d, clist len=%d" % (i, len(clist)))
        # Extend clist with a neighbour cell
        c = get_next_unfinished_cell(clist)
        if fulldebug: l.debug("Next unfinished cell is " + str(c.ii) + " after ei=%d" % ((c.ei)-1))
        if not c:
            l.debug("No cell at i=%d" % i)
        nc = get_next_neighbour_cell(c, clist)
        if fulldebug: l.debug("Next neighbour cell is " + str(nc.ii) + " after ei=%d" % nc.ei)
        if not nc:
            l.warning("Region is incomplete! Final cell count %d/%d" % (i, n))
            return clist
        clist.append(nc)
        decrease_cell_counters(nc, clist)
        if i % print_interval == 0:
            l.debug("Ordering cells, progress: %d/%d" % (i, n))

    l.debug("Sorting completed, cell count %d/%d" % (len(clist), n))
    if fulldebug:
        for i in range(len(clist)):
            l.debug(str(i) + " " + str(clist[i]))

    return clist
