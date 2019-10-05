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

# Note: Unstructured Grids assumes that PolyMesh face definitions are
# ordered in such a way that they first define internal faces, and
# then boundary faces (each boundary patch at a time). This is normal
# way for OpenFOAM meshes, since it allows easy definition of boundary
# patches by startFace and nFaces. So, even though PolyMesh
# definition allows definition of cell label "-1" for faces, it is
# never encountered in a well-ordered PolyMesh.

# Note: OpenFOAM includes renumberMesh utility, which orders cells and
# optimizes bandwidth. It is recommended to run renumberMesh to
# optimize PolyMesh.

# Note: Face normals for all boundary faces must always point from
# domain towards outside. Therefore neighbour list is always shorter
# than owner list.

# Note: Neighbour face list can contain a cell with larger cell index
# than largest cell index in owner face list. This can happen e.g. for
# meshes created with SnappyHexMesh when a cell is located inside the
# mesh refinement region. Cell traversal path can form an "inward
# spiral" structure, in which case the cell in spiral center is
# defined completely by neighbour cells (all face normals of that cell
# point inwards. If domain cell traversal ends with such a spiral,
# then neighbour list ends having a larger cell index than owner list.


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


##### IMPORT #####

class UG_OT_ImportPolyMesh(bpy.types.Operator, ImportHelper):
    '''Import OpenFOAM PolyMesh Files into Blender as UG Data'''
    bl_idname = "unstructured_grids.import_openfoam_polymesh"
    bl_label = "Import OpenFOAM PolyMesh (UG)"

    @classmethod
    def poll(cls, context):
        return context.mode in {'OBJECT','EDIT_MESH'}

    def execute(self, context):
        read_polymesh_files(self)
        return {'FINISHED'}


def read_polymesh_files(self):
    '''Reads PolyMesh files' contents from a folder into strings'''

    import os
    ug_props = bpy.context.scene.ug_props
    dirpath = os.path.dirname(self.filepath)

    # Compulsory files
    filenames = ['boundary', 'faces', 'neighbour', 'owner', 'points']
    for f in filenames:
        varname = "text_" + f
        filepath = os.path.join(dirpath, f)
        l.debug("Reading in as string: %s" % filepath)

        if not (os.path.isfile(filepath)):
            self.report({'ERROR'}, "Could not find %r" \
                        % filepath)
            return None

        with open(filepath, 'r') as infile:
            setattr(ug_props, varname, infile.read())

    # Optional files
    filenames = ['cellZones', 'faceZones']
    for f in filenames:
        varname = "text_" + f
        filepath = os.path.join(dirpath, f)
        l.debug("Reading in as string: %s" % filepath)

        if not (os.path.isfile(filepath)):
            l.debug("Could not find %r" % filepath)
            continue

        with open(filepath, 'r') as infile:
            setattr(ug_props, varname, infile.read())

    polymesh_boundary_ingroup_fix()
    polymesh_to_ugdata(self)
    return None


def polymesh_boundary_ingroup_fix():
    '''Reformats ingroup entries spanning several lines into one line,
    because otherwise multiline entry breaks regex matching logic
    in polymesh_get_boundary
    '''

    import re
    ug_props = bpy.context.scene.ug_props
    text = ug_props.text_boundary

    inside = False # boolean for marking boundary entries in text
    result = ''
    rec1 = re.compile(r'^\s*inGroups\s*$', re.M)
    rec2 = re.compile(r'\;', re.M)

    for line in text.splitlines():
        regex = rec1.search(line)
        if regex:
            # Initialize
            inside = True
            res = "        inGroups        "
        elif inside:
            res += line + ' '
            regex = rec2.search(line)
            if regex:
                # Reached end of inGroup
                inside = False
                result += res + '\n'
        else:
            result += line + '\n'

    ug_props.text_boundary = result
    return None


class UG_OT_PolyMeshToUG(bpy.types.Operator):
    '''(Re)generate UG Data from PolyMesh Data in UG Storage'''
    bl_idname = "unstructured_grids.polymesh_to_ug"
    bl_label = "(Re)generate UG Data from UG Storage"

    @classmethod
    def poll(cls, context):
        ug_props = bpy.context.scene.ug_props
        return context.mode in {'OBJECT','EDIT_MESH'} and len(ug_props.text_points) > 0

    def execute(self, context):
        polymesh_to_ugdata(self)
        return {'FINISHED'}


def polymesh_to_ugdata(self):
    '''Convert OpenFOAM polyMesh data from text files
    into UG data structures and Blender mesh
    '''

    ug_props = bpy.context.scene.ug_props
    if len(ug_props.text_points) == 0:
        return None

    ug.hide_other_objects()
    ob = ug.initialize_ug_object()
    verts = polymesh_get_verts(ug_props.text_points)
    [edges, faces] = polymesh_get_faces( \
        ug_props.text_owner, ug_props.text_neighbour, ug_props.text_faces)
    polymesh_get_boundary(ug_props.text_boundary)

    # Create vertices and faces into mesh object
    ob.data.from_pydata(verts, edges, faces)
    ob.data.validate()

    # Boundary patches
    apply_materials_to_boundaries(ob)
    ug.update_ugboundaries()

    # Cell and face zones
    polymesh_get_zone('cell', ug_props.text_cellZones)
    polymesh_get_zone('face', ug_props.text_faceZones)
    apply_vertex_groups_to_zones(ob)
    ug.update_ugzones() # not necessary after import, but yes after regenerate
    bpy.ops.object.mode_set(mode='OBJECT')


def polymesh_get_verts(text):
    '''Creates list of vertex triplets from PolyMesh points text string'''

    import re
    verts = [] # list of x, y, z point coordinate triplets

    rec = re.compile(r'^\(([dDeE\d\.\-]+)\s+([dDeE\d\.\-]+)\s+([dDeE\d\.\-]+)\)', re.M)
    i = 0
    for line in text.splitlines():
        regex = rec.search(line)
        if regex:
            x = float(regex.group(1))
            y = float(regex.group(2))
            z = float(regex.group(3))
            verts.append(tuple([x, y, z]))
            ug.UGVertex(i)
            if i % print_interval == 0:
                l.debug("... processed vertex count: %d" % i)
            i += 1

    l.debug("Number of vertices read: %d" % len(verts))
    return verts


def polymesh_get_faces(text_owner, text_neighbour, text_faces):
    '''Creates edge and face list from PolyMesh owner, neighbour and
    faces text blocks. Initializes UG faces and UG cells.
    '''

    # TODO: Separate initialization? Function is getting big.

    edges = [] # List of edge vertex index pairs, to be generated
    faces = [] # List of face vertex index lists, to be generated
    gen_edges = bpy.context.scene.ug_props.generate_internal_edges

    # Read in owner and neighbour lists
    owner = polymesh_get_intlist(text_owner)
    l.debug("Number of faces that have owner: %d" % len(owner))
    neighbour = polymesh_get_intlist(text_neighbour)
    l.debug("Number of faces that have neighbour: %d" % len(neighbour))
    l.debug("Number of owner minus neighbour: %d" % (len(owner) - len(neighbour)))
    face_verts = polymesh_get_list_intlist(text_faces)
    l.debug("Number of face definitions: %d" % len(face_verts))


    # Populate list of ugcells
    for i in range(max(max(owner), max(neighbour)) + 1):
        # Add new UGCell
        ug.UGCell()
        if i % print_interval == 0:
            l.debug("... processed cell count: %d" % i)
    l.debug("Final cell count: %d" % i)

    # Create faces at boundary and only edges for internal faces
    for i in range(len(face_verts)):
        # Add ugface and facemap entry
        f = ug.UGFace(face_verts[i])

        # Part 1. Process owner
        # Add owner cell index
        f.owner = ug.ugcells[owner[i]]
        # Add face to owner's ugfaces list
        f.owner.ugfaces.append(f)

        for v in face_verts[i]:
            # Add owner cell info for vertices
            clist = ug.ugverts[v].ugcells
            if not f.owner in clist:
                clist.append(f.owner)
            # Add vertex to cell's ugverts
            if ug.ugverts[v] not in f.owner.ugverts:
                f.owner.ugverts.append(ug.ugverts[v])

        # Part 2. Process neighbour
        if i < len(neighbour):
            # Add neighbour cell index
            f.neighbour = ug.ugcells[neighbour[i]]
            # Add face to neighbour's ugfaces list
            f.neighbour.ugfaces.append(f)

            for v in face_verts[i]:
                # Add neighbour cell info for vertices
                clist = ug.ugverts[v].ugcells
                if not f.neighbour in clist:
                    clist.append(f.neighbour)
                # Add vertex to cell's ugverts
                if ug.ugverts[v] not in f.neighbour.ugverts:
                    f.neighbour.ugverts.append(ug.ugverts[v])

            # Add edges for geometry generation if needed
            if gen_edges:
                for j in range(len(face_verts[i])):
                    edges.append(tuple([face_verts[i][j-1], face_verts[i][j]]))

        else:
            # Boundary face, add index and create face geometry for mesh object
            f.bi = len(faces)
            faces.append(tuple(face_verts[i]))
            ug.facemap[f.bi] = f # Add face to facemap

        if i % print_interval == 0:
            l.debug("... processed face count: %d" % i)
    l.debug("Final face count: %d" % i)

    l.debug("Number of edge definitions: %d" % len(edges))
    l.debug("Number of face definitions: %d" % len(faces))

    return edges, faces


def polymesh_get_intlist(text):
    '''Creates integer list from argument PolyMesh integer text block'''

    import re
    iList = [] # list of integers to be generated
    inside = False # boolean for marking integer list in text

    rec1 = re.compile(r'^\(', re.M)
    rec2 = re.compile(r'^\)', re.M)
    rec3 = re.compile(r'^(\d+)', re.M)

    for line in text.splitlines():
        # Opening of integer list by single parenthesis
        regex = rec1.search(line)
        if regex:
            inside = True

        # Closing of integer list by single parenthesis
        regex2 = rec2.search(line)
        if regex2:
            inside = False

        # Integer, at start of line
        regex3 = rec3.search(line)
        if inside and regex3:
            iList.append(int(regex3.group(1)))

    return iList


def polymesh_get_list_intlist(text):
    '''Creates list of integer lists from argument PolyMesh
    text block
    '''

    # TODO: Get rid of code duplication

    import re
    iList = [] # list of integers lists to be generated
    inside = False # boolean for marking integer list in text

    rec1 = re.compile(r'^\(', re.M)
    rec2 = re.compile(r'^\)', re.M)
    rec3 = re.compile(r'^\d+\(([\d\s]+)\)', re.M)

    for line in text.splitlines():
        # Opening of integer list by single parenthesis
        regex = rec1.search(line)
        if regex:
            inside = True

        # Closing of integer list by single parenthesis
        regex2 = rec2.search(line)
        if regex2:
            inside = False

        # List of integer list within parenthesis
        regex3 = rec3.search(line)
        if inside and regex3:
            vals = regex3.group(1).split()
            valList = []
            for val in vals:
                valList.append(int(val))
            iList.append(valList)

    return iList


def polymesh_get_boundary(text):
    '''Creates boundary objects from PolyMesh boundary text string'''

    import re
    inside = False # boolean for marking boundary entries in text

    rec1 = re.compile(r'^\(', re.M)
    rec2 = re.compile(r'^\)', re.M)
    rec3 = re.compile(r'^\s+([\w\%\:\-\.]+)$', re.M)
    rec4 = re.compile(r'^\s+type\s+(\w+)\;$', re.M)
    rec5 = re.compile(r'^\s+inGroups\s+([\w\s\(\)\<\>]+)\;\s*$', re.M)
    rec6 = re.compile(r'^\s+nFaces\s+(\d+)\;$', re.M)
    rec7 = re.compile(r'^\s+startFace\s+(\d+)\;$', re.M)

    for line in text.splitlines():
        # Opening of integer list by single parenthesis
        regex = rec1.search(line)
        if regex:
            inside = True

        # Closing of integer list by single parenthesis
        regex = rec2.search(line)
        if regex:
            inside = False
            if patch.nFaces == 0 or patch.startFace == -1:
                l.error("Boundary definition " + str(patch.patchname) \
                        + " is broken")
                return None
            # Add ugfaces to boundary
            for i in range(patch.startFace, patch.startFace + patch.nFaces):
                patch.ugfaces.append(ug.ugfaces[i])

        if not inside:
            continue

        # New entry is a word (with possibly special characters) on its own line
        regex = rec3.search(line)
        if regex:
            patchname = str(regex.group(1))
            l.debug("Reading in boundary patch %d: %s" % (len(ug.ugboundaries), patchname))
            patch = ug.UGBoundary(patchname)
            continue

        # type
        regex = rec4.search(line)
        if regex:
            patch.typename = str(regex.group(1))
            continue

        # inGroups
        regex = rec5.search(line)
        if regex:
            patch.inGroups = str(regex.group(1))
            continue

        # nFaces
        regex = rec6.search(line)
        if regex:
            patch.nFaces = int(regex.group(1))
            continue

        # startFace
        regex = rec7.search(line)
        if regex:
            patch.startFace = int(regex.group(1))
            continue

    return None


def polymesh_get_zone(zonetype, text):
    '''Creates zone objects from PolyMesh zone text string for type
    (either face or cell)
    '''

    import re
    inside1 = False # boolean for marking top level entries in text
    inside2 = False # boolean for marking bottom level entries in text
    label = '' # label definition preceding list of integers
    iList = [] # list of integers inside bottom level entry
    zone = None # UGZone

    rec1 = re.compile(r'^\s*\(', re.M)
    rec2 = re.compile(r'^\s*\)', re.M)
    rec3 = re.compile(r'^\s*([a-zA-Z][\w\%\:\-\.]*)$', re.M)
    rec4 = re.compile(r'^\s*(\w+\s+[\w\<\>]+)\s*$', re.M)
    rec5 = re.compile(r'^(\d+)', re.M)

    for line in text.splitlines():
        # Opening parenthesis
        regex = rec1.search(line)
        if regex:
            if inside1:
                inside2 = True
            else:
                inside1 = True

        # Closing parenthesis
        regex = rec2.search(line)
        if regex:
            if inside2:
                # Add indices to correct UGZone list
                if label.startswith('cellLabels'):
                    for i in iList:
                        zone.ugcells.append(ug.ugcells[i])
                    l.debug("CellZone: " + str(zonename) + ": %d cells" % len(iList))
                elif label.startswith('faceLabels'):
                    for i in iList:
                        zone.ugfaces.append(ug.ugfaces[i])
                    l.debug("FaceZone: " + str(zonename) + ": %d faces" % len(iList))
                elif label.startswith('flipMap'):
                    # face flipMap is stored as a dummy list for now
                    zone.flipMap = iList
                    l.debug("Face FlipMap: " + str(zonename) + ": %d faces" % len(iList))
                else:
                    l.error("Unknown label " + str(label))
                    return None

                # Reinitialize
                iList = []
                inside2 = False
                label = ''
            else:
                inside1 = False

        if not inside1:
            continue

        if not inside2:
            # New zone is a single word (with possibly special characters) on its own line
            regex = rec3.search(line)
            if regex:
                zonename = str(regex.group(1))
                l.debug("Reading in zone [%d]" % len(ug.ugzones))
                zone = ug.UGZone(zonetype, zonename)

            # Get label definition line
            else:
                regex = rec4.search(line)
                if regex:
                    label = str(regex.group(1))

            continue

        # Integer, at start of line
        regex = rec5.search(line)
        if regex:
            iList.append(int(regex.group(1)))

    return None


def apply_materials_to_boundaries(ob):
    '''Sets materials to faces in object ob according to boundary assignments'''

    mati = 0 # Material index
    facecount = 0

    # Clear object's material slots - not needed but this should do it:
    # for i in range(len(ob.material_slots)):
    #    ob.active_material_index = i
    #    bpy.ops.object.material_slot_select()
    #    bpy.ops.object.material_slot_remove()

    # Process each boundary
    for b in ug.ugboundaries:
        # Create new material if needed
        l.debug("Material for %s: %d" % (b.patchname, mati))
        if b.patchname in bpy.data.materials:
            mat = bpy.data.materials[b.patchname]
        else:
            mat = bpy.data.materials.new(name=b.patchname)

        # Create new material slot to object and set material
        bpy.ops.object.material_slot_add()
        ob.active_material = mat

        # Color is given to all patches but not the one named
        # 'default', because that is used for new boundary faces. User
        # probably wants to assign manually those, so we leave them
        # uncolored to help spot them.
        if b.patchname != 'default':
            ob.data.materials[mati].diffuse_color = get_face_color(mati)

        # Set material index for mesh faces
        for i in range(b.nFaces):
            ob.data.polygons[facecount].material_index = mati
            facecount += 1
        mati += 1


def get_face_color(mati):
    '''Gives a color to argument material number'''

    base_colors = [(0,0,1,1), (1,0,0,1), (0,1,0,1), \
             (0.7,0.7,0,1), (0,0.7,0.7,1), (0.7,0,0.7,1)]
    if mati < len(base_colors):
        return base_colors[mati]

    # Get random colors after base colors
    import random
    random.seed(10043 + mati)
    [r, g, b] = [random.random() for i in range(3)]
    return [r, g, b, 1.0]


def apply_vertex_groups_to_zones(ob):
    '''Set/create vertex groups for object ob according to cell and face
    zones data in ugzones
    '''

    def assign_to_new_vertex_group(ob, z):
        '''Assigns currently selected vertices to new vertex group in object
        ob for zone z
        '''

        bpy.ops.object.mode_set(mode="EDIT")
        bpy.ops.object.vertex_group_assign_new()
        # Set vertex group name
        vgname = z.zonetype + "Zone_" + z.zonename
        vg = ob.vertex_groups[-1]
        vg.name = vgname
        bpy.ops.object.mode_set(mode="OBJECT")


    mode = ob.mode # Save original mode

    for z in ug.ugzones:
        # Deselect all vertices, edges and faces
        bpy.ops.object.mode_set(mode="EDIT")
        bpy.ops.mesh.select_all(action='DESELECT')
        bpy.ops.object.mode_set(mode="OBJECT")

        if z.zonetype == 'cell':
            n = ug_op.select_vertices_from_ugcells(ob, z.ugcells)
            if n > 0:
                assign_to_new_vertex_group(ob, z)

        elif z.zonetype == 'face':
            n = ug_op.select_vertices_from_ugfaces(ob, z.ugfaces)
            if n > 0:
                assign_to_new_vertex_group(ob, z)

    bpy.ops.object.mode_set(mode=mode)


##################
##### EXPORT #####
##################


class UG_OT_ExportPolyMesh(bpy.types.Operator, ExportHelper):
    '''Export UG Data as OpenFOAM PolyMesh Files'''
    bl_idname = "unstructured_grids.export_openfoam_polymesh"
    bl_label = "Export OpenFOAM PolyMesh (UG)"

    filename_ext = ".polyMesh" # Dummy, required by ExportHelper

    @classmethod
    def poll(cls, context):
        return context.mode in {'OBJECT','EDIT_MESH'} and ug.exists_ug_state()

    def execute(self, context):
        ug.update_ug_all_from_blender(self)
        write_polymesh_files(self)
        return {'FINISHED'}


def ugdata_to_polymesh(self):
    '''Convert up-to-date UG data into PolyMesh text data strings
    and store those in ug_props.text_* string properties
    '''

    ob = bpy.data.objects[ug.obname]
    update_text_points(ob)
    owneri, neighbouri = update_ei_and_text_faces(ob)
    update_text_owner_neighbour(owneri, neighbouri)
    update_text_boundary()
    update_text_cell_zones()
    # update_text_face_zones() # TODO
    return None

def update_text_points(ob):
    '''Updates PolyMesh points string contents from Blender object vertices'''

    text_verts = '' # Text for vertex coordinates
    n = 0
    for ugv in ug.ugverts:
        if ugv.deleted:
            ugv.ei = -1
            continue
        # Update export index
        ugv.ei = n
        v = ob.data.vertices[ugv.bi]
        text_verts += "(" + "%.6g" % v.co.x + " " \
                      + "%.6g" % v.co.y + " " \
                      + "%.6g" % v.co.z + ")\n"
        n += 1

    # Generate new text
    text = of_file_header('vectorField', 'points') + "\n"
    text += str(n) + "\n(\n"
    text += text_verts + ")\n"
    bpy.context.scene.ug_props.text_points = text
    l.debug("text_points updated points: %d" % len(ob.data.vertices))
    return None


def update_ei_and_text_faces(ob):
    '''Updates export indices (ei) of ugcells, ugfaces and ugboundaries
    and generates PolyMesh faces text string contents for object ob.
    Exported face index is updated according to ugboundary ugfaces.
    Update is done in two phases:

    1. internal face pass:

    Generates export cell indices according to PolyMesh requirement that
    face normal points from lower cell index to higher cell index.
    Face normal is determined by right hand rule from vertex list.

    2. boundary face pass:

    Boundar faces are numbered according to boundary assignments.

    New owner and neighbour index lists are generated alongside face
    indexing and returned.
    '''

    def gen_line(ugverts):
        '''Construct face definition text line from vertex indices list'''

        line = str(len(ugverts)) + "("
        for j in range(len(ugverts) - 1):
            line += str(ugverts[j].ei) + " "
        line += str(ugverts[-1].ei) + ")\n"
        return line

    def reset_ei():
        '''Resets export indices'''

        for c in ug.ugcells:
            c.ei = -1
        for f in ug.ugfaces:
            f.ei = -1

    def internal_face_pass(clist):
        '''Generate face definition text for internal faces from
        argument cell list clist.
        '''

        text = ''
        fei = 0 # face export index
        cei = 0 # cell export index
        export_ugfaces = [] # list of ugfaces in export order
        owneri = [] # owner cell index list
        neighbouri = [] # neighbour cell index list

        # Go through each internal face of cells and number cells and
        # internal faces
        for c in clist:
            if c.deleted:
                continue
            c.ei = cei # Set cell index
            for f in c.ugfaces:
                if f.deleted:
                    continue
                if not f.neighbour:
                    continue
                if f.owner == c and f.neighbour.ei != -1:
                    # This leads to "Faces not in upper triangular order" errors
                    l.warning("Face normal should be mirrored for face ei %d" % f.ei)
                if f.ei != -1:
                    continue
                export_ugfaces.append(f)
                f.ei = fei # Set face index
                text += gen_line(f.ugverts) # Add definition line and proceed
                fei += 1
            cei += 1

        # Update owner and neighbour index lists
        for f in export_ugfaces:
            owneri.append(f.owner.ei)
            neighbouri.append(f.neighbour.ei)

        return text, fei, owneri, neighbouri

    def boundary_face_pass(fei, ob, owneri):
        '''Generate face definition text for boundary faces.
        Return text and number of internal faces.
        '''

        startind = fei # Face index for start of boundary faces
        text = ''

        for patch in ug.ugboundaries:
            if patch.deleted:
                continue
            # Update boundary index numbers
            patch.startFace = fei
            for f in patch.ugfaces:
                # Sanity check
                if f.ei != -1:
                    l.error("Face export index already exists for face %d" % f.bi)
                    return None
                f.ei = fei # Set face index
                owneri.append(f.owner.ei) # Append owner index
                text += gen_line(f.ugverts)
                fei += 1
        return text, fei, owneri

    # Ordering of cells and initialization
    #ordered_ugcells = ug.order_ugcells_by_internal_face_search() # Very slow
    #ordered_ugcells = ug.order_ugcells_by_BFS() # Somewhat slow
    ordered_ugcells = ug.ugcells # No ordering of cells, no guarantee of single region
    reset_ei()

    # Internal pass
    text_internal, i, owneri, neighbouri = internal_face_pass(ordered_ugcells)
    l.debug("text_faces updated internal faces: %d", i)

    # Boundary pass
    text_boundary, i, owneri = boundary_face_pass(i, ob, owneri)
    l.debug("text_faces updated total number of faces: %d", i)

    # Generate text string
    text = of_file_header('faceList', 'faces') + "\n"
    text += str(i) + "\n(\n"
    text += text_internal + text_boundary + ")\n"

    bpy.context.scene.ug_props.text_faces = text
    return owneri, neighbouri


def update_text_owner_neighbour(owneri, neighbouri):
    '''Updates PolyMesh owner and neighbour text string contents from
    argument integer lists.
    '''

    nall = len(owneri)
    nneighbour = len(neighbouri)

    # Generate text string
    text_owner = of_file_header('labelList', 'owner') + "\n"
    text_owner += str(nall) + "\n(\n"
    text_neighbour = of_file_header('labelList', 'neighbour') + "\n"
    text_neighbour += str(nneighbour) + "\n(\n"

    for i in owneri:
        text_owner += str(i) + "\n"
    for i in neighbouri:
        text_neighbour += str(i) + "\n"

    text_owner += ")\n"
    text_neighbour += ")\n"

    bpy.context.scene.ug_props.text_owner = text_owner
    bpy.context.scene.ug_props.text_neighbour = text_neighbour
    l.debug("text_owner updated faces: %d" % nall)
    l.debug("text_neighbour updated faces: %d" % nneighbour)
    return None


def update_text_boundary():
    '''Updates PolyMesh boundary text string contents from UG data'''

    btext = '' # generated boundary entries
    nboundaries = 0 # number of boundaries

    for patch in ug.ugboundaries:
        if patch.deleted:
            continue

        text = "    " + patch.patchname + "\n    {\n"
        text += "        type            " + patch.typename + ";\n"
        if patch.inGroups != '':
            text += "        inGroups        " + patch.inGroups + ";\n"
        text += "        nFaces          " + str(patch.nFaces) + ";\n"
        text += "        startFace       " + str(patch.startFace) + ";\n"
        text += "    }\n"
        btext += text
        nboundaries += 1

    # Generate new text
    text = of_file_header('polyBoundaryMesh', 'boundary') + "\n"
    text += str(nboundaries) + "\n(\n"
    text += btext + ")\n"

    bpy.context.scene.ug_props.text_boundary = text
    l.debug("text_boundary updated patches: %d" % nboundaries)
    return None


def update_text_cell_zones():
    '''Update PolyMesh cellZones text string contents from UG data'''

    ztext = '' # generated zone entries
    nzones = 0 # number of cell zones

    # Zone texts
    for z in ug.ugzones:
        if z.zonetype != 'cell':
            continue
        if z.deleted:
            continue
        text = z.zonename + "\n"
        text += "{\n"
        text += "type " + z.zonetype + "Zone;\n"
        text += "cellLabels List<label>\n"
        text += str(len(z.ugcells)) + "\n"
        text += "(\n"
        for c in z.ugcells:
            text += str(c.ei) + "\n"
        text += ");\n}\n\n"
        ztext += text
        nzones +=1

    # Generate new text
    if nzones == 0:
        bpy.context.scene.ug_props.text_cellZones = ''
        l.debug("text_cellZones is empty (no cell zones defined)")
        return None

    text = of_file_header('regIOobject', 'cellZones') + "\n"
    text += str(nzones) + "\n(\n"
    text += ztext + ")\n"

    bpy.context.scene.ug_props.text_cellZones = text
    l.debug("updated cell zones: %d" % nzones)
    return None


def write_polymesh_files(self):
    '''Write contents of data strings into PolyMesh files'''

    import os
    ug_props = bpy.context.scene.ug_props
    dirpath = os.path.dirname(self.filepath)

    # Compulsory files
    filenames = ['boundary', 'faces', 'neighbour', 'owner', 'points']
    for f in filenames:
        varname = "text_" + f
        filepath = os.path.join(dirpath, f)
        l.debug("Writing to: %s" % filepath)

        with open(filepath, 'w') as outfile:
            outfile.write(getattr(ug_props, varname))

    # Optional files
    filenames = ['cellZones', 'faceZones']
    for f in filenames:
        varname = "text_" + f
        filepath = os.path.join(dirpath, f)
        if len(getattr(ug_props, varname)) == 0:
            continue

        l.debug("Writing to: %s" % filepath)
        with open(filepath, 'w') as outfile:
            outfile.write(getattr(ug_props, varname))

    self.report({'INFO'}, "Exported PolyMesh to %r" % dirpath)
    return None


def of_file_header(class_name, object_name):
    '''Returns OpenFOAM dictionary file header using argument
    name for class and object type names
    '''

    h = "FoamFile\n{\n"
    h += "    version     2.0;\n"
    h += "    format      ascii;\n"
    h += "    class       " + class_name + ";\n"
    h += "    object      " + object_name + ";\n}\n"
    return h
