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
from bpy_extras.io_utils import (
    ImportHelper,
    ExportHelper,
    orientation_helper,
    axis_conversion,
)
from . import ug
import logging
l = logging.getLogger(__name__)

##### IMPORT #####

class UG_OT_ImportPolyMesh(bpy.types.Operator, ImportHelper):
    '''Import OpenFOAM PolyMesh as Unstructured Grid'''
    bl_idname = "unstructured_grids.import_openfoam_polymesh"
    bl_label = "Import OpenFOAM PolyMesh"

    def execute(self, context):
        read_polymesh_files(self)
        return {'FINISHED'}


def read_polymesh_files(self):
    '''Reads PolyMesh files' contents from a folder into strings'''

    import os
    ug_props = bpy.context.scene.ug_props
    dirpath = os.path.dirname(self.filepath)

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
    for line in text.splitlines():
        regex = re.search(r'^\s*inGroups\s*$', line, re.M)
        if regex:
            # Initialize
            inside = True
            res = "        inGroups        "
        elif inside:
            res += line + ' '
            regex = re.search(r'\;', line, re.M)
            if regex:
                # Reached end of inGroup
                inside = False
                result += res + '\n'
        else:
            result += line + '\n'

    ug_props.text_boundary = result
    return None


class UG_OT_PolyMeshToUG(bpy.types.Operator):
    '''Generate UG data and mesh object from OpenFOAM PolyMesh file contents'''
    bl_idname = "unstructured_grids.polymesh_to_ug"
    bl_label = "Generate UG from polyMesh texts"

    def execute(self, context):
        polymesh_to_ugdata(self)
        return {'FINISHED'}


def polymesh_to_ugdata(self):
    '''Convert OpenFOAM polyMesh data from text files
    into UG data structures and Blender mesh
    '''

    ob = initialize_ug_object()
    ug_props = bpy.context.scene.ug_props
    verts = polymesh_get_verts(ug_props.text_points)
    [edges, faces] = polymesh_get_faces( \
        ug_props.text_owner, ug_props.text_neighbour, ug_props.text_faces)
    polymesh_get_boundary(ug_props.text_boundary)
    # Create vertices and faces into mesh object
    ob.data.from_pydata(verts, edges, faces)
    ob.data.validate()
    apply_materials_to_boundaries(ob)


def initialize_ug_object():
    '''Creates and returns an initialized and empty UG mesh object and
    initializes UG data
    '''

    # Zero UG data
    ug.ugfaces = []
    ug.ugcells = []
    ug.ugboundaries = []

    # Initialize mesh object
    name = "Unstructured Grid"   
    if name in bpy.data.objects:
        l.debug("Delete existing object " + name)
        bpy.ops.object.select_all(action='DESELECT')
        bpy.data.objects[name].select_set(True)
        mesh = bpy.data.objects[name].data
        bpy.ops.object.delete()
        bpy.data.meshes.remove(mesh)

    l.debug("Create and activate new mesh object " + name)
    mesh_data = bpy.data.meshes.new(name)
    ob = bpy.data.objects.new(name, mesh_data)
    bpy.context.scene.collection.objects.link(ob)
    bpy.context.view_layer.objects.active = bpy.data.objects[name]
    ob.select_set(True)

    return ob


def polymesh_get_verts(text):
    '''Creates list of vertex triplets from PolyMesh points text string'''

    import re
    verts = [] # list of x, y, z point coordinate triplets

    for line in text.splitlines():
        regex = re.search(r'^\(([dDeE\d\.\-]+)\s+([dDeE\d\.\-]+)\s+([dDeE\d\.\-]+)\)', line, re.M)
        if regex:
            x = float(regex.group(1))
            y = float(regex.group(2))
            z = float(regex.group(3))
            verts.append(tuple([x, y, z]))

    l.debug("Number of coordinate triplets read: %d" % len(verts))
    return verts


def polymesh_get_faces(text_owner, text_neighbour, text_faces):
    '''Creates edge and face list from PolyMesh owner, neighbour and
    faces text blocks
    '''

    edges = [] # List of edge vertex index pairs, to be generated
    faces = [] # List of face vertex index lists, to be generated
    gen_edges = bpy.context.scene.ug_props.generate_internal_edges

    # Read in owner and neighbour lists
    owner = polymesh_get_intlist(text_owner)
    neighbour = polymesh_get_intlist(text_neighbour)
    face_verts = polymesh_get_list_intlist(text_faces)

    # Populate list of ugcells
    for i in range(max(owner) + 1):
        # Add new entry to list of ugcells
        ugcell = ug.UGCell(i)
        ug.ugcells.append(ugcell)

    # Create faces at boundary and only edges for internal faces
    for i in range(len(face_verts)):
        # Add to list of ugfaces
        ugface = ug.UGFace(i, face_verts[i])
        ug.ugfaces.append(ugface)
        # Add owner cell index
        ugface.owneri = owner[i]
        # Add face to owner's faces list
        ug.ugcells[owner[i]].faces.append(i)

        # Add geometry to object
        if i < len(neighbour):
            # Add neighbour cell index
            ugface.neighbouri = neighbour[i]
            # Add face to neighbour's faces list
            ug.ugcells[neighbour[i]].faces.append(i)

            # Add edges if needed
            if gen_edges:
                for j in range(len(face_verts[i])):
                    edges.append(tuple([face_verts[i][j-1], face_verts[i][j]]))
        else:
            # Boundary face, add faces
            faces.append(tuple(face_verts[i]))

    l.debug("Number of edge index pairs generated: %d" % len(edges))
    l.debug("Number of boundary face index lists generated: %d" % len(faces))

    return edges, faces


def polymesh_get_intlist(text):
    '''Creates integer list from argument PolyMesh integer text block'''

    import re
    iList = [] # list of integers to be generated
    inside = False # boolean for marking integer list in text

    for line in text.splitlines():
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

    l.debug("Number of integers read: %d" % len(iList))
    return iList


def polymesh_get_list_intlist(text):
    '''Creates list of integer lists from argument PolyMesh
    text block
    '''

    # TODO: Get rid of code duplication

    import re
    iList = [] # list of integers lists to be generated
    inside = False # boolean for marking integer list in text

    for line in text.splitlines():
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

    l.debug("Number of integer lists read: %d" % len(iList))
    return iList


def polymesh_get_boundary(text):
    '''Creates boundary objects from PolyMesh boundary text string'''

    import re
    inside = False # boolean for marking boundary entries in text

    for line in text.splitlines():
        # Opening of integer list by single parenthesis
        regex = re.search(r'^\(', line, re.M)
        if regex:
            inside = True

        # Closing of integer list by single parenthesis
        regex2 = re.search(r'^\)', line, re.M)
        if regex2:
            inside = False

        if not inside:
            continue

        # New entry is a word (with possibly special characters) on its own line
        regex = re.search(r'^\s+([\w\%\:\-]+)$', line, re.M)
        if regex:
            patchname = str(regex.group(1))
            l.debug("Reading in boundary patch %d: %s" % (len(ug.ugboundaries), patchname))
            patch = ug.UGBoundary(patchname)
            ug.ugboundaries.append(patch)
            continue

        # type
        regex = re.search(r'^\s+type\s+(\w+)\;$', line, re.M)
        if regex:
            patch.typename = str(regex.group(1))
            continue

        # inGroups
        regex = re.search(r'^\s+inGroups\s+([\w\s\(\)]+)\;\s*$', line, re.M)
        if regex:
            patch.inGroups = str(regex.group(1))
            continue

        # nFaces
        regex = re.search(r'^\s+nFaces\s+(\d+)\;$', line, re.M)
        if regex:
            patch.nFaces = int(regex.group(1))
            continue

        # startFace
        regex = re.search(r'^\s+startFace\s+(\d+)\;$', line, re.M)
        if regex:
            patch.startFace = int(regex.group(1))
            continue

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
        ob.data.materials[mati].diffuse_color = get_face_color(mati)

        # Set material index to ug.ugfaces
        for i in range(b.startFace, b.startFace + b.nFaces):
            ug.ugfaces[i].mati = mati
        # Set material index for mesh faces
        for i in range(b.nFaces):
            ob.data.polygons[facecount].material_index = mati
            facecount += 1
        mati += 1


def get_face_color(mati):
    '''Gives a color to argument material number'''

    base_colors = [(0.3,0.3,0.3,1), (0,0,1,1), (1,0,0,1), (0,1,0,1), \
             (0.7,0.7,0,1), (0,0.7,0.7,1), (0.7,0,0.7,1)]
    if mati < len(base_colors):
        return base_colors[mati]

    # Get random colors after base colors
    import random
    random.seed(10043 + mati)
    [r, g, b] = [random.random() for i in range(3)]
    return [r, g, b, 1.0]


##### EXPORT #####


class UG_OT_ExportPolyMesh(bpy.types.Operator, ExportHelper):
    '''Export OpenFOAM PolyMesh as Unstructured Grid'''
    bl_idname = "unstructured_grids.export_openfoam_polymesh"
    bl_label = "Export OpenFOAM PolyMesh"

    filename_ext = ".polyMesh" # Dummy, required by ExportHelper
    def execute(self, context):
        ugdata_to_polymesh(self)
        write_polymesh_files(self)
        return {'FINISHED'}


class UG_OT_UGToPolyMesh(bpy.types.Operator):
    '''Generate OpenFOAM PolyMesh file contents from Unstructed Grid data'''
    bl_idname = "unstructured_grids.ug_to_polymesh"
    bl_label = "Generate polyMesh texts from UG data"

    def execute(self, context):
        ugdata_to_polymesh(self)
        return {'FINISHED'}


def ugdata_to_polymesh(self):
    '''Convert UG data into polymesh text data strings'''

    obname = "Unstructured Grid"
    if not obname in bpy.data.objects:
        self.report({'ERROR'}, "No object named %r" % obname)
        return {'FINISHED'}
    ob = bpy.data.objects[obname]
    update_text_points(ob)
    update_ugdata_and_text_faces(ob)
    update_text_owner_neighbour()
    update_text_boundary()
    return None

def update_text_points(ob):
    '''Updates PolyMesh points string contents from Blender object vertices'''

    # Generate new text
    text = of_file_header('vectorField', 'points') + "\n"
    text += str(len(ob.data.vertices)) + "\n(\n"
    for v in ob.data.vertices:
        text += "(" + "%.6g" % v.co.x + " " \
                + "%.6g" % v.co.y + " " \
                + "%.6g" % v.co.z + ")\n"
    text += ")\n"
    bpy.context.scene.ug_props.text_points = text
    l.debug("text_points updated points: %d" % len(ob.data.vertices))
    return None


def update_ugdata_and_text_faces(ob):
    '''Updates UG data (properties of ugcells, ugfaces and ugboundaries)
    and generates PolyMesh faces text string contents for object ob.
    UGFace indexing is updated from current Blender face material
    assignments, but otherwise it is assumed that UG data is up-to-date.
    Update is done in two phases, at the same time as face text
    definitions are generated:

    1. internal face pass:

    Generates cell indices according to PolyMesh requirement that
    face normal points from lower cell index to higher cell index.
    Face normal is determined by right hand rule from vertex list.

    2. boundary face pass:

    Update boundary list to conform to current object material slot
    and face material assignments. Boundary faces are numbered
    accordingly to material assignments.
    '''

    def update_ugboundaryfaces():
        '''Regenerates list of boundary faces from ugfaces'''
        ug.ugboundaryfaces = []
        for f in ug.ugfaces:
            if f.deleted:
                continue
            if f.neighbouri == -1:
                ug.ugboundaryfaces.append(f)

    def gen_line(verts):
        '''Construct face definition text line from verts list'''
        line = str(len(verts)) + "("
        for j in range(len(verts) - 1):
            line += str(verts[j]) + " "
        line += str(verts[-1]) + ")\n"
        return line

    def internal_face_pass():
        '''Generate face definition text for internal faces.
        Return text and number of internal faces.
        '''

        text = ''
        ind = 0 # face index
        celli = 0 # cell index

        # Reset all cell indices
        for c in ug.ugcells:
            c.celli = -1

        faces = [f for f in ug.ugfaces if f.neighbouri != -1 and not f.deleted]
        for f in faces:
            f.facei = ind # Set face index
            text += gen_line(f.verts) # Add definition line and proceed

            # Set cell indices
            if ug.ugcells[f.owneri].celli == -1:
                ug.ugcells[f.owneri].celli = celli
                celli += 1
            if ug.ugcells[f.neighbouri].celli == -1:
                ug.ugcells[f.neighbouri].celli = celli
                celli += 1
            ind += 1
        return text, ind

    def boundary_face_pass(ind, ob):
        '''Generate face definition text for boundary faces.
        Return text and number of internal faces.
        '''

        text = ''

        # Mark all boundaries as deleted
        for b in ug.ugboundaries:
            b.deleted = True

        for sloti in range(len(ob.material_slots)):
            mat = ob.material_slots[sloti]
            l.debug("Processing patch " + mat.name)

            # Find or create UGBoundary for material in this slot
            patch = [b for b in ug.ugboundaries if b.patchname == mat.name][0]
            if not patch:
                patch = ug.UGBoundary(mat.name)
                ug.ugboundaries.append(patch)

            # Update boundary properties
            patch.deleted = False
            patch.mati = sloti
            patch.startFace = ind

            # Find faces whose material index equals slot number
            faces = [p for p in ob.data.polygons if p.material_index == sloti]
            patch.nFaces = len(faces)

            for p in faces:
                ugface = ug.get_ugface_from_polygon(p)
                # Update face properties, generate text and increase face index
                ugface.mati = sloti
                ugface.facei = ind
                text += gen_line(ugface.verts)
                ind += 1
        return text, ind

    update_ugboundaryfaces()
    text_internal, i = internal_face_pass()
    l.debug("text_faces updated internal faces: %d", i)
    text_boundary, i = boundary_face_pass(i, ob)
    l.debug("text_faces updated total number of faces: %d", i)

    # Generate text string
    text = of_file_header('faceList', 'faces') + "\n"
    text += str(i) + "\n(\n"
    text += text_internal + text_boundary + ")\n"

    bpy.context.scene.ug_props.text_faces = text
    return None


def update_text_owner_neighbour():
    '''Updates PolyMesh owner and neighbour text string contents from UG data'''

    all_faces = [f for f in ug.ugfaces if not f.deleted]
    nall = len(all_faces)
    neighbour_faces = [f for f in ug.ugfaces if f.neighbouri != -1 and not f.deleted]
    nneighbour = len(neighbour_faces)

    # Generate text string
    text_owner = of_file_header('labelList', 'owner') + "\n"
    text_owner += str(nall) + "\n(\n"
    text_neighbour = of_file_header('labelList', 'neighbour') + "\n"
    text_neighbour += str(nneighbour) + "\n(\n"

    for f in all_faces:
        text_owner += str(f.owneri) + "\n"
    for f in neighbour_faces:
        text_neighbour += str(f.neighbouri) + "\n"

    text_owner += ")\n"
    text_neighbour += ")\n"

    bpy.context.scene.ug_props.text_owner = text_owner
    bpy.context.scene.ug_props.text_neighbour = text_neighbour
    l.debug("text_owner updated faces: %d" % nall)
    l.debug("text_neighbour updated faces: %d" % nneighbour)
    return None


def update_text_boundary():
    '''Updates PolyMesh boundary text string contents from UG data'''

    mati = 0 # Boundary patch number
    btext = '' # generated boundary entries
    nboundaries = 0 # number of boundaries
    for i in range(len(ug.ugboundaries)):
        patchlist = [b for b in ug.ugboundaries if b.mati == i and not b.deleted]
        if len(patchlist) != 1:
            continue
        patch = patchlist[0]
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


def write_polymesh_files(self):
    '''Write contents of data strings into PolyMesh files'''

    import os
    ug_props = bpy.context.scene.ug_props
    dirpath = os.path.dirname(self.filepath)

    filenames = ['boundary', 'faces', 'neighbour', 'owner', 'points']
    for f in filenames:
        varname = "text_" + f
        filepath = os.path.join(dirpath, f)
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
