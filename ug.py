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

import logging
l = logging.getLogger(__name__)

class UGCell:
    '''Class for Unstructured Grid cell data objects'''

    deleted = False # boolean for marking deleted cells
    faces = [] # list of indices of faces in ugfaces that make up this cell
    celli = -1 # cell index number, is updated on export only

    def __init__(self, i):
        '''Initialize UGCell with cell index'''
        self.celli = i
        self.faces = []

ugcells = [] # global list of all UGCells


class UGFace:
    '''Class for Unstructured Grid face data objects'''

    deleted = False # boolean for marking deleted faces
    owneri = -1 # index of UGCell in ugfaces owning this face
    neighbouri = -1 # index of UGCell in ugfaces neighbouring this face
    facei = -1 # face index number, is updated on export only
    verts = [] # list of ordered vertex indices used by this face
    mati = -1 # boundary patch (material) index number

    def __init__(self, i, verts):
        '''Initialize UGFace with face index and vertex list'''
        self.facei = i
        self.verts = verts

ugfaces = [] # global list of all UGFaces
ugboundaryfaces = [] # global short list of boundary UGFaces only

class UGBoundary:
    '''Class for Unstructured Grid patch (boundary face) objects'''

    deleted = False # boolean for patches which contain zero boundary faces
    patchname = 'default' # patch name
    typename = 'default' # patch type name
    inGroups = '' # name for patch group this patch is part of
    nFaces = 0 # number of faces in patch
    startFace = -1 # index to first face for this patch
    mati = -1 # boundary patch (object material slot) index number

    def __init__(self, patchname):
        '''Initialize UGBoundary with name patchname'''
        self.patchname = patchname

ugboundaries = [] # global list of all UGBoundaries


def get_ugface_from_polygon(p):
    '''Returns UGFace matching argument mesh face polygon p.
    Assumes polygon represents a boundary face (not internal face).
    '''

    numverts = len(p.vertices)
    pv0 = p.vertices[0] # first vertex index
    vertfaces = [f for f in ugboundaryfaces if p.vertices[0] in f.verts]

    # Find match for all vertices
    for f in vertfaces:
        n = 0
        for vi in p.vertices:
            if vi in f.verts:
                n += 1
        if n == numverts:
            return f

    verts = [i for i in p.vertices]
    l.debug("No face found with verts " + str(verts))
    return None
