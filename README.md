# Unstructured Grids Add-on for Blender

<p align="left"><img src="examples/ug_title.png"></p>

## Introduction

Unstructured Grids (UG) is an add-on for [Blender
2.80](https://www.blender.org) for importing, editing and exporting
of 3D volume meshes composed of arbitrary polyhedron cells
(a.k.a [3D unstructured grids](https://en.wikipedia.org/wiki/Unstructured_grid)).
Editing includes tasks like moving of selected vertices, and assigning selected
boundary faces to named patches. The user of the add-on is assumed to know
Blender modelling and material systems on a basic level.


## Features and Limitations

- Since volume meshes are not natively supported in Blender, 
  cell and face information related to unstructured grids are kept in
  separate Python object data model and stored as text strings.
  Internal faces or edges are not shown in Blender, but vertices are
  visible in Edit Mode.

- Cell description is compatible with
  [OpenFOAM PolyMesh description](https://cfd.direct/openfoam/user-guide/mesh-description/).
  Unstructured grid is defined by lists of cells, cell faces and face
  vertices.

- Except for moving of vertices, assigning of boundary faces to materials
  (boundary patches) and assigning vertices to vertex groups (zones),
  modifications of unstructured grids rely on special operators which
  keep UG Data and Blender mesh object contents in sync.


## Use Case Examples

- Change boundary patch assignments for existing/new patches (select
  faces in Blender, then assign selection to existing/new
  material).

  <p align="left"><img src="examples/ug_boundary_patch_assign.png"></p>

- Moving of vertices. This can be applied for tasks like:
  
  - Elongation/stretching of cells (by using Proportional Editing in
    Blender)

  - Curving simulation domain, to model e.g. pipes (by applying Curve
    Modifier in Blender)

  <p align="left"><img src="examples/ug_stretch_and_bend.png"></p>

  - Scaling/moving/rotation of a selected part of the mesh


## Status

This add-on is currently at proof-of-concept stage and work in progress.
Currently implemented features include:

- Import and Export of uncompressed PolyMesh files (boundary, faces,
  neighbour, owner, points). PolyMesh Import and Export operators are located
  in File menu under Import and Export.

- Unstructured grid data is saved as text strings (UG Storage) inside Blend files.

- Undo (CTRL + Z) is not working correctly, but undo is supported indirectly:
  Use operator *Update UG Data and Storage from Blender* to sync
  changes made in Blender for safekeeping in UG Storage. To undo, use
  operator *Restore UG Data from UG Storage*.

- Boundary face patches can be changed by assigning material for faces

- Selection of cells via operators (UG Select Cells Inclusive/Exclusive)

- Cell and face zones are visualized by vertex groups. Vertex groups
  for cell zones can be edited.


## Installation

- Add-on code is available at
  https://github.com/tkeskita/unstructured_grids. To download add-on from
  Github, Select “Clone or download”, then “Download ZIP”.

- Start Blender 2.80, go to “File” –> “User Preferences” –> “Add-ons” –> “Install” –> open the add-on zip file.

- Activate the “Unstructured Grids for Blender” add-on in Preferences. Add-on is located in
  Mesh category.

- Click “Save Preferences” to autoload add-on every time Blender is started

- Note: Python Logging console messages may be useful in case of problems.
  More information in file *__init__.py*.


## Development Ideas for Future

- Add support for face zone editing (face normal direction and flipMap)

- Deletion of selected cells to carve voids into the domain

- Extrusion of new cells (one or several layers) from selected faces.
  At least initially this would be without any sense of topology (no
  detection of corners or inward edges, just create side boundary
  faces blindly).

- Subdivision of cells to produce smaller cells

- Merge selected cells

- Import/export of VTK Unstructured Grids


### OpenFOAM Trade Mark Notice

This offering is not approved or endorsed by OpenCFD Limited, producer
and distributor of the OpenFOAM software via www.openfoam.com, and
owner of the OPENFOAM® and OpenCFD® trade marks.
