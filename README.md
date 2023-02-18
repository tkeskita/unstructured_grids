# Unstructured Grids for Blender

<p align="left"><img src="docs/images/ug_title.png"></p>

## Introduction

[Unstructured Grids (UG)](https://github.com/tkeskita/unstructured_grids)
is an add-on for [Blender](https://www.blender.org)
for creating, importing, editing and exporting of
3D volume meshes composed of arbitrary polyhedron cells (a.k.a [3D
unstructured grids](https://en.wikipedia.org/wiki/Unstructured_grid)).
The motivation for this work was the lack of open source volume
mesh editors. There are numerous finite volume mesh generators,
but practically no general or visual editors.

Add-on handles mesh topology and geometry (the definition of cells).
Field data (cell data or point data) for the mesh is disregarded.
Editing includes tasks like moving of selected vertices, deletion of
existing cells, extrusion of new cells, and assignation of selected
faces and cells to named boundaries and zones. The user of the add-on
is assumed to know Blender modelling and material systems on a basic
level.

The add-on has been tested with
[Blender 2.93 LTS version](https://www.blender.org/) and
[OpenFOAM Foundation](https://openfoam.org/) version 9 of OpenFOAM.


## Documentation

Documentation (made using [Sphinx](https://www.sphinx-doc.org/en/master/))
is located in docs directory of the sources and is viewable online at
https://unstructured-grids.readthedocs.io.


## Issues and feedback

Please use
[Github issues](https://github.com/tkeskita/unstructured_grids/issues)
for bug reporting. There is also
[a discussion thread on CFD Online](https://www.cfd-online.com/Forums/openfoam-community-contributions/219493-unstructured-grids-add-blender.html)
and another one on [VTK discourse forum](https://discourse.vtk.org/t/unstructured-grids-for-blender/1959).

If you use this add-on, please star the project in GitHub!

## Development Ideas for Future

- Don't save Unstructured Grid object in blend file to save disk space
  (mesh is anyways thrown away and regenerated from UG text
  strings)

- Add Gzip option to OpenFOAM PolyMesh export

- Split Cells (Edge Division) (e.g. hexahedral cell is split to 8
  cells, tetrahedral cell is split into 4 cells)

- Dissolve Internal Faces (merge cells)

- Merge Overlapping Boundary Faces (e.g. to merge two mesh regions)

- Separate Internal Faces (convert internal face to two boundary faces)

- Split Cells (Center Slicing) (polyhedron cell face vertices are
  connected to cell center)

- Show Cells (creates temporary object with faces around
  selected cells)

- Fix Intersecting Boundaries (move boundary vertices so
  that boundary faces don't intersect each other)

- Support for OpenFOAM sets?

- Mesh Quality metrics + operator to select bad cells. Definition of
  bad cells depend on target simulation type, so there would need to
  be several metrics. This would be nice-to-have but not essential,
  since it is possible to use *checkMesh* utility of OpenFOAM, or
  e.g. Mesh Quality filter in Paraview or e.g.
  [BVTKNodes](https://github.com/tkeskita/BVtkNodes)
  to evaluate mesh quality.


### OpenFOAM Trade Mark Notice

This offering is not approved or endorsed by OpenCFD Limited, producer
and distributor of the OpenFOAM software via www.openfoam.com, and
owner of the OPENFOAM® and OpenCFD® trade marks.
