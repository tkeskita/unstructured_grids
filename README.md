# Unstructured Grids Add-on for Blender

## Introduction

Unstructured Grids (UG) is an add-on for [Blender
2.8](https://www.blender.org/2-8) for importing, editing and exporting
of 3D volume meshes composed of arbitrary polyhedron cells
(a.k.a [unstructured grids](https://en.wikipedia.org/wiki/Unstructured_grid)).
The aim is to create a tool that can be used for minor editing of
unstructured grids (like moving vertices, and assigning boundary faces
to named patches) and possibly some special meshing cases.
Since volume meshes are not natively supported in Blender, it is
currently unknown how usable Blender can become for unstructured
grids.

## Idea Description

- Cell and face information related to unstructured grids are kept in
  separate Python object data model and stored as text strings

- Cell description is compatible with
  [OpenFOAM PolyMesh description](https://cfd.direct/openfoam/user-guide/v7-mesh-description/),
  which describes unstructured grid by lists of cells, cell faces and face vertices

- Besides moving of vertices and assigning of boundary faces to materials
  corresponding to boundary patches, modifications of unstructured grids
  rely on special operator which keep data model and Blender mesh
  object contents in sync

## Status

This add-on is currently at proof-of-concept stage and work in progress.
Currently implemented features include:

- Import and Export of PolyMesh files (boundary, faces, neighbour, owner, points)

- Storage of unstructured grid data as strings inside Blend files

- Boundary face patches can be changed by assigning material for faces

### OpenFOAM Trade Mark Notice

This offering is not approved or endorsed by OpenCFD Limited, producer
and distributor of the OpenFOAM software via www.openfoam.com, and
owner of the OPENFOAM® and OpenCFD® trade marks.
