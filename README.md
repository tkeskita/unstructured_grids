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

This add-on is currently at proof-of-concept stage and work in progress.

## Idea Description

- Cell and face information related to unstructured grids are kept in
  separate Python object data model

- Besides moving of vertices, modifications of unstructured grids
  rely on special operator which keep data model and Blender mesh
  object contents in sync

- UG should support OpenFOAM PolyMesh features. Work is started by
  implementing PolyMesh import and export functions.
