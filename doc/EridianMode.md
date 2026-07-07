# Eridian Mode

Eridian Mode is an experimental *wire printing* mode for PrusaSlicer. Instead of
slicing a model into solid, stacked layers, it prints the object as a sparse 3D
**truss** — a lattice of struts drawn in the air — using far less material and time
at the cost of being a skeletal representation rather than a solid part.

It is named after the Eridians in Andy Weir's *Project Hail Mary*, who build
things as robust, over‑engineered frameworks.

## The idea in one picture

A normal print is a stack of filled 2D layers. Eridian Mode keeps the "stack of
layers" idea but replaces each solid layer with a **flat triangular net**, and
replaces the solid walls between layers with **pillars and diagonal braces**:

```
   layer k     ●───●───●───●        flat triangular net (follows the outline)
               │ ╲ │ ╲ │ ╲ │        pillars (vertical) + braces (leaning)
   layer k-1   ●───●───●───●
               │ ╲ │ ╲ │ ╲ │
   layer k-2   ●───●───●───●
                    ...
   bed         ▔▔▔▔▔▔▔▔▔▔▔▔▔▔
```

The result is a cage that follows the silhouette of the object at every height.

## Step 1 — Re-sample the model into layers

The model is sliced normally, then the cross-sections are re-sampled at a coarse
vertical step, the **pillar height** (default 2 mm). Each sample is the outline
(an `ExPolygon`, with holes) of the model at that height. These become the layers
of the truss.

## Step 2 — Build a net for each layer

Each layer's outline is turned into a triangular net that *conforms to the shape*:

1. **Boundary.** The outer contour and every hole are walked and re-sampled at the
   **lattice spacing** (default 5 mm), dropping a vertex every step. These vertices
   sit exactly on the object's outline, so the net traces the silhouette.

2. **Interior.** Points are taken from a single **global triangular lattice** that
   is shared by every layer. Only the lattice points that fall safely inside the
   outline (inset by half a spacing, and not too close to a boundary vertex) are
   kept. Because the lattice is global, an interior point sits at the same `(x, y)`
   on every layer — which is what lets pillars run perfectly vertical.

3. **Triangulation.** The boundary vertices (as edge constraints) and the interior
   points are fed to a **constrained Delaunay triangulation** (CGAL). This produces
   a clean triangular mesh that respects the outline and holes. The edges of that
   mesh are the flat net that gets printed for the layer.

## Step 3 — Connect the layers

Going up one layer at a time, every vertex of layer *k* is tied down to layer
*k − 1*:

- **Pillars.** An interior vertex exists at the same lattice coordinate on both
  layers, so it is joined straight down to itself: a vertical pillar.

- **Braces.** A boundary vertex (or any vertex without a partner directly below) is
  joined to the *nearest* vertex on the layer below. Where the object is vertical
  this is almost straight down; where the object flares out or overhangs, the
  nearest point below is set inward, so the strut leans outward — automatically
  forming the diagonal **braces** that carry overhangs. Struts longer than three
  lattice spacings are skipped to avoid wild jumps across gaps.

Each strut is a single diagonal extrusion move; the molten filament cools quickly
enough to hold its shape in the air.

## Step 4 — Emit G-code

The truss is written directly through the G-code writer as 3D moves
(`extrude_to_xyz`), reusing PrusaSlicer's normal start/end G-code, temperatures and
statistics. Flat-net lines use the *flat* flow and speed; pillars and braces use
the *pillar* flow and speed (typically slower, to give the airborne filament time
to solidify).

## Parameters

| Setting | Meaning | Default |
| --- | --- | --- |
| `eridian_mode` | Master toggle for the mode | off |
| `eridian_pillar_height` | Vertical step between layers | 2 mm |
| `eridian_lattice_spacing` | Triangle edge length of the net | 5 mm |
| `eridian_flow` | Flow multiplier for all truss lines | 100 % |
| `eridian_speed_flat` | Speed of the flat lattice lines | 5 mm/s |
| `eridian_speed_pillar` | Speed of pillars and braces | 5 mm/s |

Smaller lattice spacing and pillar height give a denser, stronger cage but a
longer print.

## Limitations (current state)

- The pillar height is quantized to the model's slicing layer height, because the
  cross-sections are re-sampled from the existing slices rather than re-sliced.
- Very thin features can make the coarsely re-sampled boundary self-intersect; the
  affected layers fall back to a sparser boundary-only net.
- The normal planar pipeline (perimeters, infill, supports) still runs before the
  truss is generated, which is wasted work.
- G-code preview and print-time estimates assume planar layers, so they may not
  render the diagonal moves correctly.
