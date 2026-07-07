#ifndef slic3r_Eridian_hpp_
#define slic3r_Eridian_hpp_

#include <functional>
#include <string>
#include <vector>

#include "libslic3r.h"
#include "ExPolygon.hpp"

namespace Slic3r {

class GCodeWriter;

// "Eridian Mode" — an alternative wire-printing algorithm (named after the Eridians in
// "Project Hail Mary"). Rather than the Cura-style ring cage, it builds a 3D truss:
//
//   * every layer is a flat triangular lattice net clipped to the cross-section,
//   * a global lattice is shared by all layers, so a vertex sits at the same (x, y) on
//     every layer and vertical "pillars" connect a vertex to itself one layer up,
//   * where a vertex appears with no pillar beneath it (an overhang), it is braced outward
//     with diagonal struts down to the active lattice neighbours on the layer below.

struct EridianParams
{
    double pillar_height    = 2.0; // vertical step between lattice layers (mm)
    double lattice_spacing  = 5.0; // triangular lattice edge length (mm)

    double mm3_per_mm_flat   = 0.; // volumetric flow for the flat net lines
    double mm3_per_mm_pillar = 0.; // volumetric flow for pillars and braces (over-extruded)

    double speed_flat   = 5.; // mm/s
    double speed_pillar = 5.; // mm/s

    double pillar_lift      = 1.; // extra travel along the pillar axis after extrusion stops (mm)
    double pillar_prime     = 0.8; // filament re-primed after each pillar lift (mm)
    double pillar_exclusion = 8.;  // keep-out radius around a printed pillar; no pillar within (mm)
};

// Emit the Eridian truss for cross-sections `layers` sliced at heights `z` (mm), in absolute
// world millimeter coordinates, through `writer`. `emit` receives ready-to-write G-code chunks.
void write_eridian_gcode(const std::vector<ExPolygons>                 &layers,
                         const std::vector<double>                     &z,
                         const EridianParams                           &params,
                         GCodeWriter                                   &writer,
                         const std::function<void(const std::string&)> &emit);

} // namespace Slic3r

#endif // slic3r_Eridian_hpp_
