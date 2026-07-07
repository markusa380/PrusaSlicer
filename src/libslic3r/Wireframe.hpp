#ifndef slic3r_Wireframe_hpp_
#define slic3r_Wireframe_hpp_

#include <functional>
#include <string>
#include <vector>

#include "libslic3r.h"
#include "Point.hpp"
#include "Polygon.hpp"

namespace Slic3r {

class Print;
class GCodeWriter;

// "Wire Printing" (a port of Cura's WirePrint / Weaver + Wireframe2gcode).
//
// The model is not sliced into solid planar layers. Instead it is re-sliced at a coarse
// "connection height" and every cross-section is printed as a ring of horizontal struts,
// consecutive rings being joined by diagonal lines extruded up into the air. Flat top/bottom
// surfaces are filled by concentric rings ("roofs" / "floors"). The result is a sparse cage.
//
// Stage 1 (Weaver): the sliced contours -> a WireFrame data model (this file + Wireframe.cpp).
// Stage 2 (Wireframe2gcode): the WireFrame -> G-code, driven through a GCodeWriter.

// All 3D geometry below is stored in unscaled millimeters (Vec3d). 2D contours are kept in the
// usual scaled integer space (Polygons) until they are turned into 3D connection segments.

enum class WeaveSegmentType {
    UP,           // diagonal move up into the air (from a lower ring to an upper node)
    DOWN,         // diagonal move back down onto the previous ring
    FLAT,         // horizontal line lying on a cross-section
    MOVE,         // travel move (no extrusion)
    DOWN_AND_FLAT // roof/floor segment that is both the end of a down move and a flat line
};

struct WeaveConnectionSegment
{
    Vec3d            to;   // destination, millimeters
    WeaveSegmentType segment_type;
    WeaveConnectionSegment(const Vec3d &to, WeaveSegmentType t) : to(to), segment_type(t) {}
};

struct PolyLine3
{
    Vec3d                               from{ Vec3d::Zero() };
    std::vector<WeaveConnectionSegment> segments;
};

struct WeaveConnectionPart
{
    PolyLine3 connection;
    int       supported_index;
    explicit WeaveConnectionPart(int top_idx) : supported_index(top_idx) {}
};

struct WeaveConnection
{
    double                           z0 = 0.; // height of the supporting polygons (mm)
    double                           z1 = 0.; // height of the supported polygons (mm)
    std::vector<WeaveConnectionPart> connections;
    Polygons                         supported; // polygons to be supported (scaled)
};

// Horizontal fill of a flat surface: concentric inset rings connected to each other.
struct WeaveRoofPart : WeaveConnection {};

struct WeaveRoof
{
    std::vector<WeaveRoofPart> roof_insets;   // connections between consecutive insets
    Polygons                   roof_outlines; // the area filled horizontally (scaled)
};

struct WeaveLayer : WeaveConnection
{
    WeaveRoof roofs; // roofs (top faces) and floors (bottom faces) filled horizontally
};

struct WireFrame
{
    WeaveRoof                bottom_infill;  // the very first (bed) layer, filled horizontally
    Polygons                 bottom_outline; // outline of the first cross-section (scaled)
    double                   z_bottom = 0.;  // height of the first cross-section (mm)
    std::vector<WeaveLayer>  layers;
};

// Geometry parameters for the Weaver stage, all in millimeters / radians.
struct WireframeParams
{
    double connection_height   = 3.0;  // wireframe_height
    double roof_inset          = 3.0;  // wireframe_roof_inset
    double nozzle_clearance    = 1.0;  // wireframe_nozzle_clearance
    double nozzle_outer_diam   = 1.0;  // wireframe_nozzle_outer_diameter
    double nozzle_expansion    = 0.785398; // wireframe_nozzle_expansion_angle (radians)

    // Minimum spacing between neighbouring nodes on a cross-section (mm).
    double node_spacing() const { return std::tan(nozzle_expansion) * connection_height + nozzle_outer_diam + nozzle_clearance; }
};

// Stage 1: build the WireFrame data model from cross-sections sliced at heights `slice_z` (mm).
WireFrame build_wireframe(const std::vector<Polygons> &slices,
                          const std::vector<double>   &slice_z,
                          const WireframeParams       &params);

// G-code generation parameters (speeds mm/s, distances mm, delays s), read from PrintConfig.
enum class WireframeStrategyKind { Compensate, Knot, Retract };

struct WireframeGCodeParams
{
    double mm3_per_mm_connection = 0.; // volumetric flow for connection lines
    double mm3_per_mm_flat       = 0.; // volumetric flow for flat lines
    double node_spacing          = 0.; // == WireframeParams::node_spacing(), mm

    double speed_bottom = 5., speed_up = 5., speed_down = 5., speed_flat = 5.; // mm/s
    double flat_delay = 0.2, bottom_delay = 0., top_delay = 0.;                // s
    double up_half_speed = 0.3;          // mm printed at half speed at start of an up move
    double top_jump = 0.6;               // knot size
    double fall_down = 0.5, drag_along = 0.6;
    double straight_before_down = 0.2;   // fraction (0..1)
    WireframeStrategyKind strategy = WireframeStrategyKind::Compensate;

    double roof_fall_down = 2., roof_drag_along = 0.8, roof_outer_delay = 0.2;
    double connection_height = 3.;
};

// Stage 2: emit the WireFrame as G-code. `emit` receives ready-to-write G-code chunks
// (the caller appends them to the output file). `writer` is used to format moves/retractions
// and must have its active extruder already set.
void write_wireframe_gcode(const WireFrame            &wire_frame,
                           GCodeWriter                &writer,
                           const WireframeGCodeParams &params,
                           const std::function<void(const std::string&)> &emit);

} // namespace Slic3r

#endif // slic3r_Wireframe_hpp_
