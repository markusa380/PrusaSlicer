#include "Wireframe.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>

#include "ClipperUtils.hpp"
#include "ExPolygon.hpp"
#include "Line.hpp"
#include "GCode/GCodeWriter.hpp"

// Port of Cura's Weaver (mesh cross-sections -> WireFrame) and Wireframe2gcode (WireFrame -> G-code).
// See Wireframe.hpp for the data model. All 3D geometry is in millimeters; 2D contours stay scaled.

namespace Slic3r {

// ---------------------------------------------------------------------------------------------
// small geometry helpers (replacing Cura's PolygonUtils)
// ---------------------------------------------------------------------------------------------

// Nearest point lying on any segment of `polys` to `p` (scaled space).
static Point nearest_point_on_polygons(const Point &p, const Polygons &polys)
{
    Point   best = p;
    double  best_d2 = std::numeric_limits<double>::max();
    for (const Polygon &poly : polys) {
        const size_t n = poly.points.size();
        if (n == 0) continue;
        for (size_t i = 0; i < n; ++i) {
            const Point &a = poly.points[i];
            const Point &b = poly.points[(i + 1) % n];
            const Vec2d  av = a.cast<double>();
            const Vec2d  bv = b.cast<double>();
            const Vec2d  pv = p.cast<double>();
            Vec2d        ab = bv - av;
            double       len2 = ab.squaredNorm();
            Vec2d        proj = av;
            if (len2 > 0) {
                double t = std::clamp((pv - av).dot(ab) / len2, 0.0, 1.0);
                proj = av + t * ab;
            }
            double d2 = (pv - proj).squaredNorm();
            if (d2 < best_d2) {
                best_d2 = d2;
                best = Point(coord_t(std::lround(proj.x())), coord_t(std::lround(proj.y())));
            }
        }
    }
    return best;
}

// Index of the polygon vertex nearest to `p`.
static size_t nearest_vertex_index(const Point &p, const Polygon &poly)
{
    size_t best = 0;
    double best_d2 = std::numeric_limits<double>::max();
    for (size_t i = 0; i < poly.points.size(); ++i) {
        double d2 = (poly.points[i] - p).cast<double>().squaredNorm();
        if (d2 < best_d2) { best_d2 = d2; best = i; }
    }
    return best;
}

// Resample every polygon of `parts` into a chain of nodes spaced `spacing` apart along the
// perimeter, starting near `start_close_to`. Mirrors Cura's Weaver::chainify_polygons.
static Polygons chainify(const Polygons &parts, Point &start_close_to, coord_t spacing)
{
    Polygons result;
    if (spacing <= 0) return result;
    for (const Polygon &poly : parts) {
        const size_t n = poly.points.size();
        if (n < 2) continue;

        // Reorder the closed loop so it begins at the vertex closest to start_close_to.
        size_t start = nearest_vertex_index(start_close_to, poly);
        std::vector<Point> loop;
        loop.reserve(n + 1);
        for (size_t k = 0; k <= n; ++k)
            loop.push_back(poly.points[(start + k) % n]);

        Polygon out;
        double  carry = 0.;
        out.points.push_back(loop.front());
        for (size_t i = 0; i + 1 < loop.size(); ++i) {
            Vec2d  a = loop[i].cast<double>();
            Vec2d  b = loop[i + 1].cast<double>();
            double seg = (b - a).norm();
            if (seg <= 0) continue;
            double pos = spacing - carry;
            while (pos <= seg) {
                Vec2d q = a + (b - a) * (pos / seg);
                out.points.emplace_back(coord_t(std::lround(q.x())), coord_t(std::lround(q.y())));
                pos += spacing;
            }
            carry = seg - (pos - spacing);
        }
        if (out.points.size() >= 2) {
            start_close_to = out.points.back();
            result.push_back(std::move(out));
        }
    }
    return result;
}

// Split a set of polygons into ExPolygons (outline + holes).
static ExPolygons split_into_parts(const Polygons &polys) { return union_ex(polys); }

// Outlines only (drop holes).
static Polygons outlines_of(const Polygons &polys)
{
    Polygons out;
    for (const ExPolygon &ex : split_into_parts(polys))
        out.push_back(ex.contour);
    return out;
}

// ---------------------------------------------------------------------------------------------
// Stage 1: Weaver
// ---------------------------------------------------------------------------------------------

namespace {

struct Weaver
{
    const WireframeParams &p;
    coord_t                spacing; // node spacing, scaled

    explicit Weaver(const WireframeParams &params)
        : p(params), spacing(coord_t(scale_(params.node_spacing()))) {}

    // For every point of `supported`, drop straight down to the nearest point of `supporting`,
    // producing a DOWN + UP zig-zag connection (Cura::connect_polygons).
    void connect_polygons(const Polygons &supporting, double z0,
                          const Polygons &supported,  double z1,
                          WeaveConnection &result)
    {
        if (supporting.empty()) return;
        result.z0 = z0;
        result.z1 = z1;
        for (size_t prt = 0; prt < supported.size(); ++prt) {
            const Polygon &upper = supported[prt];
            result.connections.emplace_back(int(prt));
            PolyLine3 &conn = result.connections.back().connection;
            bool first = true;
            for (const Point &up : upper.points) {
                Point  low = nearest_point_on_polygons(up, supporting);
                Vec3d  low3(unscale<double>(low.x()), unscale<double>(low.y()), z0);
                Vec3d  up3 (unscale<double>(up.x()),  unscale<double>(up.y()),  z1);
                if (first) conn.from = low3;
                else       conn.segments.emplace_back(low3, WeaveSegmentType::DOWN);
                conn.segments.emplace_back(up3, WeaveSegmentType::UP);
                first = false;
            }
        }
    }

    // chainify `parts1` then connect it down to `parts0` (Cura::connect).
    void connect(const Polygons &parts0, double z0, const Polygons &parts1, double z1, WeaveConnection &result)
    {
        if (parts1.empty()) return;
        Point start = parts0.empty() ? parts1.back().points.back() : parts0.back().points.back();
        result.supported = chainify(parts1, start, spacing);
        if (parts0.empty()) return;
        connect_polygons(parts0, z0, result.supported, z1, result);
    }

    // Fill a flat top (roof) with concentric rings connected inward (Cura::fillRoofs).
    // A hard cap on the number of concentric insets, to guarantee the fill loops terminate even
    // for pathological geometry (holes, degenerate offsets, roof_inset == 0).
    static constexpr int max_insets = 4096;

    static double total_area(const Polygons &ps)
    {
        double a = 0.;
        for (const Polygon &pl : ps) a += std::abs(pl.area());
        return a;
    }

    void fill_roofs(const Polygons &supporting, const Polygons &to_be_supported, int direction, double z, WeaveRoof &horizontals)
    {
        if (supporting.empty() || p.roof_inset <= 0.) return;
        Polygons roofs = diff(supporting, to_be_supported);
        roofs = offset(offset(roofs, -float(scale_(p.roof_inset))), float(scale_(p.roof_inset)));
        if (roofs.empty()) return;

        Polygons roof_outlines, roof_holes;
        for (const ExPolygon &part : split_into_parts(roofs)) {
            roof_outlines.push_back(part.contour);
            for (const Polygon &hole : part.holes) {
                roof_holes.push_back(hole);
                roof_holes.back().reverse();
            }
        }

        Polygons supporting_outlines = outlines_of(supporting);

        Polygons last_supported = supporting;
        // `inset0` advances by the raw inward offset, which shrinks monotonically to empty (this is
        // what terminates the loop). `inset1` (clipped to the roof, holes kept open) is what we print.
        Polygons inset0 = supporting_outlines;
        for (int guard = 0; ! inset0.empty() && guard < max_insets; ++ guard) {
            Polygons grown  = offset(inset0, float(direction) * float(scale_(p.roof_inset)), jtRound);
            Polygons inset1 = intersection(grown, roof_outlines);
            inset1 = union_(inset1, roof_holes);
            if (inset1.empty()) break;

            horizontals.roof_insets.emplace_back();
            connect(last_supported, z, inset1, z, horizontals.roof_insets.back());

            last_supported = horizontals.roof_insets.back().supported;
            inset0 = grown;
        }

        append(horizontals.roof_outlines, std::move(roofs));
    }

    // Fill a flat bottom (floor) with concentric rings connected outward (Cura::fillFloors).
    void fill_floors(const Polygons &supporting, const Polygons &to_be_supported, int direction, double z, WeaveRoof &horizontals)
    {
        if (supporting.empty() || to_be_supported.empty() || p.roof_inset <= 0.) return;
        Polygons floors = diff(to_be_supported, supporting);
        floors = offset(offset(floors, -float(scale_(p.roof_inset))), float(scale_(p.roof_inset)));
        if (floors.empty()) return;

        Polygons last_supported = supporting;
        // `outset0` grows outward but is clipped to the floor area; it plateaus once the floor is
        // fully covered, so we stop when the covered area no longer increases.
        Polygons outset0   = supporting;
        double   prev_area = -1.;
        for (int guard = 0; ! outset0.empty() && guard < max_insets; ++ guard) {
            Polygons grown   = offset(outset0, float(direction) * float(scale_(p.roof_inset)), jtRound);
            Polygons outset1 = intersection(grown, floors);
            if (outset1.empty()) break;
            double area = total_area(outset1);
            if (prev_area >= 0. && area <= prev_area + EPSILON) break; // converged: floor filled
            prev_area = area;

            horizontals.roof_insets.emplace_back();
            connect(last_supported, z, outset1, z, horizontals.roof_insets.back());

            last_supported = horizontals.roof_insets.back().supported;
            outset0 = grown;
        }

        append(horizontals.roof_outlines, std::move(floors));
    }

    // Build the roofs (top faces) and floors (bottom faces) of a layer (Cura::createHorizontalFill).
    void create_horizontal_fill(WeaveLayer &layer, const Polygons &layer_above)
    {
        const float bridgable = float(scale_(p.connection_height));
        { // roofs: what is here but not covered above
            Polygons to_be_supported = offset(layer_above, bridgable);
            fill_roofs(layer.supported, to_be_supported, -1, layer.z1, layer.roofs);
        }
        { // floors: what is above but not here
            Polygons to_be_supported = offset(layer_above, -bridgable);
            fill_floors(layer.supported, to_be_supported, 1, layer.z1, layer.roofs);
        }
    }
};

} // namespace

WireFrame build_wireframe(const std::vector<Polygons> &slices,
                          const std::vector<double>   &slice_z,
                          const WireframeParams       &params)
{
    Weaver    weaver(params);
    WireFrame wire_frame;

    // Find the first non-empty cross-section.
    size_t start = 0;
    while (start < slices.size() && slices[start].empty()) ++start;
    if (start >= slices.size()) return wire_frame;

    wire_frame.bottom_outline = slices[start];
    wire_frame.z_bottom       = slice_z[start];

    Point start_point = wire_frame.bottom_outline.empty()
        ? Point(0, 0)
        : (wire_frame.bottom_outline.front().bounding_box().center());

    // Chainify every subsequent cross-section into evenly spaced node rings.
    for (size_t i = start + 1; i < slices.size(); ++i) {
        Polygons chained = chainify(slices[i], start_point, weaver.spacing);
        if (chained.empty()) continue;
        wire_frame.layers.emplace_back();
        WeaveLayer &layer = wire_frame.layers.back();
        layer.z0 = slice_z[i - 1];
        layer.z1 = slice_z[i];
        layer.supported = std::move(chained);
        start_point = layer.supported.back().points.back();
    }

    // Horizontal fills (roofs + floors) for every layer.
    for (size_t i = 0; i < wire_frame.layers.size(); ++i) {
        const Polygons empty;
        const Polygons &above = (i + 1 < wire_frame.layers.size()) ? wire_frame.layers[i + 1].supported : empty;
        weaver.create_horizontal_fill(wire_frame.layers[i], above);
    }

    // Vertical connections between consecutive layers.
    {
        const Polygons *lower_top = &wire_frame.bottom_outline;
        double          last_z    = wire_frame.z_bottom;
        for (WeaveLayer &layer : wire_frame.layers) {
            weaver.connect_polygons(*lower_top, last_z, layer.supported, layer.z1, layer);
            append(layer.supported, layer.roofs.roof_outlines);
            lower_top = &layer.supported;
            last_z    = layer.z1;
        }
    }

    // Top roof and bottom floor.
    if (! wire_frame.layers.empty()) {
        WeaveLayer &top = wire_frame.layers.back();
        Polygons    none;
        weaver.fill_roofs(top.supported, none, -1, top.z1, top.roofs);
        weaver.fill_roofs(wire_frame.bottom_outline, none, -1, wire_frame.layers.front().z0, wire_frame.bottom_infill);
    }

    return wire_frame;
}

// ---------------------------------------------------------------------------------------------
// Stage 2: Wireframe2gcode
// ---------------------------------------------------------------------------------------------

namespace {

struct WireGCode
{
    GCodeWriter                                   &writer;
    const WireframeGCodeParams                    &p;
    const std::function<void(const std::string&)> &emit;
    double                                         cur_speed = -1.;

    WireGCode(GCodeWriter &w, const WireframeGCodeParams &params, const std::function<void(const std::string&)> &e)
        : writer(w), p(params), emit(e) {}

    void set_speed(double mm_s)
    {
        if (mm_s <= 0.) mm_s = 1.;
        if (std::abs(mm_s - cur_speed) < 1e-6) return;
        cur_speed = mm_s;
        emit(writer.set_speed(mm_s * 60.));
    }

    double e_for(const Vec3d &from, const Vec3d &to, double mm3_per_mm)
    {
        double len = (to - from).norm();
        return writer.extruder()->e_per_mm3() * mm3_per_mm * len;
    }

    void extrude(const Vec3d &to, double mm_s, double mm3_per_mm, const std::string_view comment = {})
    {
        set_speed(mm_s);
        Vec3d from = writer.get_position();
        emit(writer.extrude_to_xyz(to, e_for(from, to, mm3_per_mm), comment));
    }

    void travel(const Vec3d &to)
    {
        set_speed(p.speed_flat);
        emit(writer.travel_to_xyz(to));
    }

    void move_with_retract(const Vec3d &to)
    {
        double d2 = (writer.get_position() - to).squaredNorm();
        double thr = p.node_spacing * p.node_spacing * 4.;
        if (d2 >= thr) emit(writer.retract());
        travel(to);
        emit(writer.unretract());
    }

    void delay(double seconds)
    {
        if (seconds <= 0.) return;
        char buf[64];
        snprintf(buf, sizeof(buf), "G4 P%d\n", int(seconds * 1000. + 0.5));
        emit(buf);
    }

    // --- upward-move strategies -------------------------------------------------------------

    void strategy_compensate(const PolyLine3 &conn, size_t idx)
    {
        const Vec3d &from = (idx == 0) ? conn.from : conn.segments[idx - 1].to;
        const Vec3d &seg_to = conn.segments[idx].to;
        Vec3d vector = seg_to - from;
        double vlen = std::max(vector.norm(), 1e-6);
        Vec3d to = seg_to + Vec3d(0, 0, p.fall_down * vlen / p.connection_height);
        Vec3d dir = vector * (p.drag_along / vlen);

        Vec3d next_point = (idx + 1 < conn.segments.size()) ? conn.segments[idx + 1].to : conn.segments[0].to;
        Vec3d next_vec = next_point - seg_to;
        Vec2d next_dir2(next_vec.x(), next_vec.y());
        double nlen = next_dir2.norm();
        Vec3d next_dir(0, 0, 0);
        if (nlen > 0) { next_dir2 *= p.drag_along / nlen; next_dir = Vec3d(next_dir2.x(), next_dir2.y(), 0); }

        Vec3d new_top = to - next_dir + dir;
        double orr = vector.norm() + next_vec.norm() + 1e-3;
        double nl  = (new_top - from).norm() + (next_point - new_top).norm() + 1e-3;
        extrude(new_top, p.speed_up * nl / orr, p.mm3_per_mm_connection * orr / nl, "wireframe up");
    }

    void strategy_knot(const PolyLine3 &conn, size_t idx)
    {
        const Vec3d &seg_to = conn.segments[idx].to;
        extrude(seg_to, p.speed_up, p.mm3_per_mm_connection, "wireframe up");
        Vec3d next_vec = (idx + 1 < conn.segments.size()) ? conn.segments[idx + 1].to - seg_to
                                                          : conn.segments[0].to - seg_to;
        Vec2d nd(next_vec.x(), next_vec.y());
        double nl = nd.norm();
        if (nl <= 0) return;
        nd *= p.top_jump / nl;
        Vec3d pos = writer.get_position();
        Vec3d back(nd.x() / 2, nd.y() / 2, p.top_jump);
        set_speed(p.speed_up);
        emit(writer.travel_to_xyz(pos - back));
        delay(p.top_delay);
        emit(writer.travel_to_xyz(pos + Vec3d(nd.x(), nd.y(), 0)));
    }

    void strategy_retract(const PolyLine3 &conn, size_t idx)
    {
        const Vec3d &seg_to = conn.segments[idx].to;
        extrude(seg_to, p.speed_up, p.mm3_per_mm_connection, "wireframe up");
        emit(writer.retract());
        set_speed(p.speed_flat);
        emit(writer.travel_to_xyz(seg_to + Vec3d(0, 0, 1.0)));
        delay(p.top_delay);
        emit(writer.unretract());
    }

    void go_down(const PolyLine3 &conn, size_t idx)
    {
        const Vec3d &from = (idx == 0) ? conn.from : conn.segments[idx - 1].to;
        const Vec3d &to   = conn.segments[idx].to;
        if (p.straight_before_down <= 0) {
            extrude(to, p.speed_down, p.mm3_per_mm_connection, "wireframe down");
        } else {
            Vec3d vec = to - from;
            Vec3d in_between = from + vec * p.straight_before_down;
            Vec3d up(in_between.x(), in_between.y(), from.z());
            extrude(up, p.speed_down, p.mm3_per_mm_connection, "wireframe down");
            extrude(to, p.speed_down, p.mm3_per_mm_connection, "wireframe down");
        }
        delay(p.bottom_delay);
        if (p.up_half_speed > 0)
            extrude(writer.get_position() + Vec3d(0, 0, p.up_half_speed), p.speed_up / 2, p.mm3_per_mm_connection * 2, "wireframe anchor");
    }

    void handle_segment(const PolyLine3 &conn, size_t idx)
    {
        switch (conn.segments[idx].segment_type) {
        case WeaveSegmentType::MOVE: move_with_retract(conn.segments[idx].to); break;
        case WeaveSegmentType::DOWN: go_down(conn, idx); break;
        case WeaveSegmentType::UP:
            if      (p.strategy == WireframeStrategyKind::Knot)    strategy_knot(conn, idx);
            else if (p.strategy == WireframeStrategyKind::Retract) strategy_retract(conn, idx);
            else                                                    strategy_compensate(conn, idx);
            break;
        default: break;
        }
    }

    void handle_roof_segment(const PolyLine3 &conn, size_t idx)
    {
        const WeaveConnectionSegment &segment = conn.segments[idx];
        const Vec3d &from = (idx == 0) ? conn.from : conn.segments[idx - 1].to;
        const WeaveConnectionSegment *next = (idx + 1 < conn.segments.size()) ? &conn.segments[idx + 1] : nullptr;
        switch (segment.segment_type) {
        case WeaveSegmentType::MOVE:
        case WeaveSegmentType::DOWN_AND_FLAT:
            if (next && next->segment_type != WeaveSegmentType::DOWN_AND_FLAT)
                move_with_retract(segment.to);
            break;
        case WeaveSegmentType::UP: {
            Vec3d to = segment.to + Vec3d(0, 0, p.roof_fall_down);
            Vec3d vector = segment.to - from;
            if (vector.squaredNorm() == 0) return;
            Vec3d dir = vector * (p.roof_drag_along / vector.norm());
            Vec3d next_vec = next ? (next->to - segment.to) : (conn.segments[0].to - segment.to);
            Vec2d nd(next_vec.x(), next_vec.y());
            Vec3d detoured = to + dir;
            if (nd.squaredNorm() > 0) { nd *= p.roof_drag_along / nd.norm(); detoured -= Vec3d(nd.x(), nd.y(), 0); }
            extrude(detoured, p.speed_up, p.mm3_per_mm_connection, "wireframe roof");
            break;
        }
        case WeaveSegmentType::DOWN:
            extrude(segment.to, p.speed_down, p.mm3_per_mm_connection, "wireframe roof");
            delay(p.roof_outer_delay);
            break;
        default: break;
        }
    }

    // Fill a set of roof/floor insets (Cura::writeFill).
    void write_fill(const std::vector<WeaveRoofPart> &insets, const Polygons &roof_outlines,
                    bool is_bottom)
    {
        for (const WeaveRoofPart &inset : insets) {
            for (const WeaveConnectionPart &part : inset.connections) {
                const std::vector<WeaveConnectionSegment> &segs = part.connection.segments;
                if (segs.empty()) continue;
                Vec3d first_from = part.connection.from;
                size_t first = 0;
                for (; first < segs.size() && segs[first].segment_type == WeaveSegmentType::MOVE; ++first)
                    first_from = segs[first].to;
                if (first == segs.size()) continue;
                move_with_retract(first_from);
                for (size_t i = first; i < segs.size(); ++i)
                    handle_roof_segment(part.connection, i);
                // top flat lines
                for (size_t i = 0; i < segs.size(); ++i) {
                    if (segs[i].segment_type == WeaveSegmentType::DOWN) continue;
                    handle_flat(segs[i], is_bottom);
                }
            }
        }
        for (const Polygon &poly : roof_outlines) {
            if (poly.points.empty()) continue;
            Vec3d z = writer.get_position();
            move_with_retract(Vec3d(unscale<double>(poly.points.back().x()), unscale<double>(poly.points.back().y()), z.z()));
            for (const Point &pt : poly.points) {
                Vec3d to(unscale<double>(pt.x()), unscale<double>(pt.y()), z.z());
                WeaveConnectionSegment seg(to, WeaveSegmentType::FLAT);
                handle_flat(seg, is_bottom);
            }
        }
    }

    void handle_flat(const WeaveConnectionSegment &segment, bool is_bottom)
    {
        if (segment.segment_type == WeaveSegmentType::MOVE) { move_with_retract(segment.to); return; }
        if (segment.segment_type == WeaveSegmentType::DOWN_AND_FLAT) return;
        extrude(segment.to, is_bottom ? p.speed_bottom : p.speed_flat, p.mm3_per_mm_flat, "wireframe flat");
        delay(p.flat_delay);
    }
};

} // namespace

void write_wireframe_gcode(const WireFrame            &wire_frame,
                           GCodeWriter                &writer,
                           const WireframeGCodeParams &params,
                           const std::function<void(const std::string&)> &emit)
{
    WireGCode wg(writer, params, emit);

    // Bottom outline lines.
    for (const Polygon &poly : wire_frame.bottom_infill.roof_outlines) {
        if (poly.points.empty()) continue;
        wg.move_with_retract(Vec3d(unscale<double>(poly.points.back().x()), unscale<double>(poly.points.back().y()), wire_frame.z_bottom));
        for (const Point &pt : poly.points)
            wg.extrude(Vec3d(unscale<double>(pt.x()), unscale<double>(pt.y()), wire_frame.z_bottom), params.speed_bottom, params.mm3_per_mm_flat, "wireframe bottom");
    }

    // Bottom infill (floor of the first layer).
    Polygons empty_outlines;
    wg.write_fill(wire_frame.bottom_infill.roof_insets, empty_outlines, /*is_bottom*/true);

    // Layers.
    for (const WeaveLayer &layer : wire_frame.layers) {
        for (const WeaveConnectionPart &part : layer.connections) {
            const PolyLine3 &conn = part.connection;
            if (conn.segments.empty()) continue;
            wg.move_with_retract(conn.from);
            for (size_t i = 0; i < conn.segments.size(); ++i)
                wg.handle_segment(conn, i);
            // top flat lines lying on this layer's cross-section
            for (const WeaveConnectionSegment &seg : conn.segments) {
                if (seg.segment_type == WeaveSegmentType::DOWN) continue;
                if (seg.segment_type == WeaveSegmentType::MOVE) { wg.move_with_retract(seg.to); continue; }
                wg.extrude(seg.to, params.speed_flat, params.mm3_per_mm_flat, "wireframe top");
                wg.delay(params.flat_delay);
            }
        }
        // roofs of this layer
        wg.write_fill(layer.roofs.roof_insets, layer.roofs.roof_outlines, /*is_bottom*/false);
    }

    emit(writer.retract());
}

} // namespace Slic3r
