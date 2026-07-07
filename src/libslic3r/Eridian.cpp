#include "Eridian.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <set>
#include <unordered_map>

#include "BoundingBox.hpp"
#include "ClipperUtils.hpp"
#include "Triangulation.hpp"
#include "GCode/GCodeWriter.hpp"

// See Eridian.hpp. For every cross-section a triangular net is built that FOLLOWS the outer
// contour: the boundary (outer contour and holes) is resampled onto the outline, interior points
// are placed on a shared global lattice, and the whole thing is triangulated with a constrained
// Delaunay triangulation. Interior points share a lattice key across layers, giving vertical
// pillars; every other vertex is connected to the nearest vertex on the layer below, which leans
// outward into a brace under overhangs.

namespace Slic3r {

namespace {

// Interior lattice vertex identity, shared across layers -> vertical pillars.
using Key = int64_t;
static inline Key make_key(int col, int row) { return (Key(int32_t(col)) << 32) | Key(uint32_t(int32_t(row))); }

struct Lattice
{
    double s, dy;
    explicit Lattice(double spacing_scaled) : s(spacing_scaled), dy(spacing_scaled * std::sqrt(3.) / 2.) {}
    Point point(int col, int row) const
    {
        const double x = col * s + ((row & 1) ? s * 0.5 : 0.0);
        const double y = row * dy;
        return Point(coord_t(std::llround(x)), coord_t(std::llround(y)));
    }
};

struct Vert
{
    Point p;
    bool  interior = false;
    Key   key = 0; // valid when interior
};

struct Net
{
    std::vector<Vert>                verts;
    std::unordered_map<Key, int>     grid;  // interior key -> vertex index
    std::vector<std::pair<int,int>>  edges; // net edges (indices into verts)
};

// Resample a closed polygon at roughly `spacing` (scaled) along its perimeter, preserving winding.
static Points resample_closed(const Polygon &poly, double spacing)
{
    Points out;
    const Points &p = poly.points;
    const size_t  n = p.size();
    if (n < 2 || spacing <= 0) return out;

    out.push_back(p[0]);
    double carry = 0.;
    for (size_t i = 0; i < n; ++i) {
        const Vec2d a = p[i].cast<double>();
        const Vec2d b = p[(i + 1) % n].cast<double>();
        const double seg = (b - a).norm();
        if (seg <= 0) continue;
        double pos = spacing - carry;
        while (pos <= seg) {
            const Vec2d q = a + (b - a) * (pos / seg);
            out.emplace_back(coord_t(std::llround(q.x())), coord_t(std::llround(q.y())));
            pos += spacing;
        }
        carry = seg - (pos - spacing);
    }
    // Drop the last point if it collapsed back onto the first.
    if (out.size() > 1 && (out.back() - out.front()).cast<double>().norm() < spacing * 0.5)
        out.pop_back();
    return out;
}

static bool point_inside(const ExPolygons &xs, const Point &pt)
{
    for (const ExPolygon &ex : xs)
        if (ex.contains(pt))
            return true;
    return false;
}

// Build the triangular net for one cross-section.
static Net build_net(const ExPolygons &xs, const Lattice &lat, coord_t spacing)
{
    Net net;
    if (xs.empty()) return net;

    Points                     pts;
    Triangulation::HalfEdges   cons;

    auto add_ring = [&](const Polygon &poly) {
        Points rp = resample_closed(poly, double(spacing));
        if (rp.size() < 3) return;
        const uint32_t base = uint32_t(pts.size());
        const uint32_t m    = uint32_t(rp.size());
        for (const Point &q : rp) { pts.push_back(q); net.verts.push_back({ q, false, 0 }); }
        for (uint32_t i = 0; i < m; ++i)
            cons.push_back({ base + ((i + m - 1) % m), base + i }); // prev -> i, following winding
    };
    for (const ExPolygon &ex : xs) {
        add_ring(ex.contour);
        for (const Polygon &h : ex.holes)
            add_ring(h);
    }
    if (pts.size() < 3 || cons.empty())
        return Net{}; // nothing printable

    // Interior points on the shared lattice, kept clear of the boundary.
    const ExPolygons inset = shrink_ex(xs, float(spacing) * 0.55f);
    if (! inset.empty()) {
        const BoundingBox bb = get_extents(inset);
        const int row_min = int(std::floor(bb.min.y() / lat.dy)) - 1;
        const int row_max = int(std::ceil (bb.max.y() / lat.dy)) + 1;
        const int col_min = int(std::floor(bb.min.x() / lat.s)) - 1;
        const int col_max = int(std::ceil (bb.max.x() / lat.s)) + 1;
        const double thr2 = (0.5 * spacing) * (0.5 * spacing);
        for (int row = row_min; row <= row_max; ++row)
            for (int col = col_min; col <= col_max; ++col) {
                const Point gp = lat.point(col, row);
                if (! point_inside(inset, gp)) continue;
                bool too_close = false;
                for (const Point &q : pts)
                    if ((q - gp).cast<double>().squaredNorm() < thr2) { too_close = true; break; }
                if (too_close) continue;
                net.grid[make_key(col, row)] = int(pts.size());
                net.verts.push_back({ gp, true, make_key(col, row) });
                pts.push_back(gp);
            }
    }

    std::sort(cons.begin(), cons.end());
    Triangulation::Indices tris = Triangulation::triangulate(pts, cons);

    std::set<std::pair<int,int>> eset;
    auto add_edge = [&](int u, int v) { if (u > v) std::swap(u, v); eset.emplace(u, v); };
    // If the constrained boundary self-intersects (thin features resampled coarsely), CGAL inserts
    // extra vertices whose index is not one of ours; drop any triangle referencing such a vertex.
    const int nverts = int(net.verts.size());
    auto valid = [nverts](const Vec3i &t) {
        return t[0] >= 0 && t[1] >= 0 && t[2] >= 0 && t[0] < nverts && t[1] < nverts && t[2] < nverts;
    };
    int kept = 0;
    for (const Vec3i &t : tris)
        if (valid(t)) { add_edge(t[0], t[1]); add_edge(t[1], t[2]); add_edge(t[2], t[0]); ++kept; }
    if (kept == 0) {
        // Degenerate / broken triangulation — fall back to the boundary loops themselves.
        for (const auto &e : cons) add_edge(int(e.first), int(e.second));
    }
    net.edges.assign(eset.begin(), eset.end());
    return net;
}

// G-code emission helper (feedrate tracking + travel/retract between struts).
struct EriGCode
{
    GCodeWriter                                   &writer;
    const std::function<void(const std::string&)> &emit;
    double                                         cur_speed = -1.;
    double                                         retract_threshold_mm2;

    EriGCode(GCodeWriter &w, double spacing_mm, const std::function<void(const std::string&)> &e)
        : writer(w), emit(e) { const double t = spacing_mm * 2.; retract_threshold_mm2 = t * t; }

    void set_speed(double mm_s)
    {
        if (mm_s <= 0.) mm_s = 1.;
        if (std::abs(mm_s - cur_speed) < 1e-6) return;
        cur_speed = mm_s;
        emit(writer.set_speed(mm_s * 60.));
    }

    void travel_to(const Vec3d &to)
    {
        if ((writer.get_position() - to).squaredNorm() >= retract_threshold_mm2)
            emit(writer.retract());
        set_speed(cur_speed > 0. ? cur_speed : 30.);
        emit(writer.travel_to_xyz(to));
        emit(writer.unretract());
    }

    void strut(const Vec3d &a, const Vec3d &b, double mm3_per_mm, double speed)
    {
        if ((writer.get_position() - a).squaredNorm() > 1e-6)
            travel_to(a);
        const double len = (b - a).norm();
        const double dE  = writer.extruder()->e_per_mm3() * mm3_per_mm * len;
        set_speed(speed);
        emit(writer.extrude_to_xyz(b, dE, "eridian"));
    }
};

} // namespace

void write_eridian_gcode(const std::vector<ExPolygons>                 &layers,
                         const std::vector<double>                     &z,
                         const EridianParams                           &params,
                         GCodeWriter                                   &writer,
                         const std::function<void(const std::string&)> &emit)
{
    if (layers.empty())
        return;

    const coord_t spacing = coord_t(scale_(params.lattice_spacing));
    const double  spacing_d = double(spacing);
    const Lattice lat{ spacing_d };

    // Build every layer's net first.
    std::vector<Net> nets(layers.size());
    for (size_t k = 0; k < layers.size(); ++k)
        nets[k] = build_net(layers[k], lat, spacing);

    EriGCode g(writer, params.lattice_spacing, emit);

    auto vertex3 = [](const Point &p, double zz) {
        return Vec3d(unscale<double>(p.x()), unscale<double>(p.y()), zz);
    };

    const double max_d2 = std::pow(3.0 * spacing, 2); // longest allowed pillar/brace, scaled²

    for (size_t k = 0; k < layers.size(); ++k) {
        const Net &here = nets[k];
        if (here.verts.empty())
            continue;

        // Connect this layer's vertices down to the layer below.
        if (k > 0) {
            const Net   &below = nets[k - 1];
            const double z_lo  = z[k - 1];
            const double z_hi  = z[k];
            for (const Vert &v : here.verts) {
                const Vert *low = nullptr;
                if (v.interior) {
                    auto it = below.grid.find(v.key);
                    if (it != below.grid.end())
                        low = &below.verts[it->second]; // shared lattice key -> vertical pillar
                }
                if (low == nullptr) {
                    // Nearest vertex below (leans into a brace under overhangs).
                    double best = max_d2;
                    for (const Vert &b : below.verts) {
                        double d = (b.p - v.p).cast<double>().squaredNorm();
                        if (d < best) { best = d; low = &b; }
                    }
                }
                if (low != nullptr)
                    g.strut(vertex3(low->p, z_lo), vertex3(v.p, z_hi), params.mm3_per_mm_pillar, params.speed_pillar);
            }
        }

        // Draw the flat triangular net at this layer's height.
        const double zz = z[k];
        for (const auto &e : here.edges)
            g.strut(vertex3(here.verts[e.first].p, zz), vertex3(here.verts[e.second].p, zz),
                    params.mm3_per_mm_flat, params.speed_flat);
    }

    emit(writer.retract());
}

} // namespace Slic3r
