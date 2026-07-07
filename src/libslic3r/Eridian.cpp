#include "Eridian.hpp"

#include <array>
#include <cmath>
#include <cstdint>
#include <unordered_set>

#include "BoundingBox.hpp"
#include "GCode/GCodeWriter.hpp"

// See Eridian.hpp. This builds a shared triangular lattice, clips it to every cross-section, and
// emits flat nets + vertical pillars + outward overhang braces as a 3D truss.

namespace Slic3r {

namespace {

// A lattice vertex identified by integer axial indices (col, row). Encoded into a single 64-bit
// key for hashing. The same (col, row) maps to the same (x, y) on every layer, which is what makes
// pillars perfectly vertical.
using Key = int64_t;
static inline Key make_key(int col, int row)
{
    return (Key(int32_t(col)) << 32) | Key(uint32_t(int32_t(row)));
}

struct Lattice
{
    double s;   // spacing, scaled
    double dy;  // row pitch, scaled = s * sqrt(3)/2

    explicit Lattice(double spacing_scaled) : s(spacing_scaled), dy(spacing_scaled * std::sqrt(3.) / 2.) {}

    Point point(int col, int row) const
    {
        const double x = col * s + ((row & 1) ? s * 0.5 : 0.0);
        const double y = row * dy;
        return Point(coord_t(std::llround(x)), coord_t(std::llround(y)));
    }

    // The six in-plane lattice neighbours of (col, row).
    std::array<std::pair<int,int>, 6> neighbours(int col, int row) const
    {
        if (row & 1)
            return {{ {col + 1, row}, {col - 1, row},
                      {col, row + 1}, {col + 1, row + 1},
                      {col, row - 1}, {col + 1, row - 1} }};
        else
            return {{ {col + 1, row}, {col - 1, row},
                      {col - 1, row + 1}, {col, row + 1},
                      {col - 1, row - 1}, {col, row - 1} }};
    }
};

struct LayerVerts
{
    std::unordered_set<Key> active;
    std::vector<std::pair<int,int>> list;
    bool has(int col, int row) const { return active.count(make_key(col, row)) != 0; }
};

// G-code emission helper: tracks feedrate and inserts travels/retracts between disjoint segments.
struct EriGCode
{
    GCodeWriter                                   &writer;
    const EridianParams                           &p;
    const std::function<void(const std::string&)> &emit;
    double                                         cur_speed = -1.;
    double                                         retract_threshold_mm2;

    EriGCode(GCodeWriter &w, const EridianParams &params, const std::function<void(const std::string&)> &e)
        : writer(w), p(params), emit(e)
    {
        const double t = params.lattice_spacing * 2.;
        retract_threshold_mm2 = t * t;
    }

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
        set_speed(p.speed_flat);
        emit(writer.travel_to_xyz(to));
        emit(writer.unretract());
    }

    // Extrude a single strut from `a` to `b`, travelling to `a` first if needed.
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

// Is the midpoint of two lattice vertices inside the cross-section? Guards against lattice edges
// spanning across concavities or holes.
static bool midpoint_inside(const ExPolygons &layer, const Point &a, const Point &b)
{
    const Point mid((a.x() + b.x()) / 2, (a.y() + b.y()) / 2);
    for (const ExPolygon &ex : layer)
        if (ex.contains(mid))
            return true;
    return false;
}

static bool point_inside(const ExPolygons &layer, const Point &pt)
{
    for (const ExPolygon &ex : layer)
        if (ex.contains(pt))
            return true;
    return false;
}

} // namespace

void write_eridian_gcode(const std::vector<ExPolygons>                 &layers,
                         const std::vector<double>                     &z,
                         const EridianParams                           &params,
                         GCodeWriter                                   &writer,
                         const std::function<void(const std::string&)> &emit)
{
    if (layers.empty())
        return;

    const Lattice lat(scale_(params.lattice_spacing));

    // Overall lattice index range from the union bounding box of all cross-sections.
    BoundingBox bbox;
    bool have_bbox = false;
    for (const ExPolygons &layer : layers) {
        if (layer.empty()) continue;
        BoundingBox b = get_extents(layer);
        if (! have_bbox) { bbox = b; have_bbox = true; }
        else             { bbox.merge(b); }
    }
    if (! have_bbox)
        return;

    const int row_min = int(std::floor(bbox.min.y() / lat.dy)) - 1;
    const int row_max = int(std::ceil (bbox.max.y() / lat.dy)) + 1;
    const int col_min = int(std::floor(bbox.min.x() / lat.s)) - 1;
    const int col_max = int(std::ceil (bbox.max.x() / lat.s)) + 1;

    // Active lattice vertices per layer.
    std::vector<LayerVerts> verts(layers.size());
    for (size_t k = 0; k < layers.size(); ++k) {
        for (int row = row_min; row <= row_max; ++row)
            for (int col = col_min; col <= col_max; ++col) {
                const Point pt = lat.point(col, row);
                if (point_inside(layers[k], pt)) {
                    verts[k].active.insert(make_key(col, row));
                    verts[k].list.emplace_back(col, row);
                }
            }
    }

    EriGCode g(writer, params, emit);

    auto vertex3 = [&](int col, int row, double zz) {
        const Point pt = lat.point(col, row);
        return Vec3d(unscale<double>(pt.x()), unscale<double>(pt.y()), zz);
    };

    // Draw the flat triangular net of layer k at height z[k]. Each forward edge (E, NE, NW) is
    // emitted once, when both endpoints are active and the edge stays inside the cross-section.
    auto draw_flat_net = [&](size_t k) {
        const LayerVerts &lv = verts[k];
        const double      zz = z[k];
        for (const auto &v : lv.list) {
            const int c = v.first, r = v.second;
            const Point a = lat.point(c, r);
            const std::pair<int,int> fwd[3] = {
                { c + 1, r },
                (r & 1) ? std::make_pair(c + 1, r + 1) : std::make_pair(c, r + 1), // NE
                (r & 1) ? std::make_pair(c, r + 1)     : std::make_pair(c - 1, r + 1) // NW
            };
            for (const auto &nb : fwd) {
                if (! lv.has(nb.first, nb.second)) continue;
                const Point b = lat.point(nb.first, nb.second);
                if (! midpoint_inside(layers[k], a, b)) continue;
                g.strut(vertex3(c, r, zz), vertex3(nb.first, nb.second, zz), params.mm3_per_mm_flat, params.speed_flat);
            }
        }
    };

    for (size_t k = 0; k < layers.size(); ++k) {
        if (k > 0) {
            const LayerVerts &below = verts[k - 1];
            const LayerVerts &here  = verts[k];
            const double      z_lo  = z[k - 1];
            const double      z_hi  = z[k];

            // Vertical pillars: vertices present on both layers.
            for (const auto &v : here.list) {
                if (below.has(v.first, v.second))
                    g.strut(vertex3(v.first, v.second, z_lo), vertex3(v.first, v.second, z_hi),
                            params.mm3_per_mm_pillar, params.speed_pillar);
            }
            // Overhang braces: a vertex with no pillar beneath it is tied down to its active
            // lattice neighbours on the layer below, forming outward-leaning triangles.
            for (const auto &v : here.list) {
                if (below.has(v.first, v.second)) continue;
                for (const auto &nb : lat.neighbours(v.first, v.second)) {
                    if (below.has(nb.first, nb.second))
                        g.strut(vertex3(nb.first, nb.second, z_lo), vertex3(v.first, v.second, z_hi),
                                params.mm3_per_mm_pillar, params.speed_pillar);
                }
            }
        }
        draw_flat_net(k);
    }

    emit(writer.retract());
}

} // namespace Slic3r
