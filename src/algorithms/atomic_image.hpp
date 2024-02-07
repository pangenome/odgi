#pragma once

#include <fcntl.h>
#include <stdarg.h>
#include <stdbool.h>
#include <unistd.h>
#include <cstdint>
#include <vector>
#include <atomic>
#include <cmath>
#include <memory>
#include <iostream>
#include <iomanip>
#include "picosha2.h"

namespace odgi {

namespace algorithms {

struct xy_d_t {
    double x = 0;
    double y = 0;
    // project a point from a source 2D range into a target 2D range
    void into(const double& source_min_x,
              const double& source_min_y,
              const double& source_width,
              const double& source_height,
              const double& target_min_x,
              const double& target_min_y,
              const double& target_width,
              const double& target_height) {
        x = (x - source_min_x) * (target_width / source_width) + target_min_x;
        y = (y - source_min_y) * (target_height / source_height) + target_min_y;
    }
};

struct line_t {
    xy_d_t s = { 0, 0 };
    xy_d_t e = { 0, 0 };
    double d_x = 0;
    double d_y = 0;
    double a = 0;
    double b = 0;
    bool vert = false;
    line_t(const xy_d_t& _s, const xy_d_t& _e)
        : s(_s), e(_e) {
        d_y = (s.y - e.y); // delta y
        d_x = (s.x - e.x); // delta x
        a = d_y / d_x; // slope
        if (d_x == 0) {
            vert = true;
        } else {
            b = s.y - a * s.x; // use one of the points to solve for the intercept
        }
    }
    // if the point occurs above the line
    bool gt(const xy_d_t& p) {
        if (vert) {
            return p.x < s.x;
        } else {
            return p.y < a * p.x + b;
        }
    }
    // if the point occurs on or below the line
    bool lte(const xy_d_t& p) {
        return !gt(p);
    }
};

struct xy_u_t {
    uint64_t x = 0;
    uint64_t y = 0;
};

typedef struct s_rgb {
	uint8_t	r;
	uint8_t	g;
	uint8_t	b;
    uint8_t a;
} RGB;

typedef union rgb_t {
    uint32_t hex;
	RGB c;
} color_t;

std::ostream& operator<<(std::ostream& out, const color_t& c);
bool operator==(const color_t& a, const color_t& b);

color_t hash_color(const std::string& s);
color_t lighten(const color_t& c, const double& f);
color_t brighten(const color_t& a, const color_t& b, const double& f);
color_t darkest(const color_t& a, const color_t& b);
color_t layer(const color_t& a, const color_t& b, const double& f);
color_t mix(const color_t& a, const color_t& b, const double& f);

std::string to_hex(const color_t& c);
std::string to_rgba(const color_t& c);
std::string to_hexrgb(const color_t& c);

const color_t COLOR_BLACK = { 0xff000000 };
const color_t COLOR_LIGHTGRAY = { 0xffD3D3D3 };
const color_t COLOR_WHITE = { 0xffffffff };

# define COLOR_MAX  255
# define COLOR_MIN  0

/*
# define INRANGE(c)	((COLOR_MAX < (c)) \
	? COLOR_MAX : ((COLOR_MIN > (c)) \
		? COLOR_MIN : (c)))
*/

struct atomic_image_buf_t {
    std::unique_ptr<std::vector<std::atomic<uint32_t>>> image;
    uint64_t height = 0;
    uint64_t width = 0;
    double source_width = 0;
    double source_height = 0;
    double source_per_px_x = 0;
    double source_per_px_y = 0;
    double source_min_x = 0;
    double source_min_y = 0;
    atomic_image_buf_t(const uint64_t& w,
                       const uint64_t& h,
                       const double& s_w,
                       const double& s_h,
                       const double& s_m_x,
                       const double& s_m_y)
        : width(w)
        , height(h)
        , source_width(s_w)
        , source_height(s_h)
        , source_min_x(s_m_x)
        , source_min_y(s_m_y)
        {
        //std::cerr << "width x height " << w << "x" << h << std::endl;
        image = std::make_unique<std::vector<std::atomic<uint32_t>>>(height * width);
        for (uint64_t i = 0; i < image->size(); ++i) {
            //std::cerr << "coloring " << COLOR_WHITE.hex << std::endl;
            (*image)[i] = COLOR_WHITE.hex; // atomic assignment
        }
        source_per_px_x = source_width / width;
        source_per_px_y = source_height / height;
    }
    std::vector<uint8_t> to_bytes() {
        std::vector<uint8_t> bytes(4 * height * width);
        for (uint64_t i = 0; i < image->size(); ++i) {
            color_t c = {(*image)[i].load()};
            uint64_t j = i * 4;
            bytes[j    ] = c.c.r;
            bytes[j + 1] = c.c.g;
            bytes[j + 2] = c.c.b;
            bytes[j + 3] = c.c.a;
        }
        return bytes;
    }
    // ablative
    void set_pixel(const uint64_t& x,
                   const uint64_t& y,
                   const color_t& c) {
        //size_t i = width * y + x;
        //std::cerr << "setting color with intensity " << f << std::endl;
        (*image)[width * y + x] = c.hex;
    }
    // layering
    void layer_pixel(const uint64_t& x,
                     const uint64_t& y,
                     const color_t& c) {
        if (x >= width || y >= height) {
            return; // bail out
        }
        size_t i = width * y + x;
        //std::cerr << "getting i=" << i << " " << y << " " << x << " " << " in image " << height << "x" << width << std::endl;
        color_t v;
        v.hex = (*image)[i].load();
        //std::cerr << "got " << v.hex << " " << (int)v.c.r << "," << (int)v.c.g << "," << (int)v.c.b << std::endl;
        //std::cerr << "layer " << c.hex << " " << (int)c.c.r << "," << (int)c.c.g << "," << (int)c.c.b << std::endl;
        v = mix(c, v, 0.5);
        //v = layer(c, v, f);
        //std::cerr << "assigned " << v.hex << " " << (int)v.c.r << "," << (int)v.c.g << "," << (int)v.c.b << std::endl;
        // there is the possibility of a race decreasing how bright pixels get
        // but, at least we won't overflow
        (*image)[i] = v.hex; // atomic assignment
    }
};


double u_ipart(double x);

double u_round(double x);

double u_fpart(double x);

double u_rfpart(double x);

void wu_draw_line(const bool steep, const double_t gradient, double_t intery,
                  const xy_d_t pxl1, const xy_d_t pxl2,
                  const color_t& color,
                  atomic_image_buf_t& image,
                  bool top, bool bottom);

xy_d_t wu_calc_endpoint(xy_d_t xy, const double_t gradient, const bool steep,
                        const color_t& color,
                        atomic_image_buf_t& image);

void wu_calc_line(xy_d_t xy0, xy_d_t xy1,
                  const color_t& color,
                  atomic_image_buf_t& image,
                  bool top, bool bottom);

void wu_rekt(xy_d_t xy0, xy_d_t xy1,
             xy_d_t xy2, xy_d_t xy3,
             const color_t& color,
             atomic_image_buf_t& image);

void wu_calc_wide_line(xy_d_t xy0, xy_d_t xy1,
                       const color_t& color,
                       atomic_image_buf_t& image,
                       const double& width = 0);

void wu_calc_rainbow(xy_d_t xy0, xy_d_t xy1, atomic_image_buf_t& image,
                     const std::vector<color_t>& colors,
                     const double& spacing = 1,
                     const double& width = 0);

}

}
