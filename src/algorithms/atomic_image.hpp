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
              const double& target_width,
              const double& target_height) {
        x = x * (target_width / source_width) - source_min_x;
        y = y * (target_height / source_height) - source_min_y;
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
} RGB;

typedef union rgb_t {
    uint32_t hex;
	RGB c;
} color_t;

const color_t COLOR_BLACK = { 0 };
const color_t COLOR_WHITE = { 0xffffff };

struct atomic_image_buf_t {
    std::unique_ptr<std::vector<std::atomic<uint32_t>>> image;
    uint64_t height = 0;
    uint64_t width = 0;
    atomic_image_buf_t(const uint64_t& h,
                       const uint64_t& w)
        : height(h)
        , width(w) {
        image = std::make_unique<std::vector<std::atomic<uint32_t>>>(height * width);
    }
    // ablative
    void set_pixel(const uint64_t& x,
                   const uint64_t& y,
                   const color_t& c, // color
                   const double& f) { // intensity
        size_t i = width * y + x;
        color_t v;
        v.c = {
            (uint8_t)std::round((double)c.c.r * f),
            (uint8_t)std::round((double)c.c.g * f),
            (uint8_t)std::round((double)c.c.b * f)
        };
        (*image)[i] = v.hex;
    }
    // layering
    void layer_pixel(const uint64_t& x,
                     const uint64_t& y,
                     const color_t& c,
                     const double& f) {
        size_t i = width * y + x;
        color_t v = {(*image)[i].load()};
        // update until we saturate
        v.c.r += std::round((double)(255 - v.c.r) / 255.0 * std::round(f * c.c.r));
        v.c.g += std::round((double)(255 - v.c.g) / 255.0 * std::round(f * c.c.g));
        v.c.b += std::round((double)(255 - v.c.b) / 255.0 * std::round(f * c.c.b));
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
                  atomic_image_buf_t& image);

xy_d_t wu_calc_endpoint(xy_d_t xy, const double_t gradient, const bool steep,
                        atomic_image_buf_t& image);

void wu_calc_line(xy_d_t xy0, xy_d_t xy1, atomic_image_buf_t& image);

}

}
