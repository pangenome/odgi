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

color_t lighten(const color_t& c, const double& f);

const color_t COLOR_BLACK = { 0xff000000 };
const color_t COLOR_WHITE = { 0xffffffff };

struct atomic_image_buf_t {
    std::unique_ptr<std::vector<std::atomic<uint32_t>>> image;
    uint64_t height = 0;
    uint64_t width = 0;
    atomic_image_buf_t(const uint64_t& w,
                       const uint64_t& h)
        : width(w)
        , height(h) {
        //std::cerr << "width x height " << w << "x" << h << std::endl;
        image = std::make_unique<std::vector<std::atomic<uint32_t>>>(height * width);
        for (uint64_t i = 0; i < image->size(); ++i) {
            //std::cerr << "coloring " << COLOR_WHITE.hex << std::endl;
            (*image)[i] = COLOR_WHITE.hex; // atomic assignment
        }
    }
    std::vector<uint8_t> to_bytes(void) {
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
                     const color_t& c,
                     const double& f) {
        size_t i = width * y + x;
        //std::cerr << "getting i=" << i << " in image " << height << "x" << width << std::endl;
        color_t v = {(*image)[i].load()};
        // update until we saturate
        v.c.r += std::round((double)(255 - v.c.r) / 255.0 * std::round(f * c.c.r));
        v.c.g += std::round((double)(255 - v.c.g) / 255.0 * std::round(f * c.c.g));
        v.c.b += std::round((double)(255 - v.c.b) / 255.0 * std::round(f * c.c.b));
        //v.c.a += std::round((double)(255 - v.c.a) / 255.0 * std::round(f * c.c.a));
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
