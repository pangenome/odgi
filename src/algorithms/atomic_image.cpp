#include "atomic_image.hpp"

// routines for drawing raster images

// Xiaolin Wu's Line Algorithm

namespace odgi {

namespace algorithms {

color_t lighten(const color_t& c, const double& f) {
    color_t l;
    l.c = {
        (uint8_t)std::round((double)c.c.r + (double)(0xff - c.c.r) * (1.0 - f)),
        (uint8_t)std::round((double)c.c.g + (double)(0xff - c.c.g) * (1.0 - f)),
        (uint8_t)std::round((double)c.c.b + (double)(0xff - c.c.b) * (1.0 - f)),
        c.c.a
    };
    /*
    std::cerr << "lightening " << (int)c.c.r << " " << (int)c.c.g << " " << (int)c.c.b
              << " by " << f << " => "
              << (int)l.c.r << " " << (int)l.c.g << " " << (int)l.c.b << std::endl;
    */
    return l;
}

// helpers

double u_ipart(double x) { return std::floor(x); }

//double u_round(double x) { return u_ipart(x + 0.5); }
double u_round(double x) { return u_ipart(x + 0.5); }

double u_fpart(double x) { return x - std::floor(x); }

double u_rfpart(double x) { return 1.0 - u_fpart(x); }

void wu_draw_line(const bool steep, const double_t gradient, double intery,
                  const xy_d_t pxl1, const xy_d_t pxl2,
                  atomic_image_buf_t& image) {
    if (steep) {
        for (uint64_t i = pxl1.x + 1; i < pxl2.x; ++i) {
            image.set_pixel(u_ipart(intery), i, lighten(COLOR_BLACK, u_rfpart(intery)));
            image.set_pixel(u_ipart(intery) + 1, i, lighten(COLOR_BLACK, u_fpart(intery)));
            intery += gradient;
        }
    } else {
        for (uint64_t i = pxl1.x + 1; i < pxl2.x; ++i) {
            image.set_pixel(i, u_ipart(intery), lighten(COLOR_BLACK, u_rfpart(intery)));
            image.set_pixel(i, u_ipart(intery) + 1, lighten(COLOR_BLACK, u_fpart(intery)));
            intery += gradient;
        }
    }
}

xy_d_t wu_calc_endpoint(xy_d_t xy, const double_t gradient, const bool steep,
                        atomic_image_buf_t& image) {
    //std::cerr << "getting endpoint for (" << xy.x << "," << xy.y << ")" << std::endl;
    const xy_d_t end = {u_round(xy.x),
                        std::max(0.0, xy.y + gradient * (u_round(xy.x) - xy.x))};
    //std::cerr << "end is " << "(" << end.x << "," << end.y << ")" << std::endl;
    const double xgap = u_rfpart(xy.x + 0.5);
    const xy_d_t pxl = {end.x, u_ipart(end.y)};
    //std::cerr << "pxl is " << "(" << pxl.x << "," << pxl.y << ")" << std::endl;

    if (steep) {
        image.set_pixel(pxl.y, pxl.x, lighten(COLOR_BLACK, u_rfpart(end.y) * xgap));
        image.set_pixel(pxl.y + 1, pxl.x, lighten(COLOR_BLACK, u_fpart(end.y) * xgap));
    } else {
        image.set_pixel(pxl.x, pxl.y, lighten(COLOR_BLACK, u_rfpart(end.y) * xgap));
        image.set_pixel(pxl.x, pxl.y + 1, lighten(COLOR_BLACK, u_fpart(end.y) * xgap));
    }

    return pxl;
}

void wu_calc_line(xy_d_t xy0, xy_d_t xy1, atomic_image_buf_t& image) {

    /*
    std::cerr << "wu_calc_line "
              << "(" << xy0.x << "," << xy0.y << ")"
              << " --> "
              << "(" << xy1.x << "," << xy1.y << ")"
              << std::endl;
    */

    const bool steep = std::abs(xy1.y - xy0.y) > std::abs(xy1.x - xy0.x);
    double gradient;

    if (steep) {
        std::swap(xy0.x, xy0.y);
        std::swap(xy1.x, xy1.y);
    }
    if (xy0.x > xy1.x) {
        std::swap(xy0.x, xy1.x);
        std::swap(xy0.y, xy1.y);
    }

    {
        const xy_d_t d = {xy1.x - xy0.x, xy1.y - xy0.y};
        gradient = d.y / d.x;
        if (0.0 == d.x) {
            gradient = 1.0;
        }
        //std::cerr << "gradient " << gradient << std::endl;
    }

    wu_draw_line(steep, gradient,
                 (xy0.y + gradient * (u_round(xy0.x) - xy0.x)) + gradient,
                 wu_calc_endpoint(xy0, gradient, steep, image),
                 wu_calc_endpoint(xy1, gradient, steep, image), image);
}


}

}
