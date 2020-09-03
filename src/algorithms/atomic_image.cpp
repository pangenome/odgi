#include "atomic_image.hpp"

// routines for drawing raster images

// Xiaolin Wu's Line Algorithm

namespace odgi {

namespace algorithms {

// helpers

double u_ipart(double x) { return std::floor(x); }

double u_round(double x) { return u_ipart(x + 0.5); }

double u_fpart(double x) { return x - std::floor(x); }

double u_rfpart(double x) { return 1.0 - u_fpart(x); }

void wu_draw_line(const bool steep, const double_t gradient, double intery,
                  const xy_d_t pxl1, const xy_d_t pxl2,
                  atomic_image_buf_t& image) {
    for (unsigned long long i = pxl1.x + 1; i < pxl2.x - 1; i++) {
        if (steep) {
            //putpxl_plot(u_ipart(intery), i, u_rfpart(intery), pxls);
            image.layer_pixel(u_ipart(intery), i, COLOR_BLACK, u_rfpart(intery));
            //putpxl_plot(u_ipart(intery) + 1, i, u_fpart(intery), pxls);
            image.layer_pixel(u_ipart(intery) + 1, i, COLOR_BLACK, u_fpart(intery));
        } else {
            //putpxl_plot(i, u_ipart(intery), u_rfpart(intery), pxls);
            image.layer_pixel(i, u_ipart(intery), COLOR_BLACK, u_rfpart(intery));
            //putpxl_plot(i, u_ipart(intery) + 1, u_fpart(intery), pxls);
            image.layer_pixel(i, u_ipart(intery) + 1, COLOR_BLACK, u_rfpart(intery));
        }
        intery += gradient;
    }
}

xy_d_t wu_calc_endpoint(xy_d_t xy, const double_t gradient, const bool steep,
                        atomic_image_buf_t& image) {
    const xy_d_t end = {u_round(xy.x),
                        xy.y + gradient * (u_round(xy.x) - xy.x)};
    const double xgap = u_rfpart(xy.x + 0.5);
    const xy_d_t pxl = {end.x, u_ipart(end.y)};

    if (steep) {
        //putpxl_plot(pxl[1], pxl[0], u_rfpart(end[1]) * xgap, pxls);
        image.layer_pixel(pxl.y, pxl.x, COLOR_BLACK, u_rfpart(end.y) * xgap);
        //putpxl_plot(pxl[1] + 1, pxl[0], u_fpart(end[1]) * xgap, pxls);
        image.layer_pixel(pxl.y + 1, pxl.x, COLOR_BLACK, u_fpart(end.y) * xgap);
    } else {
        //putpxl_plot(pxl[0], pxl[1], u_rfpart(end[1]) * xgap, pxls);
        image.layer_pixel(pxl.x, pxl.y, COLOR_BLACK, u_rfpart(end.y) * xgap);
        //putpxl_plot(pxl[0], pxl[1] + 1, u_fpart(end[1]) * xgap, pxls);
        image.layer_pixel(pxl.x, pxl.y + 1, COLOR_BLACK, u_fpart(end.y) * xgap);
    }
    return pxl;
}

void wu_calc_line(xy_d_t xy0, xy_d_t xy1, atomic_image_buf_t& image) {

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
    }

    wu_draw_line(steep, gradient,
                 (xy0.y + gradient * (u_round(xy0.x) - xy0.x)) + gradient,
                 wu_calc_endpoint(xy0, gradient, steep, image),
                 wu_calc_endpoint(xy1, gradient, steep, image), image);
}


}

}
