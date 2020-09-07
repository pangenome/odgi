#include "atomic_image.hpp"

// routines for drawing raster images

// Xiaolin Wu's Line Algorithm

namespace odgi {

namespace algorithms {

color_t lighten(const color_t& c, const double& f) {
    color_t l;
    l.c = {
        (uint8_t)std::round((double)c.c.r + (double)(0xff - c.c.r) * f), //(1.0 - f)),
        (uint8_t)std::round((double)c.c.g + (double)(0xff - c.c.g) * f), //(1.0 - f)),
        (uint8_t)std::round((double)c.c.b + (double)(0xff - c.c.b) * f), //(1.0 - f)),
        c.c.a
    };
    /*
    std::cerr << "lightening " << (int)c.c.r << " " << (int)c.c.g << " " << (int)c.c.b
              << " by " << f << " => "
              << (int)l.c.r << " " << (int)l.c.g << " " << (int)l.c.b << std::endl;
    */
    return l;
}

color_t darken(const color_t& c, const double& f) {
    color_t l;
    l.c = {
        (uint8_t)std::round((double)c.c.r - (double)(c.c.r) * f),
        (uint8_t)std::round((double)c.c.g - (double)(c.c.g) * f),
        (uint8_t)std::round((double)c.c.b - (double)(c.c.b) * f),
        c.c.a
    };
    /*
    std::cerr << "lightening " << (int)c.c.r << " " << (int)c.c.g << " " << (int)c.c.b
              << " by " << f << " => "
              << (int)l.c.r << " " << (int)l.c.g << " " << (int)l.c.b << std::endl;
    */
    return l;
}

// layer a onto b with f intensity
color_t brighten(const color_t& a, const color_t& b, const double& f) {
	color_t out = a;

    auto inrange = [](double x) {
                       if (COLOR_MAX < x) {
                           return (uint8_t) COLOR_MAX;
                       } else if (COLOR_MIN > x) {
                           return (uint8_t) COLOR_MIN;
                       } else {
                           return(uint8_t) std::round(x);
                       }
                   };

    //std::cerr << "got a = " << a.hex << " " << (int)a.c.r << "," << (int)a.c.g << "," << (int)a.c.b << std::endl;
    //std::cerr << "got b = " << b.hex << " " << (int)b.c.r << "," << (int)b.c.g << "," << (int)b.c.b << std::endl;
    //std::cerr << "f = " << f << std::endl;
	if (.0d < f) {
		out.c = (RGB) { inrange((double)a.c.r + ((double)b.c.r * (1.0d - f))),
                        inrange((double)a.c.g + ((double)b.c.g * (1.0d - f))),
                        inrange((double)a.c.b + ((double)b.c.b * (1.0d - f))),
                        255 };
    }
    //std::cerr << "return out = " << out.hex << " " << (int)out.c.r << "," << (int)out.c.g << "," << (int)out.c.b << std::endl;
	return out;
}

color_t layer(const color_t& a, const color_t& b) {
    color_t l;
    // take the darker color
    //l.hex = b.hex;
    l.c.r = std::min(a.c.r, b.c.r);
    l.c.g = std::min(a.c.g, b.c.g);
    l.c.b = std::min(a.c.b, b.c.b);
    //l.c.r = (a.c.r > b.c.r ? b.c.r : a.c.r);
    //l.c.g = (a.c.g > b.c.g ? b.c.g : a.c.g);
    //l.c.b = (a.c.b > b.c.b ? b.c.b : a.c.b);
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
            image.layer_pixel(u_ipart(intery), i, COLOR_BLACK, u_rfpart(intery));
            image.layer_pixel(u_ipart(intery) + 1, i, COLOR_BLACK, u_fpart(intery));
            intery += gradient;
        }
    } else {
        for (uint64_t i = pxl1.x + 1; i < pxl2.x; ++i) {
            image.layer_pixel(i, u_ipart(intery), COLOR_BLACK, u_rfpart(intery));
            image.layer_pixel(i, u_ipart(intery) + 1, COLOR_BLACK, u_fpart(intery));
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
        image.layer_pixel(pxl.y, pxl.x, COLOR_BLACK, u_rfpart(end.y) * xgap);
        image.layer_pixel(pxl.y + 1, pxl.x, COLOR_BLACK, u_fpart(end.y) * xgap);
    } else {
        image.layer_pixel(pxl.x, pxl.y, COLOR_BLACK, u_rfpart(end.y) * xgap);
        image.layer_pixel(pxl.x, pxl.y + 1, COLOR_BLACK, u_fpart(end.y) * xgap);
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
        //wu_draw_line
}

void wu_calc_multiline(xy_d_t xy0, xy_d_t xy1, atomic_image_buf_t& image,
                       const double& width, const double& overlay) {

    if (width == 0) {
        wu_calc_line(xy0, xy1, image);
    }

    const xy_d_t d = { xy1.x - xy0.x, xy1.y - xy0.y };
    double gradient = d.y / d.x;
    if (d.x == 0.0) {
        gradient = d.y > 0 ? 1.0 : -1.0;
    }
    double inv_gradient = d.x / d.y;
    if (d.y == 0.0) {
        inv_gradient = d.x > 0 ? 1.0 : -1.0;
    }

    // width is given in bp space (units in the base layout)
    // we will generate a series of lines parallel to the center line
    // to simulate a line with the given width
    double width_in_px = width / image.source_per_px_y;

    xy_d_t xyA = xy0;
    xy_d_t xyB = xy1;
    // how for to get to a Y such that the length is w/2
    double move_in_x = width_in_px / ( 2 * sqrt(pow(inv_gradient, 2) + 1));
    double move_in_y = sqrt(pow(width_in_px/2,2) - pow(move_in_x, 2));

    if (gradient == 0.0) {
        move_in_x = 0.0;
        move_in_y = width_in_px / 2.0;
    } else if (inv_gradient == 0.0) {
        move_in_x = width_in_px / 2.0;
        move_in_y = 0.0;
    }

    // adjust the moves to reflect our line
    if (xyA.x > xyB.x && xyA.y > xyB.y) {
        move_in_x = -move_in_x;
        move_in_y = move_in_y;
    } else if (xyA.x > xyB.x && xyA.y <= xyB.y) {
        move_in_x = -move_in_x;
        move_in_y = -move_in_y;
    } else if (xyA.x <= xyB.x && xyA.y > xyB.y) {
        move_in_x = move_in_x;
        move_in_y = move_in_y;
    } else { //if //(xyA.x < xyB.x && xyA.y < xyB.y) {
        move_in_x = move_in_x;
        move_in_y = -move_in_y;
    }

    xyA.x -= move_in_x;
    xyA.y -= move_in_y;
    xyB.x -= move_in_x;
    xyB.y -= move_in_y;

    double pix = image.source_per_px_x / overlay;
    // make sure we always make at least two steps
    uint64_t total_steps = std::ceil(std::max((double)2, width_in_px / pix));
    double step_x = 2*move_in_x / (double)total_steps;
    double step_y = 2*move_in_y / (double)total_steps;

    for (uint64_t i = 0; i < total_steps; ++i) {
        xyA.x += step_x;
        xyA.y += step_y;
        xyB.x += step_x;
        xyB.y += step_y;
        wu_calc_line(xyA, xyB, image);
    }
}

}

}
