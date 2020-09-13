#include "atomic_image.hpp"

// routines for drawing raster images

// Xiaolin Wu's Line Algorithm

namespace odgi {

namespace algorithms {

color_t hash_color(const std::string& s) {
    // use a sha256 to get a few bytes that we'll use for a color
    picosha2::byte_t hashed[picosha2::k_digest_size];
    picosha2::hash256(s.begin(), s.end(), hashed, hashed + picosha2::k_digest_size);
    uint8_t path_r = hashed[24];
    uint8_t path_g = hashed[8];
    uint8_t path_b = hashed[16];
    double path_r_f = (double) path_r / (double) (std::numeric_limits<uint8_t>::max());
    double path_g_f = (double) path_g / (double) (std::numeric_limits<uint8_t>::max());
    double path_b_f = (double) path_b / (double) (std::numeric_limits<uint8_t>::max());
    double sum = path_r_f + path_g_f + path_b_f;
    path_r_f /= sum;
    path_g_f /= sum;
    path_b_f /= sum;
    double f = std::min(1.5, 1.0 / std::max(std::max(path_r_f, path_g_f), path_b_f));
    color_t c;
    //c.c.r = (uint8_t) std::round(255 * std::min(path_r_f * f, (double) 1.0));
    //c.c.g = (uint8_t) std::round(255 * std::min(path_g_f * f, (double) 1.0));
    //c.c.b = (uint8_t) std::round(255 * std::min(path_b_f * f, (double) 1.0));
    c.c.r = (uint8_t) std::round(255 * std::min(path_r_f, (double) 1.0));
    c.c.g = (uint8_t) std::round(255 * std::min(path_g_f, (double) 1.0));
    c.c.b = (uint8_t) std::round(255 * std::min(path_b_f, (double) 1.0));
    c.c.a = 255;
    return c;
}

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
color_t layer(const color_t& a, const color_t& b, const double& f) {
    return a;
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
		out.c = (RGB) { inrange((double)b.c.r - ((double)(255 * (1-f) - a.c.r))),
                        inrange((double)b.c.g - ((double)(255 * (1-f) - a.c.g))),
                        inrange((double)b.c.b - ((double)(255 * (1-f) - a.c.b))),
                        255 };
        /*
        out.c = (RGB) { inrange(a.c.r + ((double)255 * (1.0d - f))),
                        inrange(a.c.g + ((double)255 * (1.0d - f))),
                        inrange(a.c.b + ((double)255 * (1.0d - f))),
                        255 };
        */
    }
    //std::cerr << "return out = " << out.hex << " " << (int)out.c.r << "," << (int)out.c.g << "," << (int)out.c.b << std::endl;
	return out;
}

// helpers

double u_ipart(double x) { return std::floor(x); }

double u_round(double x) { return u_ipart(x + 0.5d); }

double u_fpart(double x) { return x - std::floor(x); }

double u_rfpart(double x) { return 1.0d - u_fpart(x); }

void wu_draw_line(const bool steep, const double_t gradient, double intery,
                  const xy_d_t pxl1, const xy_d_t pxl2, const color_t& color,
                  atomic_image_buf_t& image, bool top, bool bottom) {
    if (steep) {
        for (double i = pxl1.x + 1; i < pxl2.x; ++i) {
            image.layer_pixel(u_ipart(intery), i, color, (!bottom ? 1.0 : 1.0-u_rfpart(intery)));
            image.layer_pixel(u_ipart(intery) + 1, i, color, (!top ? 1.0 : 1.0-u_fpart(intery)));
            intery += gradient;
        }
    } else {
        for (double i = pxl1.x + 1; i < pxl2.x; ++i) {
            image.layer_pixel(i, u_ipart(intery), color, (!bottom ? 1.0 : 1.0-u_rfpart(intery)));
            image.layer_pixel(i, u_ipart(intery) + 1, color, (!top ? 1.0 : 1.0-u_fpart(intery)));
            intery += gradient;
        }
    }
}

xy_d_t wu_calc_endpoint(xy_d_t xy, const double_t gradient, const bool steep,
                        const color_t& color,
                        atomic_image_buf_t& image) {
    //std::cerr << "getting endpoint for (" << xy.x << "," << xy.y << ")" << std::endl;
    const xy_d_t end = {u_round(xy.x),
                        std::max(0.0, xy.y + gradient * (u_round(xy.x) - xy.x))};
    //std::cerr << "end is " << "(" << end.x << "," << end.y << ")" << std::endl;
    const double xgap = u_rfpart(xy.x + 0.5d);
    const xy_d_t pxl = {end.x, u_ipart(end.y)};
    //std::cerr << "pxl is " << "(" << pxl.x << "," << pxl.y << ")" << std::endl;

    if (steep) {
        image.layer_pixel(pxl.y, pxl.x, color, (1.0 - u_rfpart(end.y)) * xgap);
        image.layer_pixel(pxl.y + 1, pxl.x, color, (1.0 - u_fpart(end.y)) * xgap);
    } else {
        image.layer_pixel(pxl.x, pxl.y, color, (1.0 - u_rfpart(end.y)) * xgap);
        image.layer_pixel(pxl.x, pxl.y + 1, color, (1.0 - u_fpart(end.y)) * xgap);
    }

    return pxl;
}

void wu_calc_line(xy_d_t xy0, xy_d_t xy1,
                  const color_t& color,
                  atomic_image_buf_t& image,
                  bool top, bool bottom) {

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
        //std::swap(top, bottom);
    }
    if (xy0.x > xy1.x) {
        std::swap(xy0.x, xy1.x);
        std::swap(xy0.y, xy1.y);
        std::swap(top, bottom);
    }

    {
        const xy_d_t d = {xy1.x - xy0.x, xy1.y - xy0.y};
        gradient = d.y / d.x;
        if (0.0 == d.x) {
            gradient = 1.0;
        }
        //std::cerr << "gradient " << gradient << std::endl;
    }

    //if (steep) {
        //std::swap(top, bottom);
//}

    wu_draw_line(steep, gradient,
                 (xy0.y + gradient * (u_round(xy0.x) - xy0.x)) + gradient,
                 wu_calc_endpoint(xy0, gradient, steep, color, image),
                 wu_calc_endpoint(xy1, gradient, steep, color, image),
                 color,
                 image,
                 top, bottom);
}

void wu_rekt(xy_d_t xy0, xy_d_t xy1,
             xy_d_t xy2, xy_d_t xy3,
             const color_t& color,
             atomic_image_buf_t& image) {

    //color_t red   = { 0xff0000ff };
    //color_t green = { 0xff00ff00 };
    //color_t blue  = { 0xffff0000 };
    //color_t purp  = { 0xffff00ff };
    
    // top line
    wu_calc_line(xy0, xy1,
                 color,
                 image,
                 true, false);

    // right line
    wu_calc_line(xy1, xy3,
                 color,
                 image,
                 false, true);

    // bottom line
    wu_calc_line(xy2, xy3,
                 color,
                 image,
                 false, true);

    // left line
    wu_calc_line(xy2, xy0,
                 color,
                 image,
                 false, true);

    // now fill with a solid color
    // follow the top line
    // moving down with the gradient given by the left line

    //color_t black  = { 0xff000000 };

    // find the extents of the rectangle
    auto xs = {xy0.x, xy1.x, xy2.x, xy3.x};
    auto mx = std::minmax_element(xs.begin(), xs.end());
    double min_x = *mx.first;
    double max_x = *mx.second;
    auto ys = {xy0.y, xy1.y, xy2.y, xy3.y};
    auto my = std::minmax_element(ys.begin(), ys.end());
    double min_y = *my.first;
    double max_y = *my.second;

    // find the min y point
    // and the next two highest
    std::vector<xy_d_t> ps = { xy0, xy1, xy2, xy3 };
    std::sort(ps.begin(), ps.end(), [](const xy_d_t& a,
                                       const xy_d_t& b) {
                                        return a.y < b.y;
                                    });

    // arranged so that l0 and l1 are < l2 and l3
    line_t l0(ps[0], ps[1]);
    line_t l1(ps[0], ps[2]);
    line_t l2(ps[1], ps[3]);
    line_t l3(ps[2], ps[3]);

    auto inside =
        [&](const xy_d_t& p) {
            return l0.lte(p) && l1.lte(p) && l2.gt(p) && l3.gt(p);
        };
    
    xy_d_t l = { u_ipart(min_x), u_ipart(min_y) };
    xy_d_t h = { u_ipart(max_x), u_ipart(max_y) };
    for (double i = l.y; i < h.y; ++i) {
        for (double j = l.x; j < h.x; ++j) {
            // draw if it's in bounds
            if (inside({j, i})) {
                image.layer_pixel(j, i, color, 0.0);
            }
        }
    }
    
}

void wu_calc_wide_line(xy_d_t xy0, xy_d_t xy1,
                       const color_t& color,
                       atomic_image_buf_t& image,
                       const double& width) {

    //if (width == 0) {
    //wu_calc_line(xy0, xy1, color, image, true, true);
//}

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

    xy_d_t xyA = xy0;
    xy_d_t xyB = xy1;
    xy_d_t xyC = xy0;
    xy_d_t xyD = xy1;

    /*
    color_t blue   = { 0xffff0000 };
    color_t cyan   = { 0xffffff00 };
    color_t orange = { 0xff0099ff };
    color_t red = { 0xff0000ff };
    color_t yellow  = { 0xff00ffff };
    color_t green  = { 0xffff0000 };
    color_t black  = { 0xff000000 };
    */

    // adjust the moves to reflect our line
    if (xy0.x > xy1.x && xy0.y > xy1.y) {
        // heading NW
        move_in_x = move_in_x;
        move_in_y = -move_in_y;
        //wu_calc_line(xy0, xy1, blue, image, true, true);
    } else if (xy0.x >= xy1.x && xy0.y <= xy1.y) {
        // heading SW
        move_in_x = -move_in_x;
        move_in_y = -move_in_y;
        //wu_calc_line(xy0, xy1, red, image, true, true);
    } else if (xy0.x < xy1.x && xy0.y > xy1.y) {
        // heading NE
        move_in_x = move_in_x;
        move_in_y = move_in_y;
        //wu_calc_line(xy0, xy1, yellow, image, true, true);
    } else { //if //(xyA.x <= xyB.x && xyA.y <= xyB.y) {
        // heading SE
        move_in_x = -move_in_x;
        move_in_y = move_in_y;
        //wu_calc_line(xy0, xy1, black, image, true, true);
    }
    
    xyA.x += move_in_x;
    xyA.y += move_in_y;

    xyB.x += move_in_x;
    xyB.y += move_in_y;

    xyC.x -= move_in_x;
    xyC.y -= move_in_y;

    xyD.x -= move_in_x;
    xyD.y -= move_in_y;

    wu_rekt(xyA, xyB, xyC, xyD, color, image);

}

void wu_calc_rainbow(xy_d_t xy0, xy_d_t xy1, atomic_image_buf_t& image,
                     const std::vector<color_t>& colors,
                     const double& spacing,
                     const double& width) {

    // determine how
    double total_width = colors.size() * (width + spacing);


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
    double width_in_px = total_width / image.source_per_px_y;

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

    xy_d_t xyA = xy0;
    xy_d_t xyB = xy1;
    //xy_d_t xyC = xy0;
    //xy_d_t xyD = xy1;

    /*
    color_t blue   = { 0xffff0000 };
    color_t cyan   = { 0xffffff00 };
    color_t orange = { 0xff0099ff };
    color_t red = { 0xff0000ff };
    color_t yellow  = { 0xff00ffff };
    color_t green  = { 0xffff0000 };
    color_t black  = { 0xff000000 };
    */

    // adjust the moves to reflect our line
    if (xy0.x > xy1.x && xy0.y > xy1.y) {
        // heading NW
        move_in_x = move_in_x;
        move_in_y = -move_in_y;
        //wu_calc_line(xy0, xy1, blue, image, true, true);
    } else if (xy0.x >= xy1.x && xy0.y <= xy1.y) {
        // heading SW
        move_in_x = -move_in_x;
        move_in_y = -move_in_y;
        //wu_calc_line(xy0, xy1, red, image, true, true);
    } else if (xy0.x < xy1.x && xy0.y > xy1.y) {
        // heading NE
        move_in_x = move_in_x;
        move_in_y = move_in_y;
        //wu_calc_line(xy0, xy1, yellow, image, true, true);
    } else { //if //(xyA.x <= xyB.x && xyA.y <= xyB.y) {
        // heading SE
        move_in_x = -move_in_x;
        move_in_y = move_in_y;
        //wu_calc_line(xy0, xy1, black, image, true, true);
    }
    
    xyA.x += move_in_x;
    xyA.y += move_in_y;

    xyB.x += move_in_x;
    xyB.y += move_in_y;

    uint64_t total_steps = colors.size();
    double step_x = 2*move_in_x / (double)total_steps;
    double step_y = 2*move_in_y / (double)total_steps;
    double r_width = std::sqrt(step_x*step_x + step_y*step_y);

    for (uint64_t i = 0; i < total_steps; ++i) {
        xyA.x -= step_x;
        xyA.y -= step_y;
        xyB.x -= step_x;
        xyB.y -= step_y;
        wu_calc_wide_line(xyA, xyB,
                          colors[i],
                          image,
                          r_width);
    }
}

}

}
