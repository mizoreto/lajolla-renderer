#include "../microfacet.h"

Spectrum eval_op::operator()(const DisneySheen &bsdf) const {
    if (dot(vertex.geometric_normal, dir_in) < 0 ||
            dot(vertex.geometric_normal, dir_out) < 0) {
        // No light below the surface
        return make_zero_spectrum();
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }

    // Homework 1: implement this!
    Spectrum baseColor = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real sheen_tint = eval(bsdf.sheen_tint, vertex.uv, vertex.uv_screen_size, texture_pool);

    Vector3 half_vector = normalize(dir_in + dir_out);
    Real dot_h_out = fabs(dot(half_vector, dir_out));
    Real dot_n_out = fabs(dot(frame.n, dir_out));

    Real luminance_basecolor = luminance(baseColor);
    Spectrum c_tint;
    if (luminance_basecolor > 0) {
        c_tint = baseColor / luminance_basecolor;
    } else {
        c_tint = make_const_spectrum(1);
    }

    Spectrum c_sheen = (1 - sheen_tint) + sheen_tint * c_tint;
    Spectrum f_sheen = c_sheen * pow(1-dot_h_out, 5) * dot_n_out;

    return f_sheen;
    
}

Real pdf_sample_bsdf_op::operator()(const DisneySheen &bsdf) const {
    if (dot(vertex.geometric_normal, dir_in) < 0 ||
            dot(vertex.geometric_normal, dir_out) < 0) {
        // No light below the surface
        return 0;
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }

    // Homework 1: implement this!
    return fmax(dot(frame.n, dir_out), Real(0)) / c_PI;
}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const DisneySheen &bsdf) const {
    if (dot(vertex.geometric_normal, dir_in) < 0) {
        // No light below the surface
        return {};
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }

    // Homework 1: implement this!
    return BSDFSampleRecord{
        to_world(frame, sample_cos_hemisphere(rnd_param_uv)),
        Real(0) /* eta */, 1 /* roughness */};
}

TextureSpectrum get_texture_op::operator()(const DisneySheen &bsdf) const {
    return bsdf.base_color;
}
