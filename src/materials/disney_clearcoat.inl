#include "../microfacet.h"
inline Real smith_masking_clearcoat(const Vector3 &v_local) {
    Real alpha2 = pow(Real(0.25), 2);
    Vector3 v2 = v_local * v_local;
    Real Lambda = (-1 + sqrt(1 + (v2.x * alpha2 + v2.y * alpha2) / v2.z)) / 2;
    return 1 / (1 + Lambda);
}

inline Real GGX_clearcoat(const Real alpha_g, const Real local_h_z) {
    Real numerator = alpha_g * alpha_g - Real(1);
    Real denominator = c_PI * log(alpha_g * alpha_g) * (Real(1) + (alpha_g * alpha_g - Real(1)) * local_h_z * local_h_z);
    return numerator / denominator;
}

Spectrum eval_op::operator()(const DisneyClearcoat &bsdf) const {
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
    Real clearcoat_gloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);
    Vector3 half_vector = normalize(dir_in + dir_out);

    // if (dot(frame.n, dir_out) < 0 || dot(frame.n, half_vector) < 0) {
    //     return make_zero_spectrum();
    // } 
    
    Real dot_half_out = fabs(dot(half_vector, dir_out));
    Real dot_n_in = fabs(dot(frame.n, dir_in));

    Real alpha_g = (Real(1) - clearcoat_gloss) * Real(0.1) + clearcoat_gloss * Real(0.001);

    // F_c
    Real eta = 1.5;
    Real R_0 = pow((eta - Real(1)), 2) / pow((eta + Real(1)), 2);
    Real F_c = schlick_fresnel(R_0, dot_half_out);

    // G_c
    Real G_c_in = smith_masking_clearcoat(to_local(frame, dir_in));
    Real G_c_out = smith_masking_clearcoat(to_local(frame, dir_out));
    Real G_c = G_c_in * G_c_out;

    // D_c
    Real local_h_z = to_local(frame, half_vector).z;
    Real D_c = GGX_clearcoat(alpha_g, local_h_z);

    return make_const_spectrum(F_c * D_c * G_c / (Real(4) * dot_n_in));
}

Real pdf_sample_bsdf_op::operator()(const DisneyClearcoat &bsdf) const {
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
    Real clearcoat_gloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);
    Vector3 half_vector = normalize(dir_in + dir_out);

    if (dot(frame.n, dir_out) < 0 || dot(frame.n, half_vector) < 0) {
        return 0;
    } 

    Real dot_half_out = fabs(dot(half_vector, dir_out));
    Real dot_half_n = fabs(dot(half_vector, frame.n));

    Real alpha_g = (Real(1) - clearcoat_gloss) * Real(0.1) + clearcoat_gloss * Real(0.001);
    Real local_h_z = to_local(frame, half_vector).z;
    Real D_c = GGX_clearcoat(alpha_g, local_h_z);

    return (D_c * dot_half_n) / (4 * dot_half_out);
}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const DisneyClearcoat &bsdf) const {
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
    Real clearcoat_gloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real alpha_g = (Real(1) - clearcoat_gloss) * Real(0.1) + clearcoat_gloss * Real(0.001);

    Real h_elevation = acos(sqrt((1 - pow(alpha_g * alpha_g, 1 - rnd_param_uv.x)) / (1 - alpha_g * alpha_g)));
    Real h_azimuth = 2 * c_PI * rnd_param_uv.y;
    Real local_h_x = sin(h_elevation) * cos(h_azimuth);
    Real local_h_y = sin(h_elevation) * sin(h_azimuth);
    Real local_h_z = cos(h_elevation);

    Vector3 local_h = Vector3(local_h_x, local_h_y, local_h_z);

    // Transform the micro normal to world space
    Vector3 half_vector = to_world(frame, local_h);
    // Reflect over the world space normal
    Vector3 reflected = normalize(-dir_in + 2 * dot(dir_in, half_vector) * half_vector);
    return BSDFSampleRecord{
        reflected,
        Real(0) /* eta */, 1 /* roughness */
    };
}

TextureSpectrum get_texture_op::operator()(const DisneyClearcoat &bsdf) const {
    return make_constant_spectrum_texture(make_zero_spectrum());
}
