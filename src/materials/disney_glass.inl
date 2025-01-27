#include "../microfacet.h"
inline Real smith_masking_glass(const Vector3 &v_local, Real roughness, Real anisotropic) {
    Real aspect = sqrt(Real(1) - Real(0.9) * anisotropic);
    Real alpha_x = fmax(0.0001, roughness * roughness / aspect);
    Real alpha_y = fmax(0.0001, roughness * roughness * aspect);
    Vector3 v2 = v_local * v_local;
    Real Lambda = (-1 + sqrt(1 + (v2.x * alpha_x * alpha_x + v2.y * alpha_y * alpha_y) / v2.z)) / 2;
    return 1 / (1 + Lambda);
}

inline Real GGX_glass(const Vector3 &local_h, Real roughness, Real anisotropic) {
    Real aspect = sqrt(Real(1) - Real(0.9) * anisotropic);
    Real alpha_x = fmax(0.0001, roughness * roughness / aspect);
    Real alpha_y = fmax(0.0001, roughness * roughness * aspect);
    Real alpha2_x = alpha_x * alpha_x;
    Real alpha2_y = alpha_y * alpha_y;
    Real denom = c_PI * alpha_x * alpha_y * pow(local_h.x * local_h.x / alpha2_x + local_h.y * local_h.y / alpha2_y + local_h.z * local_h.z, 2);

    return 1/denom;
}

Spectrum eval_op::operator()(const DisneyGlass &bsdf) const {
    bool reflect = dot(vertex.geometric_normal, dir_in) *
                   dot(vertex.geometric_normal, dir_out) > 0;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    // If we are going into the surface, then we use normal eta
    // (internal/external), otherwise we use external/internal.
    Real eta = dot(vertex.geometric_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;
    Spectrum baseColor = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    // Clamp roughness to avoid numerical issues.
    roughness = std::clamp(roughness, Real(0.01), Real(1));
    Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool); 

    Vector3 half_vector;
    if (reflect) {
        half_vector = normalize(dir_in + dir_out);
    } else {
        // "Generalized half-vector" from Walter et al.
        // See "Microfacet Models for Refraction through Rough Surfaces"
        half_vector = normalize(dir_in + dir_out * eta);
    }

    // Flip half-vector if it's below surface
    if (dot(half_vector, frame.n) < 0) {
        half_vector = -half_vector;
    }

    Real dot_h_in = dot(half_vector, dir_in);
    Real dot_h_out = dot(half_vector, dir_out);
    Real dot_n_in = fabs(dot(frame.n, dir_in));
    Real dot_n_out = fabs(dot(frame.n, dir_out));

    // F_g
    Real F_g = fresnel_dielectric(dot_h_in, eta);
    
    // D_g
    Real D_g = GGX_glass(to_local(frame, half_vector), roughness, anisotropic);
    
    // G_g
    Real G_in = smith_masking_glass(to_local(frame, dir_in), roughness, anisotropic);
    Real G_out = smith_masking_glass(to_local(frame, dir_out), roughness, anisotropic);
    Real G_g = G_in * G_out;

    // Calculate bsdf
    if (reflect) {
        return baseColor * F_g * D_g * G_g / (4 * dot_n_in);
    } else {
        // Snell-Descartes law predicts that the light will contract/expand 
        // due to the different index of refraction. So the normal BSDF needs
        // to scale with 1/eta^2. However, the "adjoint" of the BSDF does not have
        // the eta term. This is due to the non-reciprocal nature of the index of refraction:
        // f(wi -> wo) / eta_o^2 = f(wo -> wi) / eta_i^2
        // thus f(wi -> wo) = f(wo -> wi) (eta_o / eta_i)^2
        // The adjoint of a BSDF is defined as swapping the parameter, and
        // this cancels out the eta term.
        // See Chapter 5 of Eric Veach's thesis "Robust Monte Carlo Methods for Light Transport Simulation"
        // for more details.
        Real eta_factor = dir == TransportDirection::TO_LIGHT ? (1 / (eta * eta)) : 1;
        Real sqrt_denom = dot_h_in + eta * dot_h_out;
        // Very complicated BSDF. See Walter et al.'s paper for more details.
        // "Microfacet Models for Refraction through Rough Surfaces"
        return sqrt(baseColor) * (eta_factor * (1 - F_g) * D_g * G_g * fabs(dot_h_out * dot_h_in)) / 
            (dot_n_in * sqrt_denom * sqrt_denom);
    }

}

Real pdf_sample_bsdf_op::operator()(const DisneyGlass &bsdf) const {
    bool reflect = dot(vertex.geometric_normal, dir_in) *
                   dot(vertex.geometric_normal, dir_out) > 0;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    // If we are going into the surface, then we use normal eta
    // (internal/external), otherwise we use external/internal.
    Real eta = dot(vertex.geometric_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;
    assert(eta > 0);

    Spectrum baseColor = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    // Clamp roughness to avoid numerical issues.
    roughness = std::clamp(roughness, Real(0.01), Real(1));
    Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool); 

    Vector3 half_vector;
    if (reflect) {
        half_vector = normalize(dir_in + dir_out);
    } else {
        // "Generalized half-vector" from Walter et al.
        // See "Microfacet Models for Refraction through Rough Surfaces"
        half_vector = normalize(dir_in + dir_out * eta);
    }

    // Flip half-vector if it's below surface
    if (dot(half_vector, frame.n) < 0) {
        half_vector = -half_vector;
    }


    Real dot_h_in = dot(half_vector, dir_in);
    Real dot_h_out = dot(half_vector, dir_out);
    Real dot_n_in = fabs(dot(frame.n, dir_in));
    Real dot_n_out = fabs(dot(frame.n, dir_out));

    // F_g
    Real F_g = fresnel_dielectric(dot_h_in, eta);

    // D_g
    Real D_g = GGX_glass(to_local(frame, half_vector), roughness, anisotropic);
    
    // G_g
    Real G_in = smith_masking_glass(to_local(frame, dir_in), roughness, anisotropic);

    // We sample the visible normals, also we use F to determine
    // whether to sample reflection or refraction
    // so PDF ~ F * D * G_in for reflection, PDF ~ (1 - F) * D * G_in for refraction.
    if (reflect) {
        return (F_g * D_g * G_in) / (4 * dot_n_in);
    } else {
        Real sqrt_denom = dot_h_in + eta * dot_h_out;
        Real dh_dout = fabs(dot_h_in * dot_h_out) / (sqrt_denom * sqrt_denom);
        return (1 - F_g) * D_g * G_in * dh_dout / dot_n_in;
    }
}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const DisneyGlass &bsdf) const {
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    Real eta = dot(vertex.geometric_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;

    Real roughness = eval(
        bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    // Clamp roughness to avoid numerical issues.
    roughness = std::clamp(roughness, Real(0.01), Real(1));
    // Sample a micro normal and transform it to world space -- this is our half-vector.
    Real alpha = roughness * roughness;
    Vector3 local_dir_in = to_local(frame, dir_in);
    Vector3 local_micro_normal =
        sample_visible_normals(local_dir_in, alpha, rnd_param_uv);

    Vector3 half_vector = to_world(frame, local_micro_normal);
    // Flip half-vector if it's below surface
    if (dot(half_vector, frame.n) < 0) {
        half_vector = -half_vector;
    }

    // F_g
    Real dot_h_in = dot(half_vector, dir_in);
    Real F_g = fresnel_dielectric(dot_h_in, eta);

    if (rnd_param_w <= F_g) {
        // Reflection
        Vector3 reflected = normalize(-dir_in + 2 * dot(dir_in, half_vector) * half_vector);
        // set eta to 0 since we are not transmitting
        return BSDFSampleRecord{reflected, Real(0) /* eta */, roughness};
    } else {
        // Refraction
        // https://en.wikipedia.org/wiki/Snell%27s_law#Vector_form
        // (note that our eta is eta2 / eta1, and l = -dir_in)
        Real h_dot_out_sq = 1 - (1 - dot_h_in * dot_h_in) / (eta * eta);
        if (h_dot_out_sq <= 0) {
            // Total internal reflection
            // This shouldn't really happen, as F will be 1 in this case.
            return {};
        }
        // flip half_vector if needed
        if (dot_h_in < 0) {
            half_vector = -half_vector;
        }
        Real h_dot_out= sqrt(h_dot_out_sq);
        Vector3 refracted = -dir_in / eta + (fabs(dot_h_in) / eta - h_dot_out) * half_vector;
        return BSDFSampleRecord{refracted, eta, roughness};
    }
}

TextureSpectrum get_texture_op::operator()(const DisneyGlass &bsdf) const {
    return bsdf.base_color;
}
