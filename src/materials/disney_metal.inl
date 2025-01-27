#include "../microfacet.h"
inline Real smith_masking_metal(const Vector3 &v_local, Real roughness, Real anisotropic) {
    Real aspect = sqrt(Real(1) - Real(0.9) * anisotropic);
    Real alpha_x = fmax(0.0001, roughness * roughness / aspect);
    Real alpha_y = fmax(0.0001, roughness * roughness * aspect);
    Vector3 v2 = v_local * v_local;
    Real Lambda = (-1 + sqrt(1 + (v2.x * alpha_x * alpha_x + v2.y * alpha_y * alpha_y) / v2.z)) / 2;
    return 1 / (1 + Lambda);
}

inline Real GGX_metal(const Vector3 &local_h, Real roughness, Real anisotropic) {
    Real aspect = sqrt(Real(1) - Real(0.9) * anisotropic);
    Real alpha_x = fmax(0.0001, roughness * roughness / aspect);
    Real alpha_y = fmax(0.0001, roughness * roughness * aspect);
    Real alpha2_x = alpha_x * alpha_x;
    Real alpha2_y = alpha_y * alpha_y;
    Real denom = c_PI * alpha_x * alpha_y * pow(local_h.x * local_h.x / alpha2_x + local_h.y * local_h.y / alpha2_y + local_h.z * local_h.z, 2);

    return 1/denom;
}

Spectrum eval_op::operator()(const DisneyMetal &bsdf) const {
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
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    // Clamp roughness to avoid numerical issues.
    roughness = std::clamp(roughness, Real(0.01), Real(1));
    Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Vector3 half_vector = normalize(dir_in + dir_out);

    // if (dot(frame.n, dir_out) < 0 || dot(frame.n, half_vector) < 0) {
    //     return make_zero_spectrum();
    // } 

    Real dot_half_out = fabs(dot(half_vector, dir_out));
    Real dot_n_in = fabs(dot(frame.n, dir_in));
    Real dot_half_n = fabs(dot(half_vector, frame.n));

    // F_m
    Spectrum F_m = schlick_fresnel(baseColor, dot_half_out);
    // D_m
    // Real D_m = GGX(dot_half_n, roughness);
    Real D_m = GGX_metal(to_local(frame, half_vector), roughness, anisotropic);
    // G_m
    Real G_in = smith_masking_metal(to_local(frame, dir_in), roughness, anisotropic);
    Real G_out = smith_masking_metal(to_local(frame, dir_out), roughness, anisotropic);
    Real G_m = G_in * G_out;

    return F_m * D_m * G_m / (Real(4) * dot_n_in);
}

Real pdf_sample_bsdf_op::operator()(const DisneyMetal &bsdf) const {
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
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    // Clamp roughness to avoid numerical issues.
    roughness = std::clamp(roughness, Real(0.01), Real(1));

    Vector3 half_vector = normalize(dir_in + dir_out);
    if (dot(frame.n, dir_out) < 0 || dot(frame.n, half_vector) < 0) {
        return 0;
    } 

    Real n_dot_h = dot(frame.n, half_vector);
    Real n_dot_in = dot(frame.n, dir_in);

    Real G = smith_masking_metal(to_local(frame, dir_in), roughness, anisotropic);
    // Real D = GTR2(n_dot_h, roughness);
    Real D = GGX_metal(to_local(frame, half_vector), roughness, anisotropic);
    
    return (G * D) / (4 * n_dot_in);
}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const DisneyMetal &bsdf) const {
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

    // Convert the incoming direction to local coordinates
    Vector3 local_dir_in = to_local(frame, dir_in);
    Real roughness = eval(
        bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    // Clamp roughness to avoid numerical issues.
    roughness = std::clamp(roughness, Real(0.01), Real(1));
    Real alpha = roughness * roughness;
    Vector3 local_micro_normal =
        sample_visible_normals(local_dir_in, alpha, rnd_param_uv);
    
    // Transform the micro normal to world space
    Vector3 half_vector = to_world(frame, local_micro_normal);
    // Reflect over the world space normal
    Vector3 reflected = normalize(-dir_in + 2 * dot(dir_in, half_vector) * half_vector);
    return BSDFSampleRecord{
        reflected,
        Real(0) /* eta */, roughness /* roughness */
    };
}

TextureSpectrum get_texture_op::operator()(const DisneyMetal &bsdf) const {
    return bsdf.base_color;
}
