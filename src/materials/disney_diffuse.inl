// Fresnel Approximation Function
inline Real fresnel_approximation(Real f, Real cos) {
    return Real(1) + (f - Real(1)) * pow(Real(1) - cos, 5);
}

Spectrum eval_op::operator()(const DisneyDiffuse &bsdf) const {
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
    Real subsurface = eval(bsdf.subsurface, vertex.uv, vertex.uv_screen_size, texture_pool);
    
    const Vector3 half_vec = normalize(dir_in + dir_out);

    // if (dot(frame.n, dir_out) < 0 || dot(frame.n, half_vec) < 0) {
    //     return make_zero_spectrum();
    // } 
    
    Real dot_n_in = fabs(dot(frame.n, dir_in));
    Real dot_n_out = fabs(dot(frame.n, dir_out));
    Real dot_half_out = fabs(dot(half_vec, dir_out));

    // baseDiffuse
    Real f_d90 = Real(0.5) + Real(2) * roughness * pow(dot_half_out, 2);
    Real f_d_in = fresnel_approximation(f_d90, dot_n_in);
    Real f_d_out = fresnel_approximation(f_d90, dot_n_out);
    Spectrum f_base_diffuse = baseColor * f_d_in * f_d_out * dot_n_out / c_PI;

    // subsurface
    Real f_ss90 = roughness * pow(dot_half_out, 2);
    Real f_ss_in = fresnel_approximation(f_ss90, dot_n_in);
    Real f_ss_out = fresnel_approximation(f_ss90, dot_n_out);
    Spectrum f_subsurface = Real(1.25) * baseColor * (f_ss_in * f_ss_out * (Real(1)/(dot_n_in + dot_n_out) - Real(0.5))  + Real(0.5)) * dot_n_out / c_PI;

    return (Real(1)-subsurface) * f_base_diffuse + subsurface * f_subsurface;
}

Real pdf_sample_bsdf_op::operator()(const DisneyDiffuse &bsdf) const {
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

std::optional<BSDFSampleRecord> sample_bsdf_op::operator()(const DisneyDiffuse &bsdf) const {
    if (dot(vertex.geometric_normal, dir_in) < 0) {
        // No light below the surface
        return {};
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }

    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    // Clamp roughness to avoid numerical issues.
    roughness = std::clamp(roughness, Real(0.01), Real(1));
    
    // Homework 1: implement this!
    return BSDFSampleRecord{
        to_world(frame, sample_cos_hemisphere(rnd_param_uv)),
        Real(0) /* eta */, roughness /* roughness */};
}

TextureSpectrum get_texture_op::operator()(const DisneyDiffuse &bsdf) const {
    return bsdf.base_color;
}
