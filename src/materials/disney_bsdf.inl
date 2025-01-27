#include "../microfacet.h"

inline Spectrum modified_eval_metal(const DisneyMetal &bsdf,  
                                const Vector3 &dir_in,
                                const Vector3 &dir_out,
                                const PathVertex &vertex,
                                const TexturePool &texture_pool,
                                const Spectrum &C0,
                                TransportDirection dir = TransportDirection::TO_LIGHT) {

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

    if (dot(frame.n, dir_out) < 0 || dot(frame.n, half_vector) < 0) {
        return make_zero_spectrum();
    } 

    Real dot_half_out = fabs(dot(half_vector, dir_out));
    Real dot_n_in = fabs(dot(frame.n, dir_in));
    Real dot_half_n = fabs(dot(half_vector, frame.n));

    // F_m
    Spectrum F_m = C0 + (1 - C0) * pow(1 - dot(half_vector, dir_out), 5);
    // D_m
    Real D_m = GGX_metal(to_local(frame, half_vector), roughness, anisotropic);
    // G_m
    Real G_in = smith_masking_metal(to_local(frame, dir_in), roughness, anisotropic);
    Real G_out = smith_masking_metal(to_local(frame, dir_out), roughness, anisotropic);
    Real G_m = G_in * G_out;

    return F_m * D_m * G_m / (Real(4) * dot_n_in);
}



Spectrum eval_op::operator()(const DisneyBSDF &bsdf) const {
    bool reflect = dot(vertex.geometric_normal, dir_in) *
                   dot(vertex.geometric_normal, dir_out) > 0;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    Spectrum base_color = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real specular_transmission = eval(bsdf.specular_transmission, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real metallic = eval(bsdf.metallic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real subsurface = eval(bsdf.subsurface, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real specular = eval(bsdf.specular, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    // Clamp roughness to avoid numerical issues.
    roughness = std::clamp(roughness, Real(0.01), Real(1));
    Real specular_tint = eval(bsdf.specular_tint, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real sheen = eval(bsdf.sheen, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real sheen_tint = eval(bsdf.sheen_tint, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real clearcoat = eval(bsdf.clearcoat, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real clearcoat_gloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real eta = dot(vertex.geometric_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;

    // create layers
    DisneyDiffuse diffuse_layer = DisneyDiffuse{bsdf.base_color, bsdf.roughness, bsdf.subsurface};
    DisneyMetal metal_layer = DisneyMetal{bsdf.base_color, bsdf.roughness, bsdf.anisotropic};
    DisneyClearcoat clearcoat_layer = DisneyClearcoat{bsdf.clearcoat_gloss};
    DisneyGlass glass_layer = DisneyGlass{bsdf.base_color, bsdf.roughness, bsdf.anisotropic, bsdf.eta};
    DisneySheen sheen_layer = DisneySheen{bsdf.base_color, bsdf.sheen_tint};

    // calculate BSDFs
    Spectrum f_diffuse = eval(diffuse_layer, dir_in, dir_out, vertex, texture_pool);
    Spectrum f_clearcoat = eval(clearcoat_layer, dir_in, dir_out, vertex, texture_pool);
    Spectrum f_glass = eval(glass_layer, dir_in, dir_out, vertex, texture_pool);
    Spectrum f_sheen = eval(sheen_layer, dir_in, dir_out, vertex, texture_pool);

    // calculate metal BSDF
    Spectrum c_tint;
    Real luminance_basecolor = luminance(base_color);
    if (luminance_basecolor > 0) {
        c_tint = base_color / luminance_basecolor;
    } else {
        c_tint = make_const_spectrum(1);
    }
    Spectrum Ks = (1 - specular_tint) + specular_tint * c_tint;
    Real R_0 = pow((eta - Real(1)), 2) / pow((eta + Real(1)), 2);
    Spectrum C0 = specular * R_0 * (1 - metallic) * Ks + metallic * base_color;
    Spectrum f_metal = modified_eval_metal(metal_layer, dir_in, dir_out, vertex, texture_pool, C0);

    Spectrum final_diffuse = (1 - specular_transmission) * (1 - metallic) * f_diffuse;
    Spectrum final_sheen = (1 - metallic) * sheen * f_sheen;
    Spectrum final_clearcoat = 0.25 * clearcoat * f_clearcoat;
    Spectrum final_glass = (1 - metallic) * specular_transmission * f_glass;
    Spectrum final_metal = (1 - specular_transmission * (1 - metallic)) * f_metal;

    if (dot(vertex.geometric_normal, dir_in) <= 0) {
        return final_glass;
    } else {
        return final_diffuse + final_sheen + final_metal + final_clearcoat + final_glass;
    }
    
}

Real pdf_sample_bsdf_op::operator()(const DisneyBSDF &bsdf) const {
    bool reflect = dot(vertex.geometric_normal, dir_in) *
                   dot(vertex.geometric_normal, dir_out) > 0;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    Spectrum base_color = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real specular_transmission = eval(bsdf.specular_transmission, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real metallic = eval(bsdf.metallic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real subsurface = eval(bsdf.subsurface, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real specular = eval(bsdf.specular, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    // Clamp roughness to avoid numerical issues.
    roughness = std::clamp(roughness, Real(0.01), Real(1));
    Real specular_tint = eval(bsdf.specular_tint, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real sheen = eval(bsdf.sheen, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real sheen_tint = eval(bsdf.sheen_tint, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real clearcoat = eval(bsdf.clearcoat, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real clearcoat_gloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real eta = dot(vertex.geometric_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;

    // create layers
    DisneyDiffuse diffuse_layer = DisneyDiffuse{bsdf.base_color, bsdf.roughness, bsdf.subsurface};
    DisneyMetal metal_layer = DisneyMetal{bsdf.base_color, bsdf.roughness, bsdf.anisotropic};
    DisneyClearcoat clearcoat_layer = DisneyClearcoat{bsdf.clearcoat_gloss};
    DisneyGlass glass_layer = DisneyGlass{bsdf.base_color, bsdf.roughness, bsdf.anisotropic, bsdf.eta};

    // pdf_sample_bsdf
    Real p_diffuse = pdf_sample_bsdf(diffuse_layer, dir_in, dir_out, vertex, texture_pool);
    // Real p_metal = pdf_sample_bsdf(metal_layer, dir_in, dir_out, vertex, texture_pool);
    Real p_clearcoat = pdf_sample_bsdf(clearcoat_layer, dir_in, dir_out, vertex, texture_pool);
    Real p_glass = pdf_sample_bsdf(glass_layer, dir_in, dir_out, vertex, texture_pool);


    // pdf_metal
    Real p_metal;
    if (dot(vertex.geometric_normal, dir_in) < 0 ||
            dot(vertex.geometric_normal, dir_out) < 0) {
        // No light below the surface
        p_metal = 0;
    } else {
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
            p_metal = 0;
        } else {
            Real n_dot_h = dot(frame.n, half_vector);
            Real n_dot_in = dot(frame.n, dir_in);

            Real G = smith_masking_metal(to_local(frame, dir_in), roughness, anisotropic);
            // Real D = GTR2(n_dot_h, roughness);
            Real D = GGX_metal(to_local(frame, half_vector), roughness, anisotropic);
            
            p_metal =  (G * D) / (4 * n_dot_in);
        }
    }
    

    // calculate weight
    Real diffuseWeight = (1 - metallic) * (1 - specular_transmission);
    Real metalWeight = (1 - specular_transmission* (1 - metallic));
    Real glassWeight = (1 - metallic) * specular_transmission;
    Real clearcoatWeight = 0.25 * clearcoat;

    Real total = diffuseWeight + metalWeight + glassWeight + clearcoatWeight;
    diffuseWeight /= total;
    metalWeight /= total;
    glassWeight /= total;
    clearcoatWeight /= total;

    if (dot(vertex.geometric_normal, dir_in) <= 0) {
        return p_glass;
    } else {
        return diffuseWeight * p_diffuse + metalWeight * p_metal + clearcoatWeight * p_clearcoat + glassWeight * p_glass;
    }
}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const DisneyBSDF &bsdf) const {
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    Spectrum base_color = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real specular_transmission = eval(bsdf.specular_transmission, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real metallic = eval(bsdf.metallic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real subsurface = eval(bsdf.subsurface, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real specular = eval(bsdf.specular, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    // Clamp roughness to avoid numerical issues.
    roughness = std::clamp(roughness, Real(0.01), Real(1));
    Real specular_tint = eval(bsdf.specular_tint, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real sheen = eval(bsdf.sheen, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real sheen_tint = eval(bsdf.sheen_tint, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real clearcoat = eval(bsdf.clearcoat, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real clearcoat_gloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real eta = dot(vertex.geometric_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;

    // create layers
    DisneyDiffuse diffuse_layer = DisneyDiffuse{bsdf.base_color, bsdf.roughness, bsdf.subsurface};
    DisneyMetal metal_layer = DisneyMetal{bsdf.base_color, bsdf.roughness, bsdf.anisotropic};
    DisneyClearcoat clearcoat_layer = DisneyClearcoat{bsdf.clearcoat_gloss};
    DisneyGlass glass_layer = DisneyGlass{bsdf.base_color, bsdf.roughness, bsdf.anisotropic, bsdf.eta};

    // calculate weight
    Real diffuseWeight = (1 - metallic) * (1 - specular_transmission);
    Real metalWeight = (1 - specular_transmission* (1 - metallic));
    Real glassWeight = (1 - metallic) * specular_transmission;
    Real clearcoatWeight = 0.25 * clearcoat;

    if (dot(vertex.geometric_normal, dir_in) <= 0) {
        return sample_bsdf(glass_layer, dir_in, vertex, texture_pool, rnd_param_uv, rnd_param_w);
    }

    Real total = diffuseWeight + metalWeight + glassWeight + clearcoatWeight;
    Real random_diffuse_metal = diffuseWeight / total;
    Real random_metal_glass = (diffuseWeight + metalWeight) / total;
    Real random_glass_clearcoat = (diffuseWeight + metalWeight + glassWeight) / total;

    if (rnd_param_w < random_diffuse_metal) {
        return sample_bsdf(diffuse_layer, dir_in, vertex, texture_pool, rnd_param_uv, rnd_param_w);
    } else if (random_diffuse_metal <= rnd_param_w  && rnd_param_w < random_metal_glass) {
        return sample_bsdf(metal_layer, dir_in, vertex, texture_pool, rnd_param_uv, rnd_param_w);
    } else if (random_metal_glass <= rnd_param_w && rnd_param_w < random_glass_clearcoat) {
        return sample_bsdf(glass_layer, dir_in, vertex, texture_pool, rnd_param_uv, rnd_param_w);
    } else {
        return sample_bsdf(clearcoat_layer, dir_in, vertex, texture_pool, rnd_param_uv, rnd_param_w);
    }
}

TextureSpectrum get_texture_op::operator()(const DisneyBSDF &bsdf) const {
    return bsdf.base_color;
}
