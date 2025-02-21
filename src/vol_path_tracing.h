#pragma once
#include "scene.h"
#include "pcg.h"
#include "phase_function.h"

// The simplest volumetric renderer: 
// single absorption only homogeneous volume
// only handle directly visible light sources
Spectrum vol_path_tracing_1(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                       (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    RayDifferential ray_diff = RayDifferential{Real(0), Real(0)}; //disable ray differentials

    std::optional<PathVertex> vertex_ = intersect(scene, ray, ray_diff);
    if (vertex_) {
        PathVertex vertex = *vertex_;
        Spectrum radiance = make_zero_spectrum();

        Spectrum sigma_a = get_sigma_a(scene.media[scene.camera.medium_id], vertex.position);
        Real t = distance(vertex.position, ray.org);
        Spectrum transmittance = exp(-sigma_a*t);
        Spectrum Le = make_zero_spectrum();
        if (is_light(scene.shapes[vertex.shape_id])) {
            Le = emission(vertex, -ray.dir, scene);
        }
        radiance = transmittance * Le;
        return radiance;
    }

    return make_zero_spectrum();
}

// The second simplest volumetric renderer: 
// single monochromatic homogeneous volume with single scattering,
// no need to handle surface lighting, only directly visible light source
Spectrum vol_path_tracing_2(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                       (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    RayDifferential ray_diff = RayDifferential{Real(0), Real(0)}; //disable ray differentials
    std::optional<PathVertex> vertex_ = intersect(scene, ray, ray_diff);

    Spectrum sigma_a = get_sigma_a(scene.media[scene.camera.medium_id], make_zero_spectrum());
    Spectrum sigma_s = get_sigma_s(scene.media[scene.camera.medium_id], make_zero_spectrum());
    Spectrum sigma_t = sigma_a + sigma_s;

    Real u = next_pcg32_real<Real>(rng);
    Spectrum t = - log(1-u) / sigma_t;

    Real t_hit = infinity<Real>();

    if (vertex_) {
        PathVertex vertex = *vertex_;
        t_hit = distance(vertex.position, ray.org);
    }

    Spectrum radiance = make_zero_spectrum();
    if (t.x < t_hit) {
        Spectrum trans_pdf = exp(-sigma_t * t) * sigma_t;
        Spectrum transmittance = exp(-sigma_t * t);

        Vector3 p_scattering = ray.org + t * ray.dir;

        Spectrum L_s1_estimate = make_zero_spectrum();
        Real L_s1_pdf = 1;

        //Sample a point on a light source
        Vector2 light_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
        Real light_w = next_pcg32_real<Real>(rng);
        Real shape_w = next_pcg32_real<Real>(rng);
        int light_id = sample_light(scene, light_w);
        const Light &light = scene.lights[light_id];
        PointAndNormal point_on_light =
            sample_point_on_light(light, p_scattering, light_uv, shape_w, scene);

        // no environment map in this hw
        Vector3 dir_light = normalize(point_on_light.position - p_scattering);
        Real G = 0;
        Ray shadow_ray{p_scattering, dir_light, 
                               get_shadow_epsilon(scene),
                               (1 - get_shadow_epsilon(scene)) *
                                   distance(point_on_light.position, p_scattering)};
        if (!occluded(scene, shadow_ray)) {
            G = max(-dot(dir_light, point_on_light.normal), Real(0)) /
                distance_squared(point_on_light.position, p_scattering);
            L_s1_pdf = light_pmf(scene, light_id) *
                pdf_point_on_light(light, point_on_light, p_scattering, scene);
            Spectrum Le = emission(light, -dir_light, Real(0), point_on_light, scene);
            PhaseFunction phaseFunction = get_phase_function(scene.media[scene.camera.medium_id]);
            Spectrum ro = eval(phaseFunction, -ray.dir, dir_light);
            L_s1_estimate = ro * Le * exp(-sigma_t * distance(point_on_light.position, p_scattering)) * G;
        }

        return (transmittance / trans_pdf) * sigma_s * (L_s1_estimate / L_s1_pdf);

    } else {
        Spectrum trans_pdf = exp(-sigma_t * t_hit);
        Spectrum transmittance = exp(-sigma_t * t_hit);
        Spectrum Le = make_zero_spectrum();
        if (is_light(scene.shapes[(*vertex_).shape_id])) {
            Le = emission(*vertex_, -ray.dir, scene);
        }
        radiance = transmittance * Le / trans_pdf;
        return radiance;
    }

    return make_zero_spectrum();
}

//update meium id
int update_medium(Ray ray, PathVertex isect, int medium) {
    if (isect.interior_medium_id != isect.exterior_medium_id) {
        // at medium transition. update medium.
        if (dot(ray.dir, isect.geometric_normal) > 0) {
            medium = isect.exterior_medium_id;
        } else {
            medium = isect.interior_medium_id;
        }
    }
    return medium;
}

// The third volumetric renderer (not so simple anymore): 
// multiple monochromatic homogeneous volumes with multiple scattering
// no need to handle surface lighting, only directly visible light source
Spectrum vol_path_tracing_3(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                       (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    RayDifferential ray_diff = RayDifferential{Real(0), Real(0)}; //disable ray differentials

    int current_medium_id = scene.camera.medium_id;
    Spectrum current_path_throughput = fromRGB(Vector3{1, 1, 1});
    Spectrum radiance = make_zero_spectrum();
    int bounces = 0;
    int max_depth = scene.options.max_depth;

    while (true) {
        bool scatter = false;
        std::optional<PathVertex> vertex_ = intersect(scene, ray, ray_diff);
        // !!set the ray.tnear to this after the ray initially came from the camera; no need to update the ray like line 212 each time
        ray.tnear = get_intersection_epsilon(scene);
        Spectrum transmittance = make_const_spectrum(1);
        Spectrum trans_pdf = make_const_spectrum(1);

        Real t_hit = infinity<Real>();

        if (vertex_) {
            PathVertex vertex = *vertex_;
            t_hit = distance(vertex.position, ray.org);
        }

        // if pass through medium
        if (current_medium_id != -1) {
            Spectrum sigma_a = get_sigma_a(scene.media[current_medium_id], make_zero_spectrum());
            Spectrum sigma_s = get_sigma_s(scene.media[current_medium_id], make_zero_spectrum());
            Spectrum sigma_t = sigma_a + sigma_s;

            Real u = next_pcg32_real<Real>(rng);
            Spectrum t = - log(1-u) / sigma_t;

            // the travel distance less that the distance to hit a surface => scattering happen
            if (t.x < t_hit) {
                scatter = true;
                trans_pdf = exp(-sigma_t * t) * sigma_t;
                transmittance = exp(-sigma_t * t);

                ray.org = ray.org + t * ray.dir;
            } else {
                trans_pdf = exp(-sigma_t * t_hit);
                transmittance = exp(-sigma_t * t_hit);

                ray.org = ray.org + t_hit * ray.dir;
            }
        } else if (t_hit != infinity<Real>()){
            ray.org = (*vertex_).position;
        }

        current_path_throughput *= transmittance / trans_pdf;

        // hit emissive surface & no scattering
        if (!scatter && vertex_) {
            Spectrum Le = make_zero_spectrum();
            if (is_light(scene.shapes[(*vertex_).shape_id])) {
                Le = emission(*vertex_, -ray.dir, scene);
            }
            radiance += current_path_throughput * Le;
        }

        // reach maximum bounces
        if (bounces == max_depth - 1 && max_depth != -1) {
            break;
        }

        // hit index-matching surface && no scattering
        if (!scatter && vertex_) {
            if ((*vertex_).material_id == -1) {
                // !! avoid self-intersection of the ray; avoid it from keep hitting the same point and not going forward; see line 155
                ray = {ray.org, ray.dir, get_intersection_epsilon(scene), infinity<Real>()};
                //ray = {ray.org, ray.dir, ray.tnear, ray.tfar};
                current_medium_id = update_medium(ray, *vertex_, current_medium_id);
                bounces += 1;
                continue;
            }
        }

        if (scatter) {
            // scattering
            Spectrum sigma_s = get_sigma_s(scene.media[current_medium_id], make_zero_spectrum());
            PhaseFunction phaseFunction = get_phase_function(scene.media[current_medium_id]);
            Vector2 phase_rnd_param_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
            std::optional<Vector3> next_dir_ = sample_phase_function(phaseFunction, -ray.dir, phase_rnd_param_uv);

            if (next_dir_) {
                Vector3 next_dir = *next_dir_;
                current_path_throughput *= (eval(phaseFunction, -ray.dir, next_dir) / pdf_sample_phase(phaseFunction, -ray.dir, next_dir)) * sigma_s;
                ray.dir = next_dir;
            } else {
                break;
            }
            
        } else {
            break;
        }

        // Russian roulette
        Real rr_prob = 1;
        if (bounces >= scene.options.rr_depth) {
            rr_prob = min(max(current_path_throughput), Real(0.95));
            if (next_pcg32_real<Real>(rng) > rr_prob) {
                // Terminate the path
                break;
            } else {
                current_path_throughput /= rr_prob;
            }
        }
        bounces += 1;
    }

    return radiance;
}

Spectrum next_event_estimation(const Scene &scene, const Ray &ray, Vector3 p, int current_medium_id, int& bounces, bool isScatter, const std::optional<PathVertex> vertex_, pcg32_state &rng) {
    // sample a point on a light source
    Vector2 light_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
    Real light_w = next_pcg32_real<Real>(rng);
    Real shape_w = next_pcg32_real<Real>(rng);
    int light_id = sample_light(scene, light_w);
    const Light &light = scene.lights[light_id];
    PointAndNormal point_on_light =
        sample_point_on_light(light, p, light_uv, shape_w, scene);
    Vector3 p_prime = point_on_light.position;

    Spectrum T_light = make_const_spectrum(1);
    int shadow_medium_id = current_medium_id;
    int shadow_bounces = 0;
    Spectrum p_trans_dir = make_const_spectrum(1);
    RayDifferential ray_diff = RayDifferential{Real(0), Real(0)}; //disable ray differentials
    int max_depth = scene.options.max_depth;

    //for NEE
    Vector3 position = p;

    while (true) {
        Vector3 dir_light = normalize(p_prime - p);
        Ray shadow_ray{p, dir_light, 
                               get_shadow_epsilon(scene),
                               (1 - get_shadow_epsilon(scene)) *
                                   distance(p_prime, p)};
        std::optional<PathVertex> shadow_vertex = intersect(scene, shadow_ray, ray_diff);
        Real next_t = distance(p, p_prime);
        if (shadow_vertex) {
            next_t = distance(p, (*shadow_vertex).position);
        }

        if (shadow_medium_id != -1) {
            Spectrum sigma_a = get_sigma_a(scene.media[shadow_medium_id], make_zero_spectrum());
            Spectrum sigma_s = get_sigma_s(scene.media[shadow_medium_id], make_zero_spectrum());
            Spectrum sigma_t = sigma_a + sigma_s;
            T_light *= exp(-sigma_t * next_t);
            p_trans_dir *= exp(-sigma_t * next_t);
        }

        // nothing is blocking, directly reach point on light
        if (!shadow_vertex) {
            break;
        } 
        // something is blocking
        else {
            //is it an opaque surface?
            if ((*shadow_vertex).material_id >= 0) {
                // we're blocked
                return make_zero_spectrum();
            }
            // it's index-matching surface
            shadow_bounces += 1;
            if (max_depth != -1 && (bounces + shadow_bounces + 1) >= max_depth) {
                //reach the max number of vertices
                return make_zero_spectrum();
            }
            shadow_medium_id = update_medium(shadow_ray, *shadow_vertex, shadow_medium_id);
            p = p + next_t * dir_light;
        }
    }

    if (T_light.x > 0) {
        Vector3 light_dir = normalize(p_prime - position);
        Real G = max(-dot(light_dir, point_on_light.normal), Real(0)) /
                        distance_squared(p_prime, position);
        Spectrum L = emission(light, -light_dir, Real(0), point_on_light, scene);
        Real pdf_nee = light_pmf(scene, light_id) *
                    pdf_point_on_light(light, point_on_light, position, scene);
        
        // scattering => medium => get phase function sampling pdf
        if (isScatter) {
            PhaseFunction phaseFunction = get_phase_function(scene.media[current_medium_id]);
            Spectrum rho = eval(phaseFunction, -ray.dir, light_dir);

            Spectrum contrib = T_light * G * rho * L / pdf_nee;
            Spectrum pdf_phase = pdf_sample_phase(phaseFunction, -ray.dir, light_dir) * G * p_trans_dir;
            Spectrum w = (pdf_nee * pdf_nee) / (pdf_nee * pdf_nee + pdf_phase * pdf_phase);

            return w * contrib;
        } 
        // hit something => surface => get bsdf sampling pdf
        else {
            const Material &mat = scene.materials[vertex_->material_id];
            Spectrum bsdf_f = eval(mat, -ray.dir, light_dir, *vertex_, scene.texture_pool);
            
            Spectrum contrib = T_light * G * bsdf_f * L / pdf_nee;
            Spectrum bsdf_pdf = pdf_sample_bsdf(mat, -ray.dir, light_dir, *vertex_, scene.texture_pool) * G * p_trans_dir;
            Spectrum w = (pdf_nee * pdf_nee) / (pdf_nee * pdf_nee + bsdf_pdf * bsdf_pdf);

            return w * contrib;
        }
        
    }

    return make_zero_spectrum();
}

// The fourth volumetric renderer: 
// multiple monochromatic homogeneous volumes with multiple scattering
// with MIS between next event estimation and phase function sampling
// still no surface lighting
Spectrum vol_path_tracing_4(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                       (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    RayDifferential ray_diff = RayDifferential{Real(0), Real(0)}; //disable ray differentials

    int current_medium_id = scene.camera.medium_id;
    Spectrum current_path_throughput = fromRGB(Vector3{1, 1, 1});
    Spectrum radiance = make_zero_spectrum();
    int bounces = 0;
    int max_depth = scene.options.max_depth;

    Real dir_pdf = 0;
    Vector3 nee_p_cache;
    Spectrum multi_trans_pdf = make_const_spectrum(1);
    bool never_scatter = true;

    while (true) {
        bool scatter = false;
        std::optional<PathVertex> vertex_ = intersect(scene, ray, ray_diff);
        Spectrum transmittance = make_const_spectrum(1);
        Spectrum trans_pdf = make_const_spectrum(1);

        Real t_hit = infinity<Real>();

        if (vertex_) {
            PathVertex vertex = *vertex_;
            t_hit = distance(vertex.position, ray.org);
        }

        // if pass through medium
        if (current_medium_id != -1) {
            Spectrum sigma_a = get_sigma_a(scene.media[current_medium_id], make_zero_spectrum());
            Spectrum sigma_s = get_sigma_s(scene.media[current_medium_id], make_zero_spectrum());
            Spectrum sigma_t = sigma_a + sigma_s;

            Real u = next_pcg32_real<Real>(rng);
            Spectrum t = - log(1-u) / sigma_t;

            // the travel distance less that the distance to hit a surface => scattering happen
            if (t.x < t_hit) {
                scatter = true;
                
                trans_pdf = exp(-sigma_t * t) * sigma_t;
                transmittance = exp(-sigma_t * t);

                ray.org = ray.org + t * ray.dir;

            } else {
                trans_pdf = exp(-sigma_t * t_hit);
                transmittance = exp(-sigma_t * t_hit);

                ray.org = ray.org + t_hit * ray.dir;
            }
        } else if (t_hit != infinity<Real>()){
            ray.org = (*vertex_).position;
        }

        multi_trans_pdf *= trans_pdf;
        current_path_throughput *= transmittance / trans_pdf;

        // hit emissive surface
        if (!scatter && vertex_) {
            if (never_scatter) {
                Spectrum Le = make_zero_spectrum();
                if (is_light(scene.shapes[(*vertex_).shape_id])) {
                    Le = emission(*vertex_, -ray.dir, scene);
                }
                radiance += current_path_throughput * Le;
            } else {
                if (is_light(scene.shapes[(*vertex_).shape_id])) {
                    PointAndNormal point_on_light{(*vertex_).position, (*vertex_).geometric_normal};
                    Real light_id = get_area_light_id(scene.shapes[(*vertex_).shape_id]);
                    // assert((vertex_)->shape_id >= 0 && (vertex_)->shape_id < scene.shapes.size());
                    // std::cout << light_id << std::endl;
                    // assert(light_id >= 0);
                    Light light = scene.lights[light_id];
                    Real pdf_nee = pdf_point_on_light(light, point_on_light, nee_p_cache, scene);

                    Real G = max(-dot(ray.dir, point_on_light.normal), Real(0)) /
                            distance_squared(point_on_light.position, nee_p_cache);

                    Spectrum dir_pdf_ = dir_pdf * multi_trans_pdf * G;
                    Spectrum w = (dir_pdf_ * dir_pdf_) / (pdf_nee * pdf_nee + dir_pdf_ * dir_pdf_);

                    radiance += current_path_throughput * w * emission(*vertex_, -ray.dir, scene);
                }
                
            }

            // Spectrum Le = make_zero_spectrum();
            // if (is_light(scene.shapes[(*vertex_).shape_id])) {
            //     Le = emission(*vertex_, -ray.dir, scene);
            // }
            // radiance += current_path_throughput * Le;
        
        }

        // reach maximum bounces
        if (bounces == max_depth - 1 && max_depth != -1) {
            break;
        }

        // hit index-matching surface
        if (!scatter && vertex_) {
            if ((*vertex_).material_id == -1) {
                ray = {ray.org, ray.dir, get_intersection_epsilon(scene), infinity<Real>()};
                current_medium_id = update_medium(ray, *vertex_, current_medium_id);
                bounces += 1;
                continue;
            }
        }

        if (scatter) {
            never_scatter = false;
            // scattering
            Spectrum sigma_s = get_sigma_s(scene.media[current_medium_id], make_zero_spectrum());
            // NEE
            nee_p_cache = ray.org;
            Spectrum nee = next_event_estimation(scene, ray, ray.org, current_medium_id, bounces, true, vertex_, rng);
            radiance += current_path_throughput * nee * sigma_s;
            // phase function sampling
            PhaseFunction phaseFunction = get_phase_function(scene.media[current_medium_id]);
            Vector2 phase_rnd_param_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
            std::optional<Vector3> next_dir_ = sample_phase_function(phaseFunction, -ray.dir, phase_rnd_param_uv);
            
            if (next_dir_) {
                Vector3 next_dir = *next_dir_;
                dir_pdf = pdf_sample_phase(phaseFunction, -ray.dir, next_dir);
                current_path_throughput *= (eval(phaseFunction, -ray.dir, next_dir) / dir_pdf) * sigma_s;
                ray.dir = next_dir;
                multi_trans_pdf = make_const_spectrum(1);
            } else {
                break;
            }
            
        } else {
            break;
            
        }

        // Russian roulette
        Real rr_prob = 1;
        if (bounces >= scene.options.rr_depth) {
            rr_prob = min(max(current_path_throughput), Real(0.95));
            if (next_pcg32_real<Real>(rng) > rr_prob) {
                // Terminate the path
                break;
            } else {
                current_path_throughput /= rr_prob;
            }
        }
        bounces += 1;
    }

    return radiance;
}

// The fifth volumetric renderer: 
// multiple monochromatic homogeneous volumes with multiple scattering
// with MIS between next event estimation and phase function sampling
// with surface lighting
Spectrum vol_path_tracing_5(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                       (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    RayDifferential ray_diff = RayDifferential{Real(0), Real(0)}; //disable ray differentials

    int current_medium_id = scene.camera.medium_id;
    Spectrum current_path_throughput = fromRGB(Vector3{1, 1, 1});
    Spectrum radiance = make_zero_spectrum();
    int bounces = 0;
    int max_depth = scene.options.max_depth;

    Real dir_pdf = 0;
    Vector3 nee_p_cache;
    Spectrum multi_trans_pdf = make_const_spectrum(1);
    bool never_scatter = true;

    Real eta_scale = Real(1);

    while (true) {
        bool scatter = false;
        std::optional<PathVertex> vertex_ = intersect(scene, ray, ray_diff);
        Spectrum transmittance = make_const_spectrum(1);
        Spectrum trans_pdf = make_const_spectrum(1);

        Real t_hit = infinity<Real>();

        if (vertex_) {
            PathVertex vertex = *vertex_;
            t_hit = distance(vertex.position, ray.org);
        }

        // if pass through medium
        if (current_medium_id != -1) {
            Spectrum sigma_a = get_sigma_a(scene.media[current_medium_id], make_zero_spectrum());
            Spectrum sigma_s = get_sigma_s(scene.media[current_medium_id], make_zero_spectrum());
            Spectrum sigma_t = sigma_a + sigma_s;

            Real u = next_pcg32_real<Real>(rng);
            Spectrum t = - log(1-u) / sigma_t;

            // the travel distance less that the distance to hit a surface => scattering happen
            if (t.x < t_hit) {
                scatter = true;
                
                trans_pdf = exp(-sigma_t * t) * sigma_t;
                transmittance = exp(-sigma_t * t);

                ray.org = ray.org + t * ray.dir;
                ray = {ray.org, ray.dir, get_intersection_epsilon(scene), infinity<Real>()};
            } else {
                trans_pdf = exp(-sigma_t * t_hit);
                transmittance = exp(-sigma_t * t_hit);

                ray.org = ray.org + t_hit * ray.dir;
                ray = {ray.org, ray.dir, get_intersection_epsilon(scene), infinity<Real>()};
            }
        } else if (t_hit != infinity<Real>()){
            // we need to remember to change the origin of the ray when the ray pass through vacuum and hit something!
            ray.org = (*vertex_).position;
            ray = {ray.org, ray.dir, get_intersection_epsilon(scene), infinity<Real>()};
        }

        multi_trans_pdf *= trans_pdf;
        current_path_throughput *= transmittance / trans_pdf;

        // hit emissive surface
        if (!scatter && vertex_) {
            if (never_scatter) {
                Spectrum Le = make_zero_spectrum();
                if (is_light(scene.shapes[(*vertex_).shape_id])) {
                    Le = emission(*vertex_, -ray.dir, scene);
                }
                radiance += current_path_throughput * Le;
            } else {
                if (is_light(scene.shapes[(*vertex_).shape_id])) {
                    PointAndNormal point_on_light{(*vertex_).position, (*vertex_).geometric_normal};
                    Real light_id = get_area_light_id(scene.shapes[(*vertex_).shape_id]);
                    // assert((vertex_)->shape_id >= 0 && (vertex_)->shape_id < scene.shapes.size());
                    // std::cout << light_id << std::endl;
                    // assert(light_id >= 0);
                    Light light = scene.lights[light_id];
                    Real pdf_nee = pdf_point_on_light(light, point_on_light, nee_p_cache, scene);

                    Real G = max(-dot(ray.dir, point_on_light.normal), Real(0)) /
                            distance_squared(point_on_light.position, nee_p_cache);

                    Spectrum dir_pdf_ = dir_pdf * multi_trans_pdf * G;
                    Spectrum w = (dir_pdf_ * dir_pdf_) / (pdf_nee * pdf_nee + dir_pdf_ * dir_pdf_);

                    radiance += current_path_throughput * w * emission(*vertex_, -ray.dir, scene);
                }
                
            }

            // Spectrum Le = make_zero_spectrum();
            // if (is_light(scene.shapes[(*vertex_).shape_id])) {
            //     Le = emission(*vertex_, -ray.dir, scene);
            // }
            // radiance += current_path_throughput * Le;
        
        }

        // reach maximum bounces
        if (bounces >= max_depth - 1 && max_depth != -1) {
            break;
        }

        // hit index-matching surface
        if (!scatter && vertex_) {
            if ((*vertex_).material_id == -1) {
                ray = {ray.org, ray.dir, get_intersection_epsilon(scene), infinity<Real>()};
                current_medium_id = update_medium(ray, *vertex_, current_medium_id);
                bounces += 1;
                continue;
            }
        }

        if (scatter) {
            never_scatter = false;
            // scattering
            Spectrum sigma_s = get_sigma_s(scene.media[current_medium_id], make_zero_spectrum());
            // NEE
            nee_p_cache = ray.org;
            Spectrum nee = next_event_estimation(scene, ray, ray.org, current_medium_id, bounces, true, vertex_, rng);
            radiance += current_path_throughput * nee * sigma_s;
            // phase function sampling
            PhaseFunction phaseFunction = get_phase_function(scene.media[current_medium_id]);
            Vector2 phase_rnd_param_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
            std::optional<Vector3> next_dir_ = sample_phase_function(phaseFunction, -ray.dir, phase_rnd_param_uv);
            
            if (next_dir_) {
                Vector3 next_dir = *next_dir_;
                dir_pdf = pdf_sample_phase(phaseFunction, -ray.dir, next_dir);
                current_path_throughput *= (eval(phaseFunction, -ray.dir, next_dir) / dir_pdf) * sigma_s;
                ray.dir = next_dir;
                multi_trans_pdf = make_const_spectrum(1);
                ray = {ray.org, ray.dir, get_intersection_epsilon(scene), infinity<Real>()};
            } else {
                break;
            }
            
        } else if(vertex_) {
            // nee
            nee_p_cache = ray.org;
            never_scatter = false;
            Spectrum nee = next_event_estimation(scene, ray, ray.org, current_medium_id, bounces, false, vertex_, rng);
            radiance += current_path_throughput * nee;
            
            // bsdf sampling
            Vector3 dir_view = -ray.dir;
            Vector2 bsdf_rnd_param_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
            Real bsdf_rnd_param_w = next_pcg32_real<Real>(rng);
            const Material &mat = scene.materials[vertex_->material_id];
            std::optional<BSDFSampleRecord> bsdf_sample_ =
            sample_bsdf(mat,
                        dir_view,
                        *vertex_,
                        scene.texture_pool,
                        bsdf_rnd_param_uv,
                        bsdf_rnd_param_w);
            if (!bsdf_sample_) {
                // BSDF sampling failed. Abort the loop.
                break;
            }
            const BSDFSampleRecord &bsdf_sample = *bsdf_sample_;
            Vector3 dir_bsdf = bsdf_sample.dir_out;
            // Update ray differentials & eta_scale
            if (bsdf_sample.eta == 0) {
                ray_diff.spread = reflect(ray_diff, vertex_->mean_curvature, bsdf_sample.roughness);
            } else {
                ray_diff.spread = refract(ray_diff, vertex_->mean_curvature, bsdf_sample.eta, bsdf_sample.roughness);
                eta_scale /= (bsdf_sample.eta * bsdf_sample.eta);
                current_medium_id = update_medium(ray, *vertex_, current_medium_id);
            }

            dir_pdf = pdf_sample_bsdf(mat, dir_view, dir_bsdf, (*vertex_), scene.texture_pool);
            Spectrum bsdf_f = eval(mat, dir_view, dir_bsdf, (*vertex_), scene.texture_pool);

            if (dir_pdf <= 0) {
                // Numerical issue -- we generated some invalid rays.
                break;
            }
            current_path_throughput *= bsdf_f / dir_pdf;
            ray.dir = dir_bsdf;
            ray = Ray{ray.org, ray.dir, get_intersection_epsilon(scene), infinity<Real>()};
            multi_trans_pdf = make_const_spectrum(1);
        }

        // Russian roulette
        Real rr_prob = 1;
        if (bounces >= scene.options.rr_depth) {
            // rr_prob = min(max(current_path_throughput), Real(0.95));
            rr_prob = min(max((1 / eta_scale) * current_path_throughput), Real(0.95));
            if (next_pcg32_real<Real>(rng) > rr_prob) {
                // Terminate the path
                break;
            } else {
                current_path_throughput /= rr_prob;
            }
        }
        bounces += 1;
    }

    return radiance;
}

Spectrum final_next_event_estimation(const Scene &scene, const Ray &ray, Vector3 p, int current_medium_id, int& bounces, bool isScatter, const std::optional<PathVertex> vertex_, pcg32_state &rng) {
    // sample a point on a light source
    Vector2 light_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
    Real light_w = next_pcg32_real<Real>(rng);
    Real shape_w = next_pcg32_real<Real>(rng);
    int light_id = sample_light(scene, light_w);
    const Light &light = scene.lights[light_id];
    PointAndNormal point_on_light =
        sample_point_on_light(light, p, light_uv, shape_w, scene);
    Vector3 p_prime = point_on_light.position;

    Spectrum T_light = make_const_spectrum(1);
    int shadow_medium_id = current_medium_id;
    int shadow_bounces = 0;
    Spectrum p_trans_dir = make_const_spectrum(1);
    Spectrum p_trans_nee = make_const_spectrum(1);
    RayDifferential ray_diff = RayDifferential{Real(0), Real(0)}; //disable ray differentials
    int max_depth = scene.options.max_depth;

    //for NEE
    Vector3 position = p;
    
    while (true) {
        Vector3 dir_light = normalize(p_prime - p);
        Ray shadow_ray{p, dir_light, 
                               get_shadow_epsilon(scene),
                               (1 - get_shadow_epsilon(scene)) *
                                   distance(p_prime, p)};
        std::optional<PathVertex> shadow_vertex = intersect(scene, shadow_ray, ray_diff);
        Real next_t = distance(p, p_prime);
        if (shadow_vertex) {
            next_t = distance(p, (*shadow_vertex).position);
        }

        if (shadow_medium_id != -1) {
            Real u = next_pcg32_real<Real>(rng);
            Real channel = std::clamp(int(u * 3), 0, 2);
            Real iteration = 0;
            Real accum_t = 0;

            while (true) {
                Spectrum majorant = get_majorant(scene.media[shadow_medium_id], shadow_ray);
                if (majorant[channel] <= 0) {
                    break;
                }

                if (iteration >= scene.options.max_null_collisions) {
                    break;
                }
                Real t = -log(1 - next_pcg32_real<Real>(rng)) / majorant[channel];
                Real dt = next_t - accum_t;
                // update accumulated distance
                accum_t = min(accum_t + t, next_t);

                Vector3 temp_p = shadow_ray.org + accum_t * shadow_ray.dir;
                Spectrum sigma_a = get_sigma_a(scene.media[shadow_medium_id], temp_p);
                Spectrum sigma_s = get_sigma_s(scene.media[shadow_medium_id], temp_p);
                Spectrum sigma_t = sigma_a + sigma_s;
                Spectrum sigma_n = majorant - sigma_t;

                if (t < dt) {
                    T_light *= exp(-majorant * t) * sigma_n / max(majorant);
                    p_trans_nee *= exp(-majorant * t) * majorant / max(majorant);
                    Spectrum real_prob = sigma_t / majorant;
                    p_trans_dir *= exp(-majorant * t) * majorant * (1 - real_prob) / max(majorant);

                    if (max(T_light) <= 0) {
                        return make_zero_spectrum();
                    }
                } else {
                    T_light *= exp(-majorant * dt);
                    p_trans_nee *= exp(-majorant * dt);
                    p_trans_dir *= exp(-majorant * dt);
                    iteration +=1;
                    break;
                }
                iteration += 1;
            }
        }

        // nothing is blocking, directly reach point on light
        if (!shadow_vertex) {
            break;
        } 
        // something is blocking
        else {
            //is it an opaque surface?
            if ((*shadow_vertex).material_id >= 0) {
                // we're blocked
                return make_zero_spectrum();
            }
            // it's index-matching surface
            shadow_bounces += 1;
            if (max_depth != -1 && (bounces + shadow_bounces + 1) >= max_depth) {
                //reach the max number of vertices
                return make_zero_spectrum();
            }
            shadow_medium_id = update_medium(shadow_ray, *shadow_vertex, shadow_medium_id);
            p = p + next_t * dir_light;
        }
    }

    if (max(T_light) > 0) {
        Vector3 light_dir = normalize(p_prime - position);
        Real G = max(-dot(light_dir, point_on_light.normal), Real(0)) /
                        distance_squared(p_prime, position);
        Spectrum L = emission(light, -light_dir, Real(0), point_on_light, scene);
        Spectrum pdf_nee = p_trans_nee * light_pmf(scene, light_id) *
                    pdf_point_on_light(light, point_on_light, position, scene);
        
        // scattering => medium => get phase function sampling pdf
        if (isScatter && current_medium_id != -1) {
            PhaseFunction phaseFunction = get_phase_function(scene.media[current_medium_id]);
            Spectrum rho = eval(phaseFunction, -ray.dir, light_dir);

            Spectrum contrib = T_light * G * rho * L / avg(pdf_nee);
            Spectrum pdf_phase = pdf_sample_phase(phaseFunction, -ray.dir, light_dir) * G * p_trans_dir;
            Spectrum w = (pdf_nee * pdf_nee) / (pdf_nee * pdf_nee + pdf_phase * pdf_phase);

            return w * contrib;
        } 
        // hit something => surface => get bsdf sampling pdf
        else {
            const Material &mat = scene.materials[vertex_->material_id];
            Spectrum bsdf_f = eval(mat, -ray.dir, light_dir, *vertex_, scene.texture_pool);
            
            Spectrum contrib = T_light * G * bsdf_f * L / avg(pdf_nee);
            Spectrum bsdf_pdf = pdf_sample_bsdf(mat, -ray.dir, light_dir, *vertex_, scene.texture_pool) * G * p_trans_dir;
            Spectrum w = (pdf_nee * pdf_nee) / (pdf_nee * pdf_nee + bsdf_pdf * bsdf_pdf);

            return w * contrib;
        }
        
    }

    return make_zero_spectrum();
}

// The final volumetric renderer: 
// multiple chromatic heterogeneous volumes with multiple scattering
// with MIS between next event estimation and phase function sampling
// with surface lighting
Spectrum vol_path_tracing(const Scene &scene,
                          int x, int y, /* pixel coordinates */
                          pcg32_state &rng) {
    // Homework 2: implememt this!
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                       (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    RayDifferential ray_diff = RayDifferential{Real(0), Real(0)}; //disable ray differentials

    int current_medium_id = scene.camera.medium_id;
    Spectrum current_path_throughput = fromRGB(Vector3{1, 1, 1});
    Spectrum radiance = make_zero_spectrum();
    int bounces = 0;
    int max_depth = scene.options.max_depth;
    bool never_scatter = true;

    Real dir_pdf = 0;
    Vector3 nee_p_cache;
    Spectrum multi_trans_pdf = make_const_spectrum(1);

    Real eta_scale = Real(1);

    while (true) {
        bool scatter = false;
        std::optional<PathVertex> vertex_ = intersect(scene, ray, ray_diff);
        Spectrum transmittance = make_const_spectrum(1);
        Spectrum trans_dir_pdf = make_const_spectrum(1);
        Spectrum trans_nee_pdf = make_const_spectrum(1);

        Real t_hit = infinity<Real>();

        if (vertex_) {
            PathVertex vertex = *vertex_;
            t_hit = distance(vertex.position, ray.org);
        }

        // if pass through medium
        if (current_medium_id != -1) {
            Spectrum majorant = get_majorant(scene.media[current_medium_id], ray);

            Real u = next_pcg32_real<Real>(rng);
            Real channel = std::clamp(int(u * 3), 0, 2);
            Real accum_t = 0;
            Real iteration = 0;

            while (true) {
                if (majorant[channel] <= 0) {
                    break;
                }

                if (iteration >= scene.options.max_null_collisions) {
                    break;
                }
                Real t = -log(1 - next_pcg32_real<Real>(rng)) / majorant[channel];
                Real dt = t_hit - accum_t;
                // update accumulated distance
                accum_t = min(accum_t + t, t_hit);
                Vector3 temp_p = ray.org + accum_t * ray.dir;
                Spectrum sigma_a = get_sigma_a(scene.media[current_medium_id], temp_p);
                Spectrum sigma_s = get_sigma_s(scene.media[current_medium_id], temp_p);
                Spectrum sigma_t = sigma_a + sigma_s;
                Spectrum sigma_n = majorant - sigma_t;

                if (t < dt) {   // haven't reached the surface
                    //sample from real/fake particle events
                    Spectrum real_prob = sigma_t / majorant;
                    if (next_pcg32_real<Real>(rng) < real_prob[channel]) {
                        // hit a real particle
                        scatter = true;
                        never_scatter = false;
                        transmittance *= exp(-majorant * t) / max(majorant);
                        trans_dir_pdf *= exp(-majorant * t) * majorant * real_prob / max(majorant);
                        ray.org = ray.org + accum_t * ray.dir;
                        ray = {ray.org, ray.dir, get_intersection_epsilon(scene), infinity<Real>()};
                        break;
                    } else {
                        // hit a fake particle
                        transmittance *= exp(-majorant * t) * sigma_n / max(majorant);
                        trans_dir_pdf *= exp(-majorant * t) * majorant * (1 - real_prob) / max(majorant);
                        trans_nee_pdf *= exp(-majorant * t) * majorant / max(majorant);
                    }
                } else {
                    // reach the surface
                    transmittance *= exp(-majorant * dt);
                    trans_dir_pdf *= exp(-majorant * dt);
                    trans_nee_pdf *= exp(-majorant * dt);
                    ray.org = ray.org + t_hit * ray.dir;
                    ray = {ray.org, ray.dir, get_intersection_epsilon(scene), infinity<Real>()};
                    break;
                }
                iteration += 1;
            }
            
            // Spectrum t = - log(1-u) / sigma_t;

            // // the travel distance less that the distance to hit a surface => scattering happen
            // if (t.x < t_hit) {
            //     scatter = true;
                
            //     trans_pdf = exp(-sigma_t * t) * sigma_t;
            //     transmittance = exp(-sigma_t * t);

            //     ray.org = ray.org + t * ray.dir;
            //     ray = {ray.org, ray.dir, get_intersection_epsilon(scene), infinity<Real>()};
            // } 
            // else {
            //     trans_pdf = exp(-sigma_t * t_hit);
            //     transmittance = exp(-sigma_t * t_hit);

            //     ray.org = ray.org + t_hit * ray.dir;
            //     ray = {ray.org, ray.dir, get_intersection_epsilon(scene), infinity<Real>()};
            // }
        } else if (t_hit != infinity<Real>()){
            // we need to remember to change the origin of the ray when the ray pass through vacuum and hit something!
            ray.org = (*vertex_).position;
            ray = {ray.org, ray.dir, get_intersection_epsilon(scene), infinity<Real>()};
        }

        multi_trans_pdf *= trans_dir_pdf;
        current_path_throughput *= transmittance / avg(trans_dir_pdf);

        // hit emissive surface
        if (!scatter && vertex_) {
            if (never_scatter) {
                Spectrum Le = make_zero_spectrum();
                if (is_light(scene.shapes[(*vertex_).shape_id])) {
                    Le = emission(*vertex_, -ray.dir, scene);
                }
                radiance += current_path_throughput * Le;
            } else {
                if (is_light(scene.shapes[(*vertex_).shape_id])) {
                    PointAndNormal point_on_light{(*vertex_).position, (*vertex_).geometric_normal};
                    Real light_id = get_area_light_id(scene.shapes[(*vertex_).shape_id]);
                    // assert((vertex_)->shape_id >= 0 && (vertex_)->shape_id < scene.shapes.size());
                    // std::cout << light_id << std::endl;
                    // assert(light_id >= 0);
                    Light light = scene.lights[light_id];
                    Spectrum pdf_nee = trans_nee_pdf * pdf_point_on_light(light, point_on_light, nee_p_cache, scene);

                    Real G = max(-dot(ray.dir, point_on_light.normal), Real(0)) /
                            distance_squared(point_on_light.position, nee_p_cache);

                    Spectrum dir_pdf_ = dir_pdf * multi_trans_pdf * G;
                    Spectrum w = (dir_pdf_ * dir_pdf_) / (pdf_nee * pdf_nee + dir_pdf_ * dir_pdf_);

                    radiance += current_path_throughput * w * emission(*vertex_, -ray.dir, scene);
                }
                
            }

            // Spectrum Le = make_zero_spectrum();
            // if (is_light(scene.shapes[(*vertex_).shape_id])) {
            //     Le = emission(*vertex_, -ray.dir, scene);
            // }
            // radiance += current_path_throughput * Le;
        
        }

        // reach maximum bounces
        if (bounces >= max_depth - 1 && max_depth != -1) {
            break;
        }

        // hit index-matching surface
        if (!scatter && vertex_) {
            if ((*vertex_).material_id == -1) {
                ray = {(*vertex_).position, ray.dir, get_intersection_epsilon(scene), infinity<Real>()};
                current_medium_id = update_medium(ray, *vertex_, current_medium_id);
                bounces += 1;
                continue;
            }
        }

        if (scatter) {
            never_scatter = false;
            // scattering
            Spectrum sigma_s = get_sigma_s(scene.media[current_medium_id], ray.org);
            // NEE
            nee_p_cache = ray.org;
            Spectrum nee = final_next_event_estimation(scene, ray, ray.org, current_medium_id, bounces, true, vertex_, rng);
            radiance += current_path_throughput * nee * sigma_s;
            // phase function sampling
            never_scatter = false;//?????????
            PhaseFunction phaseFunction = get_phase_function(scene.media[current_medium_id]);
            Vector2 phase_rnd_param_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
            std::optional<Vector3> next_dir_ = sample_phase_function(phaseFunction, -ray.dir, phase_rnd_param_uv);
            
            if (next_dir_) {
                Vector3 next_dir = *next_dir_;
                dir_pdf = pdf_sample_phase(phaseFunction, -ray.dir, next_dir);
                current_path_throughput *= (eval(phaseFunction, -ray.dir, next_dir) / dir_pdf) * sigma_s;
                ray.dir = next_dir;
                multi_trans_pdf = make_const_spectrum(1);
                ray = {ray.org, ray.dir, get_intersection_epsilon(scene), infinity<Real>()};
            } else {
                break;
            }
            
        } else if(vertex_) {
            // nee
            nee_p_cache = ray.org;
            never_scatter = false;
            Spectrum nee = final_next_event_estimation(scene, ray, ray.org, current_medium_id, bounces, false, vertex_, rng);
            radiance += current_path_throughput * nee;
            
            // bsdf sampling
            Vector3 dir_view = -ray.dir;
            Vector2 bsdf_rnd_param_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
            Real bsdf_rnd_param_w = next_pcg32_real<Real>(rng);
            const Material &mat = scene.materials[vertex_->material_id];
            std::optional<BSDFSampleRecord> bsdf_sample_ =
            sample_bsdf(mat,
                        dir_view,
                        *vertex_,
                        scene.texture_pool,
                        bsdf_rnd_param_uv,
                        bsdf_rnd_param_w);
            if (!bsdf_sample_) {
                // BSDF sampling failed. Abort the loop.
                break;
            }
            const BSDFSampleRecord &bsdf_sample = *bsdf_sample_;
            Vector3 dir_bsdf = bsdf_sample.dir_out;
            // Update ray differentials & eta_scale
            if (bsdf_sample.eta == 0) {
                ray_diff.spread = reflect(ray_diff, vertex_->mean_curvature, bsdf_sample.roughness);
            } else {
                ray_diff.spread = refract(ray_diff, vertex_->mean_curvature, bsdf_sample.eta, bsdf_sample.roughness);
                eta_scale /= (bsdf_sample.eta * bsdf_sample.eta);
                current_medium_id = update_medium(ray, *vertex_, current_medium_id);
            }

            dir_pdf = pdf_sample_bsdf(mat, dir_view, dir_bsdf, (*vertex_), scene.texture_pool);
            Spectrum bsdf_f = eval(mat, dir_view, dir_bsdf, (*vertex_), scene.texture_pool);

            if (dir_pdf <= 0) {
                // Numerical issue -- we generated some invalid rays.
                break;
            }
            current_path_throughput *= bsdf_f / dir_pdf;
            ray.dir = dir_bsdf;
            ray = Ray{ray.org, ray.dir, get_intersection_epsilon(scene), infinity<Real>()};
            multi_trans_pdf = make_const_spectrum(1);
        }

        // Russian roulette
        Real rr_prob = 1;
        if (bounces >= scene.options.rr_depth) {
            // rr_prob = min(max(current_path_throughput), Real(0.95));
            rr_prob = min(max((1 / eta_scale) * current_path_throughput), Real(0.95));
            if (next_pcg32_real<Real>(rng) > rr_prob) {
                // Terminate the path
                break;
            } else {
                current_path_throughput /= rr_prob;
            }
        }
        bounces += 1;
    }

    return radiance;
}
