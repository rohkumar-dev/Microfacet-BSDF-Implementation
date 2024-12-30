#include "material.h"
#include "muni/camera.h"
#include "muni/common.h"
#include "muni/image.h"
#include "muni/material.h"
#include "muni/math_helpers.h"
#include "muni/obj_loader.h"
#include "muni/ray_tracer.h"
#include "muni/sampler.h"
#include "muni/scenes/box.h"
#include "muni/triangle.h"
#include "ray_tracer.h"
#include "spdlog/spdlog.h"
#include "triangle.h"
#include <cmath>
#include <iostream>
#include <thread>
#include <vector>
#include <string>

#define PI 3.141592653589793
#define sqrt(a) (pow(a, 0.5))
#define norm(a) (sqrt(a.x * a.x + a.y * a.y + a.z * a.z))
#define normSquared(a) (a.x * a.x + a.y * a.y + a.z * a.z)

using namespace muni;
using std::vector;
using std::get;
using std::cout;
using std::endl;

RayTracer::Octree octree{};

/** Offset the ray origin to avoid self-intersection.
    \param[in] ray_pos The original ray origin.
    \param[in] normal The normal of the surface at the hit point.
    \return The offset ray origin.
*/
Vec3f offset_ray_origin(Vec3f ray_pos, Vec3f normal) {
    return ray_pos + EPS * normal;
}

/** Check if the triangle is an emitter.
    \param[in] tri The triangle to check
    \return True if the triangle is an emitter, false otherwise.
*/
bool is_emitter(const Triangle &tri) { return tri.emission != Vec3f{0.0f}; }

/** Evaluate the radiance of the area light. We **do not** check whether the hit
 point is on the light source, so make sure
 *  the hit point is on the light source before calling this function.
    \param[in] light_dir The **outgoing** direction from the light source to the
 scene. \return The radiance of the light source.
*/
Vec3f eval_area_light(const Vec3f light_dir) {
    if (dot(light_dir, BoxScene::light_normal) > 0.0f)
        return BoxScene::light_color;
    return Vec3f{0.0f};
}

/** Sample a point on the area light with a uniform distribution.
    \param[in] samples A 2D uniform random sample.
    \return A tuple containing the sampled position, the normal of the light
 source, and the PDF value.
*/
std::tuple<Vec3f, Vec3f, float> sample_area_light(Vec2f samples) {
    Vec3f pos = {BoxScene::light_x + samples[0] * BoxScene::light_len_x, BoxScene::light_y + samples[1] * BoxScene::light_len_y, BoxScene::light_z};
    Vec3f normal = BoxScene::light_normal;
    float pdf = BoxScene::inv_light_area;
    return {pos, normal, pdf};
}

Vec3f shade_with_light_sampling(Triangle tri, Vec3f p, Vec3f wo) {
    Vec3f L_dir{0.0f}, L_ind{0.0f};

    // Contribution from the light source
    const auto [light_pos, light_normal, pdf_light] = sample_area_light(UniformSampler::next2d());
    Vec3f wi1 = light_pos - p;
    float dist_to_light_squared = normSquared(wi1);
    wi1 = normalize(wi1);
    const auto [hit1, t1, nearest_tri1] = RayTracer::closest_hit(p, wi1, octree, BoxScene::triangles);

    bool tri_contains_lambertian = std::holds_alternative<Lambertian>(BoxScene::materials[tri.material_id]);

    if (is_emitter(nearest_tri1)) {
        Vec3f Li = eval_area_light(-wi1);
        float cos = std::max(0.0f, dot(normalize(tri.face_normal), wi1));
        float cos_prime = std::max(dot(-wi1, normalize(light_normal)), 0.0f);
        
        if (tri_contains_lambertian) {
            Vec3f fr = get<Lambertian>(BoxScene::materials[tri.material_id]).eval();
            L_dir = Li * fr * cos / (pdf_light * dist_to_light_squared / cos_prime);
        }
        
    }

    const float p_rr = 0.8f;
    if (UniformSampler::next1d() > p_rr) return L_dir;

    if (tri_contains_lambertian) {
        Lambertian material = std::get<Lambertian>(BoxScene::materials[tri.material_id]);
        const auto [wi2, pdf_wi] = material.sample(tri.face_normal, UniformSampler::next2d());
        const auto [hit2, t2, nearest_tri2] = RayTracer::closest_hit(p, wi2, octree, BoxScene::triangles);

        Vec3f fr = material.eval();

        if (hit2 && !is_emitter(nearest_tri2)) {
            float cos = std::max(dot(normalize(tri.face_normal), normalize(wi2)), 0.0f);
            Vec3f q = offset_ray_origin(p + normalize(wi2) * t2, nearest_tri2.face_normal);

            L_ind = shade_with_light_sampling(nearest_tri2, q, -wi2) * fr * cos / p_rr / pdf_wi;
            // spdlog::info("Lambertian: {}", L_ind);
        }
    } else {
        Dielectric material = get<Dielectric>(BoxScene::materials[tri.material_id]);
        const auto [wi2, pdf_wi] = material.sample(wo, tri.face_normal, UniformSampler::next3d());
        const auto [hit2, t2, nearest_tri2] = RayTracer::closest_hit(p, wi2, octree, BoxScene::triangles);

        float fr = material.eval(wo, wi2, tri.face_normal);

        if (hit2 && !is_emitter(nearest_tri2) && pdf_wi > 0.0f) {
            float cos = abs(dot(normalize(tri.face_normal), normalize(wi2)));
            Vec3f q = offset_ray_origin(p + normalize(wi2) * t2, nearest_tri2.face_normal);
            L_ind = shade_with_light_sampling(nearest_tri2, q, -wi2) * fr * cos / p_rr / pdf_wi;
            // spdlog::info("Sampled Direction: {}", fr);
            // spdlog::info("L_ind: {}, rec: {}, fr: {}, cos: {}, pdf_wi: {}", L_ind, rec, fr, abs(cos), pdf_wi);

        }
    }

    return L_dir + L_ind;
}

Vec3f path_tracing_with_light_sampling(Vec3f ray_pos, Vec3f ray_dir) {
    const auto [is_ray_hit, t_min, nearest_tri] =
        RayTracer::closest_hit(ray_pos, ray_dir, octree, BoxScene::triangles);
    if (!is_ray_hit) return Vec3f{0.0f};
    const Vec3f hit_position = ray_pos + t_min * ray_dir;
    if (is_emitter(nearest_tri)) return eval_area_light(-ray_dir);

    return shade_with_light_sampling(nearest_tri, hit_position, -ray_dir);
}

// Define the function to be executed in each thread
void renderRow(int row, int max_spp, const Camera& camera, Image& image) {
    for (int x = 0; x < image.width; x++) {
        image(x, row) = Vec3f{0.0f};
        for (int sample = 0; sample < max_spp; sample++) {
            const float u = (x + UniformSampler::next1d()) / image.width;
            const float v = (row + UniformSampler::next1d()) / image.height;
            Vec3f ray_direction = camera.generate_ray(u, (1.0f - v));
            image(x, row) += clamp(path_tracing_with_light_sampling(camera.position, ray_direction), Vec3f(0.0f), Vec3f(50.0f));
        }
        image(x, row) /= static_cast<float>(max_spp);
    }

    if (row % 25 == 0)
        spdlog::info("Finished row {}", row);
}

int main(int argc, char **argv) {


    spdlog::info("\n"
                 "----------------------------------------------\n"
                 "Welcome to CS 190I Assignment 4: Microfacet Materials\n"
                 "----------------------------------------------");
    // const unsigned int max_spp = 64;
    const unsigned int image_width = 1080;
    const unsigned int image_height = 1080;
    // Some prepereations
    Image image{.width = image_width,
                .height = image_height,
                .pixels = std::vector<Vec3f>(image_width * image_height)};
    Camera camera{.vertical_field_of_view = 38.6f,
                  .aspect = static_cast<float>(image_width) / image_height,
                  .focal_distance = 0.8f,
                  .position = Vec3f{0.278f, 0.8f, 0.2744f},
                  .view_direction = Vec3f{0.0f, -1.0f, 0.0f},
                  .up_direction = Vec3f{0.0f, 0.0f, 1.0f},
                  .right_direction = Vec3f{-1.0f, 0.0f, 0.0f}};
    camera.init();
    UniformSampler::init(190);


    // =============================================================================================
    // Change the material ID after you have implemented the Microfacet BRDF
    // Diffuse
    // const int bunny_material_id = 0;
    // Glass
    const int bunny_material_id = 5;
    
    // Load the scene
    // If program can't find the bunny.obj file, use xmake run -w . or move the bunny.obj file to the 
    // same directory as the executable file.
    const std::string obj_path = "./bunny.obj";
    std::vector<Triangle> obj_triangles = load_obj(obj_path, bunny_material_id);
    BoxScene::triangles.insert(BoxScene::triangles.end(),
                               std::make_move_iterator(obj_triangles.begin()),
                               std::make_move_iterator(obj_triangles.end()));
    octree.build_octree(BoxScene::triangles);

    int num_threads = std::thread::hardware_concurrency();
    std::vector<std::thread> threads;
    spdlog::info("Found {} threads", num_threads);

    // =============================================================================================
    // Path Tracing with light sampling
    std::vector<int> max_spps{512};
    for (int max_spp : max_spps) {
        spdlog::info("Path Tracing with light sampling: rendering started!");
        threads.clear();
        for (int y = 0; y < image.height; y++)
            threads.emplace_back(renderRow, y, max_spp, std::ref(camera), std::ref(image));

        for (auto& thread : threads)
            thread.join();

        spdlog::info("Path Tracing with light sampling: Rendering finished!");
        image.save_with_tonemapping("./path_tracing_with_light_sampling" + std::to_string(max_spp) + ".png");
    }

    // =============================================================================================
    return 0;
}
