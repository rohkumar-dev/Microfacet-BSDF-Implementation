#pragma once
#include "common.h"
#include "math_helpers.h"
#include <algorithm>
#include <cmath>
#include <iostream>
// #include <winuser.h>
using std::tuple;
using std::clamp;
using std::cout;
using std::endl;
using std::swap;

#define PI 3.141592653589
#define INV_PI 0.31830989
#define EXP 2.718281828459045

namespace muni {
struct Lambertian {
    Vec3f albedo;

    /** Evaluates the BRDF for the Lambertian material.
      \return The BRDF (fr) value.
  */
    Vec3f eval() const {
        return albedo * (float) INV_PI; 
     }

    /** Samples the BRDF for the Lambertian material.
      \param[in] normal The normal of the surface.
      \param[in] u A random number in (0,1)^2.
      \return A tuple containing the sampled direction in world space and the PDF.
    */
    std::tuple<Vec3f, float> sample(Vec3f normal, Vec2f u) const {
      float r = sqrt(std::max(0.0f, 1.0f - u[0] * u[0]));
      float phi = 2 * PI * u[1];

      Vec3f v = normalize(Vec3f(r * std::cos(phi), r * std::sin(phi), u[0]));
      Vec3f dir = from_local(v, normal);

      return {dir, 1.0/(2.0 * PI)};
    }
    /** Computes the PDF for the Lambertian material.
      \param[in] wo The outgoing direction in world space.
      \param[in] wi The light incident direction in world space.
      \return The PDF value.
    */
    float pdf(Vec3f wo_world, Vec3f wi_world, Vec3f normal) const {
      if (dot(wo_world, normal) < 0.0f) return 0;
      return INV_PI / 4.0f;
    }

};


struct Dielectric {
  float eta, roughness;

  float Fresnel(float cos_thetaI) const {
    cos_thetaI = std::clamp(cos_thetaI, -1.0f, 1.0f);
    bool entering = cos_thetaI > 0.0f;
    float ei = eta, et = 1.0f;
    if (!entering)
      std::swap(ei, et);

    float sin_thetaT = ei / et * sqrt(std::max(0.0f, 1.0f - cos_thetaI * cos_thetaI));
    if (sin_thetaT >= 1.0f)
      return 1.0f;

    cos_thetaI = abs(cos_thetaI);
    float cos_thetaT = sqrt(std::max(0.0f, 1.0f - sin_thetaT * sin_thetaT));

    float Rs = (et * cos_thetaI - ei * cos_thetaT) / (et * cos_thetaI + ei * cos_thetaT);
    float Rp = (ei * cos_thetaI - et * cos_thetaT) / (ei * cos_thetaI + et * cos_thetaT);

    return 0.5f * (Rs * Rs + Rp * Rp);
  }

  float D(Vec3f h) const {   
    if (h.z < 0.0f) return 0.0f;

    float cos_thetaH_sq = h.z * h.z;
    float tan_thetaH_sq = (1.0f - cos_thetaH_sq) / cos_thetaH_sq;
    float alpha_sq = roughness * roughness;

    float sqrt_denom = alpha_sq + tan_thetaH_sq;

    return alpha_sq * INV_PI / (cos_thetaH_sq * cos_thetaH_sq * sqrt_denom * sqrt_denom);
  }

  float G(Vec3f wo, Vec3f wi, Vec3f h) const {
    return G1(wo, h) * G1(wi, h);
  }

  float G1(Vec3f w, Vec3f h) const {
    if (dot(w, h) * w.z <= 0.0f)
      return 0.0f;

    float tan_thetaW_sq = (1.0f - w.z * w.z) / (w.z * w.z);
    if (tan_thetaW_sq <= 0.0f)
      return 1.0f;

    float tan_thetaW = sqrt(tan_thetaW_sq);
    float root = roughness * tan_thetaW;

    return 2.0f / (1.0f + sqrt(1.0f + root * root));
  }

  float pdf(Vec3f h) const {
    return D(h) * h.z;
  }

  float eval(Vec3f wo_world, Vec3f wi_world, Vec3f n) const {
    Vec3f wo = normalize(to_local(wo_world, n)), wi = normalize(to_local(wi_world, n));
    bool reflect = (wi.z * wo.z) > 0.0f;
    
    float etaT = (wo.z < 0.0f) ? 1.0f : eta;
    float etaI = (wo.z < 0.0f) ? eta : 1.0f;

    float reflect_coeff = (wi.z > 0.0f) ? 1.0f : -1.0f;
    Vec3f h = (reflect) ? normalize(wi + wo) * reflect_coeff : normalize(wi * etaI + wo * etaT);

    if (wo.z < 0.0f) wo *= -1.0f;
    if (wi.z < 0.0f) wi *= -1.0f;

    float Dr = D(h);
    float Fr = Fresnel(dot(wi, h));
    float Gr = G(wo, wi, h);

    if (reflect) {
      float res = Fr * Dr * Gr / (4.0f * wi.z * wo.z);
      // spdlog::info("h: {}, Fr: {}, Dr: {}, Gr: {}, eval: {}", h, Fr, Dr, Gr, res);
      return res;
    }

    float sqrt_denom = etaI * dot(wi, h) + etaT * dot(wo, h);

    float res = abs((1.0f - Fr) * Dr * Gr * etaT * etaT * dot(wi, h) * dot(wo, h) / (wi.z * wo.z * sqrt_denom * sqrt_denom)); 
    // spdlog::info("h: {}, Fr: {}, Dr: {}, Gr: {}, eval: {}", h, Fr, Dr, Gr, res);
    return res;
  }

  tuple<Vec3f, float> sample_normal(Vec2f u) const {
    float thetaM = atan(roughness * sqrt(u.x) / sqrt(1.0f - u.x));
    float phiM = 2 * PI * u.y;

    Vec3f m = normalize(Vec3f{sin(thetaM) * sin(phiM), sin(thetaM) * cos(phiM), cos(thetaM)});
    float pdf_m = pdf(m);

    return {m, pdf_m};
  }


  tuple<Vec3f, float> sample(Vec3f wo_world, Vec3f n, Vec3f u) const {
    Vec3f wo = normalize(to_local(wo_world, n));

    auto [m, pdf_m] = sample_normal(Vec2f{u.y, u.z});
    // spdlog::info("m: {}, pdf_m: {}", m, pdf_m);
  
    float cos_thetaO = dot(wo, m);
    float F = Fresnel(cos_thetaO);
    // spdlog::info("Fresnel_m: {}, Fresnel: {}", F, Fresnel(wo.z));

    if (u.x < F) {
      Vec3f wi = mirror_reflect(-wo, m);
      Vec3f wi_world = from_local(normalize(wi), n);
      return {wi_world, F};
    }

    bool entering = wo.z > 0.0f;
    float ei = eta, et = 1.0f;
    if (!entering) swap(ei, et);
    
    float sin_thetaO2 = 1.0f - cos_thetaO * cos_thetaO;
    float etaEff = ei / et;

    float sin_thetaT2 = etaEff * etaEff * sin_thetaO2;

    if (sin_thetaT2 >= 1.0f)
      return {Vec3f{}, 0.0f};

    float cos_thetaT = sqrt(1.0f - sin_thetaT2);

    if (entering)
      cos_thetaT *= -1;

    Vec3f wt = {etaEff * -wo.x, etaEff * -wo.y, cos_thetaT};

    // spdlog::info("wo: {}, wt: {}", wo, wt);

    Vec3f wt_world = from_local(normalize(wt), n);

    return {wt_world, (1.0f - F)};
  }



};


}  // namespace muni
