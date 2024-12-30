# Dielectric Microfacet BSDF Model

## Overview

This project is an implementation of a **Dielectric Microfacet Bidirectional Scattering Distribution Function (BSDF)** using the **GGX normal distribution function** for physically-based rendering. The BSDF accounts for both **ideal reflection** and **refraction** using **Fresnel equations** and supports varying roughness parameters for realistic rendering of dielectric materials such as glass. 

<p align="center">
  <img src="https://github.com/user-attachments/assets/88f16cfa-d1f0-4b27-8f43-962a0afbe64f" alt="Roughness 0.005" width="45%" />
  <img src="https://github.com/user-attachments/assets/255f38b1-e14e-40f6-8085-871ca854db16" alt="Roughness 0.5" width="45%" />
</p>

<p align="center">
  <b>Figure 1:</b> Rendered images with Roughness values - Left: 0.005 (smooth surface), Right: 0.5 (rough surface).
</p>

## Features

- **Ideal Dielectric Model**: Implements perfect reflection and transmission based on Fresnel equations.  
- **GGX Microfacet Distribution**: Samples microfacet normals to simulate rough surfaces.  
- **Roughness Parameterization**: Adjusts the model to simulate surfaces ranging from smooth to rough.  
- **Physically Accurate Fresnel Effect**: Computes Fresnel reflectance and transmittance for varying incident angles.  
- **Efficient Sampling**: Includes optimized sampling for reflection and refraction directions.  



## Implementation Details

### Scene Setup

- **Modified Cornell Box**: Includes a glass bunny with:
  - Green back wall and red floor to distinguish reflection and transmission effects.  
- **Index of Refraction**: Fixed at **1.5**, representative of glass.  

### Key Components

#### 1. Ideal Reflection and Transmission
- **Reflection**: Light reflects at the same angle as incidence.  
- **Transmission**: Handles refraction with Snell’s law and accounts for orientation changes.  
- **Sampling**: Fresnel probabilities determine reflection vs transmission.  

```cpp
tuple<Vec3f, float> sample(Vec3f wo_world, Vec3f n, Vec3f u) const {
    Vec3f wo = normalize(to_local(wo_world, n));
    Vec3f wi = mirror_reflect(-wo, m);
    Vec3f wi_world = from_local(normalize(wi), n);
    return {wi_world, 1.0f};
}
```
<p align="center">
  <img src="https://github.com/user-attachments/assets/f9929dfd-4f5c-4933-8c77-4b17459b3058" alt="Reflection" width="45%" />
  <img src="https://github.com/user-attachments/assets/d254fc65-36bd-4eae-8ff7-040bcbaa3420" alt="Transmission" width="45%" />
</p>

<p align="center">
  <b>Figure 2:</b> Ideal Reflection (left) and Ideal Transmission (right) results using the BSDF model.
</p>

#### 2. GGX Microfacet Distribution
- Models the microfacet distribution using a roughness parameter `g`.  
- Implements **Shadow-Masking** functions to ensure energy conservation.  

```cpp
float D(Vec3f m) const {
    if (m.z < 0.0f) return 0.0f;
    float cos_thetaM2 = m.z * m.z;
    float tan_thetaM2 = (1.0f - cos_thetaM2) / cos_thetaM2;
    float alpha2 = roughness * roughness;
    float sqrt_denom = alpha2 + tan_thetaM2;
    return alpha2 * INV_PI / 
           (cos_thetaM2 * cos_thetaM2 * sqrt_denom * sqrt_denom);
}
```

The GGX Distribution Function determines the probability of microfacet orientations based on the roughness parameter, producing smoother or rougher surfaces depending on its value.

#### 3. Fresnel Equations
- Computes reflectance based on the angle of incidence and relative indices of refraction.
- Determines probability of reflection vs transmission.

```cpp
float Fresnel(float cos_thetaI) const {
    cos_thetaI = clamp(cos_thetaI, -1.0f, 1.0f);
    bool entering = cos_thetaI > 0.0f;
    float ei = eta, et = 1.0f;
    if (!entering) swap(ei, et);

    float sin_thetaT = ei / et * sqrt(max(0.0f, 1.0f - cos_thetaI * cos_thetaI));
    if (sin_thetaT >= 1.0f) return 1.0f;

    float cos_thetaT = sqrt(max(0.0f, 1.0f - sin_thetaT * sin_thetaT));

    float Rs = (et * cos_thetaI - ei * cos_thetaT) / 
               (et * cos_thetaI + ei * cos_thetaT);
    float Rp = (ei * cos_thetaI - et * cos_thetaT) / 
               (ei * cos_thetaI + et * cos_thetaT);

    return 0.5f * (Rs * Rs + Rp * Rp);
}
```

The Fresnel equations calculate the reflectance based on the angle of incidence, producing realistic behavior where reflection increases at grazing angles and transmission dominates at normal incidence.

#### 4. Shadow-Masking Functions
- Ensures that microfacets are not shadowed or masked, preserving energy conservation in the model.

```cpp
float G(Vec3f wo, Vec3f wi, Vec3f m) const {
    return G1(wo, m) * G1(wi, m);
}

float G1(Vec3f w, Vec3f m) const {
    if (dot(w, m) * w.z <= 0.0f) return 0.0f;

    float tan_thetaW2 = (1.0f - w.z * w.z) / (w.z * w.z);
    if (tan_thetaW2 <= 0.0f) return 1.0f;
    float tan_thetaW = sqrt(tan_thetaW2);
    float root = roughness * tan_thetaW;
    return 2.0f / (1.0f + sqrt(1.0f + root * root));
}
```
These functions account for geometric occlusion, which occurs when facets block light from reaching others.

#### 5. Sampling the BSDF
- The sampling function generates microfacet normals using the GGX distribution and evaluates reflection or transmission based on Fresnel probabilities.

```cpp
tuple<Vec3f, float> sample(Vec3f wo_world, Vec3f n, Vec3f u) const {
    Vec3f wo = normalize(to_local(wo_world, n));

    const auto [m, pdf_m] = sample_normal(Vec2f{u[1], u[2]});
    float cos_thetaO = dot(wo, m);
    float F = Fresnel(cos_thetaO);

    if (u[0] < F) {
        Vec3f wi = mirror_reflect(-wo, m);
        Vec3f wi_world = from_local(normalize(wi), n);
        return {wi_world, F * pdf_m};
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
    if (entering) cos_thetaT *= -1;

    Vec3f wt = {etaEff * -wo.x, etaEff * -wo.y, cos_thetaT};
    Vec3f wt_world = from_local(normalize(wt), n);

    return {wt_world, (1.0f - F) * pdf_m};
}
```

This function calculates probabilities for both reflection and transmission and determines which behavior occurs based on random sampling and Fresnel equations.



## References

- **PBRT Textbook** - Concepts for BSDF and Fresnel equations.  
- **Course Material** - Provided project guidelines and sample implementations.  



## Acknowledgments

- **Professor's Guidance** - Inspired by the CS190I: Offline Rendering Final Project requirements.  
- **PBRT Textbook** - Key references for BSDF theory and implementations.  



## License

This project is licensed under the **MIT License**.  
See the [LICENSE](LICENSE) file for details.  
