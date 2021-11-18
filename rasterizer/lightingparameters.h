#ifndef LIGHTINGPARAMETERS_H
#define LIGHTINGPARAMETERS_H

struct LightingParameters
{
    double lightDir[3];  // The direction of the light source
    double Ka;           // The coefficient for ambient lighting
    double Kd;           // The coefficient for diffuse lighting
    double Ks;           // The coefficient for specular lighting
    double alpha;        // The exponent term for specular lighting

    // Initialize the values of a light's values
    LightingParameters(void)
    {
         lightDir[0] = -0.6;
         lightDir[1] = 0;
         lightDir[2] = -0.8;
         Ka = 0.3;
         Kd = 0.7;
         Ks = 2.3;
         alpha = 2.5;
    };
};

#endif
