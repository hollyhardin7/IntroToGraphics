#ifndef CAMERA_H
#define CAMERA_H

extern double SineParameterize(int imageIndex, int maxImages, int ramp);

class Camera
{ 
  public:
    double   near, far;
    double   angle; 
    double   position[3]; 
    double   focus[3];
    double   up[3];

    Camera(int i, int max);
};

// Provided code
Camera::Camera(int imageIndex, int maxImages)
{
    double t = SineParameterize(imageIndex, maxImages, maxImages/10);
    near = 5;
    far = 200;
    angle = M_PI/6;
    position[0] = 40*sin(2*M_PI*t);
    position[1] = 40*cos(2*M_PI*t);
    position[2] = 40;
    focus[0] = 0;
    focus[1] = 0;
    focus[2] = 0;
    up[0] = 0;
    up[1] = 1;
    up[2] = 0;
}

#endif
