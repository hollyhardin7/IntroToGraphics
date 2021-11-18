#ifndef SCREEN_H
#define SCREEN_H

#include <cmath>

extern int width;
extern int height;
extern double GetCeiling(double f);

class Screen
{ 
  public:
      unsigned char *buffer;
      double *depth;
      int width, height;
      
      Screen(unsigned char *buffer);
      void SetPixel(unsigned char *buffer, int row, int col, double r, double g, double b, double s);
};

Screen::Screen(unsigned char *b)
{
    // Initialize the width, height, depth, number of pixels, and buffer of the image
    width = ::width;
    height = ::height;
    int numPrimaryColors = 3;
    int numPixels = width*height*numPrimaryColors;
    depth = new double[numPixels];
    for (int i = 0; i < numPixels; i++)
       depth[i] = -1;
    buffer = new unsigned char[numPixels];
    for (int i = 0; i < numPixels; i++)
       buffer[i] = 0;
    buffer = b;
}

void Screen::SetPixel(unsigned char *buffer, int row, int col, double r, double g, double b, double shading)
{
    // If the color value multipied by the shading is too high, then set it to the max value
    if (1 < r*shading)
       buffer[3*(row*width + col) + 0] = (unsigned char) 255;
    else
       buffer[3*(row*width + col) + 0] = (unsigned char) GetCeiling(255*r*shading);
    if (1 < g*shading)
       buffer[3*(row*width + col) + 1] = (unsigned char) 255;
    else
       buffer[3*(row*width + col) + 1] = (unsigned char) GetCeiling(255*g*shading);
    if (1 < b*shading)
       buffer[3*(row*width + col) + 2] = (unsigned char) 255;
    else
       buffer[3*(row*width + col) + 2] = (unsigned char) GetCeiling(255*b*shading);
}

#endif
