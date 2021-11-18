#ifndef IMAGE_H
#define IMAGE_H

#include <cmath>

extern int width;
extern int height;
extern int numPixels;
extern double GetCeiling(double f);

class Image
{ 
  public:
      unsigned char *buffer;
      double *depth;
      int width, height;
      
      Image();
      void SetPixel(unsigned char *buffer, int row, int col, double r, double g, double b, double s);
};

Image::Image()
{
    // Initialize the width, height, depth of the image
    width = ::width;
    height = ::height;
    depth = new double[numPixels];
    for (int i = 0; i < numPixels; i++)
       depth[i] = -1;

    // Initialize the RGB values for every pixel in the image
    int numPrimaryColors = 3;
    int total = numPrimaryColors*numPixels;
    buffer = new unsigned char[total];
    for (int j = 0; j < total; j++)
       buffer[j] = 0;
}

void Image::SetPixel(unsigned char *buffer, int row, int col, double r, double g, double b, double shading)
{
    // If the color value multipied by the shading is too high, then set it to the max value
    if (1 < r*shading)
    {
       buffer[3*(row*width + col) + 0] = (unsigned char) 255;
    }
    else
    {
       buffer[3*(row*width + col) + 0] = (unsigned char) GetCeiling(255*r*shading);
    }
    if (1 < g*shading)
    {
       buffer[3*(row*width + col) + 1] = (unsigned char) 255;
    }
    else
    {
       buffer[3*(row*width + col) + 1] = (unsigned char) GetCeiling(255*g*shading);
    }
    if (1 < b*shading)
    {
       buffer[3*(row*width + col) + 2] = (unsigned char) 255;
    }
    else
    {
       buffer[3*(row*width + col) + 2] = (unsigned char) GetCeiling(255*b*shading);
    }
}

#endif
