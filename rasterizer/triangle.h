#ifndef TRIANGLE_H
#define TRIANGLE_H

#include <cmath>
#include <algorithm>
#include "lightingparameters.h"
#include "matrix.h"
#include "screen.h"
#include "camera.h"

enum direction { UP, DOWN, RIGHT, ARBITRARY };

extern void SetLightReflectionProportion(double *r, double *l, double *n);
extern void SwapValues(double &x1, double &x2);
extern int WithinBounds(int r, int c);

class Triangle
{
  public:
     // Variables
     direction   d; // Direction of triangle
     double      X[3]; // X value of coordinates
     double      Y[3]; // Y value of coordinates
     double      Z[3]; // Z value of coordinates
     double      slope[2]; // Slopes of triangle (m in y=mx+b)
     double      offset[2]; // Offset of triangle (b in y=mx+b)
     double      colors[3][3]; // RBG values for each coordinate
     double      normals[3][3]; // Normal for each coordinate's color (line that is perpendicular of the triangle)
     double      shading[3]; // Shading intensity for each coordinate
     int         topIndex; // Index of top vertex
     int         middleIndex; // Index of middle vertex 
     int         bottomIndex; // Index of bottom Vertex

     // Getters
     double   GetSlope(double x, double y); // Get the slope of the triangle
     double   GetOffset(double y, double m, double x); // Get the offset of the triangle
     double   GetEndpoint(double x, int i); // Get the endpoint (x in y=mx+b)
     double   GetShading(double cameraPosition[3], int i); // Get the shading of the triangle

     // Setters
     void   SetDirection(); // Set the direction of the triangle
     void   SetCoordinates(int i, double coordinates[3]); // Set the coordinates for a vertex 
     void   SetValuesInSlopeInterceptForm(double x1, double y1, double y2, int index); // Set the values in y=mx+b
     void   SetSlopes(); // Set the slopes
     void   SetOffsets(); // Set the offsets
     void   SetColors(int i, double r, double g, double b); // Set the colors for a vertex
     void   SetNormals(int i, float *normals, vtkIdType *ptIds); // Set the normals for a vertex
     void   SetShading(double cameraPosition[3]); // Set the shading intensity
     void   InterpolateTriangle(Triangle original, double *color, int order[3]); // Method for iterpolating the shading, depth, and colors of a triangle
     void   InterpolateColor(double A, double B, double X, double *FA, double *FB, double *FX); // Method for interpolating the RGB colors
     void   Interpolate(double A, double B, double X, double FA, double FB, double &FX); // Method for interpolating values
     void   SetIndicesGoingUp(Triangle original); // Set the indices for each vertex in a triangle pointing up
     void   SetIndicesGoingDown(Triangle original); // Set the indices for each vertex in a triangle pointing down
     void   UpdateTriangle (Triangle original, int originalIndex, int triangleIndex); // Update the coordinates, shading, and color of the triangle
     void   Rasterize(unsigned char *b, Screen &s); // Taking an imaged described in a vector graphics format and create a raster image
     void   TransformTriangle(Matrix matrix); // Transform the triangle according to the matrix
};

double Triangle::GetSlope(double x, double y)
{
   // If slopes are the same, return infinity
   // Otherwise, we can get slope by finding the difference of rise over run
   if (X[0] == x)
      return INFINITY;
   else
      return ((y - Y[0]) / (x - X[0]));
}

double Triangle::GetOffset(double y, double m, double x)
{
  // If the slope is a horizontal line, return infinity
  // Otherwise, we can get slope by solving b in equation y = mx + b
  if (isinf(m))
    return 0.0;
  else
    return (y - (m * x));
}

double Triangle::GetEndpoint(double row, int i)
{
   // If the slope is a vertical line, return the next value for X
   // Otherwise, we can solve for endpoint by solving x in equation y = mx + b
   if (isinf(slope[i]))
      return X[i+1];
   else
      return ((row - offset[i]) / slope[i]);
}

double Triangle::GetShading(double cameraPosition[3], int i)
{
   // Determine shading according to definition of equation
   LightingParameters lp;
   double ambient = lp.Ka;
   double diffuse = fabs(lp.Kd * GetDotProduct(lp.lightDir, normals[i]));
   double specular;
   double *r = new double[3];
   SetLightReflectionProportion(r, lp.lightDir, normals[i]);
   double v[3] = { cameraPosition[0] - X[i], cameraPosition[1] - Y[i], cameraPosition[2] - Z[i] };
   Normalize(v);
   Normalize(r);
   specular = lp.Ks * pow(GetDotProduct(r, v), lp.alpha);
   if (isnan(specular))
      specular = 0;
   delete [] r;
   return ambient + diffuse + specular;
}

void Triangle::SetDirection()
{
   // Determine the direction the triangle is pointing
   if ((Y[0] == Y[1] && Y[1] < Y[2]) || (Y[0] == Y[2] && Y[2] < Y[1]) || (Y[1] == Y[2] && Y[2] < Y[0]))
      d = UP;
   else if ((Y[0] == Y[1] && Y[2] < Y[1]) || (Y[0] == Y[2] && Y[1] < Y[2]) || (Y[1] == Y[2] && Y[0] < Y[2]))
      d = DOWN;
   else if ((X[0] == X[1] && X[1] < X[2]) || (X[0] == X[2] && X[2] < X[1]) || (X[1] == X[2] && X[2] < X[0]))
      d = RIGHT;
   else
      d = ARBITRARY;
}

void Triangle::SetCoordinates(int i, double coordinates[3])
{
    // Set the coordinates of the triangle
    X[i] = coordinates[0];
    Y[i] = coordinates[1];
    Z[i] = coordinates[2];
}

void Triangle::SetValuesInSlopeInterceptForm(double x1, double y1, double y2, int index)
{
   // Determine the slope and offset given x and y provided by definition
   double x = x1;
   double y = y1;
   double m = GetSlope(x, y);
   double b = GetOffset(y, m, x);
   if (isinf(m))
      x = x1;
   else
      x = ((y2 - b) / m);
   X[index] = x;
   Y[index] = y2;
}

void Triangle::SetSlopes()
{
  // Set the slope of the triangle
  for (int i = 0; i < 2; i++)
     slope[i] = GetSlope(X[i+1], Y[i+1]);
}

void Triangle::SetOffsets()
{
  // Set the offset of the triangle
  for (int i = 0; i < 2; i++)
     offset[i] = GetOffset(Y[i+1], slope[i], X[i+1]);
}

void Triangle::SetColors(int i, double r, double g, double b)
{   
    // Set the RGB values of the triangle
    colors[i][0] = r;
    colors[i][1] = g;
    colors[i][2] = b;
}

void Triangle::SetNormals(int i, float *normals, vtkIdType *ptIds)
{
   // Set the normals for each color
   for (int j = 0; j < 3; j++)
      (*this).normals[i][j] = normals[3*ptIds[i]+j];
}

void Triangle::SetShading(double cameraPosition[3])
{
   // Set the shading for each triangle
   for (int i = 0; i < 3; i++)
      shading[i] = GetShading(cameraPosition, i);
}

void Triangle::InterpolateTriangle(Triangle original, double *color, int order[3])
{  
   // Interpolate the shading, depth, and colors
   Interpolate(original.Y[order[0]], original.Y[order[1]], original.Y[order[2]], original.shading[order[0]], original.shading[order[1]], shading[1]);
   Interpolate(original.Y[order[0]], original.Y[order[1]], original.Y[order[2]], original.Z[order[0]], original.Z[order[1]], Z[1]);
   InterpolateColor(original.Y[order[0]], original.Y[order[1]], original.Y[order[2]], original.colors[order[0]], original.colors[order[1]], color);
}

void Triangle::InterpolateColor(double A, double B, double X, double FA[3], double FB[3], double FX[3])
{
   // Interpolate each color
   for (int i = 0; i < 3; i++)
      Interpolate(A, B, X, FA[i], FB[i], FX[i]);
}

void Triangle::Interpolate(double A, double B, double X, double FA, double FB, double &FX)
{
   // Interpolate the value by definition
   double t;
   if (A != B)
      t = ((X - A) / (B - A));
   else
      t = NAN;
   if (isnan(t))
      FX = FA;
   else
      FX = (FA + (t * (FB - FA)));
}

void Triangle::SetIndicesGoingUp(Triangle original)
{
    // Find the index of the top vertex
    topIndex = 0;
    int possibilities[2] = { 1, 2 };
    if (original.Y[topIndex] < original.Y[1])
    {
        topIndex = 1;
        possibilities[0] = 0;
        possibilities[1] = 2;
    }
    if (original.Y[topIndex] < original.Y[2])
    {
        topIndex = 2;
        possibilities[0] = 0;
        possibilities[1] = 1;
    }

    // Find the indices of the bottom and middle vertices
    bottomIndex = possibilities[0];
    middleIndex = possibilities[1];
    if (original.Y[bottomIndex] < original.Y[middleIndex])
    {
        bottomIndex = possibilities[0];
        middleIndex = possibilities[1];
    }
    else
    {
        bottomIndex = possibilities[1];
        middleIndex = possibilities[0];
    }

   // Update the triangles coordinates, colors, and shading
   UpdateTriangle(original, topIndex, 0);
   if (original.Y[bottomIndex] == original.Y[middleIndex])
   {
      UpdateTriangle(original, middleIndex, 1);
      UpdateTriangle(original, bottomIndex, 2);
   }
   else
   {
      UpdateTriangle(original, middleIndex, 2);

      SetValuesInSlopeInterceptForm(original.X[bottomIndex], original.Y[bottomIndex], original.Y[middleIndex], 1);
      double *color = new double[3];
      if (original.Y[topIndex] > original.Y[middleIndex])
      {
         int order[3] = { bottomIndex, topIndex, middleIndex };
         InterpolateTriangle(original, color, order);
      }
      else
      {
         int order[3] = { bottomIndex, middleIndex, topIndex };
         InterpolateTriangle(original, color, order);
      }
      SetColors(1, color[0], color[1], color[2]);
      delete [] color;
   }

   // Switch x values if bottom-left is bigger than bottom-right
   if (X[1] > X[2])
   {
      SwapValues(X[1], X[2]);
      SwapValues(Z[1], Z[2]);
      SwapValues(shading[1], shading[2]);
      for (int i = 0; i < 3; i++)
         SwapValues(colors[1][i], colors[2][i]);
   }
}

void Triangle::SetIndicesGoingDown(Triangle original)
{
    // Find the index of the bottom vertex
    bottomIndex = 0;
    int possibilities[2] = { 1, 2 };
    if (original.Y[bottomIndex] > original.Y[1])
    {
        bottomIndex = 1;
        possibilities[0] = 0;
        possibilities[1] = 2;
    }
    if (original.Y[bottomIndex] > original.Y[2])
    {
        bottomIndex = 2;
        possibilities[0] = 0;
        possibilities[1] = 1;
    }

    // Find the indices of the top and middle vertices
    topIndex = possibilities[0];
    middleIndex = possibilities[1];
    if (original.Y[topIndex] < original.Y[middleIndex])
    {
        topIndex = possibilities[0];
        middleIndex = possibilities[1];
    }
    else
    {
        topIndex = possibilities[1];
        middleIndex = possibilities[0];
    }

    // Save the vertices as top = 0, top-left = 1, and top-right = 2
    UpdateTriangle(original, bottomIndex, 0);

    if (original.Y[topIndex] == original.Y[middleIndex])
    {
       UpdateTriangle(original, topIndex, 1);
       UpdateTriangle(original, middleIndex, 2);
    }
    else
    {
       UpdateTriangle(original, topIndex, 1);

       SetValuesInSlopeInterceptForm(original.X[middleIndex], original.Y[middleIndex], original.Y[topIndex], 2);

       double color[3] = { 0, 0, 0 };
        if (original.Y[topIndex] > original.Y[middleIndex])
        {
           int order[3] = { bottomIndex, topIndex, middleIndex };
           InterpolateTriangle(original, color, order);
        }
        else
        {
           int order[3] = { bottomIndex, middleIndex, topIndex };
           InterpolateTriangle(original, color, order);
        }
        SetColors(2, color[0], color[1], color[2]);
    }

    // Switch x, depth, and color values if top-left is bigger than top-right
    if (X[1] > X[2])
    {
      SwapValues(X[1], X[2]);
      SwapValues(Z[1], Z[2]);
      SwapValues(shading[1], shading[2]);
      for (int i = 0; i < 3; i++)
         SwapValues(colors[1][i], colors[2][i]);
    }
}

void Triangle::UpdateTriangle(Triangle original, int originalIndex, int triangleIndex)
{
   // Update the coordinates, shading, and color of the triangle
   double coordinates[3] = { original.X[originalIndex], original.Y[originalIndex], original.Z[originalIndex] };
   SetCoordinates(triangleIndex, coordinates);
   shading[triangleIndex] = original.shading[originalIndex];
   SetColors(triangleIndex, original.colors[originalIndex][0], original.colors[originalIndex][1], original.colors[originalIndex][2]);
}

void Triangle::Rasterize(unsigned char *b, Screen &s)
{
   // Set the slopes and offsets for the triangle to be rasterized
   SetSlopes();
   SetOffsets();

   double rowMin = ceil((*std::min_element(Y, Y+3)));
   double rowMax = floor((*std::max_element(Y, Y+3)));
   for (int row = rowMin; row <= rowMax; row++)
   {
       double *leftColor = new double[3];
       double *rightColor = new double[3];
       double leftDepth;
       double rightDepth;
       double leftShade;
       double rightShade;

       InterpolateColor(Y[0], Y[1], row, colors[0], colors[1], leftColor);
       InterpolateColor(Y[0], Y[2], row, colors[0], colors[2], rightColor);
       Interpolate(Y[0], Y[1], row, Z[0], Z[1], leftDepth);
       Interpolate(Y[0], Y[2], row, Z[0], Z[2], rightDepth);
       Interpolate(Y[0], Y[1], row, shading[0], shading[1], leftShade);
       Interpolate(Y[0], Y[1], row, shading[0], shading[2], rightShade);

       double leftEnd = GetEndpoint(row, 0);
       double rightEnd = GetEndpoint(row, 1);

       for (int col = ceil(leftEnd); col <= floor(rightEnd); col++)
       {
          if (WithinBounds(row, col))
          {
             double *color = new double[3];
             InterpolateColor(leftEnd, rightEnd, col, leftColor, rightColor, color);

             double depth;
             Interpolate(leftEnd, rightEnd, col, leftDepth, rightDepth, depth);

             double shade;
             Interpolate(leftEnd, rightEnd, col, leftShade, rightShade, shade);

             if (depth > s.depth[row*width + col])
             {
                s.SetPixel(b, row, col, color[0], color[1], color[2], shade);
                s.depth[row*width + col] = depth;
             }

             delete [] color;
          }
       }

       delete [] leftColor;
       delete [] rightColor;
   }
}

void Triangle::TransformTriangle(Matrix matrix)
{
   // Apply the transformation matrix for every triangle coordinate
   for (int i = 0; i < 3; i++)
   {
      double points_in[4] = { X[i], Y[i], Z[i], 1 };
      double points_out[4] = { 0, 0, 0, 0 };
      matrix.TransformPoint(points_in, points_out);

          double w = points_out[3];
          
          for (int k = 0; k < 3; k++)
             points_out[k] = points_out[k] / w;

      X[i] = points_out[0] / points_out[3];
      Y[i] = points_out[1] / points_out[3];
      Z[i] = points_out[2] / points_out[3];
   }
}

#endif
