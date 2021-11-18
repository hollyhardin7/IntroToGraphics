#ifndef SPACE_H
#define SPACE_H

#include "matrix.h"
#include "camera.h"
#include <cmath>

extern void Normalize(double *v);
extern void SetCrossProduct(double *a, double *u, double *w);
extern double GetDotProduct(double *u, double *v);

class Space
{
  public:
     Matrix   GetWorldSpace();
     Matrix   GetCameraSpace(Camera c);
     Matrix   GetImageSpace(Matrix cartesian, Matrix camera);
     Matrix   GetScreenSpace(Camera c);
     Matrix   GetDeviceSpace(Matrix cartesian);
};

Matrix Space::GetWorldSpace()
{
   // Identity matrix that will be used to find x, y, and z coordinates
   Matrix m;
   for (int i = 0; i < 4; i++)
      m.A[i][i] = 1;
   return m;
}

Matrix Space::GetCameraSpace(Camera c)
{
   // The vector for the origin's coordinates
   double *o = new double[3];
   for (int i = 0; i < 3; i++)
      o[i] = c.position[i];

   // The basis vector for the Z coordinate
   // W = O-focus by definition
   double *w = new double[3];
   for (int i = 0; i < 3; i++)
      w[i] = o[i] - c.focus[i];
   Normalize(w);

   // The basis vector for the X coordinate
   // U = up x (O-focus) by definition
   double *u = new double[3];
   SetCrossProduct(u, c.up, w);
   Normalize (u);

   // The basis vector for the Y coordinate
   // V = (O-focus) x u by definition
   double *v = new double[3];
   SetCrossProduct(v, w, u);
   Normalize(v);

   // Create and return the matrix of combined vectors
   Matrix m;
   m.PutVectorsInMatrix(u, v, w, o);
   return m;
}

Matrix Space::GetImageSpace(Matrix cartesian, Matrix camera)
{
   // Matrix that will hold camera's transformed coordinates
   Matrix m;

   double temp[4] = { 1, 0, 0, 0 };
   double *temp2 = new double[4];
   double u[4] = { camera.A[0][0], camera.A[1][0], camera.A[2][0], camera.A[3][0] };
   double v[4] = { camera.A[0][1], camera.A[1][1], camera.A[2][1], camera.A[3][1] };
   double w[4] = { camera.A[0][2], camera.A[1][2], camera.A[2][2], camera.A[3][2] };
   double t[4] = { 0-camera.A[0][3], 0-camera.A[1][3], 0-camera.A[2][3], 0-camera.A[3][3] };

   double *d = new double[4];
   SetCrossProduct(d, u, v);
   double denominator = GetDotProduct(d, w);

   double numerator;
   for (int i = 0; i < 3; i++)
   {
       if (i != 0)
       {
           temp[i] = 1;
           temp[i-1] = 0;
       }

       for (int j = 0; j < 3; j++)
       {
           if (j == 0)
           {
               SetCrossProduct(temp2, temp, v);
               numerator = GetDotProduct(temp2, w);
           }
           else if (j == 1)
           {
               SetCrossProduct(temp2, u, temp);
               numerator = GetDotProduct(temp2, w);
           }
           else if (j == 2)
           {
               SetCrossProduct(temp2, u, v);
               numerator = GetDotProduct(temp2, temp);
           }
           m.A[i][j] = numerator / denominator;
       }
   }

   SetCrossProduct(temp2, t, v);
   numerator = GetDotProduct(temp2, w);
   m.A[3][0] = numerator / denominator;

   SetCrossProduct(temp2, u, t);
   numerator = GetDotProduct(temp2, w);
   m.A[3][1] = numerator / denominator;

   SetCrossProduct(temp2, u, v);
   numerator = GetDotProduct(temp2, t);
   m.A[3][2] = numerator / denominator;

   m.A[3][3] = 1;

   Matrix transform = Matrix::ComposeMatrices(cartesian, m);

   return transform;
}


Matrix Space::GetScreenSpace(Camera c)
{
   // Create a screen space by defintion
   Matrix m;
   m.A[0][0] = 1/tan(c.angle/2);
   m.A[1][1] = 1/tan(c.angle/2);
   m.A[2][2] = (c.far+c.near)/(c.far-c.near);
   m.A[2][3] = -1;
   m.A[3][2] = (2*c.far*c.near)/(c.far-c.near);

   // Return the screen space 
   return m;
}

Matrix Space::GetDeviceSpace(Matrix cartesian)
{
   // Create a matrix and apply manipulation by definition
   Matrix m;
   m.A[0][0] = width / 2;
   m.A[1][1] = height / 2;
   m.A[2][2] = 1;
   m.A[3][0] = width / 2;
   m.A[3][1] = height / 2;
   m.A[3][3] = 1;

   // Create the output matrix by multiping it by the cartesian matrix 
   Matrix m2 = Matrix::ComposeMatrices(cartesian, m);

   // Return the device image
   return m2;
}

#endif
