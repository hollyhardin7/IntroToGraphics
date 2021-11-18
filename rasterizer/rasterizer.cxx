#include <iostream>
#include <vtkDataSet.h>
#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <math.h>
#include <algorithm>
#include <cmath>
#include <stdlib.h>

#include "lightingparameters.h"
#include "screen.h"
#include "matrix.h"
#include "camera.h"
#include "space.h"
#include "triangle.h"

#define NORMALS

using namespace std;

int width = 1000;
int height = 1000;

std::vector<Triangle> GetTriangles(void)
{
   // Provided by instructor
   vtkPolyDataReader *rdr = vtkPolyDataReader::New();
   rdr->SetFileName("vtk_geometry.vtk");
   cerr << "Reading" << endl;
   rdr->Update();
   cerr << "Done reading" << endl;
   if (rdr->GetOutput()->GetNumberOfCells() == 0)
   {
       cerr << "Unable to open file!!" << endl;
       exit(EXIT_FAILURE);
   }
   vtkPolyData *pd = rdr->GetOutput();
   int numTris = pd->GetNumberOfCells();
   vtkPoints *pts = pd->GetPoints();
   vtkCellArray *cells = pd->GetPolys();
   vtkDoubleArray *var = (vtkDoubleArray *) pd->GetPointData()->GetArray("hardyglobal");
   double *color_ptr = var->GetPointer(0);
   vtkFloatArray *n = (vtkFloatArray *) pd->GetPointData()->GetNormals();
   float *normals = n->GetPointer(0);
   std::vector<Triangle> tris(numTris);
   vtkIdType npts;
   vtkIdType *ptIds;
   int idx;
   for (idx = 0, cells->InitTraversal() ; cells->GetNextCell(npts, ptIds) ; idx++)
   {
       if (npts != 3)
       {
           cerr << "Non-triangles!! ???" << endl;
           exit(EXIT_FAILURE);
       }
       double *pt = NULL;
       pt = pts->GetPoint(ptIds[0]);
       tris[idx].X[0] = pt[0];
       tris[idx].Y[0] = pt[1];
       tris[idx].Z[0] = pt[2];
#ifdef NORMALS
       tris[idx].normals[0][0] = normals[3*ptIds[0]+0];
       tris[idx].normals[0][1] = normals[3*ptIds[0]+1];
       tris[idx].normals[0][2] = normals[3*ptIds[0]+2];
#endif
       pt = pts->GetPoint(ptIds[1]);
       tris[idx].X[1] = pt[0];
       tris[idx].Y[1] = pt[1];
       tris[idx].Z[1] = pt[2];
#ifdef NORMALS
       tris[idx].normals[1][0] = normals[3*ptIds[1]+0];
       tris[idx].normals[1][1] = normals[3*ptIds[1]+1];
       tris[idx].normals[1][2] = normals[3*ptIds[1]+2];
#endif
       pt = pts->GetPoint(ptIds[2]);
       tris[idx].X[2] = pt[0];
       tris[idx].Y[2] = pt[1];
       tris[idx].Z[2] = pt[2];
#ifdef NORMALS
       tris[idx].normals[2][0] = normals[3*ptIds[2]+0];
       tris[idx].normals[2][1] = normals[3*ptIds[2]+1];
       tris[idx].normals[2][2] = normals[3*ptIds[2]+2];
#endif
       // 1->2 interpolate between light blue, dark blue
       // 2->2.5 interpolate between dark blue, cyan
       // 2.5->3 interpolate between cyan, green
       // 3->3.5 interpolate between green, yellow
       // 3.5->4 interpolate between yellow, orange
       // 4->5 interpolate between orange, brick
       // 5->6 interpolate between brick, salmon
       double mins[7] = { 1, 2, 2.5, 3, 3.5, 4, 5 };
       double maxs[7] = { 2, 2.5, 3, 3.5, 4, 5, 6 };
       unsigned char RGB[8][3] = { { 71, 71, 219 }, 
                                   { 0, 0, 91 },
                                   { 0, 255, 255 },
                                   { 0, 128, 0 },
                                   { 255, 255, 0 },
                                   { 255, 96, 0 },
                                   { 107, 0, 0 },
                                   { 224, 76, 76 } 
                                 };
       for (int j = 0 ; j < 3 ; j++)
       {
           float val = color_ptr[ptIds[j]];
           int r;
           for (r = 0 ; r < 7 ; r++)
           {
               if (mins[r] <= val && val < maxs[r])
                   break;
           }
           if (r == 7)
           {
               cerr << "Could not interpolate color for " << val << endl;
               exit(EXIT_FAILURE);
           }
           double proportion = (val-mins[r]) / (maxs[r]-mins[r]);
           tris[idx].colors[j][0] = (RGB[r][0]+proportion*(RGB[r+1][0]-RGB[r][0]))/255.0;
           tris[idx].colors[j][1] = (RGB[r][1]+proportion*(RGB[r+1][1]-RGB[r][1]))/255.0;
           tris[idx].colors[j][2] = (RGB[r][2]+proportion*(RGB[r+1][2]-RGB[r][2]))/255.0;
       }
   }

   return tris;
}

vtkImageData *NewImage(int width, int height)
{
    vtkImageData *image = vtkImageData::New(); // Create a new image from VTK library
    image->SetDimensions(width, height, 1); // Set the dimensions of the image
    image->AllocateScalars(VTK_UNSIGNED_CHAR, 3); // Allocate the scalars of the image
    return image; // Return the image
}

void WriteImage(vtkImageData *image, const char *filename)
{
   string full_filename = filename; // Create a string for the filename
   full_filename += ".png"; // Append .png extension to the filename
   vtkPNGWriter *writer = vtkPNGWriter::New(); // Create a VTK image writer
   writer->SetInputData(image); // Give the VTK write the image data
   writer->SetFileName(full_filename.c_str()); // Set the name of the VTK writer filename
   writer->Write(); // Write the data to the VTK writer
   writer->Delete(); // Delete the VTK writer
}

double GetCeiling(double n)
{
   return ceil(n-0.00001); // Return the ceiling of a number
}

int WithinBounds(int row, int column)
{
   // Check if the row and column are within bounds of the height and width
   if (0 <= row && row < height &&  0 <= column && column < width)
      return 1;
   else
      return 0;
}

void SwapValues(double &x1, double &x2)
{
   // Swap the values of numbers
   double temp = x1;
   x1 = x2;
   x2 = temp;
}

double SineParameterize(int curFrame, int nFrames, int ramp)
{
   // Provided by instructor
   int nNonRamp = nFrames-2*ramp;
   double height = 1./(nNonRamp + 4*ramp/M_PI);
   if (curFrame < ramp)
   {
       double factor = 2*height*ramp/M_PI;
       double eval = cos(M_PI/2*((double)curFrame)/ramp);
       return (1.-eval)*factor;
   }
   else if (curFrame > nFrames-ramp)
   {
       int amount_left = nFrames-curFrame;
       double factor = 2*height*ramp/M_PI;
       double eval =cos(M_PI/2*((double)amount_left/ramp));
       return 1. - (1-eval)*factor;
   }
   double amount_in_quad = ((double)curFrame-ramp);
   double quad_part = amount_in_quad*height;
   double curve_part = height*(2*ramp)/M_PI;
   return quad_part+curve_part;
}

double GetDotProduct(double *a, double *b)
{
    return ((a[0] * b[0]) + (a[1] * b[1]) + (a[2] * b[2])); // Calculate the dot product by definition and return it
}

void SetLightReflectionProportion(double *r, double *l, double *n)
{
   // Set the value of the light reflection proportion by definition
   double c = GetDotProduct(l, n);
   for (int i = 0; i < 3; i++)
      r[i] = 2*c*n[i] - l[i];
}

void Normalize(double *v)
{
   // Normalize the vector by definition
   double norm = sqrt((v[0]*v[0])+(v[1]*v[1])+(v[2]*v[2]));
   v[0] = v[0] / norm;
   v[1] = v[1] / norm;
   v[2] = v[2] / norm;
}

void CalculateShading(Camera c, Triangle &t, int i)
{
   // Provided by instructor
   LightingParameters lp;
   double ambient = lp.Ka;
   double diffuse = fabs(lp.Kd * GetDotProduct(lp.lightDir, t.normals[i]));
   double specular;
   double *r = new double[3];
   SetLightReflectionProportion(r, lp.lightDir, t.normals[i]);
   double v[3] = { c.position[0] - t.X[i], c.position[1] - t.Y[i], c.position[2] - t.Z[i] };
   Normalize(v);
   Normalize(r);
   specular = lp.Ks * pow(GetDotProduct(r, v), lp.alpha);
   if (isnan(specular))
      specular = 0;
   t.shading[i] = ambient + diffuse + specular;
}

void ScanLine(vector<Triangle> triangles, unsigned char *b, Screen s)
{
   // Find the direction the triangle is pointing, then the indices and rasterize accordingly
   for (int i = 0; i < triangles.size(); i++)
   {
       triangles[i].SetDirection();
       if (triangles[i].d == UP)
       {
           Triangle down = triangles[i];
           down.SetIndicesGoingUp(triangles[i]);
           down.Rasterize(b, s);
       }
       else if (triangles[i].d == DOWN)
       { 
           Triangle down = triangles[i];
           down.SetIndicesGoingDown(triangles[i]);
           down.Rasterize(b, s);
       }
       else
       {
           Triangle up = triangles[i];
           Triangle down = triangles[i];

           up.SetIndicesGoingUp(triangles[i]);
           down.SetIndicesGoingDown(triangles[i]);

           up.Rasterize(b, s);
           down.Rasterize(b, s);
       }
   }
}

void SetCrossProduct(double *result, double *u, double *w)
{
   // Set the cross product by definition
   result[0] = (u[1]*w[2]) - (u[2]*w[1]);
   result[1] = (u[2]*w[0]) - (u[0]*w[2]);
   result[2] = (u[0]*w[1]) - (u[1]*w[0]);
}

Matrix GetTransformationMatrix(Camera c)
{
    Space space; // Create a space that will get all of the matrices of the data
    Matrix worldSpace = space.GetWorldSpace(); // Create an identity matrix to determine transformed x, y, and z coordinates
    Matrix cameraSpace = space.GetCameraSpace(c); // Create a basis matrix that holds x, y, and z coordinates of the camera's location
    Matrix imageSpace = space.GetImageSpace(worldSpace, cameraSpace); // Create a camera space that holds the actual x, y, and z coordinates of a camera's location
    Matrix screenSpace = space.GetScreenSpace(c); // Create an image space that holds the pixels according to the camera's location
    Matrix deviceSpace = space.GetDeviceSpace(worldSpace); // Create a device space that is reflected from the world space
    Matrix result = Matrix::ComposeMatrices(imageSpace, screenSpace); // Create a temporary matrix to hold the output image of the screen from camera's location
    Matrix total_transform = Matrix::ComposeMatrices(result, deviceSpace); // Create the matrix that holds the output image from the camera's location
    return total_transform; // Return the matrix of the output matrix
}

void TransformTrianglesToDevicesSpace(vector<Triangle> triangles, Camera c, vtkImageData *image, unsigned char *buffer, Screen screen, int f)
{
   // For every triangle, set the shading, transform the matrix, then write the image data to the file
   for (int i = 0; i < triangles.size(); i++)
   {
      triangles[i].SetShading(c.position);

      Matrix total_transform = GetTransformationMatrix(c);

      for (int j = 0; j < 3; j++)
      {
          double points_in[4] = { triangles[i].X[j], triangles[i].Y[j], triangles[i].Z[j], 1 };
          double points_out[4] = { 0, 0, 0, 0 };
          total_transform.TransformPoint(points_in, points_out);
          double w = points_out[3];

          for (int k = 0; k < 3; k++)
             points_out[k] = points_out[k] / w;

          triangles[i].X[j] = points_out[0];
          triangles[i].Y[j] = points_out[1];
          triangles[i].Z[j] = points_out[2];
      }
   }

   ScanLine(triangles, buffer, screen); // Rasterize the image

   char temp[32]; // Buffer for filename
   sprintf(temp, "frame%03d", f); // Write to the filename buffer
   WriteImage(image, temp); // Write the image data
}

int main()
{
   int numImages = 100; // Total number of images
   vtkImageData *image = NewImage(width, height); // Create a VTK object
   unsigned char *buffer = (unsigned char *) image->GetScalarPointer(0,0,0); // The the buffer as the image's origin
   int numPrimaryColors = 3;
   int numPixels = width*height*numPrimaryColors; // Total number of pixels

   vector<Triangle> triangles = GetTriangles(); // Get the triangles from the VTK geomerty file

   // For every triangle, initialize the buffer and create an image
   for (int i = 0; i < numImages; i++)
   {
       // Initialize the buffer
       for (int j = 0; j < numPixels; j++)
          buffer[j] = 0;

       Screen screen(buffer); // Create a screen object
       Camera c(i, numImages); // Create an camera object
       TransformTrianglesToDevicesSpace(triangles, c, image, buffer, screen, i);
   }

   // Free the image data
   free(image);
   free(buffer);
}
