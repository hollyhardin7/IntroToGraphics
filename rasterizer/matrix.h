#ifndef MATRIX_H
#define MATRIX_H

class Matrix
{
  public:
    double          A[4][4];

    Matrix();
    void            PutVectorsInMatrix(double *u, double *v, double *w, double *o);
    void            TransformPoint(const double *ptIn, double *ptOut);
    static Matrix   ComposeMatrices(const Matrix m1, const Matrix m2);
};

Matrix::Matrix()
{
   // Initialize matrix's values to zero
   for (int i = 0; i < 4; i++)
       for (int j = 0; j < 4; j++)
           A[i][j] = 0.0;
}

void Matrix::PutVectorsInMatrix(double *u, double *v, double *w, double *o)
{
   // Combine the vectors into matrix form
   for (int i = 0; i < 4; i++)
      A[i][0] = u[i];
   for (int i = 0; i < 4; i++)
      A[i][1] = v[i];
   for (int i = 0; i < 4; i++)
      A[i][2] = w[i];
   for (int i = 0; i < 4; i++)
      A[i][3] = o[i];
}

void Matrix::TransformPoint(const double *ptIn, double *ptOut)
{
   // Apply the transformation matrix to the original points
   for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++)
         ptOut[i] += ptIn[j]*A[j][i];
}

Matrix Matrix::ComposeMatrices(const Matrix M1, const Matrix M2)
{
    // Multiply the matrices
    Matrix result;
    for (int i = 0 ; i < 4 ; i++)
        for (int j = 0 ; j < 4 ; j++)
        {
            result.A[i][j] = 0;
            for (int k = 0 ; k < 4 ; k++)
                result.A[i][j] += M1.A[i][k]*M2.A[k][j];
        }

    return result;
}

#endif
