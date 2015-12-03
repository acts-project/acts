#ifndef EVENTPRIMITIVES_SYMMETRICMATRIXHELPERS_H
#define EVENTPRIMITIVES_SYMMETRICMATRIXHELPERS_H

#include <cmath>

    //////////////////////////////////////////////////////////////////////
    // Inverse calculation for symmetric matrix (5x5) with Gauss-Jordan
    //
    // Input parameters  : V(0/24) - only the elements of the lower
    //                               triangle of the matrix are accessed
    //                               0
    //                               1  6
    //                               2  7  12
    //                               3  8  13 18
    //                               4  9  14 19 24
    //
    // Output parameters : Vn(0/24) symmetric ////////////////////////////


    bool inverseSym5(Matrix<Scalar, 5,5>& out) const {

     if(this == &out || this->cols() != 5 || this->rows() != 5){
         return false;
     }
     const double * a = this->data();


     double* b = out.data();
     if   (a[ 0] == 0. ) {
         return false;
     }
     double x1    = 1./a[ 0],
            x2    =-a[ 1]*x1,
            x3    =-a[ 2]*x1,
            x4    =-a[ 3]*x1,
            x5    =-a[ 4]*x1;
     if   ((b[ 0] = a[ 6]+a[ 1]*x2)==0.){
         return false;
     }
            b[ 1] = a[ 7]+a[ 2]*x2;
            b[ 6] = a[12]+a[ 2]*x3;
            b[ 2] = a[ 8]+a[ 3]*x2;
            b[ 7] = a[13]+a[ 3]*x3;
            b[12] = a[18]+a[ 3]*x4;
            b[ 3] = a[ 9]+a[ 4]*x2;
            b[ 8] = a[14]+a[ 4]*x3;
            b[13] = a[19]+a[ 4]*x4;
            b[18] = a[24]+a[ 4]*x5;
     double y1    = 1./b[ 0],
            y2    =-b[ 1]*y1,
            y3    =-b[ 2]*y1,
            y4    =-b[ 3]*y1,
            y5    = x2   *y1;
     if   ((b[ 0] = b[ 6]+b[ 1]*y2)==0.) {
         return false;
     }
            b[ 1] = b[ 7]+b[ 2]*y2;
            b[ 6] = b[12]+b[ 2]*y3;
            b[ 2] = b[ 8]+b[ 3]*y2;
            b[ 7] = b[13]+b[ 3]*y3;
            b[12] = b[18]+b[ 3]*y4;
            b[ 3] = x3   +x2   *y2;
            b[ 8] = x4   +x2   *y3;
            b[13] = x5   +x2   *y4;
            b[18] = x1   +x2   *y5;
            x2    =-b[ 1]*(x1 = 1./b[ 0]);
            x3    =-b[ 2]* x1;
            x4    = b[ 3]* x1;
            x5    = y2   * x1;
     if   ((b[ 0] = b[ 6]+b[ 1]*x2)==0.) {
         return false;
     }
            b[ 1] = b[ 7]+b[ 2]*x2;
            b[ 6] = b[12]+b[ 2]*x3;
            b[ 2] = b[ 8]+b[ 3]*x2;
            b[ 7] = b[13]+b[ 3]*x3;
            b[12] = b[18]+b[ 3]*x4;
            b[ 3] = y3   +y2   *x2;
            b[ 8] = y4   +y2   *x3;
            b[13] = y5   +y2   *x4;
            b[18] = y1   +y2   *x5;
            y2    =-b[ 1]*(y1 = 1./b[ 0]);
            y3    = b[ 2]* y1;
            y4    = b[ 3]* y1;
            y5    = x2   * y1;
     if   ((b[ 0] = b[ 6]+b[ 1]*y2)==0.) {
         return false;
     }
            b[ 1] = b[ 7]+b[ 2]*y2;
            b[ 6] = b[12]+b[ 2]*y3;
            b[ 2] = b[ 8]+b[ 3]*y2;
            b[ 7] = b[13]+b[ 3]*y3;
            b[12] = b[18]+b[ 3]*y4;
            b[ 3] = x3   +x2   *y2;
            b[ 8] = x4   +x2   *y3;
            b[13] = x5   +x2   *y4;
            b[18] = x1   +x2   *y5;
            b[ 4] = b[ 1]*(b[24] = 1./b[ 0]);
            b[ 9] = b[ 2]* b[24];
            b[14] = b[ 3]* b[24];
            b[19] = y2   * b[24];
            b[ 0] = b[ 6]+b[ 1]*b[ 4];
            b[ 1] = b[ 7]+b[ 2]*b[ 4];
            b[ 6] = b[12]+b[ 2]*b[ 9];
            b[ 2] = b[ 8]+b[ 3]*b[ 4];
            b[ 7] = b[13]+b[ 3]*b[ 9];
            b[12] = b[18]+b[ 3]*b[14];
            b[ 3] = y3   +y2   *b[ 4];
            b[ 8] = y4   +y2   *b[ 9];
            b[13] = y5   +y2   *b[14];
            b[18] = y1   +y2   *b[19];

            b[ 5] = b[ 1];
            b[10] = b[ 2];
            b[15] = b[ 3];
            b[20] = b[ 4];

            b[11] = b[ 7];
            b[16] = b[ 8];
            b[21] = b[ 9];

            b[17] = b[13];
            b[22] = b[14];

            b[23] = b[19];

            return true;

    }

    ////////////////////////////////////////////////////////////////////////////////////////////
    // Similarity calculation for symmetric matrix (5x5)
    //
    // Input parameters  : V(0/24) - only the elements of the lower
    //                               triangle of the matrix are accessed
    //                               0
    //                               1  6
    //                               2  7  12
    //                               3  8  13 18
    //                               4  9  14 19 24
    //                      Jac(0/24) -jacobiam of transformation
    //                               0  5 10 15 20
    //                               1  6 11 16 21
    //                               2  7 12 17 22
    //                               3  8 13 18 23
    //                               4  9 14 19 24
    //
    // Output parameters : Vn(0/24) = Jac*Vn*Jact symmetric /////////

    Matrix<Scalar,RowsAtCompileTime,ColsAtCompileTime> similaritySym5(Matrix<Scalar,RowsAtCompileTime,ColsAtCompileTime>& JacMat) const {


     if(this->rows() != 5 || this->cols() != 5 || JacMat.rows() != 5 || JacMat.cols() !=5){
         return (*this).similarity(JacMat);
     }

     const double * V = this->data();
     const double * Jac = JacMat.data();
     double * Vn;
     Matrix<Scalar, 5, 5> out;
     Vn = out.data();

     double a11 = (Jac[ 0]*V[ 0]+Jac[ 5]*V[ 1]+Jac[10]*V[ 2])+(Jac[15]*V[ 3]+Jac[20]*V[ 4]);
     double a12 = (Jac[ 0]*V[ 1]+Jac[ 5]*V[ 6]+Jac[10]*V[ 7])+(Jac[15]*V[ 8]+Jac[20]*V[ 9]);
     double a13 = (Jac[ 0]*V[ 2]+Jac[ 5]*V[ 7]+Jac[10]*V[12])+(Jac[15]*V[13]+Jac[20]*V[14]);
     double a14 = (Jac[ 0]*V[ 3]+Jac[ 5]*V[ 8]+Jac[10]*V[13])+(Jac[15]*V[18]+Jac[20]*V[19]);
     double a15 = (Jac[ 0]*V[ 4]+Jac[ 5]*V[ 9]+Jac[10]*V[14])+(Jac[15]*V[19]+Jac[20]*V[24]);

     Vn[ 0] = (a11*Jac[ 0]+a12*Jac[ 5]+a13*Jac[10])+(a14*Jac[15]+a15*Jac[20]);

     double a21 = (Jac[ 1]*V[ 0]+Jac[ 6]*V[ 1]+Jac[11]*V[ 2])+(Jac[16]*V[ 3]+Jac[21]*V[ 4]);
     double a22 = (Jac[ 1]*V[ 1]+Jac[ 6]*V[ 6]+Jac[11]*V[ 7])+(Jac[16]*V[ 8]+Jac[21]*V[ 9]);
     double a23 = (Jac[ 1]*V[ 2]+Jac[ 6]*V[ 7]+Jac[11]*V[12])+(Jac[16]*V[13]+Jac[21]*V[14]);
     double a24 = (Jac[ 1]*V[ 3]+Jac[ 6]*V[ 8]+Jac[11]*V[13])+(Jac[16]*V[18]+Jac[21]*V[19]);
     double a25 = (Jac[ 1]*V[ 4]+Jac[ 6]*V[ 9]+Jac[11]*V[14])+(Jac[16]*V[19]+Jac[21]*V[24]);

     Vn[ 1] = (a21*Jac[ 0]+a22*Jac[ 5]+a23*Jac[10])+(a24*Jac[15]+a25*Jac[20]);
     Vn[ 6] = (a21*Jac[ 1]+a22*Jac[ 6]+a23*Jac[11])+(a24*Jac[16]+a25*Jac[21]);

     double a31 = (Jac[ 2]*V[ 0]+Jac[ 7]*V[ 1]+Jac[12]*V[ 2])+(Jac[17]*V[ 3]+Jac[22]*V[ 4]);
     double a32 = (Jac[ 2]*V[ 1]+Jac[ 7]*V[ 6]+Jac[12]*V[ 7])+(Jac[17]*V[ 8]+Jac[22]*V[ 9]);
     double a33 = (Jac[ 2]*V[ 2]+Jac[ 7]*V[ 7]+Jac[12]*V[12])+(Jac[17]*V[13]+Jac[22]*V[14]);
     double a34 = (Jac[ 2]*V[ 3]+Jac[ 7]*V[ 8]+Jac[12]*V[13])+(Jac[17]*V[18]+Jac[22]*V[19]);
     double a35 = (Jac[ 2]*V[ 4]+Jac[ 7]*V[ 9]+Jac[12]*V[14])+(Jac[17]*V[19]+Jac[22]*V[24]);

     Vn[ 2] = (a31*Jac[ 0]+a32*Jac[ 5]+a33*Jac[10])+(a34*Jac[15]+a35*Jac[20]);
     Vn[ 7] = (a31*Jac[ 1]+a32*Jac[ 6]+a33*Jac[11])+(a34*Jac[16]+a35*Jac[21]);
     Vn[12] = (a31*Jac[ 2]+a32*Jac[ 7]+a33*Jac[12])+(a34*Jac[17]+a35*Jac[22]);

     double a41 = (Jac[ 3]*V[ 0]+Jac[ 8]*V[ 1]+Jac[13]*V[ 2])+(Jac[18]*V[ 3]+Jac[23]*V[ 4]);
     double a42 = (Jac[ 3]*V[ 1]+Jac[ 8]*V[ 6]+Jac[13]*V[ 7])+(Jac[18]*V[ 8]+Jac[23]*V[ 9]);
     double a43 = (Jac[ 3]*V[ 2]+Jac[ 8]*V[ 7]+Jac[13]*V[12])+(Jac[18]*V[13]+Jac[23]*V[14]);
     double a44 = (Jac[ 3]*V[ 3]+Jac[ 8]*V[ 8]+Jac[13]*V[13])+(Jac[18]*V[18]+Jac[23]*V[19]);
     double a45 = (Jac[ 3]*V[ 4]+Jac[ 8]*V[ 9]+Jac[13]*V[14])+(Jac[18]*V[19]+Jac[23]*V[24]);

     Vn[ 3] = (a41*Jac[ 0]+a42*Jac[ 5]+a43*Jac[10])+(a44*Jac[15]+a45*Jac[20]);
     Vn[ 8] = (a41*Jac[ 1]+a42*Jac[ 6]+a43*Jac[11])+(a44*Jac[16]+a45*Jac[21]);
     Vn[13] = (a41*Jac[ 2]+a42*Jac[ 7]+a43*Jac[12])+(a44*Jac[17]+a45*Jac[22]);
     Vn[18] = (a41*Jac[ 3]+a42*Jac[ 8]+a43*Jac[13])+(a44*Jac[18]+a45*Jac[23]);

     double a51 = (Jac[ 4]*V[ 0]+Jac[ 9]*V[ 1]+Jac[14]*V[ 2])+(Jac[19]*V[ 3]+Jac[24]*V[ 4]);
     double a52 = (Jac[ 4]*V[ 1]+Jac[ 9]*V[ 6]+Jac[14]*V[ 7])+(Jac[19]*V[ 8]+Jac[24]*V[ 9]);
     double a53 = (Jac[ 4]*V[ 2]+Jac[ 9]*V[ 7]+Jac[14]*V[12])+(Jac[19]*V[13]+Jac[24]*V[14]);
     double a54 = (Jac[ 4]*V[ 3]+Jac[ 9]*V[ 8]+Jac[14]*V[13])+(Jac[19]*V[18]+Jac[24]*V[19]);
     double a55 = (Jac[ 4]*V[ 4]+Jac[ 9]*V[ 9]+Jac[14]*V[14])+(Jac[19]*V[19]+Jac[24]*V[24]);

     Vn[ 4] = (a51*Jac[ 0]+a52*Jac[ 5]+a53*Jac[10])+(a54*Jac[15]+a55*Jac[20]);
     Vn[ 9] = (a51*Jac[ 1]+a52*Jac[ 6]+a53*Jac[11])+(a54*Jac[16]+a55*Jac[21]);
     Vn[14] = (a51*Jac[ 2]+a52*Jac[ 7]+a53*Jac[12])+(a54*Jac[17]+a55*Jac[22]);
     Vn[19] = (a51*Jac[ 3]+a52*Jac[ 8]+a53*Jac[13])+(a54*Jac[18]+a55*Jac[23]);
     Vn[24] = (a51*Jac[ 4]+a52*Jac[ 9]+a53*Jac[14])+(a54*Jac[19]+a55*Jac[24]);

     Vn[ 5] = Vn[ 1];
     Vn[10] = Vn[ 2];
     Vn[15] = Vn[ 3];
     Vn[20] = Vn[ 4];

     Vn[11] = Vn[ 7];
     Vn[16] = Vn[ 8];
     Vn[21] = Vn[ 9];

     Vn[17] = Vn[13];
     Vn[22] = Vn[14];

     Vn[23] = Vn[19];

     return out;
    }

#endif
