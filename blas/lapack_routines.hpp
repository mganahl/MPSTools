#ifndef lapack_routines_h_
#define lapack_routines_h_

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <iostream>
#include <string>
#include "../boost_typedefs.hpp"
//#include "BoostToolsLapack/LapackSVD.h"
#include "boost/multi_array.hpp"
#include "helperfunctions.hpp"

using namespace BoostTypedefs;
using namespace std;

int tridiag(double *diag,int length_diag, double* sub_diag,double*eigvec,char JOBZ);
int tridiag(Complex *diag,int length_diag, Complex* sub_diag,double*eigvec,char JOBZ);

int svd(double *matrix, double *U,double*VT,double *sing_val, int dim1,int dim2);
int rsvd(double *matrix, double *U,double*VT,double *sing_val, int dim1,int dim2);


//tested and working
int BoostSVD(const CMatrixType & matrix, CMatrixType &U,MatrixType &D,CMatrixType &V_dagger);
//tested and working
int BoostSVD(const CMatrixType & matrix, CMatrixType &U,VectorType &D,CMatrixType &V_dagger);
//tested and working
int BoostSVD(const MatrixType & matrix, MatrixType &U,MatrixType &D,MatrixType &V_dagger);
//tested and working
int BoostSVD(const MatrixType & matrix, MatrixType &U,VectorType &D,MatrixType &V_dagger);

//tested and working
int BoostSVD(const CMatrixTypeColMaj & matrix, CMatrixTypeColMaj &U,MatrixTypeColMaj &D,CMatrixTypeColMaj &V_dagger);
//tested and working
int BoostSVD(const CMatrixTypeColMaj & matrix, CMatrixTypeColMaj &U,VectorType &D,CMatrixTypeColMaj &V_dagger);


//do not pass as reference!! matrix is changed in lapack svd routine!!!

//tested and working
int BoostSVD(MatrixTypeColMaj  matrix, MatrixTypeColMaj &U,MatrixTypeColMaj &D,MatrixTypeColMaj &V_dagger);
//tested and working
int BoostSVD(MatrixTypeColMaj  matrix, MatrixTypeColMaj &U,VectorType &D,MatrixTypeColMaj &V_dagger);

int svd_MultiArray(const boost::multi_array<Real,2> &mat,boost::multi_array<Real,2>&U,boost::multi_array<Real,2> &D,boost::multi_array<Real,2> &V_dagger);
int svd_CMultiArray(const boost::multi_array<Complex,2> &mat,boost::multi_array<Complex,2>&U,boost::multi_array<Complex,2> &D,boost::multi_array<Complex,2> &V_dagger);

int diag_hermitian(double *matrix, double *eigen, int dim);
void BoostDiag(CMatrixTypeColMaj &M,VectorType &E);
void BoostDiag(MatrixTypeColMaj &M,VectorType &E);

#endif
