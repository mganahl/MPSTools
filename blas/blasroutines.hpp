/*
 * blasroutines.hpp
 *
 *  Created on: 07.09.2010
 *      Author: martinjg
 */

#ifndef BLASROUTINES_HPP_
#define BLASROUTINES_HPP_
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "helperfunctions.hpp"
#include "../boost_typedefs.hpp"
#include "boost/multi_array.hpp"

//using namespace Boost;
using namespace BoostTypedefs;
Real dot_prod(int, double*, double*);
Complex dot_prod(int, Complex*, Complex*);

void cmat_mat_prod(Complex alpha,double*matA ,int dim1A,int dim2A, double* matB,int dim1B, int dim2B,double* mat_out,char tra,char trb);
void mat_mat_prod(double alpha,double*matA ,int dim1A,int dim2A, double* matB, int dim1B, int dim2B,double* mat_out,char tra,char trb);

void mat_mat_prod_2(double alpha,double*matA ,int dim1A,int dim2A, double* matB, int dim1B, int dim2B,double* mat_out,char tra,char trb);

//tested and working
void MatMatProd(const Complex alpha,const CMatrixTypeColMaj  & Ar1,const CMatrixTypeColMaj &Ar2,CMatrixTypeColMaj & out, char tra,char trb,bool test=false);
//tested and working
void MatMatProd(const Complex alpha,const VectorType  &Ar1,const CMatrixTypeColMaj &Ar2,CMatrixTypeColMaj & out, char tra,char trb);
//tested and working
void MatMatProd(const Complex alpha,const CMatrixTypeColMaj  & Ar1,const VectorType &Ar2,CMatrixTypeColMaj & out, char tra,char trb);
//tested and working
void MatMatProd(const double alpha,const MatrixTypeColMaj  & Ar1,const MatrixTypeColMaj &Ar2,MatrixTypeColMaj & out, char tra,char trb,bool test=false);
//tested and working
void MatMatProd(const double alpha,const VectorType  & Ar1,const MatrixTypeColMaj &Ar2,MatrixTypeColMaj & out, char tra,char trb);
//tested and working
void MatMatProd(const double alpha,const MatrixTypeColMaj & Ar1,const VectorType &Ar2,MatrixTypeColMaj & out, char tra,char trb);

void MatMatProd_2(const double alpha,const MatrixTypeColMaj & Ar1,const MatrixTypeColMaj &Ar2,MatrixTypeColMaj & out, char tra,char trb,bool test=false);
//TODO: check the rwo routines!
void MatMatProd(const double alpha,const MatrixType& Ar1,const VectorType&Ar2,MatrixType& out, char tra,char trb);
void MatMatProd(const Complex alpha,const CMatrixType& Ar1,const VectorType&Ar2,CMatrixType& out, char tra,char trb);


void MatMatProd(const double alpha,const VectorType& Ar1,const MatrixType&Ar2,MatrixType& out, char tra,char trb);
void MatMatProd(const Complex alpha,const VectorType& Ar1,const CMatrixType&Ar2,CMatrixType& out, char tra,char trb);

void MatMatProd(const double alpha,const MatrixType& Ar1,const MatrixType&Ar2,MatrixType& out, char tra,char trb);
void MatMatProd(const Complex alpha,const CMatrixType& Ar1,const CMatrixType&Ar2,CMatrixType& out, char tra,char trb);



void cmat_vec_prod(double*mat,int M,int N,double*vec_in,double*vec_out,char tr);
void mat_vec_prod(double*mat,int M,int N,double*vec_in,double*vec_out,char tr);

int tridiag(double *diag,int length_diag, double* sub_diag,double*eigvec,char JOBZ = 'N');
#endif /* BLASROUTINES_HPP_ */
