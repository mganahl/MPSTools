/*
 * blasroutines.cpp
 *
 *  Created on: 07.09.2010
 *      Author: martinjg
 */


#include "blasroutines.hpp"



extern "C"
{
/*
    SUBROUTINE ZGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
*     .. Scalar Arguments ..
    DOUBLE COMPLEX ALPHA,BETA
    INTEGER K,LDA,LDB,LDC,M,N
    CHARACTER TRANSA,TRANSB
*     ..
*     .. Array Arguments ..
    DOUBLE COMPLEX A(LDA,*),B(LDB,*),C(LDC,*)
*     ..
*
*  Purpose
*  =======
*
*  ZGEMM  performs one of the matrix-matrix operations
*
*     C := alpha*op( A )*op( B ) + beta*C,
*
*  where  op( X ) is one of
*
*     op( X ) = X   or   op( X ) = X'   or   op( X ) = conjg( X' ),
*
*  alpha and beta are scalars, and A, B and C are matrices, with op( A )
*  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
*
*  Arguments
*  ==========
*
*  TRANSA - CHARACTER*1.
*           On entry, TRANSA specifies the form of op( A ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSA = 'N' or 'n',  op( A ) = A.
*
*              TRANSA = 'T' or 't',  op( A ) = A'.
*
*              TRANSA = 'C' or 'c',  op( A ) = conjg( A' ).
*
*           Unchanged on exit.
*
*  TRANSB - CHARACTER*1.
*           On entry, TRANSB specifies the form of op( B ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSB = 'N' or 'n',  op( B ) = B.
*
*              TRANSB = 'T' or 't',  op( B ) = B'.
*
*              TRANSB = 'C' or 'c',  op( B ) = conjg( B' ).
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry,  M  specifies  the number  of rows  of the  matrix
*           op( A )  and of the  matrix  C.  M  must  be at least  zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry,  N  specifies the number  of columns of the matrix
*           op( B ) and the number of columns of the matrix C. N must be
*           at least zero.
*           Unchanged on exit.
*
*  K      - INTEGER.
*           On entry,  K  specifies  the number of columns of the matrix
*           op( A ) and the number of rows of the matrix op( B ). K must
*           be at least  zero.
*           Unchanged on exit.
*
*  ALPHA  - COMPLEX*16      .
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - COMPLEX*16       array of DIMENSION ( LDA, ka ), where ka is
*           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
*           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
*           part of the array  A  must contain the matrix  A,  otherwise
*           the leading  k by m  part of the array  A  must contain  the
*           matrix A.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
*           LDA must be at least  max( 1, m ), otherwise  LDA must be at
*           least  max( 1, k ).
*           Unchanged on exit.
*
*  B      - COMPLEX*16       array of DIMENSION ( LDB, kb ), where kb is
*           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
*           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
*           part of the array  B  must contain the matrix  B,  otherwise
*           the leading  n by k  part of the array  B  must contain  the
*           matrix B.
*           Unchanged on exit.
*
*  LDB    - INTEGER.
*           On entry, LDB specifies the first dimension of B as declared
*           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
*           LDB must be at least  max( 1, k ), otherwise  LDB must be at
*           least  max( 1, n ).
*           Unchanged on exit.
*
*  BETA   - COMPLEX*16      .
*           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
*           supplied as zero then C need not be set on input.
*           Unchanged on exit.
*
*  C      - COMPLEX*16       array of DIMENSION ( LDC, n ).
*           Before entry, the leading  m by n  part of the array  C must
*           contain the matrix  C,  except when  beta  is zero, in which
*           case C need not be set on entry.
*           On exit, the array  C  is overwritten by the  m by n  matrix
*           ( alpha*op( A )*op( B ) + beta*C ).
*
*  LDC    - INTEGER.
*           On entry, LDC specifies the first dimension of C as declared
*           in  the  calling  (sub)  program.   LDC  must  be  at  least
*           max( 1, m ).
*           Unchanged on exit.
*/


  void zgemm_(char* TRANSA,char *TRANSB,int* M,int* N,int * K,double* ALPHA,double *A,int* LDA,double *B,int *LDB,double* BETA,double *C,int* LDC);

}
//no check if dimensions of matrices match!!!! mat_out has to be allocated before it is passed to cmat_mat_prod! note that you have to allocate space for
//real and imaginary parts of the matrices
void cmat_mat_prod(Complex alpha,double*matA ,int dim1A,int dim2A, double* matB,int dim1B, int dim2B,double* mat_out,char tra,char trb)
{
	char Transa = tra,Transb=trb;
	int M,N,Ka,Kb;
	int lda = dim1A,ldb=dim1B,ldc;
	if(tra=='n'||tra=='N')
	{
		M = dim1A;
		Ka = dim2A;
	}
	else
	{
		M = dim2A;
		Ka = dim1A;
	}
	if(trb =='n'||trb=='N') {
		N = dim2B;
		Kb = dim1B;
	}
	else{
		N = dim1B;
		Kb = dim2B;
	}
	ldc = M;
	assert(Ka==Kb);


	double alpha_[2] ={alpha.real(),alpha.imag()};
	double beta[2] = {0,0};

	zgemm_(&Transa,&Transb,&M,&N,&Ka,alpha_,matA,&lda,matB,&ldb,beta,mat_out,&ldc);

}


extern "C"
{
	extern Real ddot_(int*,double*,int*,double*,int*);
}

Real dot_prod(int size,double* v1,double*v2)
{
	int inc1 = 1,inc2 = 1;;
	double out= ddot_(&size,v1,&inc1,v2,&inc2);
	return out;
}

Complex dot_prod(int size,  Complex* v1, Complex*v2)
{
	double *in1 = new double[2*size];
	double *in2 = new double[2*size];
	double* out = new double [2];
	for(int i=0;i<size;i++)
	{
		in1[2*i] = v1[i].real();
		in1[2*i+1] = -v1[i].imag();//this is the complex conjugation
		in2[2*i] = v2[i].real();
		in2[2*i+1] = v2[i].imag();
	}
	cmat_mat_prod(Complex (1.0),in1,1,size,in2,size,1,out,'n','n');
	Complex res;
	res.real() = out[0];
	res.imag() = out[1];
	delete [] in1;
	delete [] in2;
	delete [] out;
	return res;
}



extern "C"
{
	/*
	SUBROUTINE DGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
	*     .. Scalar Arguments ..
	      DOUBLE PRECISION ALPHA,BETA
	      INTEGER INCX,INCY,LDA,M,N
	      CHARACTER TRANS
	*     ..
	*     .. Array Arguments ..
	      DOUBLE PRECISION A(LDA,*),X(*),Y(*)
	*     ..
	*
	*  Purpose
	*  =======
	*
	*  DGEMV  performs one of the matrix-vector operations
	*
	*     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
	*
	*  where alpha and beta are scalars, x and y are vectors and A is an
	*  m by n matrix.
	*
	*  Arguments
	*  ==========
	*
	*  TRANS  - CHARACTER*1.
	*           On entry, TRANS specifies the operation to be performed as
	*           follows:
	*
	*              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
	*
	*              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
	*
	*              TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y.
	*
	*           Unchanged on exit.
	*
	*  M      - INTEGER.
	*           On entry, M specifies the number of rows of the matrix A.
	*           M must be at least zero.
	*           Unchanged on exit.
	*
	*  N      - INTEGER.
	*           On entry, N specifies the number of columns of the matrix A.
	*           N must be at least zero.
	*           Unchanged on exit.
	*
	*  ALPHA  - DOUBLE PRECISION.
	*           On entry, ALPHA specifies the scalar alpha.
	*           Unchanged on exit.
	*
	*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
	*           Before entry, the leading m by n part of the array A must
	*           contain the matrix of coefficients.
	*           Unchanged on exit.
	*
	*  LDA    - INTEGER.
	*           On entry, LDA specifies the first dimension of A as declared
	*           in the calling (sub) program. LDA must be at least
	*           max( 1, m ).
	*           Unchanged on exit.
	*
	*  X      - DOUBLE PRECISION array of DIMENSION at least
	*           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
	*           and at least
	*           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
	*           Before entry, the incremented array X must contain the
	*           vector x.
	*           Unchanged on exit.
	*
	*  INCX   - INTEGER.
	*           On entry, INCX specifies the increment for the elements of
	*           X. INCX must not be zero.
	*           Unchanged on exit.
	*
	*  BETA   - DOUBLE PRECISION.
	*           On entry, BETA specifies the scalar beta. When BETA is
	*           supplied as zero then Y need not be set on input.
	*           Unchanged on exit.
	*
	*  Y      - DOUBLE PRECISION array of DIMENSION at least
	*           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
	*           and at least
	*           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
	*           Before entry with BETA non-zero, the incremented array Y
	*           must contain the vector y. On exit, Y is overwritten by the
	*           updated vector y.
	*
	*  INCY   - INTEGER.
	*           On entry, INCY specifies the increment for the elements of
	*           Y. INCY must not be zero.
	*           Unchanged on exit.
	*
*/

  void dgemv_(char* TRANS,int* M,int* N,double* ALPHA,double *A,int* LDA,double *X,int* INCX,double* BETA,double *Y,int* INCY);

}

//vec_out is not checked for dimension consistency!! dimension mismatch may cause unpredictable output!
void mat_vec_prod(double alpha,double*mat,int M,int N,double*vec_in,double*vec_out,char tr)
{
  char Trans = tr;
  //  double alpha = 1;
  int lda = M;
  int incx = 1;
  double beta = 0;
  //double y[M];
  int incy =1;

  dgemv_(&Trans,&M,&N,&alpha,mat,&lda,vec_in,&incx,&beta,vec_out,&incy);

}




extern "C"
{/*
	SUBROUTINE ZGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
	*     .. Scalar Arguments ..
	      DOUBLE COMPLEX ALPHA,BETA
	      INTEGER INCX,INCY,LDA,M,N
	      CHARACTER TRANS
	*     ..
	*     .. Array Arguments ..
	      DOUBLE COMPLEX A(LDA,*),X(*),Y(*)
	*     ..
	*
	*  Purpose
	*  =======
	*
	*  ZGEMV  performs one of the matrix-vector operations
	*
	*     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,   or
	*
	*     y := alpha*conjg( A' )*x + beta*y,
	*
	*  where alpha and beta are scalars, x and y are vectors and A is an
	*  m by n matrix.
	*
	*  Arguments
	*  ==========
	*
	*  TRANS  - CHARACTER*1.
	*           On entry, TRANS specifies the operation to be performed as
	*           follows:
	*
	*              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
	*
	*              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
	*
	*              TRANS = 'C' or 'c'   y := alpha*conjg( A' )*x + beta*y.
	*
	*           Unchanged on exit.
	*
	*  M      - INTEGER.
	*           On entry, M specifies the number of rows of the matrix A.
	*           M must be at least zero.
	*           Unchanged on exit.
	*
	*  N      - INTEGER.
	*           On entry, N specifies the number of columns of the matrix A.
	*           N must be at least zero.
	*           Unchanged on exit.
	*
	*  ALPHA  - COMPLEX*16      .
	*           On entry, ALPHA specifies the scalar alpha.
	*           Unchanged on exit.
	*
	*  A      - COMPLEX*16       array of DIMENSION ( LDA, n ).
	*           Before entry, the leading m by n part of the array A must
	*           contain the matrix of coefficients.
	*           Unchanged on exit.
	*
	*  LDA    - INTEGER.
	*           On entry, LDA specifies the first dimension of A as declared
	*           in the calling (sub) program. LDA must be at least
	*           max( 1, m ).
	*           Unchanged on exit.
	*
	*  X      - COMPLEX*16       array of DIMENSION at least
	*           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
	*           and at least
	*           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
	*           Before entry, the incremented array X must contain the
	*           vector x.
	*           Unchanged on exit.
	*
	*  INCX   - INTEGER.
	*           On entry, INCX specifies the increment for the elements of
	*           X. INCX must not be zero.
	*           Unchanged on exit.
	*
	*  BETA   - COMPLEX*16      .
	*           On entry, BETA specifies the scalar beta. When BETA is
	*           supplied as zero then Y need not be set on input.
	*           Unchanged on exit.
	*
	*  Y      - COMPLEX*16       array of DIMENSION at least
	*           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
	*           and at least
	*           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
	*           Before entry with BETA non-zero, the incremented array Y
	*           must contain the vector y. On exit, Y is overwritten by the
	*           updated vector y.
	*
	*  INCY   - INTEGER.
	*           On entry, INCY specifies the increment for the elements of
	*           Y. INCY must not be zero.
	*           Unchanged on exit.
	*
*/


  void zgemv_(char* TRANS,int* M,int* N,double* ALPHA,double *A,int* LDA,double *X,int* INCX,double* BETA,double *Y,int* INCY);

}

//vec_out is not checked for dimension consistency!! dimension mismatch may cause unpredictable output!
void cmat_vec_prod(double*mat,int M,int N,double*vec_in,double*vec_out,char tr)
{

	char Trans = tr;
	double alpha[2] = {1,0};
	int lda = M;
	int incx = 1;
	double beta[2] = {0,0};
	//double y[M];
	int incy =1;

	zgemv_(&Trans,&M,&N,alpha,mat,&lda,vec_in,&incx,beta,vec_out,&incy);

}


//=?===========================================================Complex Matrix Products========================================================
//tested and working
void MatMatProd(const Complex alpha,const VectorType &Ar1,const CMatrixTypeColMaj &Ar2,CMatrixTypeColMaj & out, char tra,char trb)
{
  double* mat1 = new double[Ar1.size() * Ar1.size()*2];
  double* mat2 = new double[Ar2.size1() * Ar2.size2()*2 ];
  //std::cout<<"after alloc"<<std::endl;
  //std::cout<<((tra=='N'||tra=='n')?Ar1.shape()[0]:Ar1.shape()[1]) * 2 * ((trb=='N'||trb=='n')?Ar2.shape()[1]:Ar2.shape()[0])<<std::endl;
  double* res = new double[((tra=='N'||tra=='n')?Ar1.size():Ar1.size()) * 2 * ((trb=='N'||trb=='n')?Ar2.size2():Ar2.size1())];
  uint dimx = Ar1.size();
  uint dimy = Ar1.size();
  for (uint y = 0; y < dimy; ++y) {
    for(uint x = 0;x<dimx;++x)
      {
	if(x==y)
	  {
	    mat1[2*y*dimx+2*x] = Ar1(x);
	    mat1[2*y*dimx+2*x+1] = 0;
	  }
	else{
	  mat1[2*y*dimx+2*x] = 0;
	  mat1[2*y*dimx+2*x+1] = 0;
	}
      }

  }
  dimx = Ar2.size1();
  dimy = Ar2.size2();
  for (uint cnt = 0; cnt < Ar2.size1() * Ar2.size2(); ++cnt) {
    mat2[2 * cnt] = Ar2.data()[cnt].real();
    mat2[2 * cnt + 1] = Ar2.data()[cnt].imag();
  }
  //std::cout<<"before fortran routine"<<std::endl;
  cmat_mat_prod(alpha,mat1, Ar1.size(), Ar1.size(), mat2,
		Ar2.size1(), Ar2.size2(), res, tra, trb);
  //std::cout<<"after fortran routine"<<std::endl;
  out.resize((tra=='N'||tra=='n')?Ar1.size():Ar1.size(),(trb=='N'||trb=='n')?Ar2.size2():Ar2.size1());
  dimx = out.size1();
  dimy = out.size2();
  for (uint y = 0; y < dimy; y++){
    for(uint x = 0;x<dimx;x++){
      out(x,y).real()=res[2*y*dimx+2*x];
      out(x,y).imag()=res[2*y*dimx+2*x+1];
    }
  }
  delete[] mat1;
  delete[] mat2;
  delete[] res;
}

//tested and working
void MatMatProd(const Complex alpha,const CMatrixTypeColMaj & Ar1,const VectorType & Ar2,CMatrixTypeColMaj & out, char tra,char trb)
{
	double* mat1 = new double[Ar1.size1() * Ar1.size2()*2];
	double* mat2 = new double[Ar2.size() * Ar2.size()* 2];


	//std::cout<<"after alloc"<<std::endl;
	//std::cout<<((tra=='N'||tra=='n')?Ar1.shape()[0]:Ar1.shape()[1]) * 2 * ((trb=='N'||trb=='n')?Ar2.shape()[1]:Ar2.shape()[0])<<std::endl;
	double* res = new double[((tra=='N'||tra=='n')?Ar1.size1():Ar1.size2()) * 2 * ((trb=='N'||trb=='n')?Ar2.size():Ar2.size())];

	uint dimx = Ar1.size1();
	uint dimy = Ar1.size2();
	for (uint cnt = 0; cnt < Ar1.size1() * Ar1.size2(); cnt++) {
			mat1[2 * cnt] = Ar1.data()[cnt].real();
			mat1[2 * cnt + 1] = Ar1.data()[cnt].imag();
	}
/*
	for (uint y = 0; y < dimy; y++) {
		for(uint x = 0;x<dimx;x++)
		{
			mat1[2*y*dimx+2*x] = Ar1(x,y).real();
			mat1[2*y*dimx+2*x+1] = Ar1(x,y).imag();
		}

	}
	*/

	dimx = Ar2.size();
	dimy = Ar2.size();
	for (uint y = 0; y < dimy; y++) {
		for(uint x = 0;x<dimx;x++)
		{
			if(x==y)
			{
				mat2[2*y*dimx+2*x] = Ar2(x);
				mat2[2*y*dimx+2*x+1] = 0;
			}
			else{
				mat2[2*y*dimx+2*x] = 0;
				mat2[2*y*dimx+2*x+1] = 0;
			}
		}

	}


	//std::cout<<"before fortran routine"<<std::endl;
	cmat_mat_prod(alpha,mat1, Ar1.size1(), Ar1.size2(), mat2,
			Ar2.size(), Ar2.size(), res, tra, trb);
	//std::cout<<"after fortran routine"<<std::endl;

	out.resize((tra=='N'||tra=='n')?Ar1.size1():Ar1.size2(),(trb=='N'||trb=='n')?Ar2.size():Ar2.size());


	dimx = out.size1();
	dimy = out.size2();

	for (uint y = 0; y < dimy; y++) {
		for(uint x = 0;x<dimx;x++)
		{
			out(x,y).real()=res[2*y*dimx+2*x];
			out(x,y).imag()=res[2*y*dimx+2*x+1];
		}

	}


	delete[] mat1;
	delete[] mat2;
	delete[] res;


}

//tested and working
void MatMatProd(const Complex alpha,const CMatrixTypeColMaj & Ar1,const CMatrixTypeColMaj &Ar2,CMatrixTypeColMaj & out, char tra,char trb,bool test)
{
	double* mat1 = new double[Ar1.size1() * Ar1.size2()*2];
	double* mat2 = new double[Ar2.size1() * Ar2.size2()* 2];


	//std::cout<<"after alloc"<<std::endl;
	//std::cout<<((tra=='N'||tra=='n')?Ar1.shape()[0]:Ar1.shape()[1]) * 2 * ((trb=='N'||trb=='n')?Ar2.shape()[1]:Ar2.shape()[0])<<std::endl;
	double* res = new double[((tra=='N'||tra=='n')?Ar1.size1():Ar1.size2()) * 2 * ((trb=='N'||trb=='n')?Ar2.size2():Ar2.size1())];


	for (uint cnt = 0; cnt < Ar1.size1() * Ar1.size2(); cnt++) {
		mat1[2 * cnt] = Ar1.data()[cnt].real();
		mat1[2 * cnt + 1] = Ar1.data()[cnt].imag();
	}
	for (uint cnt = 0; cnt < Ar2.size1() * Ar2.size2(); cnt++) {
			mat2[2 * cnt] = Ar2.data()[cnt].real();
			mat2[2 * cnt + 1] = Ar2.data()[cnt].imag();
	}
	//std::cout<<"before fortran routine"<<std::endl;
	cmat_mat_prod(alpha,mat1, Ar1.size1(), Ar1.size2(), mat2,
			Ar2.size1(), Ar2.size2(), res, tra, trb);
	//std::cout<<"after fortran routine"<<std::endl;

	out.resize((tra=='N'||tra=='n')?Ar1.size1():Ar1.size2(),(trb=='N'||trb=='n')?Ar2.size2():Ar2.size1());
	if(test){
		std::cout<<"start copying elements"<<std::endl;
		//my_pause();
	}
	for (uint cnt = 0; cnt < out.size1()*out.size2(); cnt++) {
		out.data()[cnt].real() = res[2*cnt];
		out.data()[cnt].imag() = res[2*cnt+1];

	}
	if(test){
		std::cout<<"end copying elements"<<std::endl;
		//my_pause();
	}

	delete[] mat1;
	delete[] mat2;
	delete[] res;


}

void MatMatProd(const Complex alpha,const VectorType& Ar1,const CMatrixType&Ar2,CMatrixType& out, char tra,char trb)
{
	double* mat1 = new double[Ar1.size() * Ar1.size()*2];
	double* mat2 = new double[Ar2.size1() * Ar2.size2()*2 ];


	//std::cout<<"after alloc"<<std::endl;
	//std::cout<<((tra=='N'||tra=='n')?Ar1.shape()[0]:Ar1.shape()[1]) * 2 * ((trb=='N'||trb=='n')?Ar2.shape()[1]:Ar2.shape()[0])<<std::endl;
	double* res = new double[((tra=='N'||tra=='n')?Ar1.size():Ar1.size()) * 2 * ((trb=='N'||trb=='n')?Ar2.size2():Ar2.size1())];

	uint dimx = Ar1.size();
	uint dimy = Ar1.size();

	for (uint y = 0; y < dimy; y++) {
		for(uint x = 0;x<dimx;x++)
		{
			if(x==y)
			{
				mat1[2*y*dimx+2*x] = Ar1(x);
				mat1[2*y*dimx+2*x+1] = 0;
			}
			else{
				mat1[2*y*dimx+2*x] = 0;
				mat1[2*y*dimx+2*x+1] = 0;
			}
		}

	}

	dimx = Ar2.size1();
	dimy = Ar2.size2();

	for (uint y = 0; y < dimy; y++) {
		for(uint x = 0;x<dimx;x++)
		{


				mat2[2*y*dimx+2*x] = Ar2(x,y).real();
				mat2[2*y*dimx+2*x+1] =Ar2(x,y).imag();


		}

	}


	//std::cout<<"before fortran routine"<<std::endl;
	cmat_mat_prod(alpha,mat1, Ar1.size(), Ar1.size(), mat2,
			Ar2.size1(), Ar2.size2(), res, tra, trb);
	//std::cout<<"after fortran routine"<<std::endl;

	out.resize((tra=='N'||tra=='n')?Ar1.size():Ar1.size(),(trb=='N'||trb=='n')?Ar2.size2():Ar2.size1());


	dimx = out.size1();
	dimy = out.size2();

	for (uint y = 0; y < dimy; y++) {
		for(uint x = 0;x<dimx;x++)
		{
			out(x,y).real()=res[2*y*dimx+2*x];
			out(x,y).imag()=res[2*y*dimx+2*x+1];
		}

	}


	delete[] mat1;
	delete[] mat2;
	delete[] res;

}


void MatMatProd(const Complex alpha,const CMatrixType& Ar1,const VectorType&Ar2,CMatrixType& out, char tra,char trb)
{
	double* mat1 = new double[Ar1.size1() * Ar1.size2()*2];
	double* mat2 = new double[Ar2.size() * Ar2.size()* 2];


	//std::cout<<"after alloc"<<std::endl;
	//std::cout<<((tra=='N'||tra=='n')?Ar1.shape()[0]:Ar1.shape()[1]) * 2 * ((trb=='N'||trb=='n')?Ar2.shape()[1]:Ar2.shape()[0])<<std::endl;
	double* res = new double[((tra=='N'||tra=='n')?Ar1.size1():Ar1.size2()) * 2 * ((trb=='N'||trb=='n')?Ar2.size():Ar2.size())];

	uint dimx = Ar1.size1();
	uint dimy = Ar1.size2();

	for (uint y = 0; y < dimy; y++) {
		for(uint x = 0;x<dimx;x++)
		{
			mat1[2*y*dimx+2*x] = Ar1(x,y).real();
			mat1[2*y*dimx+2*x+1] = Ar1(x,y).imag();
		}

	}
	dimx = Ar2.size();
	dimy = Ar2.size();

	for (uint y = 0; y < dimy; y++) {
		for(uint x = 0;x<dimx;x++)
		{
			if(x==y)
			{
				mat2[2*y*dimx+2*x] = Ar2(x);
				mat2[2*y*dimx+2*x+1] = 0;
			}
			else{
				mat2[2*y*dimx+2*x] = 0;
				mat2[2*y*dimx+2*x+1] = 0;
			}
		}

	}


	//std::cout<<"before fortran routine"<<std::endl;
	cmat_mat_prod(alpha,mat1, Ar1.size1(), Ar1.size2(), mat2,
			Ar2.size(), Ar2.size(), res, tra, trb);
	//std::cout<<"after fortran routine"<<std::endl;

	out.resize((tra=='N'||tra=='n')?Ar1.size1():Ar1.size2(),(trb=='N'||trb=='n')?Ar2.size():Ar2.size());


	dimx = out.size1();
	dimy = out.size2();

	for (uint y = 0; y < dimy; y++) {
		for(uint x = 0;x<dimx;x++)
		{
			out(x,y).real()=res[2*y*dimx+2*x];
			out(x,y).imag()=res[2*y*dimx+2*x+1];
		}

	}


	delete[] mat1;
	delete[] mat2;
	delete[] res;

}

void MatMatProd(const Complex alpha,const CMatrixType & Ar1,const CMatrixType &Ar2,CMatrixType & out, char tra,char trb)
{
	double* mat1 = new double[Ar1.size1() * Ar1.size2()*2];
	double* mat2 = new double[Ar2.size1() * Ar2.size2()* 2];


	//std::cout<<"after alloc"<<std::endl;
	//std::cout<<((tra=='N'||tra=='n')?Ar1.shape()[0]:Ar1.shape()[1]) * 2 * ((trb=='N'||trb=='n')?Ar2.shape()[1]:Ar2.shape()[0])<<std::endl;
	double* res = new double[((tra=='N'||tra=='n')?Ar1.size1():Ar1.size2()) * 2 * ((trb=='N'||trb=='n')?Ar2.size2():Ar2.size1())];

	uint dimx = Ar1.size1();
	uint dimy = Ar1.size2();

	for (uint y = 0; y < dimy; y++) {
		for(uint x = 0;x<dimx;x++)
		{
			mat1[2*y*dimx+2*x] = Ar1(x,y).real();
			mat1[2*y*dimx+2*x+1] = Ar1(x,y).imag();
		}

	}


	dimx = Ar2.size1();
	dimy = Ar2.size2();

	for (uint y = 0; y < dimy; y++) {
		for(uint x = 0;x<dimx;x++)
		{
			mat2[2*y*dimx+2*x] = Ar2(x,y).real();
			mat2[2*y*dimx+2*x+1] = Ar2(x,y).imag();
		}

	}
	//std::cout<<"before fortran routine"<<std::endl;
	cmat_mat_prod(alpha,mat1, Ar1.size1(), Ar1.size2(), mat2,
			Ar2.size1(), Ar2.size2(), res, tra, trb);
	//std::cout<<"after fortran routine"<<std::endl;

	out.resize((tra=='N'||tra=='n')?Ar1.size1():Ar1.size2(),(trb=='N'||trb=='n')?Ar2.size2():Ar2.size1());


	dimx = out.size1();
	dimy = out.size2();

	for (uint y = 0; y < dimy; y++) {
		for(uint x = 0;x<dimx;x++)
		{
			out(x,y).real()=res[2*y*dimx+2*x];
			out(x,y).imag()=res[2*y*dimx+2*x+1];
		}

	}


	delete[] mat1;
	delete[] mat2;
	delete[] res;


}






//=================================================================Real matrix products===============================================



//this is tested and working!!!!!
void MatMatProd(const double alpha,const MatrixTypeColMaj & Ar1,const MatrixTypeColMaj &Ar2,MatrixTypeColMaj & out, char tra,char trb,bool test)
{
	double *mat1,*mat2;
	double *res;
	mat1 = (double*)Ar1.data().begin();
	mat2 = (double*)Ar2.data().begin();

	out.resize((tra=='N'||tra=='n')?Ar1.size1():Ar1.size2(),(trb=='N'||trb=='n')?Ar2.size2():Ar2.size1());
	res = out.data().begin();

	mat_mat_prod(alpha,mat1, Ar1.size1(), Ar1.size2(), mat2,
			Ar2.size1(), Ar2.size2(), res, tra, trb);
}
void MatMatProd_2(const double alpha,const MatrixTypeColMaj & Ar1,const MatrixTypeColMaj &Ar2,MatrixTypeColMaj & out, char tra,char trb,bool test)
{
	double *mat1,*mat2;
	double *res;
	mat1 = (double*)Ar1.data().begin();
	mat2 = (double*)Ar2.data().begin();

	out.resize((tra=='N'||tra=='n')?Ar1.size1():Ar1.size2(),(trb=='N'||trb=='n')?Ar2.size2():Ar2.size1());
	res = out.data().begin();

	mat_mat_prod_2(alpha,mat1, Ar1.size1(), Ar1.size2(), mat2,
			Ar2.size1(), Ar2.size2(), res, tra, trb);
}

void MatMatProd(const double alpha,const VectorType & Ar1,const MatrixTypeColMaj &Ar2,MatrixTypeColMaj & out, char tra,char trb)
{
	int dim1 = Ar1.size();
	double *mat1,*mat2;
	double *res;
	mat1 = new double[dim1*dim1];
	for (uint y = 0; y < dim1; y++) {
		for(uint x = 0;x<dim1;x++)
		{
			if(x==y)
			{
				mat1[y*dim1+x] = Ar1(x);
			}
			else{
				mat1[y*dim1+x] = 0;
			}
		}

	}
	mat2 = (double*)Ar2.data().begin();
	out.resize(Ar1.size(),(trb=='N'||trb=='n')?Ar2.size2():Ar2.size1());
	boost_setTo(out,0.0);
	res = out.data().begin();

	mat_mat_prod(alpha,mat1,dim1, dim1, mat2,
			Ar2.size1(), Ar2.size2(), res, tra, trb);
	delete[] mat1;
}
void MatMatProd(const double alpha,const MatrixTypeColMaj & Ar1,const VectorType &Ar2,MatrixTypeColMaj & out, char tra,char trb)
{
	double *mat1;
	uint dim2 = Ar2.size();
	double *res;
	mat1 = (double*)Ar1.data().begin();
	double *mat2 = new double[dim2*dim2];
	for (uint y = 0; y < dim2; y++) {
		for(uint x = 0;x<dim2;x++)
		{
			if(x==y)
			{
				mat2[y*dim2+x] = Ar2(x);
			}
			else{
				mat2[y*dim2+x] = 0;
			}
		}

	}

	out.resize((tra=='N'||tra=='n')?Ar1.size1():Ar1.size2(),dim2);
	boost_setTo(out,0.0);
	res = out.data().begin();

	mat_mat_prod(alpha,mat1, Ar1.size1(), Ar1.size2(), mat2,
			Ar2.size(), Ar2.size(), res, tra, trb);
	delete[] mat2;
}
void MatMatProd(const double alpha,const VectorType& Ar1,const MatrixType&Ar2,MatrixType& out, char tra,char trb)
{
	double* mat1 = new double[Ar1.size() * Ar1.size()];
	double* mat2 = new double[Ar2.size1() * Ar2.size2()];


	//std::cout<<"after alloc"<<std::endl;
	//std::cout<<((tra=='N'||tra=='n')?Ar1.shape()[0]:Ar1.shape()[1]) * ((trb=='N'||trb=='n')?Ar2.shape()[1]:Ar2.shape()[0])<<std::endl;
	double* res = new double[((tra=='N'||tra=='n')?Ar1.size():Ar1.size()) * ((trb=='N'||trb=='n')?Ar2.size2():Ar2.size1())];

	uint dimx = Ar1.size();
	uint dimy = Ar1.size();

	for (uint y = 0; y < dimy; y++) {
		for(uint x = 0;x<dimx;x++)
		{
			if(x==y)
			{
				mat1[y*dimx+x] = Ar1(x);
			}
			else{
				mat1[y*dimx+x] = 0;
			}
		}

	}

	dimx = Ar2.size1();
	dimy = Ar2.size2();

	for (uint y = 0; y < dimy; y++) {
		for(uint x = 0;x<dimx;x++)
		{


				mat2[y*dimx+x] = Ar2(x,y);


		}

	}


	//std::cout<<"before fortran routine"<<std::endl;
	mat_mat_prod(alpha,mat1, Ar1.size(), Ar1.size(), mat2,
			Ar2.size1(), Ar2.size2(), res, tra, trb);
	//std::cout<<"after fortran routine"<<std::endl;

	out.resize((tra=='N'||tra=='n')?Ar1.size():Ar1.size(),(trb=='N'||trb=='n')?Ar2.size2():Ar2.size1());


	dimx = out.size1();
	dimy = out.size2();

	for (uint y = 0; y < dimy; y++) {
		for(uint x = 0;x<dimx;x++)
		{
			out(x,y)=res[y*dimx+x];
		}

	}


	delete[] mat1;
	delete[] mat2;
	delete[] res;

}


void MatMatProd(const double alpha,const MatrixType& Ar1,const VectorType&Ar2,MatrixType& out, char tra,char trb)
{

	double* mat1 = new double[Ar1.size1() * Ar1.size2()];
		double* mat2 = new double[Ar2.size() * Ar2.size()];



		//std::cout<<"after alloc"<<std::endl;
		//std::cout<<((tra=='N'||tra=='n')?Ar1.shape()[0]:Ar1.shape()[1]) * 2 * ((trb=='N'||trb=='n')?Ar2.shape()[1]:Ar2.shape()[0])<<std::endl;
		double* res = new double[((tra=='N'||tra=='n')?Ar1.size1():Ar1.size2()) *((trb=='N'||trb=='n')?Ar2.size():Ar2.size())];

		uint dimx = Ar1.size1();
		uint dimy = Ar1.size2();

		for (uint y = 0; y < dimy; y++) {
			for(uint x = 0;x<dimx;x++)
			{
				mat1[y*dimx+x] = Ar1(x,y);
			}

		}


		dimx = Ar2.size();
		dimy = Ar2.size();

		for (uint y = 0; y < dimy; y++) {
			for(uint x = 0;x<dimx;x++)
			{
				if(x==y){mat2[y*dimx+x] = Ar2(x);}
				else{mat2[y*dimx+x]=0;}
			}

		}
		//std::cout<<"before fortran routine"<<std::endl;
		mat_mat_prod(alpha,mat1, Ar1.size1(), Ar1.size2(), mat2,
				Ar2.size(), Ar2.size(), res, tra, trb);
		//std::cout<<"after fortran routine"<<std::endl;

		out.resize((tra=='N'||tra=='n')?Ar1.size1():Ar1.size2(),(trb=='N'||trb=='n')?Ar2.size():Ar2.size());


		dimx = out.size1();
		dimy = out.size2();

		for (uint y = 0; y < dimy; y++) {
			for(uint x = 0;x<dimx;x++)
			{
				out(x,y)=res[y*dimx+x];
			}

		}


		delete[] mat1;
		delete[] mat2;
		delete[] res;




}

void MatMatProd(const double alpha,const MatrixType & Ar1,const MatrixType &Ar2,MatrixType & out, char tra,char trb)
{
	double* mat1 = new double[Ar1.size1() * Ar1.size2()];
	double* mat2 = new double[Ar2.size1() * Ar2.size2()];



	//std::cout<<"after alloc"<<std::endl;
	//std::cout<<((tra=='N'||tra=='n')?Ar1.shape()[0]:Ar1.shape()[1]) * 2 * ((trb=='N'||trb=='n')?Ar2.shape()[1]:Ar2.shape()[0])<<std::endl;
	double* res = new double[((tra=='N'||tra=='n')?Ar1.size1():Ar1.size2()) *((trb=='N'||trb=='n')?Ar2.size2():Ar2.size1())];

	uint dimx = Ar1.size1();
	uint dimy = Ar1.size2();

	for (uint y = 0; y < dimy; y++) {
		for(uint x = 0;x<dimx;x++)
		{
			mat1[y*dimx+x] = Ar1(x,y);
		}

	}


	dimx = Ar2.size1();
	dimy = Ar2.size2();

	for (uint y = 0; y < dimy; y++) {
		for(uint x = 0;x<dimx;x++)
		{
			mat2[y*dimx+x] = Ar2(x,y);
		}

	}
	//std::cout<<"before fortran routine"<<std::endl;
	mat_mat_prod(alpha,mat1, Ar1.size1(), Ar1.size2(), mat2,
			Ar2.size1(), Ar2.size2(), res, tra, trb);
	//std::cout<<"after fortran routine"<<std::endl;

	out.resize((tra=='N'||tra=='n')?Ar1.size1():Ar1.size2(),(trb=='N'||trb=='n')?Ar2.size2():Ar2.size1());


	dimx = out.size1();
	dimy = out.size2();

	for (uint y = 0; y < dimy; y++) {
		for(uint x = 0;x<dimx;x++)
		{
			out(x,y)=res[y*dimx+x];
		}

	}


	delete[] mat1;
	delete[] mat2;
	delete[] res;


}
extern "C"
{
	/*
  SUBROUTINE DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
*     .. Scalar Arguments ..
      DOUBLE PRECISION ALPHA,BETA
      INTEGER K,LDA,LDB,LDC,M,N
      CHARACTER TRANSA,TRANSB
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*),B(LDB,*),C(LDC,*)
*     ..
*
*  Purpose
*  =======
*
*  DGEMM  performs one of the matrix-matrix operations
*
*     C := alpha*op( A )*op( B ) + beta*C,
*
*  where  op( X ) is one of
*
*     op( X ) = X   or   op( X ) = X',
*
*  alpha and beta are scalars, and A, B and C are matrices, with op( A )
*  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
*
*  Arguments
*  ==========
*
*  TRANSA - CHARACTER*1.
*           On entry, TRANSA specifies the form of op( A ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSA = 'N' or 'n',  op( A ) = A.
*
*              TRANSA = 'T' or 't',  op( A ) = A'.
*
*              TRANSA = 'C' or 'c',  op( A ) = A'.
*
*           Unchanged on exit.
*
*  TRANSB - CHARACTER*1.
*           On entry, TRANSB specifies the form of op( B ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSB = 'N' or 'n',  op( B ) = B.
*
*              TRANSB = 'T' or 't',  op( B ) = B'.
*
*              TRANSB = 'C' or 'c',  op( B ) = B'.
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry,  M  specifies  the number  of rows  of the  matrix
*           op( A )  and of the  matrix  C.  M  must  be at least  zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry,  N  specifies the number  of columns of the matrix
*           op( B ) and the number of columns of the matrix C. N must be
*           at least zero.
*           Unchanged on exit.
*
*  K      - INTEGER.
*           On entry,  K  specifies  the number of columns of the matrix
*           op( A ) and the number of rows of the matrix op( B ). K must
*           be at least  zero.
*           Unchanged on exit.
*
*  ALPHA  - DOUBLE PRECISION.
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
*           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
*           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
*           part of the array  A  must contain the matrix  A,  otherwise
*           the leading  k by m  part of the array  A  must contain  the
*           matrix A.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
*           LDA must be at least  max( 1, m ), otherwise  LDA must be at
*           least  max( 1, k ).
*           Unchanged on exit.
*
*  B      - DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is
*           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
*           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
*           part of the array  B  must contain the matrix  B,  otherwise
*           the leading  n by k  part of the array  B  must contain  the
*           matrix B.
*           Unchanged on exit.
*
*  LDB    - INTEGER.
*           On entry, LDB specifies the first dimension of B as declared
*           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
*           LDB must be at least  max( 1, k ), otherwise  LDB must be at
*           least  max( 1, n ).
*           Unchanged on exit.
*
*  BETA   - DOUBLE PRECISION.
*           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
*           supplied as zero then C need not be set on input.
*           Unchanged on exit.
*
*  C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).
*           Before entry, the leading  m by n  part of the array  C must
*           contain the matrix  C,  except when  beta  is zero, in which
*           case C need not be set on entry.
*           On exit, the array  C  is overwritten by the  m by n  matrix
*           ( alpha*op( A )*op( B ) + beta*C ).
*
*  LDC    - INTEGER.
*           On entry, LDC specifies the first dimension of C as declared
*           in  the  calling  (sub)  program.   LDC  must  be  at  least
*           max( 1, m ).
*           Unchanged on exit.
*/
  void dgemm_(char* TRANSA,char *TRANSB,int* M,int* N,int * K,double* ALPHA,const double *A,int* LDA,const double *B,int *LDB,double* BETA,double *C,int* LDC);

}
//no check if dimensions of matrices match!!!! mat_out has to be allocated before it is passed to cmat_mat_prod!
void mat_mat_prod(double alpha,double*matA ,int dim1A,int dim2A, double* matB, int dim1B, int dim2B,double* mat_out,char tra,char trb)
{
	char Transa = tra,Transb=trb;

	int M,N,Ka,Kb;

	//double alpha = 1;
	int lda = dim1A,ldb=dim1B,ldc;
	if(tra=='n'||tra=='N')
	{
		M = dim1A;
		Ka = dim2A;
	}
	else
	{
		M = dim2A;
		Ka = dim1A;
	}
	if(trb =='n'||trb=='N') {
		N = dim2B;
		Kb = dim1B;
	}
	else{
		N = dim1B;
		Kb = dim2B;
	}
	ldc = M;
	assert(Ka==Kb);


	double beta = 0;

	dgemm_(&Transa,&Transb,&M,&N,&Ka,&alpha,matA,&lda,matB,&ldb,&beta,mat_out,&ldc);
}


//this just incorporates an additonal number alpha by which the product may be multiplied
void mat_mat_prod_2(double alpha,double*matA ,int dim1A,int dim2A, double* matB, int dim1B, int dim2B,double* mat_out,char tra,char trb)
{
	char Transa = tra,Transb=trb;

	int M,N,Ka,Kb;

	//	double alpha = 1;
	int lda = dim1A,ldb=dim1B,ldc;
	if(tra=='n'||tra=='N')
	{
		M = dim1A;
		Ka = dim2A;
	}
	else
	{
		M = dim2A;
		Ka = dim1A;
	}
	if(trb =='n'||trb=='N') {
		N = dim2B;
		Kb = dim1B;
	}
	else{
		N = dim1B;
		Kb = dim2B;
	}
	ldc = M;
	assert(Ka==Kb);


	double beta = 0;

	dgemm_(&Transa,&Transb,&M,&N,&Ka,&alpha,matA,&lda,matB,&ldb,&beta,mat_out,&ldc);
}
