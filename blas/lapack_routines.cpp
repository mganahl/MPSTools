#include"lapack_routines.hpp"

//Produces SVD of array matrix. Return value is return variable info of zgesvd_. if info = 0, zgesvd_ has
//succesfully finished the job.
//U contains left singular vectors of array matrix, VT contains right singular vectors of matrix,
//D contains singular values of array matrix in decreasing order.
//Note: matrix = U*D*VT
//If matrix is m by n array with m<n, then U is m by m, D is m by m, VT is m by n.
//If matrix is m by n with m>n, then U is m by n, D is n by n, VT is n by n
//routine produces the same matrices as matlab routine svd(mat,'econ').
extern "C"
{

  /*SUBROUTINE ZGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT,
    $                   WORK, LWORK, RWORK, INFO )
    *
    *  -- LAPACK driver routine (version 3.2) --
    *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    *     November 2006
    *
    *     .. Scalar Arguments ..
    CHARACTER          JOBU, JOBVT
    INTEGER            INFO, LDA, LDU, LDVT, LWORK, M, N
    *     ..
    *     .. Array Arguments ..
    DOUBLE PRECISION   RWORK( * ), S( * )
    COMPLEX*16         A( LDA, * ), U( LDU, * ), VT( LDVT, * ),
    $                   WORK( * )
    *     ..
    *
    *  Purpose
    *  =======
    *
    *  ZGESVD computes the singular value decomposition (SVD) of a complex
    *  M-by-N matrix A, optionally computing the left and/or right singular
    *  vectors. The SVD is written
    *
    *       A = U * SIGMA * conjugate-transpose(V)
    *
    *  where SIGMA is an M-by-N matrix which is zero except for its
    *  min(m,n) diagonal elements, U is an M-by-M unitary matrix, and
    *  V is an N-by-N unitary matrix.  The diagonal elements of SIGMA
    *  are the singular values of A; they are real and non-negative, and
    *  are returned in descending order.  The first min(m,n) columns of
    *  U and V are the left and right singular vectors of A.
    *
    *  Note that the routine returns V**H, not V.
    *
    *  Arguments
    *  =========
    *
    *  JOBU    (input) CHARACTER*1
    *          Specifies options for computing all or part of the matrix U:
    *          = 'A':  all M columns of U are returned in array U:
    *          = 'S':  the first min(m,n) columns of U (the left singular
    *                  vectors) are returned in the array U;
    *          = 'O':  the first min(m,n) columns of U (the left singular
    *                  vectors) are overwritten on the array A;
    *          = 'N':  no columns of U (no left singular vectors) are
    *                  computed.
    *
    *  JOBVT   (input) CHARACTER*1
    *          Specifies options for computing all or part of the matrix
    *          V**H:
    *          = 'A':  all N rows of V**H are returned in the array VT;
    *          = 'S':  the first min(m,n) rows of V**H (the right singular
    *                  vectors) are returned in the array VT;
    *          = 'O':  the first min(m,n) rows of V**H (the right singular
    *                  vectors) are overwritten on the array A;
    *          = 'N':  no rows of V**H (no right singular vectors) are
    *                  computed.
    *
    *          JOBVT and JOBU cannot both be 'O'.
    *
    *  M       (input) INTEGER
    *          The number of rows of the input matrix A.  M >= 0.
    *
    *  N       (input) INTEGER
    *          The number of columns of the input matrix A.  N >= 0.
    *
    *  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
    *          On entry, the M-by-N matrix A.
    *          On exit,
    *          if JOBU = 'O',  A is overwritten with the first min(m,n)
    *                          columns of U (the left singular vectors,
    *                          stored columnwise);
    *          if JOBVT = 'O', A is overwritten with the first min(m,n)
    *                          rows of V**H (the right singular vectors,
    *                          stored rowwise);
    *          if JOBU .ne. 'O' and JOBVT .ne. 'O', the contents of A
    *                          are destroyed.
    *
    *  LDA     (input) INTEGER
    *          The leading dimension of the array A.  LDA >= max(1,M).
    *
    *  S       (output) DOUBLE PRECISION array, dimension (min(M,N))
    *          The singular values of A, sorted so that S(i) >= S(i+1).
    *
    *  U       (output) COMPLEX*16 array, dimension (LDU,UCOL)
    *          (LDU,M) if JOBU = 'A' or (LDU,min(M,N)) if JOBU = 'S'.
    *          If JOBU = 'A', U contains the M-by-M unitary matrix U;
    *          if JOBU = 'S', U contains the first min(m,n) columns of U
    *          (the left singular vectors, stored columnwise);
    *          if JOBU = 'N' or 'O', U is not referenced.
    *
    *  LDU     (input) INTEGER
    *          The leading dimension of the array U.  LDU >= 1; if
    *          JOBU = 'S' or 'A', LDU >= M.
    *
    *  VT      (output) COMPLEX*16 array, dimension (LDVT,N)
    *          If JOBVT = 'A', VT contains the N-by-N unitary matrix
    *          V**H;
    *          if JOBVT = 'S', VT contains the first min(m,n) rows of
    *          V**H (the right singular vectors, stored rowwise);
    *          if JOBVT = 'N' or 'O', VT is not referenced.
    *
    *  LDVT    (input) INTEGER
    *          The leading dimension of the array VT.  LDVT >= 1; if
    *          JOBVT = 'A', LDVT >= N; if JOBVT = 'S', LDVT >= min(M,N).
    *
    *  WORK    (workspace/output) COMPLEX*16 array, dimension (MAX(1,LWORK))
    *          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
    *
    *  LWORK   (input) INTEGER
    *          The dimension of the array WORK.
    *          LWORK >=  MAX(1,2*MIN(M,N)+MAX(M,N)).OMPLEX*16 array,
    *          For good performance, LWORK should generally be larger.
    *
    *          If LWORK = -1, then a workspace query is assumed; the routine
    *          only calculates the optimal size of the WORK array, returns
    *          this value as the first entry of the WORK array, and no error
    *          message related to LWORK is issued by XERBLA.
    *
    *  RWORK   (workspace) DOUBLE PRECISION array, dimension (5*min(M,N))
    *          On exit, if INFO > 0, RWORK(1:MIN(M,N)-1) contains the
    *          unconverged superdiagonal elements of an upper bidiagonal
    *          matrix B whose diagonal is in S (not necessarily sorted).
    *          B satisfies A = U * B * VT, so it has the same singular
    *          values as A, and singular vectors related by U and VT.
    *
    *  INFO    (output) INTEGER
    *          = 0:  successful exit.
    *          < 0:  if INFO = -i, the i-th argument had an illegal value.
    *          > 0:  if ZBDSQR did not converge, INFO specifies how many
    *                superdiagonals of an intermediate bidiagonal form B
    *                did not converge to zero. See the description of RWORK
    *                above for details.
    *
    *  =====================================================================
    */

  void zgesvd_( char* JOBU, char* JOBVT, int* M, int *N, double *A, int *LDA, double * S, double * U, int* LDU, double *VT, int *LDVT, double *WORK, int *LWORK, double *RWORK,int * INFO );

}

//memory for all inputs and outputs has to be allocated outside svd!!! dimensions must match!!
int svd(double *matrix, double *U,double*VT,double *D, int dim1,int dim2)
{

  //int lwork_extention = 10;
  //extents the minimum amount of memory needed by zgesvd_ by a factor 2. See Lapack documentation for more information
  char JOBU = 'S';
  char JOBVT = 'S';

  int lda = dim1,ldu = dim1;
  int ldvt = (dim1<=dim2?dim1:dim2);

  int info;

  //compute optimal lwork
  int lwork = -1;
  double *rwork = new double [5 * (dim1<=dim2?dim1:dim2)];
  double *temp = new double[2];
  zgesvd_(&JOBU,&JOBVT, &dim1, &dim2, matrix, &lda, D,U, &ldu, VT, &ldvt, temp, &lwork, rwork,&info );
  lwork =static_cast<int>(temp[0]);
  delete[] temp;


  double *work = new double[2*lwork];
  zgesvd_( &JOBU,&JOBVT, &dim1, &dim2, matrix, &lda, D,U, &ldu, VT, &ldvt, work, &lwork, rwork,&info );

  delete[] work;
  delete[] rwork;

  return info;
}

extern "C"
{
  /*
    SUBROUTINE DGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT,
    $                   WORK, LWORK, INFO )
    *
    *  -- LAPACK driver routine (version 3.2) --
    *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    *     November 2006
    *
    *     .. Scalar Arguments ..
    CHARACTER          JOBU, JOBVT
    INTEGER            INFO, LDA, LDU, LDVT, LWORK, M, N
    *     ..
    *     .. Array Arguments ..
    DOUBLE PRECISION   A( LDA, * ), S( * ), U( LDU, * ),
    $                   VT( LDVT, * ), WORK( * )
    *     ..
    *
    *  Purpose
    *  =======
    *
    *  DGESVD computes the singular value decomposition (SVD) of a real
    *  M-by-N matrix A, optionally computing the left and/or right singular
    *  vectors. The SVD is written
    *
    *       A = U * SIGMA * transpose(V)
    *
    *  where SIGMA is an M-by-N matrix which is zero except for its
    *  min(m,n) diagonal elements, U is an M-by-M orthogonal matrix, and
    *  V is an N-by-N orthogonal matrix.  The diagonal elements of SIGMA
    *  are the singular values of A; they are real and non-negative, and
    *  are returned in descending order.  The first min(m,n) columns of
    *  U and V are the left and right singular vectors of A.
    *
    *  Note that the routine returns V**T, not V.
    *
    *  Arguments
    *  =========
    *
    *  JOBU    (input) CHARACTER*1
    *          Specifies options for computing all or part of the matrix U:
    *          = 'A':  all M columns of U are returned in array U:
    *          = 'S':  the first min(m,n) columns of U (the left singular
    *                  vectors) are returned in the array U;
    *          = 'O':  the first min(m,n) columns of U (the left singular
    *                  vectors) are overwritten on the array A;
    *          = 'N':  no columns of U (no left singular vectors) are
    *                  computed.
    *
    *  JOBVT   (input) CHARACTER*1
    *          Specifies options for computing all or part of the matrix
    *          V**T:
    *          = 'A':  all N rows of V**T are returned in the array VT;
    *          = 'S':  the first min(m,n) rows of V**T (the right singular
    *                  vectors) are returned in the array VT;
    *          = 'O':  the first min(m,n) rows of V**T (the right singular
    *                  vectors) are overwritten on the array A;
    *          = 'N':  no rows of V**T (no right singular vectors) are
    *                  computed.
    *
    *          JOBVT and JOBU cannot both be 'O'.
    *
    *  M       (input) INTEGER
    *          The number of rows of the input matrix A.  M >= 0.
    *
    *  N       (input) INTEGER
    *          The number of columns of the input matrix A.  N >= 0.
    *
    *  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
    *          On entry, the M-by-N matrix A.
    *          On exit,
    *          if JOBU = 'O',  A is overwritten with the first min(m,n)
    *                          columns of U (the left singular vectors,
    *                          stored columnwise);
    *          if JOBVT = 'O', A is overwritten with the first min(m,n)
    *                          rows of V**T (the right singular vectors,
    *                          stored rowwise);
    *          if JOBU .ne. 'O' and JOBVT .ne. 'O', the contents of A
    *                          are destroyed.
    *
    *  LDA     (input) INTEGER
    *          The leading dimension of the array A.  LDA >= max(1,M).
    *
    *  S       (output) DOUBLE PRECISION array, dimension (min(M,N))
    *          The singular values of A, sorted so that S(i) >= S(i+1).
    *
    *  U       (output) DOUBLE PRECISION array, dimension (LDU,UCOL)
    *          (LDU,M) if JOBU = 'A' or (LDU,min(M,N)) if JOBU = 'S'.
    *          If JOBU = 'A', U contains the M-by-M orthogonal matrix U;
    *          if JOBU = 'S', U contains the first min(m,n) columns of U
    *          (the left singular vectors, stored columnwise);
    *          if JOBU = 'N' or 'O', U is not referenced.
    *
    *  LDU     (input) INTEGER
    *          The leading dimension of the array U.  LDU >= 1; if
    *          JOBU = 'S' or 'A', LDU >= M.
    *
    *  VT      (output) DOUBLE PRECISION array, dimension (LDVT,N)
    *          If JOBVT = 'A', VT contains the N-by-N orthogonal matrix
    *          V**T;
    *          if JOBVT = 'S', VT contains the first min(m,n) rows of
    *          V**T (the right singular vectors, stored rowwise);
    *          if JOBVT = 'N' or 'O', VT is not referenced.
    *
    *  LDVT    (input) INTEGER
    *          The leading dimension of the array VT.  LDVT >= 1; if
    *          JOBVT = 'A', LDVT >= N; if JOBVT = 'S', LDVT >= min(M,N).
    *
    *  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
    *          On exit, if INFO = 0, WORK(1) returns the optimal LWORK;
    *          if INFO > 0, WORK(2:MIN(M,N)) contains the unconverged
    *          superdiagonal elements of an upper bidiagonal matrix B
    *          whose diagonal is in S (not necessarily sorted). B
    *          satisfies A = U * B * VT, so it has the same singular values
    *          as A, and singular vectors related by U and VT.
    *
    *  LWORK   (input) INTEGER
    *          The dimension of the array WORK.
    *          LWORK >= MAX(1,3*MIN(M,N)+MAX(M,N),5*MIN(M,N)).
    *          For good performance, LWORK should generally be larger.
    *
    *          If LWORK = -1, then a workspace query is assumed; the routine
    *          only calculates the optimal size of the WORK array, returns
    *          this value as the first entry of the WORK array, and no error
    *          message related to LWORK is issued by XERBLA.
    *
    *  INFO    (output) INTEGER
    *          = 0:  successful exit.
    *          < 0:  if INFO = -i, the i-th argument had an illegal value.
    *          > 0:  if DBDSQR did not converge, INFO specifies how many
    *                superdiagonals of an intermediate bidiagonal form B
    *                did not converge to zero. See the description of WORK
    *                above for details.
    *
    *  =====================================================================


    */

  void dgesvd_(char* JOBU,char* JOBVT,int* M,int*  N,double* A,int* LDA,double* S,double* U,int* LDU,double* VT,int* LDVT,double* WORK,int*  LWORK,int*  INFO );

}
//memory for all inputs and outputs has to be allocated outside rsvd!!! dimensions must match!!

int rsvd(double *matrix, double *U,double*VT,double *sing_val, int dim1,int dim2)
{

  int lwork_extention = 1;
  //extents the minimum amount of memory needed by zgesvd_ by a factor 2. See Lapack documentation for more information
  char JOBU = 'S';
  char JOBVT = 'S';
  int lda = dim1;
  int ldvt = dim1>dim2?dim2:dim1;
  int info;

  int lwork = 1*lwork_extention;
  if((3*(dim1<dim2?dim1:dim2)+(dim1<dim2?dim2:dim1))< 5*(dim1<dim2?dim1:dim2))lwork=5*(dim1<dim2?dim1:dim2);
  if((3*(dim1<dim2?dim1:dim2)+(dim1<dim2?dim2:dim1))>= 5*(dim1<dim2?dim1:dim2))lwork=3*(dim1<dim2?dim1:dim2)+(dim1<dim2?dim2:dim1);


  // workspace variables

  double *work = new double[lwork];



  dgesvd_( &JOBU,&JOBVT, &dim1, &dim2, matrix, &lda, sing_val,U, &lda, VT, &ldvt, work, &lwork,&info );

  delete[] work;




  return info;
}


//Argumente der fortran routine:
/*SUBROUTINE DSTEV( JOBZ, N, D, E, Z, LDZ, WORK, INFO )
 *
 *  -- LAPACK driver routine (version 3.1) --
 *     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
 *     November 2006
 *
 *     .. Scalar Arguments ..
 CHARACTER          JOBZ
 INTEGER            INFO, LDZ, N
 *     ..
 *     .. Array Arguments ..
 DOUBLE PRECISION   D( * ), E( * ), WORK( * ), Z( LDZ, * )
 *     ..
 *
 *  Purpose
 *  =======
 *
 *  DSTEV computes all eigenvalues and, optionally, eigenvectors of a
 *  real symmetric tridiagonal matrix A.
 *
 *  Arguments
 *  =========
 *
 *  JOBZ    (input) CHARACTER*1
 *          = 'N':  Compute eigenvalues only;
 *          = 'V':  Compute eigenvalues and eigenvectors.
 *
 *  N       (input) INTEGER
 *          The order of the matrix.  N >= 0.
 *
 *  D       (input/output) DOUBLE PRECISION array, dimension (N)
 *          On entry, the n diagonal elements of the tridiagonal matrix
 *          A.
 *          On exit, if INFO = 0, the eigenvalues in ascending order.
 *
 *  E       (input/output) DOUBLE PRECISION array, dimension (N-1)
 *          On entry, the (n-1) subdiagonal elements of the tridiagonal
 *          matrix A, stored in elements 1 to N-1 of E.
 *          On exit, the contents of E are destroyed.
 *
 *  Z       (output) DOUBLE PRECISION array, dimension (LDZ, N)
 *          If JOBZ = 'V', then if INFO = 0, Z contains the orthonormal
 *          eigenvectors of the matrix A, with the i-th column of Z
 *          holding the eigenvector associated with D(i).
 *          If JOBZ = 'N', then Z is not referenced.
 *
 *  LDZ     (input) INTEGER
 *          The leading dimension of the array Z.  LDZ >= 1, and if
 *          JOBZ = 'V', LDZ >= max(1,N).
 *
 *  WORK    (workspace) DOUBLE PRECISION array, dimension (max(1,2*N-2))
 *          If JOBZ = 'N', WORK is not referenced.
 *
 *  INFO    (output) INTEGER
 *          = 0:  successful exit
 *          < 0:  if INFO = -i, the i-th argument had an illegal value
 *          > 0:  if INFO = i, the algorithm failed to converge; i
 *                off-diagonal elements of E did not converge to zero.
 */

extern "C" {
extern void dstev_(char *JOBZ, int *N, double *D, double *E, double *Z, int *LDZ,
		double *WORK, int *INFO);

}
//please make sure that eigvec is allocated if you want the eigenvectors; else this causes crash
int tridiag(double *diag,int length_diag, double* sub_diag,double*eigvec,char JOBZ)
{
	//char JOBZ = 'V';
	assert(JOBZ=='V'||JOBZ=='N');
	int LDZ = length_diag;
	double* WORK = new double[2*length_diag-2];
	int info;
	dstev_(&JOBZ,&length_diag,diag,sub_diag,eigvec,&LDZ,WORK,&info);
	delete [] WORK;
	return info;

}
//please make sure that eigvec is allocated if you want the eigenvectors; else this causes crash
int tridiag(Complex *diag,int length_diag, Complex* sub_diag,double*eigvec,char JOBZ)
{
	//char JOBZ = 'V';
	assert(JOBZ=='V'||JOBZ=='N');
	int LDZ = length_diag;
	double* WORK = new double[2*length_diag-2];
	int info;

	double*ndiag = new double[length_diag];
	double*nsdiag = new double[length_diag-1];
	for(int i = 0;i<length_diag;i++)
	{
		ndiag[i] = diag[i].real();
		assert(abs(diag[i].imag())<1e-14);

	}
	for(int i = 0;i<length_diag-1;i++){
		nsdiag[i] = sub_diag[i].real();
		assert(abs(sub_diag[i].imag())<1e-14);
	}

	dstev_(&JOBZ,&length_diag,ndiag,nsdiag,eigvec,&LDZ,WORK,&info);

	for(uint i = 0;i<length_diag;i++){diag[i].real() = ndiag[i];diag[i].imag()=0.0;}

	delete [] nsdiag;
	delete [] ndiag;
	delete [] WORK;

	return info;

}


int BoostSVD(const MatrixType & mat, MatrixType &U,VectorType &D,MatrixType &V_dagger)
{
	int dim1 = mat.size1(),dim2 = mat.size2();
	U.resize(dim1,(dim2<dim1?dim2:dim1));
	boost_setTo(U,0.0);
	V_dagger.resize(dim1<dim2?dim1:dim2,dim2);
	boost_setTo(V_dagger,0.0);
	D.resize(dim1<dim2?dim1:dim2);
	boost_setTo(D,0.0);

	double *matrix = new double[dim1*dim2];

	//dimension of u_in has to be chosen according to zegsvd (see manual):dim(u_in) = (dim1,Min(dim1,dim2)));
	double *u_in = new double[dim1*(dim1<dim2?dim1:dim2)];

	//dimension of vt_in: dim(vt_in) = (Min(dim1,dim2),dim2)
	double *vt_in = new double[(dim1<dim2?dim1:dim2)*dim2];

	//dimension of sing_val: this is a vector of doubles (singular values): dim(sing_val) = min(dim1,dim2)
	double *sing_val = new double[(dim1<dim2?dim1:dim2)];

	for(int cnt2 = 0;cnt2<dim2;cnt2++)
	{
		for(int cnt1 = 0;cnt1<dim1;cnt1++)
		{
			matrix[cnt2*dim1+cnt1] = mat(cnt1,cnt2);
			//matrix[cnt2*2*dim1+2*cnt1+1] = mat(cnt1,cnt2).imag();
		}
	}

	int val= rsvd(matrix, u_in,vt_in,sing_val,dim1,dim2);
/*	for (uint cnt = 0;cnt<2*dim1*2*dim2;cnt++)
	{
		std::cout<<vt_in[cnt]<<" "<<std::endl;
	}
std::cout<<std::endl;
*/
	for(int cnt2 = 0;cnt2<(dim1<dim2?dim1:dim2);cnt2++)
	{
		for(int cnt1 = 0;cnt1<dim1;cnt1++)
		{
			U(cnt1,cnt2)=u_in[cnt2*dim1+cnt1];
		}
	}
	for(int cnt2 = 0;cnt2<dim2;cnt2++)
	{
		for(int cnt1 = 0;cnt1<(dim1<dim2?dim1:dim2);cnt1++)
		{
			V_dagger(cnt1,cnt2)=vt_in[cnt2*(dim1<dim2?dim1:dim2)+cnt1];

		}
	}

		for(int cnt2 = 0;cnt2<(dim1<dim2?dim1:dim2);cnt2++)
		{
			D(cnt2)=sing_val[cnt2];
		}
		 delete[] u_in;
		 delete[] vt_in;
		 delete[] sing_val;
		 delete[] matrix;


	  return val;
}

int BoostSVD(const MatrixType & mat, MatrixType &U,MatrixType &D,MatrixType &V_dagger)
{
	int dim1 = mat.size1(),dim2 = mat.size2();
	U.resize(dim1,(dim2<dim1?dim2:dim1));
	boost_setTo(U,0.0);
	V_dagger.resize(dim1<dim2?dim1:dim2,dim2);
	boost_setTo(V_dagger,0.0);
	D.resize(dim1<dim2?dim1:dim2,dim1<dim2?dim1:dim2);
	boost_setTo(D,0.0);

	double *matrix = new double[dim1*dim2];

	//dimension of u_in has to be chosen according to zegsvd (see manual):dim(u_in) = (dim1,Min(dim1,dim2)));
	double *u_in = new double[dim1*(dim1<dim2?dim1:dim2)];

	//dimension of vt_in: dim(vt_in) = (Min(dim1,dim2),dim2)
	double *vt_in = new double[(dim1<dim2?dim1:dim2)*dim2];

	//dimension of sing_val: this is a vector of doubles (singular values): dim(sing_val) = min(dim1,dim2)
	double *sing_val = new double[(dim1<dim2?dim1:dim2)];

	for(int cnt2 = 0;cnt2<dim2;cnt2++)
	{
		for(int cnt1 = 0;cnt1<dim1;cnt1++)
		{
			matrix[cnt2*dim1+cnt1] = mat(cnt1,cnt2);
			//matrix[cnt2*2*dim1+2*cnt1+1] = mat(cnt1,cnt2).imag();
		}
	}

	int val= rsvd(matrix, u_in,vt_in,sing_val,dim1,dim2);
/*	for (uint cnt = 0;cnt<2*dim1*2*dim2;cnt++)
	{
		std::cout<<vt_in[cnt]<<" "<<std::endl;
	}
std::cout<<std::endl;
*/
	for(int cnt2 = 0;cnt2<(dim1<dim2?dim1:dim2);cnt2++)
	{
		for(int cnt1 = 0;cnt1<dim1;cnt1++)
		{
			U(cnt1,cnt2)=u_in[cnt2*dim1+cnt1];
		}
	}
	for(int cnt2 = 0;cnt2<dim2;cnt2++)
	{
		for(int cnt1 = 0;cnt1<(dim1<dim2?dim1:dim2);cnt1++)
		{
			V_dagger(cnt1,cnt2)=vt_in[cnt2*(dim1<dim2?dim1:dim2)+cnt1];

		}
	}

		for(int cnt2 = 0;cnt2<(dim1<dim2?dim1:dim2);cnt2++)
		{
			D(cnt2,cnt2)=sing_val[cnt2];
		}
		FixPhase(U,V_dagger);
		 delete[] u_in;
		 delete[] vt_in;
		 delete[] sing_val;
		 delete[] matrix;


	  return val;
}
int BoostSVD(MatrixTypeColMaj  mat, MatrixTypeColMaj &U,VectorType &D,MatrixTypeColMaj &V_dagger)
{
	int dim1 = mat.size1(),dim2 = mat.size2();
	U.resize(dim1,(dim2<dim1?dim2:dim1));
	boost_setTo(U,0.0);
	V_dagger.resize(dim1<dim2?dim1:dim2,dim2);
	boost_setTo(V_dagger,0.0);
	D.resize(dim1<dim2?dim1:dim2,dim1<dim2?dim1:dim2);
	boost_setTo(D,0.0);

	double *matrix = (double*)mat.data().begin();
	double *u_in = U.data().begin();
	double *vt_in = V_dagger.data().begin();

	double *sing_val = D.data().begin();
	int val= rsvd(matrix, u_in,vt_in,sing_val,dim1,dim2);
	FixPhase(U,V_dagger);
	  return val;


}

int BoostSVD(MatrixTypeColMaj mat,MatrixTypeColMaj &U,MatrixTypeColMaj &D,MatrixTypeColMaj &V_dagger)
{
	int dim1 = mat.size1(),dim2 = mat.size2();
	U.resize(dim1,(dim2<dim1?dim2:dim1));
	boost_setTo(U,0.0);
	V_dagger.resize(dim1<dim2?dim1:dim2,dim2);
	boost_setTo(V_dagger,0.0);
	D.resize(dim1<dim2?dim1:dim2,dim1<dim2?dim1:dim2);
	boost_setTo(D,0.0);

	double *matrix = (double*)mat.data().begin();
	double *u_in = U.data().begin();
	double *vt_in = V_dagger.data().begin();
	double *sing_val = new double[(dim1<dim2?dim1:dim2)];
	int val= rsvd(matrix, u_in,vt_in,sing_val,dim1,dim2);
		for(int cnt2 = 0;cnt2<(dim1<dim2?dim1:dim2);cnt2++)
		{
			D(cnt2,cnt2)=sing_val[cnt2];
		}
		 delete[] sing_val;
		 FixPhase(U,V_dagger);
	  return val;
}

int BoostSVD(const CMatrixType & mat, CMatrixType &U,VectorType &D,CMatrixType &V_dagger)
{
	int dim1 = mat.size1(),dim2 = mat.size2();
	U.resize(dim1,(dim2<dim1?dim2:dim1));
	boost_setTo(U,Complex(0));
	V_dagger.resize(dim1<dim2?dim1:dim2,dim2);
	boost_setTo(V_dagger,Complex(0));
	D.resize(dim1<dim2?dim1:dim2);
	boost_setTo(D,0.0);

	double *matrix = new double[dim1*2*dim2];

	//dimension of u_in has to be chosen according to zegsvd (see manual):dim(u_in) = (dim1,Min(dim1,dim2)));
	double *u_in = new double[dim1*2*(dim1<dim2?dim1:dim2)];

	//dimension of vt_in: dim(vt_in) = (Min(dim1,dim2),dim2)
	double *vt_in = new double[2*(dim1<dim2?dim1:dim2)*dim2];

	//dimension of sing_val: this is a vector of doubles (singular values): dim(sing_val) = min(dim1,dim2)
	double *sing_val = new double[(dim1<dim2?dim1:dim2)];

	for(int cnt2 = 0;cnt2<dim2;cnt2++)
	{
		for(int cnt1 = 0;cnt1<dim1;cnt1++)
		{
			matrix[cnt2*2*dim1+2*cnt1] = mat(cnt1,cnt2).real();
			matrix[cnt2*2*dim1+2*cnt1+1] = mat(cnt1,cnt2).imag();
		}
	}
	int val= svd(matrix, u_in,vt_in,sing_val,dim1,dim2);

/*	for (uint cnt = 0;cnt<2*dim1*2*dim2;cnt++)
	{
		std::cout<<vt_in[cnt]<<" "<<std::endl;
	}
std::cout<<std::endl;
*/
	for(int cnt2 = 0;cnt2<(dim1<dim2?dim1:dim2);cnt2++)
	{
		for(int cnt1 = 0;cnt1<dim1;cnt1++)
		{
			U(cnt1,cnt2)=Complex(u_in[cnt2*2*dim1+2*cnt1],u_in[cnt2*2*dim1+2*cnt1+1]);
		}
	}
	for(int cnt2 = 0;cnt2<dim2;cnt2++)
	{
		for(int cnt1 = 0;cnt1<(dim1<dim2?dim1:dim2);cnt1++)
		{
			V_dagger(cnt1,cnt2)=Complex(vt_in[cnt2*2*(dim1<dim2?dim1:dim2)+2*cnt1],vt_in[cnt2*2*(dim1<dim2?dim1:dim2)+2*cnt1+1]);

		}
	}

		for(int cnt2 = 0;cnt2<(dim1<dim2?dim1:dim2);cnt2++)
		{
			D(cnt2)=sing_val[cnt2];
		}
		 delete[] u_in;
		 delete[] vt_in;
		 delete[] sing_val;
		 delete[] matrix;

		 FixPhase(U,V_dagger);
	  return val;
}

int BoostSVD(const CMatrixType & mat, CMatrixType &U,MatrixType &D,CMatrixType &V_dagger)
{
	int dim1 = mat.size1(),dim2 = mat.size2();
	U.resize(dim1,(dim2<dim1?dim2:dim1));
	boost_setTo(U,Complex(0));
	V_dagger.resize(dim1<dim2?dim1:dim2,dim2);
	boost_setTo(V_dagger,Complex(0));
	D.resize(dim1<dim2?dim1:dim2,dim1<dim2?dim1:dim2);
	boost_setTo(D,0.0);

	double *matrix = new double[dim1*2*dim2];

	//dimension of u_in has to be chosen according to zegsvd (see manual):dim(u_in) = (dim1,Min(dim1,dim2)));
	double *u_in = new double[dim1*2*(dim1<dim2?dim1:dim2)];

	//dimension of vt_in: dim(vt_in) = (Min(dim1,dim2),dim2)
	double *vt_in = new double[2*(dim1<dim2?dim1:dim2)*dim2];

	//dimension of sing_val: this is a vector of doubles (singular values): dim(sing_val) = min(dim1,dim2)
	double *sing_val = new double[(dim1<dim2?dim1:dim2)];

	for(int cnt2 = 0;cnt2<dim2;cnt2++)
	{
		for(int cnt1 = 0;cnt1<dim1;cnt1++)
		{
			matrix[cnt2*2*dim1+2*cnt1] = mat(cnt1,cnt2).real();
			matrix[cnt2*2*dim1+2*cnt1+1] = mat(cnt1,cnt2).imag();
		}
	}

	int val= svd(matrix, u_in,vt_in,sing_val,dim1,dim2);
/*	for (uint cnt = 0;cnt<2*dim1*2*dim2;cnt++)
	{
		std::cout<<vt_in[cnt]<<" "<<std::endl;
	}
std::cout<<std::endl;
*/
	for(int cnt2 = 0;cnt2<(dim1<dim2?dim1:dim2);cnt2++)
	{
		for(int cnt1 = 0;cnt1<dim1;cnt1++)
		{
			U(cnt1,cnt2)=Complex(u_in[cnt2*2*dim1+2*cnt1],u_in[cnt2*2*dim1+2*cnt1+1]);
		}
	}
	for(int cnt2 = 0;cnt2<dim2;cnt2++)
	{
		for(int cnt1 = 0;cnt1<(dim1<dim2?dim1:dim2);cnt1++)
		{
			V_dagger(cnt1,cnt2)=Complex(vt_in[cnt2*2*(dim1<dim2?dim1:dim2)+2*cnt1],vt_in[cnt2*2*(dim1<dim2?dim1:dim2)+2*cnt1+1]);

		}
	}

		for(int cnt2 = 0;cnt2<(dim1<dim2?dim1:dim2);cnt2++)
		{
			D(cnt2,cnt2)=sing_val[cnt2];
		}
		 delete[] u_in;
		 delete[] vt_in;
		 delete[] sing_val;
		 delete[] matrix;

		 FixPhase(U,V_dagger);
	  return val;
}

int BoostSVD(const CMatrixTypeColMaj & mat,CMatrixTypeColMaj &U,MatrixTypeColMaj &D,CMatrixTypeColMaj &V_dagger)
{
  int dim1 = mat.size1(),dim2 = mat.size2();
  U.resize(dim1,(dim2<dim1?dim2:dim1));
  boost_setTo(U,Complex(0.0));
  V_dagger.resize(dim1<dim2?dim1:dim2,dim2);
  boost_setTo(V_dagger,Complex(0.0));
  D.resize(dim1<dim2?dim1:dim2,dim1<dim2?dim1:dim2);
  boost_setTo(D,0.0);

  double *matrix = new double[2*dim1*dim2];

  //dimension of u_in has to be chosen according to zegsvd (see manual):dim(u_in) = (dim1,Min(dim1,dim2)));
  double *u_in = new double[dim1*2*(dim1<dim2?dim1:dim2)];

  //dimension of vt_in: dim(vt_in) = (Min(dim1,dim2),dim2)
  double *vt_in = new double[(dim1<dim2?dim1:dim2)*2*dim2];

  //dimension of sing_val: this is a vector of doubles (singular values): dim(sing_val) = min(dim1,dim2)
  double *sing_val = new double[(dim1<dim2?dim1:dim2)];


  for(int cnt = 0;cnt<dim1*dim2;cnt++)
    {
      matrix[2*cnt] = mat.data()[cnt].real();
      matrix[2*cnt+1] = mat.data()[cnt].imag();


    }
  int val= svd(matrix, u_in,vt_in,sing_val,dim1,dim2);
  /*	for (uint cnt = 0;cnt<2*dim1*2*dim2;cnt++)
	{
	std::cout<<vt_in[cnt]<<" "<<std::endl;
	}
	std::cout<<std::endl;
  */
  for(uint cnt = 0;cnt<U.size1()*U.size2();cnt++)
    {
      U.data()[cnt].real() =u_in[2*cnt];
      U.data()[cnt].imag() =u_in[2*cnt+1];


    }

  for(uint cnt = 0;cnt<V_dagger.size1()*V_dagger.size2();cnt++)
    {
      V_dagger.data()[cnt].real()=vt_in[2*cnt];
      V_dagger.data()[cnt].imag()=vt_in[2*cnt+1];


    }


  for(int cnt2 = 0;cnt2<(dim1<dim2?dim1:dim2);cnt2++)
    {
      D(cnt2,cnt2)=sing_val[cnt2];
    }
  delete[] u_in;
  delete[] vt_in;
  delete[] sing_val;
  delete[] matrix;
  FixPhase(U,V_dagger);

  return val;
}
int BoostSVD(const CMatrixTypeColMaj & mat,CMatrixTypeColMaj &U,VectorType &D,CMatrixTypeColMaj &V_dagger)
{
	int dim1 = mat.size1(),dim2 = mat.size2();
	U.resize(dim1,(dim2<dim1?dim2:dim1));
	boost_setTo(U,Complex(0.0));
	V_dagger.resize(dim1<dim2?dim1:dim2,dim2);
	boost_setTo(V_dagger,Complex(0.0));
	D.resize(dim1<dim2?dim1:dim2);
	boost_setTo(D,0.0);

	double *matrix = new double[2*dim1*dim2];

	//dimension of u_in has to be chosen according to zegsvd (see manual):dim(u_in) = (dim1,Min(dim1,dim2)));
	double *u_in = new double[dim1*2*(dim1<dim2?dim1:dim2)];

	//dimension of vt_in: dim(vt_in) = (Min(dim1,dim2),dim2)
	double *vt_in = new double[(dim1<dim2?dim1:dim2)*2*dim2];

	//dimension of sing_val: this is a vector of doubles (singular values): dim(sing_val) = min(dim1,dim2)
	double *sing_val = D.data().begin();


	for(int cnt = 0;cnt<dim1*dim2;cnt++)
	{
		matrix[2*cnt] = mat.data()[cnt].real();
		matrix[2*cnt+1] = mat.data()[cnt].imag();


	}
	int val= svd(matrix, u_in,vt_in,sing_val,dim1,dim2);
/*	for (uint cnt = 0;cnt<2*dim1*2*dim2;cnt++)
	{
		std::cout<<vt_in[cnt]<<" "<<std::endl;
	}
std::cout<<std::endl;
*/
	for(uint cnt = 0;cnt<U.size1()*U.size2();cnt++)
	{
			U.data()[cnt].real() =u_in[2*cnt];
			U.data()[cnt].imag() =u_in[2*cnt+1];


	}

		for(uint cnt = 0;cnt<V_dagger.size1()*V_dagger.size2();cnt++)
		{
			V_dagger.data()[cnt].real()=vt_in[2*cnt];
			V_dagger.data()[cnt].imag()=vt_in[2*cnt+1];


		}

		 delete[] u_in;
		 delete[] vt_in;
		 delete[] matrix;

		 FixPhase(U,V_dagger);
	  return val;
}



// Diagonalisation of Hermitian, dense matrices using the LAPACK-routine
// ZHEEVD (divide-and-conquer-algorithm)

// SUBROUTINE ZHEEVD( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK,
//                 LRWORK, IWORK, LIWORK, INFO )

/*  JOBZ    (input) CHARACTER*1
*          = 'N':  Compute eigenvalues only;
*          = 'V':  Compute eigenvalues and eigenvectors.
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  Upper triangle of A is stored;
*          = 'L':  Lower triangle of A is stored.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input/output) COMPLEX*16 array, dimension (LDA, N)
*          On entry, the Hermitian matrix A.  If UPLO = 'U', the
*          leading N-by-N upper triangular part of A contains the
*          upper triangular part of the matrix A.  If UPLO = 'L',
*          the leading N-by-N lower triangular part of A contains
*          the lower triangular part of the matrix A.
*          On exit, if JOBZ = 'V', then if INFO = 0, A contains the
*          orthonormal eigenvectors of the matrix A.
*          If JOBZ = 'N', then on exit the lower triangle (if UPLO='L')
*          or the upper triangle (if UPLO='U') of A, including the
*          diagonal, is destroyed.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  W       (output) DOUBLE PRECISION array, dimension (N)
*          If INFO = 0, the eigenvalues in ascending order.
*
*  WORK    (workspace/output) COMPLEX*16 array, dimension (MAX(1,LWORK))
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The length of the array WORK.
*          If N <= 1,                LWORK must be at least 1.
*          If JOBZ  = 'N' and N > 1, LWORK must be at least N + 1.
*          If JOBZ  = 'V' and N > 1, LWORK must be at least 2*N + N**2.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal sizes of the WORK, RWORK and
*          IWORK arrays, returns these values as the first entries of
*          the WORK, RWORK and IWORK arrays, and no error message
*          related to LWORK or LRWORK or LIWORK is issued by XERBLA.
*
*  RWORK   (workspace/output) DOUBLE PRECISION array,
*                                         dimension (LRWORK)
*          On exit, if INFO = 0, RWORK(1) returns the optimal LRWORK.
*
*  LRWORK  (input) INTEGER
*          The dimension of the array RWORK.
*          If N <= 1,                LRWORK must be at least 1.
*          If JOBZ  = 'N' and N > 1, LRWORK must be at least N.
*          If JOBZ  = 'V' and N > 1, LRWORK must be at least
*                         1 + 5*N + 2*N**2.
*
*          If LRWORK = -1, then a workspace query is assumed; the
*          routine only calculates the optimal sizes of the WORK, RWORK
*          and IWORK arrays, returns these values as the first entries
*          of the WORK, RWORK and IWORK arrays, and no error message
*          related to LWORK or LRWORK or LIWORK is issued by XERBLA.
*
*  IWORK   (workspace/output) INTEGER array, dimension (MAX(1,LIWORK))
*          On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK.
*
*  LIWORK  (input) INTEGER
*          The dimension of the array IWORK.
*          If N <= 1,                LIWORK must be at least 1.
*          If JOBZ  = 'N' and N > 1, LIWORK must be at least 1.
*          If JOBZ  = 'V' and N > 1, LIWORK must be at least 3 + 5*N.
*
*          If LIWORK = -1, then a workspace query is assumed; the
*          routine only calculates the optimal sizes of the WORK, RWORK
*          and IWORK arrays, returns these values as the first entries
*          of the WORK, RWORK and IWORK arrays, and no error message
*          related to LWORK or LRWORK or LIWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i and JOBZ = 'N', then the algorithm failed
*                to converge; i off-diagonal elements of an intermediate
*                tridiagonal form did not converge to zero;
*                if INFO = i and JOBZ = 'V', then the algorithm failed
*                to compute an eigenvalue while working on the submatrix
*                lying in rows and columns INFO/(N+1) through
*                mod(INFO,N+1). */




extern "C"
{
  void zheevd_ (char *jobz, char *uplo, int *n, double *a, int *lda, double *w,
           double *work, int *lwork, double *rwork, int *lrwork, int *iwork, int *liwork, int *info);

}

int diag_hermitian(double *matrix, double *eigen, int dim)
{
  char jobz = 'V';
  char uplo = 'U';
  // n = dim, a = matrix, w = eigen
  int lda = dim;
  int info;

  // workspace variables
  int lwork = 2 * dim + dim * dim;
  int lrwork = 1 + 5 * dim + 2 * dim * dim;
  int liwork = 3 + 5 * dim;

  double *work = new double[lwork * 2];
  double *rwork = new double[lrwork];
  int *iwork = new int[liwork];

  zheevd_ (&jobz, &uplo, &dim, matrix, &lda, eigen,
             work, &lwork, rwork, &lrwork, iwork, &liwork, &info);

  delete[] iwork;
  delete[] rwork;
  delete[] work;

  return info;
}




void BoostDiag(CMatrixTypeColMaj &M,VectorType &E)
{

	double* matrix = new double[2*M.size1()*M.size2()];
	double* eigen = new double[2*M.size1()*M.size2()];


	for(uint ind = 0;ind<M.size1()*M.size2();ind++)
	{
		matrix[2*ind] = M.data()[ind].real();
		matrix[2*ind+1] = M.data()[ind].imag();
	}


	diag_hermitian(matrix, eigen,M.size1());
	E.resize(M.size1());
	for(uint ind = 0;ind<M.size1()*M.size2();ind++)
	{
		M.data()[ind].real()=matrix[2*ind] ;
		M.data()[ind].imag()=matrix[2*ind+1];
	}
	for(uint ind = 0;ind<M.size1();ind++)
	{
		E.data()[ind] = eigen[ind];
	}
	delete [] eigen;
	delete [] matrix;
}
void BoostDiag(MatrixTypeColMaj &M,VectorType &E)
{

	double* matrix = new double[2*M.size1()*M.size2()];
	double* eigen = new double[M.size1()];


	for(uint ind = 0;ind<M.size1()*M.size2();ind++)
	{
		matrix[2*ind] = M.data()[ind];
		matrix[2*ind+1] =0.0;
	}
	E.resize(M.size1());
	diag_hermitian(matrix, eigen,M.size1());
	for(uint ind = 0;ind<M.size1()*M.size2();ind++)
	{
		M.data()[ind] = matrix[2*ind] ;
	}
	for(uint ind = 0;ind<M.size1();ind++)
	{
		E.data()[ind] = eigen[ind];
	}
	delete [] eigen;
	delete [] matrix;

}
















