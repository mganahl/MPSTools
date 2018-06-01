#ifndef TYPEDEFS_HPP_
#define TYPEDEFS_HPP_

#define BEGIN_TIME(var) static double var=0; var-=clock();
#define END_TIME(var)  var+=clock(); cout<<#var<<": "<<var/CLOCKS_PER_SEC<<endl;


#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/vector_sparse.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/operation_sparse.hpp>
#include <boost/lexical_cast.hpp>
#include <complex>

namespace BoostTypedefs
{

typedef boost::numeric::ublas::range Range;

//basic numeric types
typedef double Real;
typedef std::complex<double> Complex;

//basic types
typedef boost::numeric::ublas::vector<Real> VectorType;
typedef boost::numeric::ublas::matrix<Real> MatrixType;
typedef boost::numeric::ublas::matrix<Real,boost::numeric::ublas::column_major> MatrixTypeColMaj;
typedef boost::numeric::ublas::matrix_range<boost::numeric::ublas::matrix<Real, boost::numeric::ublas::column_major> >MatrixRangeTypeColMaj;
typedef boost::numeric::ublas::compressed_vector<Real> SparseVectorType;
typedef boost::numeric::ublas::compressed_matrix<Real> SparseMatrixType;

//basic complex types
typedef boost::numeric::ublas::vector<Complex> CVectorType;
typedef boost::numeric::ublas::matrix<Complex> CMatrixType;
typedef boost::numeric::ublas::matrix<Complex,boost::numeric::ublas::column_major> CMatrixTypeColMaj;
typedef boost::numeric::ublas::matrix_range<boost::numeric::ublas::matrix<Complex, boost::numeric::ublas::column_major> >CMatrixRangeTypeColMaj;
typedef boost::numeric::ublas::compressed_vector<Complex> SparseCVectorType;
typedef boost::numeric::ublas::compressed_matrix<Complex> SparseCMatrixType;
typedef boost::numeric::ublas::matrix_range<CMatrixType> CMatrixTypeRange;
//integer types
typedef boost::numeric::ublas::vector<int> IntVectorType;
typedef boost::numeric::ublas::matrix<int> IntMatrixType;
typedef boost::numeric::ublas::vector<uint> UintVectorType;
typedef boost::numeric::ublas::vector<uint> CoordVectorType;

//matrix proxies
typedef boost::numeric::ublas::matrix_column<MatrixType> MatrixColumnType;
typedef boost::numeric::ublas::matrix_row<MatrixType> MatrixRowType;
typedef boost::numeric::ublas::matrix_range<MatrixType> MatrixRangeType;
typedef boost::numeric::ublas::matrix_range<SparseMatrixType> SparseMatrixRangeType;

//complex matrix proxies
typedef boost::numeric::ublas::matrix_column<CMatrixType> CMatrixColumnType;
typedef boost::numeric::ublas::matrix_row<CMatrixType> CMatrixRowType;
typedef boost::numeric::ublas::matrix_range<CMatrixType> CMatrixRangeType;

//const matrix proxies
typedef boost::numeric::ublas::matrix_column<const MatrixType> ConstMatrixColumnType;
typedef boost::numeric::ublas::matrix_range<const MatrixType> ConstMatrixRangeType;
typedef boost::numeric::ublas::matrix_range<const SparseMatrixType> ConstSparseMatrixRangeType;

//vector proxies
typedef boost::numeric::ublas::vector_range<VectorType> VectorRangeType;
typedef boost::numeric::ublas::vector_range<CoordVectorType> CoordVectorRangeType;

//complex vector proxies
typedef boost::numeric::ublas::vector_range<CVectorType> CVectorRangeType;

//const vector proxies
typedef boost::numeric::ublas::vector_range<const VectorType> ConstVectorRangeType;
typedef boost::numeric::ublas::vector_range<const CVectorType> ConstCVectorRangeType;

//special matrices/vectors
typedef boost::numeric::ublas::identity_matrix<Real> IdentityMatrixType;
typedef boost::numeric::ublas::zero_matrix<Real> ZeroMatrixType;
typedef boost::numeric::ublas::zero_vector<Real> ZeroVectorType;
typedef boost::numeric::ublas::scalar_vector<Real> ScalarVectorType;

//complex special matrices/vectors
typedef boost::numeric::ublas::identity_matrix<Complex> IdentityCMatrixType;
typedef boost::numeric::ublas::zero_matrix<Complex> ZeroCMatrixType;
typedef boost::numeric::ublas::zero_vector<Complex> ZeroCVectorType;
typedef boost::numeric::ublas::scalar_vector<Complex> ScalarCVectorType;

//range
typedef boost::numeric::ublas::range RangeType;

}
#endif /*TYPEDEFS_HPP_*/
