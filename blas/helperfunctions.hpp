#ifndef helperfunctions_h_
#define helperfunctions_h_

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <iostream>
#include <string>
#include "../boost_typedefs.hpp"
//#include  "BoostToolsLapack/LapackSVD.h"
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include "boost/multi_array.hpp"
#include<complex>
using namespace BoostTypedefs;
using namespace std;






template<class T>
uint VecAbsMaxInd(boost::numeric::ublas::vector<T>& vec)
{
    uint ind=0;
    Real val=0;
    for (uint i=0; i<vec.size(); i++)
    {
        if(abs(vec(i))>val)
        {
            val=abs(vec(i));
            ind=i;
        }
    }
    return ind;
}
template<typename StorageOrder>
void FixPhase(boost::numeric::ublas::matrix<Real,StorageOrder>& U, boost::numeric::ublas::matrix<Real,StorageOrder>& VT)
{
  using namespace boost::numeric::ublas;
  Real phi;
  boost::numeric::ublas::vector<Real> vec;
  uint dim=U.size2();
  assert(dim==VT.size1());
  for (uint i=0; i<dim; i++)
    {
      vec=VectorType(row(VT,i));
      phi=(vec(VecAbsMaxInd(vec)))>0?1:-1;
      row(VT,i) *=phi;
      column(U,i) *= phi;
    }
}
template<typename StorageOrder>
void FixPhase(boost::numeric::ublas::matrix<Complex,StorageOrder>& U, boost::numeric::ublas::matrix<Complex,StorageOrder>& VT)
{
	using namespace boost::numeric::ublas;
    Real phi;
    boost::numeric::ublas::vector<Complex> vec;
    uint dim=U.size2();
    assert(dim==VT.size1());
    for (uint i=0; i<dim; i++)
    {
        vec=CVectorType(row(VT,i));
        phi=arg(vec(VecAbsMaxInd(vec)));
        row(VT,i) *=Complex(exp(-Complex(0,1)*phi));
        column(U,i) *=Complex(exp(Complex(0,1)*phi));
    }
}


template<typename T>
boost::numeric::ublas::vector<T> cut(boost::numeric::ublas::vector<T> vec,T delta)
{
	typename boost::numeric::ublas::vector<T>::iterator it=vec.begin();
	boost::numeric::ublas::vector<T> out;
	uint l = 0;//
	//first test if end is reached, then if value is larger than delta; other wise program
	//crashes
	while(it!=vec.end()&&fabs(*it)>delta)
	{
		l +=1;
		it++;
	}
	out.resize(l);
	for(uint ind = 0;ind<l;ind++)
	{
		out(ind) = vec(ind);
	}
	return out;
}
//reduces size of matrix to newsize1,newsize2
template<typename T>
void cut(boost::numeric::ublas::matrix<T,boost::numeric::ublas::column_major> &mat,boost::numeric::ublas::matrix<T,boost::numeric::ublas::column_major> &matout, uint newsize1,uint newsize2)
{
	if (mat.size1()<newsize1){cout<<"in cut: newsize1<=size1() failed"<<endl;abort();}
	if (mat.size2()<newsize2){cout<<"in cut: newsize2<=size2() failed"<<endl;abort();}
	matout.resize(newsize1,newsize2);

	//boost::numeric::ublas::matrix_range<boost::numeric::ublas::matrix<T,boost::numeric::ublas::column_major> > Mat_range(Matrix,nummapl.find(it->first.first)->second,nummapr.find(it->first.third)->second);
	matout=boost::numeric::ublas::matrix_range<boost::numeric::ublas::matrix<T,boost::numeric::ublas::column_major> >(mat,Range(0,newsize1),Range(0,newsize2));
}



bool isequal(UintVectorType,UintVectorType);

void my_pause();
std::string itoa(int num,int base,int precision);
/*
void boost_setTo(CMatrixType &mat, Complex value);
void boost_setTo(MatrixType &mat, Real value);
void boost_setTo(IntMatrixType &mat, int value);
void boost_setTo(MatrixType &mat, int value);
void boost_setTo(VectorType &vec, Real value);
void boost_setTo(UintVectorType &vec,uint value);
void boost_setTo(CVectorType &vec, Complex value);
*/

void generate_next_occupation(std::vector<uint> &occ,std::vector<uint> &max_occ);
uint compute_single_index(std::vector<uint> occ,std::vector<uint> max_occ);
void compute_multi_index(uint x,std::vector<uint> max_occ,std::vector<uint> &occ);

void printboostCmatrix(const CMatrixType&ar);
void printboostCmatrix(const CMatrixTypeColMaj&ar);
void print2dboostCarray(const boost::multi_array<Complex,2> &ar);

Real Anorm(const boost::multi_array<Complex,3>&ar);
Real Anorm(const boost::multi_array<Complex,2>&ar);

//please ensure that both multi_arrays have the same storage order!!
Complex AScalarProd(const boost::multi_array<Complex,3>&ar1,boost::multi_array<Complex,3>&ar2);
Complex AScalarProd(const boost::multi_array<Complex,2>&ar1,boost::multi_array<Complex,2>&ar2);

bool VectorCompareLEQ(UintVectorType vec1,UintVectorType vec2);
bool VectorCompareLESS( UintVectorType vec1, UintVectorType vec2);
bool VectorCompareLEQ( VectorType vec1, VectorType vec2);

Real myFabs(const Real num);
Real myFabs(const Complex num);


template<class T>
bool isequal(const boost::numeric::ublas::matrix<T,boost::numeric::ublas::column_major> mat1,const boost::numeric::ublas::matrix<T,boost::numeric::ublas::column_major> mat2)
{
	assert(mat1.size1()==mat2.size1());
	assert(mat1.size2()==mat2.size2());

	bool comp=true;
	uint cnt = 0;
	while(cnt!=mat1.size1()*mat1.size2())
	{
		if(mat1.data()[cnt]!=mat2.data()[cnt]){comp = false;break;}
		else{cnt++;}
	}
	return comp;
}
template<class T>
bool isequal(const boost::numeric::ublas::matrix<T> mat1,const boost::numeric::ublas::matrix<T> mat2)
{
	assert(mat1.size1()==mat2.size1());
	assert(mat1.size2()==mat2.size2());

	bool comp=true;
	uint cnt = 0;
	while(cnt!=mat1.size1()*mat1.size2())
	{
		if(mat1.data()[cnt]!=mat2.data()[cnt]){comp = false;break;}
		else{cnt++;}
	}
	return comp;
}

template <typename T>
void boost_setTo(boost::numeric::ublas::matrix<T> &mat, T value)
{
     uint size1 = mat.size1(),size2 = mat.size2();
     for(uint x = 0;x<size1;x++) for(uint y=0;y<size2;y++)mat(x,y)=value;
}

template <typename T>
void boost_setTo(boost::numeric::ublas::matrix<T,boost::numeric::ublas::column_major> &mat, T value)
{
     uint size1 = mat.size1(),size2 = mat.size2();
     for(uint x = 0;x<size1;x++) for(uint y=0;y<size2;y++)mat(x,y)=value;
}
template <typename T>
void boost_setTo(boost::numeric::ublas::vector<T> &vec, T value)
{
	uint size = vec.size();
	for(uint y=0;y<size;y++){vec(y)=value; }
}

template<typename T>
void printboostmatrix(const boost::numeric::ublas::matrix<T> &mat)
{
  cout.precision(4);
	uint size1 = mat.size1(),size2 = mat.size2();
	cout<<"Ar(:,:) = ["<<endl;
	for(uint x = 0;x<size1;x++){for(uint y=0;y<size2;y++){std::cout<<mat(x,y)<<" ";}std::cout<<std::endl;}
	cout<<"]"<<endl;
}
template<typename T>
void printboostmatrix(const boost::numeric::ublas::matrix<T,boost::numeric::ublas::column_major,boost::numeric::ublas::unbounded_array<T> > &mat)
{
 cout.precision(4);
	uint size1 = mat.size1(),size2 = mat.size2();
	for(uint x = 0;x<size1;x++){for(uint y=0;y<size2;y++){std::cout<<mat(x,y)<<" ";}std::cout<<std::endl;}
}
template<typename T>
void print4dboostarray(const boost::multi_array<T,4> &ar)
{using std::endl;
using std::cout;
	uint size1 = ar.shape()[0],size2 = ar.shape()[1],size3 =ar.shape()[2],size4 =ar.shape()[3];;
	for(uint w=0;w<size1;w++){for(uint z=0;z<size2;z++){cout<<"Ar("<<w+1<<","<<z+1<<",:,:)=["<<endl<<endl;for(uint x = 0;x<size3;x++){for(uint y=0;y<size4;y++){std::cout<<ar[w][z][x][y].real()<<"+1i*"<<ar[w][z][x][y].imag()<<" ";}cout<<endl;}cout<<"]"<<endl;}}

}

template<typename T>
void print3dboostarray(const boost::multi_array<T,3> &ar)
{using std::endl;
using std::cout;
	uint size1 = ar.shape()[0],size2 = ar.shape()[1],size3 =ar.shape()[2];
	for(uint z=0;z<size3;z++){cout<<"Ar(:,:,"<<z+1<<")=["<<endl<<endl;for(uint x = 0;x<size1;x++){for(uint y=0;y<size2;y++){std::cout<<ar[x][y][z].real()<<"+1i*"<<ar[x][y][z].imag()<<" ";}cout<<endl;}cout<<"]"<<endl;}

}

template<typename T>
void print2dboostarray(const boost::multi_array<T,2> &ar)
{using std::endl;
using std::cout;
	uint size1 = ar.shape()[0],size2 = ar.shape()[1];
	cout<<"Ar(:,:)="<<endl<<endl;for(uint x = 0;x<size1;x++){for(uint y=0;y<size2;y++){std::cout<<ar[x][y].real()<<"+1i*"<<ar[x][y].imag()<<" ";}cout<<endl;}cout<<endl;
}





#endif

