#ifndef __TENSOR_IZTOK_H__
#define __TENSOR_IZTOK_H__
/****************************************************************************
 * tensor.h                                                                 *
 * DESCRIPTION:                                                             *
 *  implementation of tensor classes.                                       *
 * DEPENDENCIES:                                                            *
 *  std::algorithm                                                          *
 *                                                                          *
 *                                                                          *
 * AUTHOR:                                                                  *
 *  iztok pizorn, universitaet wien, 2009                                   *
 *                                                                          *
 ****************************************************************************/

#include <cassert>
#include <cstring>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <set>
#include "valarray.hpp"
#include "keyTypes.hpp"
//#define VECTORCLASS cVector
//#define VECTORCLASS ivalarray

template<typename T>
T arrayproduct(const unsigned int n, const T* a)
{
    T res=1;
    for(int i=0;i<n;++i) res *= a[i];
    return res;
}


namespace std
{
  inline double conj(double x) { return x; }
  inline double real(double x) { return x; }
  inline double imag(double x) { return 0.; }
  inline double norm(double x) { return x*x; }
}

template<size_t N> int arrprod(const int* A) { return A[0] * arrprod<N-1>(A+1); }
template<> int arrprod<1>(const int* A) { return *A; }


template<size_t N> unsigned int arrprod(const unsigned int* A) { return A[0] * arrprod<N-1>(A+1); }
template<> unsigned int arrprod<1>(const unsigned int* A) { return *A; }

template<size_t N> unsigned long arrprod(const unsigned long* A) { return A[0] * arrprod<N-1>(A+1); }
template<> unsigned long arrprod<1>(const unsigned long* A) { return *A; }



template<typename T>
void _reccopy(T* dst, T* src, const unsigned int n, const unsigned int* newm, const unsigned int* oldm);

template<size_t Nl, typename T> class tensor: public ivalarray<T>
{
  unsigned int _mm[Nl];
public:
        // constructor
        tensor();
        tensor(const int tsize);
        tensor(const unsigned int* newm);
        tensor(const ukey<Nl>& newm);
        // link-constructor
        tensor(const int tsize, T* space);
        tensor(const unsigned int* newm, T* src) __attribute__((deprecated));
        
        // copy constructor
        tensor(const tensor<Nl,T>& );
        
        // read-only
        void getM(int* m) const __attribute__((deprecated)) { for(int i=0;i<Nl;++i) m[i]=_mm[i]; }
        void getM(unsigned int* m) const __attribute__((deprecated)) { memcpy(m, _mm, sizeof(_mm)); }
        const unsigned int* M() const { return _mm; }
        unsigned int M(const int which) const { return _mm[which]; }
        const ukey<Nl> shape() const { return ukey<Nl>(_mm); }
        unsigned int effsize() const { return arrprod<Nl>(_mm); }

        void resize(const unsigned int newsize, const bool copycontents=true) { ivalarray<T>::resize(newsize,copycontents); }
        void resize(const unsigned int* newM, const bool copycontents=true);
        void resize(const ukey<Nl>& newM, const bool copycontents=true);

        void reshape_in_place(const ukey<Nl>& newM);
        //void reshape(const ukey<Nl>& newM) __attribute__((deprecated));
        
        template<int Ml> tensor<Ml,T> reshape(const ukey<Ml>& newM) const;
        template<int Ml> tensor<Ml,T> releg(const ukey<Ml>& order) const;
        tensor<Nl,T> permute(const ukey<Nl>& P) const;
        
        // ACTION
        //void setM(const int* newm) { memcpy(_mm, newm, sizeof(int)*Nl); }
        void setM(const unsigned int* newm) { ivalarray<T>::resize( arrprod<Nl>(newm), true); memcpy(_mm, newm, sizeof(unsigned int)*Nl); }
        void setM(const int which, const int val) { _mm[which] = val; ivalarray<T>::resize( arrprod<Nl>(_mm), true); }
        void truncate(const unsigned int* newm);
        void truncate(const ukey<Nl>& shape) { return truncate(shape.raw() ); }
        void one() { this->resize(1, false); for(int i=0;i<Nl;++i) _mm[i]=1; (*this)[0]=1.0; }
        void copy(tensor<Nl,T>* , const unsigned int* mnew = NULL);
        tensor<Nl,T> copy() const { tensor<Nl,T> X(*this); return X; }
        double renorm();
        double norm() const;
        double maxval() const;
        
        // access
        T at(int i0,int i1,int i2,int i3,int i4,int i5) const {assert( i0<_mm[0] && i1 < _mm[1] && i2<_mm[2] && i3<_mm[3] && i4<_mm[4]&&i5<_mm[5]);  return (*this)[ i0 + _mm[0]*(i1 + _mm[1]*(i2+_mm[2]*(i3+_mm[3]*(i4+_mm[4]*i5))))]; }
        T& at(int i0,int i1,int i2,int i3,int i4,int i5) {assert( i0<_mm[0] && i1 < _mm[1] && i2<_mm[2] && i3<_mm[3] && i4<_mm[4]&&i5<_mm[5]);  return (*this)[ i0 + _mm[0]*(i1 + _mm[1]*(i2+_mm[2]*(i3+_mm[3]*(i4+_mm[4]*i5))))]; }
        //
        T at(int i0,int i1,int i2,int i3,int i4) const {assert( i0<_mm[0] && i1 < _mm[1] && i2<_mm[2] && i3<_mm[3] && i4<_mm[4]);  return (*this)[ i0 + _mm[0]*(i1 + _mm[1]*(i2+_mm[2]*(i3+_mm[3]*i4)))]; }
        T& at(int i0,int i1,int i2,int i3,int i4) {assert( i0<_mm[0] && i1 < _mm[1] && i2<_mm[2] && i3<_mm[3] && i4<_mm[4]);  return (*this)[ i0 + _mm[0]*(i1 + _mm[1]*(i2+_mm[2]*(i3+_mm[3]*i4)))]; }
        T at(int i0, int i1, int i2, int i3) const {assert( i0<_mm[0] && i1 < _mm[1] && i2<_mm[2] && i3<_mm[3]); return (*this)[ i0 + _mm[0]*(i1 + _mm[1]*(i2+_mm[2]*i3))]; }
        T& at(int i0, int i1, int i2, int i3) {assert( i0<_mm[0] && i1 < _mm[1] && i2<_mm[2] && i3<_mm[3]); return (*this)[ i0 + _mm[0]*(i1 + _mm[1]*(i2+_mm[2]*i3))]; }
        T at(int i0, int i1, int i2) const {assert( i0<_mm[0] && i1 < _mm[1] && i2<_mm[2]);  return (*this)[ i0 + _mm[0]*(i1 + _mm[1]*i2)]; }
        T& at(int i0, int i1, int i2) {assert( i0<_mm[0] && i1 < _mm[1] && i2<_mm[2]);  return (*this)[ i0 + _mm[0]*(i1 + _mm[1]*i2)]; }
        T at(int i0, int i1) const {assert( i0<_mm[0] && i1 < _mm[1]);  return (*this)[ i0 + _mm[0]*i1]; }
        T& at(int i0, int i1) {assert( i0<_mm[0] && i1 < _mm[1]);  return (*this)[ i0 + _mm[0]*i1 ]; }
        const T& at(int i0) const { return ivalarray<T>::operator[](i0); }
        T& at(int i0) { return ivalarray<T>::operator[](i0); }
        T at(const unsigned int* aridx) const { return operator[](arr2idx(aridx)); }
        T& at(const unsigned int* aridx) { return ivalarray<T>::operator[](arr2idx(aridx)); }
        

        T& operator()(const unsigned int i0,const unsigned int i1,const unsigned int i2,const unsigned int i3,const unsigned int i4, const unsigned int i5) {assert( i0<_mm[0] && i1 < _mm[1] && i2<_mm[2] && i3<_mm[3] && i4<_mm[4] && i5<_mm[5]);  return (*this)[ i0 + _mm[0]*(i1 + _mm[1]*(i2+_mm[2]*(i3+_mm[3]*(i4+_mm[4]*i5))))]; }
        const T& operator()(const unsigned int i0,const unsigned int i1,const unsigned int i2,const unsigned int i3,const unsigned int i4, const unsigned int i5) const {assert( i0<_mm[0] && i1 < _mm[1] && i2<_mm[2] && i3<_mm[3] && i4<_mm[4] && i5<_mm[5]); return (*this)[ i0 + _mm[0]*(i1 + _mm[1]*(i2+_mm[2]*(i3+_mm[3]*(i4+_mm[4]*i5))))]; }
        
        T& operator()(const unsigned int i0,const unsigned int i1,const unsigned int i2,const unsigned int i3,const unsigned int i4) {assert( i0<_mm[0] && i1 < _mm[1] && i2<_mm[2] && i3<_mm[3] && i4<_mm[4]);  return (*this)[ i0 + _mm[0]*(i1 + _mm[1]*(i2+_mm[2]*(i3+_mm[3]*i4)))]; }
        const T& operator()(const unsigned int i0,const unsigned int i1,const unsigned int i2,const unsigned int i3,const unsigned int i4) const {assert( i0<_mm[0] && i1 < _mm[1] && i2<_mm[2] && i3<_mm[3] && i4<_mm[4]); return (*this)[ i0 + _mm[0]*(i1 + _mm[1]*(i2+_mm[2]*(i3+_mm[3]*i4)))]; }
        //
        const T& operator()(int i0, int i1, int i2, int i3) const {assert( i0<_mm[0] && i1 < _mm[1] && i2<_mm[2] && i3<_mm[3]); return (*this)[ i0 + _mm[0]*(i1 + _mm[1]*(i2+_mm[2]*i3))]; }
        T& operator()(int i0, int i1, int i2, int i3) {assert( i0<_mm[0] && i1 < _mm[1] && i2<_mm[2] && i3<_mm[3]);  return (*this)[ i0 + _mm[0]*(i1 + _mm[1]*(i2+_mm[2]*i3))]; }
        //
        const T& operator()(int i0, int i1, int i2) const {assert( i0<_mm[0] && i1 < _mm[1] && i2<_mm[2]);  return (*this)[ i0 + _mm[0]*(i1 + _mm[1]*i2)]; }
        T& operator()(int i0, int i1, int i2) {assert( i0<_mm[0] && i1 < _mm[1] && i2<_mm[2]);  return (*this)[ i0 + _mm[0]*(i1 + _mm[1]*i2)]; }
        //
        const T& operator()(int i0, int i1) const { assert( i0<_mm[0] && i1 < _mm[1]); return (*this)[ i0 + _mm[0]*i1]; }
        T& operator()(int i0, int i1) { assert( i0<_mm[0] && i1 < _mm[1]); return (*this)[ i0 + _mm[0]*i1 ]; }
        //
        const T& operator()(const unsigned int* aridx) const { return at(arr2idx(aridx)); }
        T& operator()(const unsigned int* aridx) { return at(arr2idx(aridx)); }

        
        unsigned int arr2idx(const unsigned int* idx) const;
        void idx2arr(const unsigned int idx, unsigned int* arr) const;


        tensor<Nl,T>& operator+=(const tensor<Nl,T>& );
        tensor<Nl,T>& operator=(const tensor<Nl,T>& );
        tensor<Nl,T>& operator-=(const tensor<Nl,T>& );
        

        void tobuffer(std::ofstream& buf) const;
        void write(std::ofstream& buf) const { tobuffer(buf); }
        void frombuffer(std::ifstream& buf);
        void read(std::ifstream& buf) { frombuffer(buf); }
    //protected:
        //unsigned int& M(const int which) { return _mm[which]; }
        
        tensor<Nl,T>& operator*=(const T&);
        
};



template<size_t Nl, typename T>
tensor<Nl,T>& tensor<Nl,T>::operator*=(const T& x)
{
    const unsigned int n = effsize();
    scal(n, x, this->raw(), 1);
    return (*this);    
}


template<size_t Nl, typename T>
tensor<Nl,T> operator*(const T& x, const tensor<Nl,T>& A)
{
    tensor<Nl,T> B(A);
    B *= x;
    return B;
}


template<size_t Nl, typename T>
tensor<Nl,T> operator*(const tensor<Nl,T>& A, const T& x)
{
    tensor<Nl,T> B(A);
    B *= x;
    return B;
}



// constructors
template<size_t Nl, typename T>
tensor<Nl,T>::tensor(): ivalarray<T>()
{
    memset(_mm,0,sizeof(_mm));
}


template<size_t Nl, typename T>
tensor<Nl,T>::tensor(const int tsize): ivalarray<T>(tsize)
{
    memset(_mm,0,sizeof(_mm));
}

template<size_t Nl, typename T>
tensor<Nl,T>::tensor(const unsigned int* newm): ivalarray<T>(arrprod<Nl>(newm))
{
    setM(newm);
}

template<size_t Nl, typename T>
tensor<Nl,T>::tensor(const ukey<Nl>& newm): ivalarray<T>(arrprod<Nl>(newm.raw()))
{
    this->setM(newm.raw());
}


template<size_t Nl, typename T>
tensor<Nl,T>::tensor(const int tsize, T* space): ivalarray<T>(tsize, space)
{
    memset(_mm,0,sizeof(_mm));
}


template<size_t Nl, typename T>
tensor<Nl,T>::tensor(const unsigned int* newm, T* src): ivalarray<T>(arrprod<Nl>(newm), src)
{
    this->setM(newm);
}

template<size_t Nl,typename T>
tensor<Nl,T>::tensor(const tensor<Nl,T>& tnsr)
    : ivalarray<T>( (const ivalarray<T>&)tnsr)
{
    for(int i=0;i<Nl;++i) _mm[i]=tnsr._mm[i];
}


template<size_t Nl,typename T>
tensor<Nl,T>& tensor<Nl,T>::operator=(const tensor<Nl,T>& tnsr)
{
    ivalarray<T>::operator=( (const ivalarray<T>&) tnsr );
    for(int i=0;i<Nl;++i) _mm[i]=tnsr._mm[i];
    return (*this);
}



template<size_t Nl,typename T>
tensor<Nl,T>& tensor<Nl,T>::operator+=(const tensor<Nl,T>& tnsr)
{
    assert( tnsr.shape() == this->shape() );
    //for(int j=0;j<Nl;++j) assert( tnsr.M(j)==this->M(j));
    ivalarray<T>::operator+=( (const ivalarray<T>&) tnsr );
    return (*this);
}



template<size_t Nl,typename T>
tensor<Nl,T>& tensor<Nl,T>::operator-=(const tensor<Nl,T>& tnsr)
{
    for(int j=0;j<Nl;++j) assert( tnsr.M(j)==this->M(j));
    ivalarray<T>::operator-=( (const ivalarray<T>&) tnsr );
    return (*this);
}



template<size_t Nl,typename T>
tensor<Nl,T> operator-(const tensor<Nl,T>& A, const tensor<Nl,T>& B)
{
    assert( A.shape() == B.shape() );
    tensor<Nl,T> C(A);
    C -= B;
    return C;
}


template<size_t Nl, typename T>
unsigned int tensor<Nl,T>::arr2idx(const unsigned int* idx) const
{
    unsigned int ari=0, fac=1;
    #pragma unroll(Nl)
    for(unsigned int i=0;i!=Nl;++i) {
        ari += idx[i]*fac;
        fac *= _mm[i];
    }
    return ari;
}

// this function is pretty expensive!!
template<size_t Nl, typename T>
void tensor<Nl,T>::idx2arr(const unsigned int idx, unsigned int* arr) const
{
    unsigned int ii = idx;
    unsigned int i=0;
    
    #pragma unroll(Nl)
    for(unsigned int i=0;i!=Nl;++i) {
        arr[i] = ii % _mm[i];
        unsigned int j  =  ii / _mm[i];
        ii = j;
    }

}


template<size_t Nl, typename T>
void tensor<Nl,T>::reshape_in_place(const ukey<Nl>& newm)
{
    if (arrprod<Nl>(_mm) != arrprod<Nl>(newm.raw()) ) {
        std::cerr << "shape mismatch.\n";
        throw "gbfafds";
    }
    memcpy(_mm, newm.raw(), sizeof(_mm) );
}


 
template<size_t Nl, typename T>
void tensor<Nl,T>::truncate(const unsigned int* newm)
{
    for(int leg=Nl-1;leg>=0;--leg) {
        if (newm[leg]==_mm[leg]) continue;
        
        assert( newm[leg] < _mm[leg] );
        
        const unsigned int lda = arrayproduct<unsigned int>(leg, _mm);
        const unsigned int grt = arrayproduct<unsigned int>(Nl-leg-1, newm+leg+1);

        for(int j=1;j<grt;++j) {
            memmove( this->raw()+j*lda*newm[leg], this->raw()+j*lda*_mm[leg], lda*newm[leg]*sizeof(T) );
        }
    }
    for(int i=0;i<Nl;++i) _mm[i]=newm[i];
}


template<size_t Nl,typename T>
void tensor<Nl,T>::copy(tensor<Nl,T>* tnsr, const unsigned int* newm)
{ 
    if (!newm) {
        memcpy(_mm, tnsr->_mm,sizeof(_mm)); 
        memcppy(this->raw(), tnsr->raw(), sizeof(T)*std::min(this->size(), tnsr->size() ));
        return;
    }
    memcpy(_mm, newm, sizeof(_mm) );
    _reccopy(this->raw(), tnsr->raw(), Nl, newm, tnsr->_mm);
}




template<size_t Nl,typename T>
void tensor<Nl,T>::resize(const ukey<Nl>& newm, const bool copycontents)
{
    resize(newm.raw(), copycontents);
}

template<size_t Nl,typename T>
void tensor<Nl,T>::resize(const unsigned int* newm, const bool copycontents)
{
    ivalarray<T>::resize(arrprod<Nl>(newm), copycontents);
    setM(newm);
}

template<size_t Nl,typename T>
double tensor<Nl,T>::maxval() const
{
    double maxv=0;
    const int nsize = this->effsize();
    for(int i=0;i<nsize;++i)
        if (std::abs( (*this)[i] ) > maxv ) 
            maxv = std::abs( (*this)[i] );
    return maxv;
}

template<size_t Nl,typename T>
double tensor<Nl,T>::renorm()
{
    assert( this->size() == this->effsize() );
    const double maxv = this->norm();
    //double maxv = this->maxval();
    ivalarray<T>::operator*=(1.0/maxv);
    //this->multiply(1.0/maxv );
    return maxv;
}


template<size_t Nl,typename T>
double tensor<Nl,T>::norm() const
{
    double nrr = 0.0;
    const unsigned int nsize = this->effsize() * (sizeof(T)/sizeof(double) );
    const double* x = (const double*)this->raw();
    for(unsigned int i=0;i<nsize;++i)
        nrr += x[i]*x[i];
    return sqrt(nrr);
}


template<size_t Nl,typename T>
void tensor<Nl,T>::tobuffer(std::ofstream& buf) const
{
    unsigned int nl = Nl;
    buf.write((char*)&nl, sizeof(unsigned int));
    buf.write((char*)_mm, sizeof(_mm));
    buf.write((char*)this->raw(), sizeof(T)*this->effsize());
}

template<size_t Nl,typename T>
void tensor<Nl,T>::frombuffer(std::ifstream& buf)
{
    unsigned int nl;
    buf.read((char*)&nl, sizeof(unsigned int));
    assert( nl == Nl );

    unsigned int mm[Nl];
    buf.read((char*)mm, sizeof(mm));
    
    this->resize(mm, false);
    buf.read((char*)this->raw(), arrprod<Nl>(mm)*sizeof(T));
}





template<typename T> void makehermitian(tensor<2,T>* A)
{
    int n = A->M(0);
    assert( A->M(1) == A->M(0) );
     
    for(int i=0;i<n;++i) {
        A->at(i,i) = std::real(A->at(i,i));
        for(int j=i+1;j<n;++j) {
            A->at(i,j) = 0.5*(A->at(i,j) + std::conj(A->at(j,i)));
            A->at(j,i) = std::conj(A->at(i,j));
        }
    }

}


template<typename T>
void _reccopy(T* dst, T* src, const unsigned int n, const unsigned int* newm, const unsigned int* oldm)
{
    if (n==1) {
        memcpy(dst, src, (*newm)*sizeof(T) );
        return;
    }
    int rnew=1, rold=1;
    for(int i=0;i<n-1;++i) {
        rnew *= newm[i];
        rold *= oldm[i];
    }
    for(int i=0;i<newm[n-1];++i) {
        _reccopy(dst + i*rnew, src + i*rold, n-1, newm, oldm);
    }
}


template<size_t Nl, typename T>
std::ostream& operator<<(std::ostream& stream, const tensor<Nl,T>& ego)
{
    unsigned int arr[Nl];
    for(unsigned int i=0;i<ego.effsize();++i) {
        T value = ego[i];
        if (value == (T)0) continue;
        index2array(Nl, i, arr, ego.M() );
        stream<<"[";
        for(int j=0;j<Nl;++j)
            stream << arr[j]<<(j<Nl-1 ? "|" : "]");
        stream << " --> " << ego[i] << std::endl;
    }
    return stream;
}

template<typename T>
std::ostream& operator<<(std::ostream& stream, const tensor<2,T>& ego)
{
    // it's just a matrix
    for(int i=0;i<ego.M(0);++i) {
        for(int j=0;j<ego.M(1);++j)
            stream << std::setw(8) << std::setprecision(3) << ego(i,j)<<" "<< std::setprecision(-1);

        stream << std::endl;
    }


    return stream;
}

#endif
