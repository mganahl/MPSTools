/*
This file contains all hamiltonians to be used int the dmrg or any other engine. the hamiltonians
are all in Matrix product operator (MPO) representation.

 */
#ifndef HAMILTONIANS_HPP__
#define HAMILTONIANS_HPP__

#include "MPOTypes.hpp"
#include "keyTypes.hpp"
#include <sstream>
#include <algorithm>
#include "utilities.hpp"


template<typename T>
class IdentityMPO:public MPO<T>
{
public:
  template<typename KeyType>
  IdentityMPO(const MPS<KeyType,T>&mps){
    MPOMat<1,T> M;
    for(uint site = 0;site<mps.size();++site){
      SparseLocalOperator<T> id(mps.getPSpace(site),mps.getPSpace(site));
      for(uint n=0;n<mps.getPSpace(site);++n)
	id[OpKeyType<1>(n,n)]=T(1.0);

      MPOMat<1,T> M;
      M.clear();
      M[ukey<2>(0,0)] = id;
      M.Dl=1,M.Dr=1;
      M.squeeze(1e-15);
      (*this)[site]=M;
    }
  }
  ~IdentityMPO(){};
};


template<typename T>
class HeisenbergMPO:public MPO<T>
{
public:
  HeisenbergMPO(VectorType Jz,VectorType Jxy,VectorType B);
  ~HeisenbergMPO(){};
  void shiftandrescale(const T E0,const T a, const T W);
};



template<typename T>
void HeisenbergMPO<T>::shiftandrescale(const T E0,const T a, const T W){
  SparseLocalOperator<T> id(2,2);
  const int N=this->size();
  const T E=E0*1.0/N;
  id[OpKeyType<1>(0,0)]=T(1.);id[OpKeyType<1>(1,1)]=T(1.);
  id*=E;
  id.setSite(0);
  (*this)[0][ukey<2>(0,0)]-=id ;
  for(uint site=1;site<N;++site){
    const uint lend=(*this)[site].Dl;
    id.setSite(site);
    (*this)[site][ukey<2>(lend-1,0)]-=id ;
  }

  typename MPOMat<1,T>::iterator o;
  for(o=(*this)[0].begin();o!=(*this)[0].end();++o)
    o->second/=a;

  if(W!=0.0){
    id[OpKeyType<1>(0,0)]=T(1.);id[OpKeyType<1>(1,1)]=T(1.);
    const T _w=W*1.0/N;
    id*=_w;
    id.setSite(0);
    (*this)[0][ukey<2>(0,0)]-=id;
      for(uint site=1;site<N;++site){
	const uint lend=(*this)[site].Dl;
	id.setSite(site);
	(*this)[site][ukey<2>(lend-1,0)]-=id ;
      }
  }
  //that's it ...
}

template<typename T>
HeisenbergMPO<T>::  HeisenbergMPO(VectorType Jz,VectorType Jxy,VectorType B){
  uint N = B.size();
  MPOMat<1,T> M;
  //  this->resize(N);
  SparseLocalOperator<T> Sz(2,2),Sp(2,2),Sm(2,2),id(2,2);//all statistics are set to bosonic (_statistics=1)
  Sz[OpKeyType<1>(0,0)]=T(-0.5);Sz[OpKeyType<1>(1,1)]=T(0.5);
  Sp[OpKeyType<1>(1,0)]=T(1);Sm[OpKeyType<1>(0,1)]=T(1);//Sp.setSign(-1);
  id[OpKeyType<1>(0,0)]=T(1);id[OpKeyType<1>(1,1)]=T(1);//Sm.setSign(-1);

  Sz.setSite(0);Sp.setSite(0);Sm.setSite(0);id.setSite(0);
  M[ukey<2>(0,0)] = -B(0)*Sz;
  M[ukey<2>(0,1)] = Jxy(0)/2.*Sm;
  M[ukey<2>(0,2)] = Jxy(0)/2.*Sp;
  M[ukey<2>(0,3)] = Jz(0)*Sz;
  M[ukey<2>(0,4)] = id;
  M.Dl=1,M.Dr=5;
  M.squeeze(1e-16);
  (*this)[0]=M;

  M.clear();
  Sz.setSite(N-1);Sp.setSite(N-1);Sm.setSite(N-1);id.setSite(N-1);
  M[ukey<2>(0,0)] = id;
  M[ukey<2>(1,0)] = Sp;
  M[ukey<2>(2,0)] = Sm;
  M[ukey<2>(3,0)] = Sz;
  M[ukey<2>(4,0)] = -B(N-1)*Sz;
  M.Dl=5,M.Dr=1;
  M.squeeze(1e-16);
  (*this)[N-1]=M;

  for(uint site = 1;site<(N-1);++site){
    MPOMat<1,T> Mb;
    Sz.setSite(site);Sp.setSite(site);Sm.setSite(site);id.setSite(site);
    Mb[ukey<2>(0,0)] = id;
    Mb[ukey<2>(1,0)] = Sp;
    Mb[ukey<2>(2,0)] = Sm;
    Mb[ukey<2>(3,0)] = Sz;
    Mb[ukey<2>(4,0)] = -B(site)*Sz;

    Mb[ukey<2>(4,1)] = Jxy(site)/2.*Sm;
    Mb[ukey<2>(4,2)] = Jxy(site)/2.*Sp;
    Mb[ukey<2>(4,3)] = Jz(site)*Sz;
    Mb[ukey<2>(4,4)] = id;
    Mb.Dl=5,Mb.Dr=5;
    Mb.squeeze(1e-16);
    (*this)[site]=Mb;
  }
}

//template<typename T>
//class HeisenbergImpMPO:public MPO<T>
//{
//public:
//  HeisenbergImpMPO(VectorType Jz,VectorType Jxy,VectorType B,uint impsite,Real Jimp,Real Bimp);
//  ~HeisenbergImpMPO(){};
//  void shiftandrescale(const T E0,const T a, const T W){};
//};
//
//
//template<typename T>
//HeisenbergImpMPO<T>::  HeisenbergImpMPO(VectorType Jz,VectorType Jxy,VectorType B,Real Jimp,Real Bimp){
//  uint N = B.size();
//  MPOMat<1,T> M;
//  //  this->resize(N);
//  SparseLocalOperator<T> Sz(2,2),Sp(2,2),Sm(2,2),id(2,2),SzI(4,4),SpI(4,4),SmI(4,4),SzC(4,4),SpC(4,4),SmC(4,4),ID(4,4);
//  Sz[OpKeyType<1>(0,0)]=T(-0.5);Sz[OpKeyType<1>(1,1)]=T(0.5);
//  Sp[OpKeyType<1>(1,0)]=T(1);Sm[OpKeyType<1>(0,1)]=T(1);//Sp.setSign(-1);
//  id[OpKeyType<1>(0,0)]=T(1);id[OpKeyType<1>(1,1)]=T(1);//Sm.setSign(-1);
//    
//  SzI=kronecker(Sz,id);
//  SpI=kronecker(Sp,id);
//  SmI=kronecker(Sm,id);
//
//  SzC=kronecker(id,Sz);
//  SpC=kronecker(id,Sp);
//  SmC=kronecker(id,Sm);
//
//  ID_=kronecker(id,id);
//
//  Sz.setSite(0);Sp.setSite(0);Sm.setSite(0);id.setSite(0),SzC.setSite(0),SpC.setSite(0),SmC.setSite(0),SzI.setSite(0),SpI.setSite(0),SmI.setSite(0);
//  if(0!=impsite){
//    M[ukey<2>(0,0)] = -B(0)*Sz;
//    M[ukey<2>(0,1)] = Jxy(0)/2.*Sm;
//    M[ukey<2>(0,2)] = Jxy(0)/2.*Sp;
//    M[ukey<2>(0,3)] = Jz(0)*Sz;
//    M[ukey<2>(0,4)] = id;
//    M.Dl=1,M.Dr=5;
//    M.squeeze(1e-16);
//    (*this)[0]=M;
//  }else if(1==impsite){
//    M[ukey<2>(0,0)] = -B(0)*SzC-Bimp*SzI+(Jimp/2.0)*(SpC*SmI+SpI*SmC);
//    M[ukey<2>(0,1)] = Jxy(0)/2.*SmC;
//    M[ukey<2>(0,2)] = Jxy(0)/2.*SpC;
//    M[ukey<2>(0,3)] = Jz(0)*SzC;
//    M[ukey<2>(0,4)] = ID_;
//    M.Dl=1,M.Dr=5;
//    M.squeeze(1e-16);
//    (*this)[0]=M;
//  }
//
//  Sz.setSite(N-1);Sp.setSite(N-1);Sm.setSite(N-1);id.setSite(N-1),SzC.setSite(N-1),SpC.setSite(N-1),SmC.setSite(N-1),SzI.setSite(N-1),SpI.setSite(N-1),SmI.setSite(N-1);
//  if(N-1!=impsite){
//    M.clear();
//    M[ukey<2>(0,0)] = id;
//    M[ukey<2>(1,0)] = Sp;
//    M[ukey<2>(2,0)] = Sm;
//    M[ukey<2>(3,0)] = Sz;
//    M[ukey<2>(4,0)] = -B(N-1)*Sz;
//    M.Dl=5,M.Dr=1;
//    M.squeeze(1e-16);
//    (*this)[N-1]=M;
//  }else if(N-1==impsite){
//    M.clear();
//    M[ukey<2>(0,0)] = id;
//    M[ukey<2>(1,0)] = Sp;
//    M[ukey<2>(2,0)] = Sm;
//    M[ukey<2>(3,0)] = Sz;
//    M[ukey<2>(4,0)] = -B(N-1)*Sz;
//    M.Dl=5,M.Dr=1;
//    M.squeeze(1e-16);
//    (*this)[N-1]=M;
//
//  }
//
//  for(uint site = 1;site<(N-1);++site){
//    MPOMat<1,T> Mb;
//    Sz.setSite(site);Sp.setSite(site);Sm.setSite(site);id.setSite(site);
//    Mb[ukey<2>(0,0)] = id;
//    Mb[ukey<2>(1,0)] = Sp;
//    Mb[ukey<2>(2,0)] = Sm;
//    Mb[ukey<2>(3,0)] = Sz;
//    Mb[ukey<2>(4,0)] = -B(site)*Sz;
//
//    Mb[ukey<2>(4,1)] = Jxy(site)/2.*Sm;
//    Mb[ukey<2>(4,2)] = Jxy(site)/2.*Sp;
//    Mb[ukey<2>(4,3)] = Jz(site)*Sz;
//    Mb[ukey<2>(4,4)] = id;
//    Mb.Dl=5,Mb.Dr=5;
//    Mb.squeeze(1e-16);
//    (*this)[site]=Mb;
//  }
//}
//

template<typename T>
class IntegrableNNNXXZMPO:public MPO<T>
{
public:
  IntegrableNNNXXZMPO(VectorType Jz,VectorType Jxy,VectorType B,Real alpha);
  ~IntegrableNNNXXZMPO(){};
  void shiftandrescale(const T E0,const T a, const T W);
};



template<typename T>
void IntegrableNNNXXZMPO<T>::shiftandrescale(const T E0,const T a, const T W){
  SparseLocalOperator<T> id(2,2);
  const int N=this->size();
  const T E=E0*1.0/N;
  //  cout<<"that's E: "<<E<<endl;
  id[OpKeyType<1>(0,0)]=T(1.);id[OpKeyType<1>(1,1)]=T(1.);
  id*=E;
  id.setSite(0);
  (*this)[0][ukey<2>(0,0)]-=id ;
  for(uint site=1;site<N;++site){
    const uint lend=(*this)[site].Dl;
    id.setSite(site);
    (*this)[site][ukey<2>(lend-1,0)]-=id ;
  }

  typename MPOMat<1,T>::iterator o;
  for(o=(*this)[0].begin();o!=(*this)[0].end();++o)
    o->second/=a;

  if(W!=0.0){
    id[OpKeyType<1>(0,0)]=T(1.);id[OpKeyType<1>(1,1)]=T(1.);
    const T _w=W*1.0/N;
    id*=_w;
    id.setSite(0);
    (*this)[0][ukey<2>(0,0)]-=id;
      for(uint site=1;site<N;++site){
	const uint lend=(*this)[site].Dl;
	id.setSite(site);
	(*this)[site][ukey<2>(lend-1,0)]-=id ;
      }
  }
  //that's it ...
}


template<typename T>
IntegrableNNNXXZMPO<T>::IntegrableNNNXXZMPO(VectorType Jz,VectorType Jxy,VectorType B,Real alpha){
  uint N = B.size();
  MPOMat<1,T> M;
  //  this->resize(N);
  SparseLocalOperator<T> Sz(2,2),Sp(2,2),Sm(2,2),id(2,2);//all statistics are set to bosonic (_statistics=1)
  Sz[OpKeyType<1>(0,0)]=T(-0.5);Sz[OpKeyType<1>(1,1)]=T(0.5);
  Sp[OpKeyType<1>(1,0)]=T(1);Sm[OpKeyType<1>(0,1)]=T(1);//Sp.setSign(-1);
  id[OpKeyType<1>(0,0)]=T(1);id[OpKeyType<1>(1,1)]=T(1);//Sm.setSign(-1);

  Sz.setSite(0);Sp.setSite(0);Sm.setSite(0);id.setSite(0);
  M[ukey<2>(0,0)] = -B(0)*Sz;
  M[ukey<2>(0,1)] = Jxy(0)/2.*Sm;
  M[ukey<2>(0,2)] = Jxy(0)/2.*Sp;
  M[ukey<2>(0,3)] = Jz(0)*Sz;
  M[ukey<2>(0,4)] = Complex(0,alpha*Jz(0)*Jxy(0)/2.0)*Sm;
  M[ukey<2>(0,5)] = Complex(0,-alpha*Jz(0)*Jxy(0)/2.0)*Sp;
  M[ukey<2>(0,6)] = Complex(0,-alpha*Jz(0)*Jxy(0)/2.0)*Sz;
  M[ukey<2>(0,7)] = Complex(0,alpha*Jz(0)*Jxy(0)/2.0)*Sz;
  M[ukey<2>(0,8)] = Complex(0,-alpha*Jxy(0)*Jxy(0)/2.0)*Sm;
  M[ukey<2>(0,9)] = Complex(0,alpha*Jxy(0)*Jxy(0)/2.0)*Sp;
  M[ukey<2>(0,10)] = id;
  M.Dl=1,M.Dr=11;
  M.squeeze(1e-16);
  (*this)[0]=M;

  M.clear();
  Sz.setSite(N-1);Sp.setSite(N-1);Sm.setSite(N-1);id.setSite(N-1);
  M[ukey<2>(0,0)] = id;
  M[ukey<2>(1,0)] = Sp;
  M[ukey<2>(2,0)] = Sm;
  M[ukey<2>(3,0)] = Sz;
  M[ukey<2>(10,0)] = -B(N-1)*Sz;
  M.Dl=11,M.Dr=1;
  M.squeeze(1e-16);
  (*this)[N-1]=M;

  for(uint site = 1;site<(N-1);++site){
    MPOMat<1,T> Mb;
    Sz.setSite(site);Sp.setSite(site);Sm.setSite(site);id.setSite(site);
    Mb[ukey<2>(0,0)] = id;
    Mb[ukey<2>(1,0)] = Sp;
    Mb[ukey<2>(2,0)] = Sm;
    Mb[ukey<2>(3,0)] = Sz;
    Mb[ukey<2>(4,3)] = Sp;
    Mb[ukey<2>(5,3)] = Sm;
    Mb[ukey<2>(6,2)] = Sp;
    Mb[ukey<2>(7,1)] = Sm;
    Mb[ukey<2>(8,1)] = Sz;
    Mb[ukey<2>(9,2)] = Sz;
    Mb[ukey<2>(10,0)] = -B(site)*Sz;


    Mb[ukey<2>(10,1)] = Jxy(site)/2.*Sm;
    Mb[ukey<2>(10,2)] = Jxy(site)/2.*Sp;
    Mb[ukey<2>(10,3)] = Jz(site)*Sz;
    Mb[ukey<2>(10,4)] = Complex(0,alpha*Jz(site)*Jxy(site)/2.0)*Sm;
    Mb[ukey<2>(10,5)] = Complex(0,-alpha*Jz(site)*Jxy(site)/2.0)*Sp;
    Mb[ukey<2>(10,6)] = Complex(0,-alpha*Jz(site)*Jxy(site)/2.0)*Sz;
    Mb[ukey<2>(10,7)] = Complex(0,alpha*Jz(site)*Jxy(site)/2.0)*Sz;
    Mb[ukey<2>(10,8)] = Complex(0,-alpha*Jxy(site)*Jxy(site)/2.0)*Sm;
    Mb[ukey<2>(10,9)] = Complex(0,alpha*Jxy(site)*Jxy(site)/2.0)*Sp;
    Mb[ukey<2>(10,10)] = id;


    Mb.Dl=11,Mb.Dr=11;
    Mb.squeeze(1e-16);
    (*this)[site]=Mb;
  }
}


template<typename T>
class HeisenbergNNNIntMPO:public MPO<T>{
public:
  HeisenbergNNNIntMPO(VectorType Jz,VectorType Jz2,VectorType Jxy,VectorType Jxy2,VectorType B);
  ~HeisenbergNNNIntMPO(){};
  void shiftandrescale(const T E0,const T a, const T W);
};


template<typename T>
void HeisenbergNNNIntMPO<T>::shiftandrescale(const T E0,const T a, const T W){
  SparseLocalOperator<T> id(2,2);
  id[OpKeyType<1>(0,0)]=T(1.);id[OpKeyType<1>(1,1)]=T(1.);
  id*=E0;
  cout<<(*this)[0]<<endl;
  (*this)[0][ukey<2>(0,0)]-=id;
  typename MPOMat<1,T>::iterator o;
  for(o=(*this)[0].begin();o!=(*this)[0].end();++o)
    o->second/=a;

  id[OpKeyType<1>(0,0)]=T(1.);id[OpKeyType<1>(1,1)]=T(1.);
  id*=W;
  (*this)[0][ukey<2>(0,0)]-=id;
  cout<<(*this)[0]<<endl;
  //that's it ...
}

template<typename T>
HeisenbergNNNIntMPO<T>::HeisenbergNNNIntMPO(VectorType Jz,VectorType Jz2,VectorType Jxy,VectorType Jxy2,VectorType B) {
  uint N = B.size();
  MPOMat<1,T> M;
  SparseLocalOperator<T> Sz(2,2),Sp(2,2),Sm(2,2),id(2,2);
  Sz[OpKeyType<1>(0,0)]=T(-0.5);Sz[OpKeyType<1>(1,1)]=T(0.5);
  Sp[OpKeyType<1>(1,0)]=T(1);Sm[OpKeyType<1>(0,1)]=T(1);
  id[OpKeyType<1>(0,0)]=T(1);id[OpKeyType<1>(1,1)]=T(1);

  Sz.setSite(0);Sp.setSite(0);Sm.setSite(0);id.setSite(0);
  M[ukey<2>(0,0)] = -B(0)*Sz;
  M[ukey<2>(0,1)] = Jxy(0)/2.*Sm;
  M[ukey<2>(0,2)] = Jxy(0)/2.*Sp;
  M[ukey<2>(0,3)] = Jz(0)*Sz;
  M[ukey<2>(0,4)] = Jxy2(0)/2.*Sm;
  M[ukey<2>(0,5)] = Jxy2(0)/2.*Sp;
  M[ukey<2>(0,6)] = Jz2(0)*Sz;
  M[ukey<2>(0,7)] = id;
  M.Dl=1,M.Dr=8;
  M.squeeze(1e-16);
  (*this)[0]=M;

  M.clear();
  Sz.setSite(N-1);Sp.setSite(N-1);Sm.setSite(N-1);id.setSite(N-1);
  M[ukey<2>(0,0)] = id;
  M[ukey<2>(1,0)] = Sp;
  M[ukey<2>(2,0)] = Sm;
  M[ukey<2>(3,0)] = Sz;
  M[ukey<2>(7,0)] = -B(N-1)*Sz;
  M.Dl=8,M.Dr=1;
  M.squeeze(1e-16);
  (*this)[N-1]=M;
  for(uint site = 1;site<(N-1);++site){
    MPOMat<1,T> Mb;
    Sz.setSite(site);Sp.setSite(site);Sm.setSite(site);id.setSite(site);
    Mb[ukey<2>(0,0)] = id;
    Mb[ukey<2>(1,0)] = Sp;
    Mb[ukey<2>(2,0)] = Sm;
    Mb[ukey<2>(3,0)] = Sz;
    Mb[ukey<2>(4,1)] = id;
    Mb[ukey<2>(5,2)] = id;
    Mb[ukey<2>(6,3)] = id;

    Mb[ukey<2>(7,0)] = -B(site)*Sz;
    Mb[ukey<2>(7,1)] = Jxy(site)/2.*Sm;
    Mb[ukey<2>(7,2)] = Jxy(site)/2.*Sp;
    Mb[ukey<2>(7,3)] = Jz(site)*Sz;
    Mb[ukey<2>(7,4)] = Jxy2(site)/2.*Sm;
    Mb[ukey<2>(7,5)] = Jxy2(site)/2.*Sp;
    Mb[ukey<2>(7,6)] = Jz2(site)*Sz;
    Mb[ukey<2>(7,7)] = id;

    Mb.Dl=8,Mb.Dr=8;
    Mb.squeeze(1e-16);
    (*this)[site]=Mb;
  }
}



//H=sum_i (-t)(cd_i c_(i+1) +h.c.) +Vn_in_(i+1)
template<typename T>
class SpinlessFermionsMPO:public MPO<T>
{
public:
  SpinlessFermionsMPO(){};
  SpinlessFermionsMPO(VectorType V,VectorType t,VectorType mu);
  ~SpinlessFermionsMPO(){};
  void shiftandrescale(const Complex E0,const Complex a, const Complex W);
  void shiftandrescale(const Real E0,const Real a, const Real W);

  //virtual void construct(const int step){};
};

template<typename T>
void SpinlessFermionsMPO<T>::shiftandrescale(const Real E0,const Real a, const Real W){
  SparseLocalOperator<Real> id(2,2);
  const int N=this->size();
  const Real E=E0*1.0/N;
  //cout<<"that's E: "<<E<<endl;
  id[OpKeyType<1>(0,0)]=E;id[OpKeyType<1>(1,1)]=E;
  id.setSite(0);
  (*this)[0][ukey<2>(0,0)]-=id;
  for(uint site=1;site<N;++site){
    const uint lend=(*this)[site].Dl;
    id.setSite(site);
    (*this)[site][ukey<2>(lend-1,0)]-=id;
  }
  if(abs(a)>1e-14){
    typename MPOMat<1,Real>::iterator o;
    for(o=(*this)[0].begin();o!=(*this)[0].end();++o)
      o->second/=a;
  }

  if(abs(W)>1e-14){
    id[OpKeyType<1>(0,0)]=T(1.);id[OpKeyType<1>(1,1)]=T(1.);
    const Real _w=W*1.0/N;
    id*=_w;
    id.setSite(0);
    (*this)[0][ukey<2>(0,0)]-=id;
      for(uint site=1;site<N;++site){
	const uint lend=(*this)[site].Dl;
	id.setSite(site);
	(*this)[site][ukey<2>(lend-1,0)]-=id ;
      }
  }
  //that's it ...
}
template<typename T>
void SpinlessFermionsMPO<T>::shiftandrescale(const Complex E0,const Complex a, const Complex W){
  SparseLocalOperator<Complex> id(2,2);
  const int N=this->size();
  const Complex E=Complex(E0.real()*1.0/N,E0.imag()*1.0/N);

  id[OpKeyType<1>(0,0)]=E;id[OpKeyType<1>(1,1)]=E;
  id.setSite(0);

  (*this)[0][ukey<2>(0,0)]-=id;
  for(uint site=1;site<N;++site){
    const uint lend=(*this)[site].Dl;
    id.setSite(site);
    (*this)[site][ukey<2>(lend-1,0)]-=id;
  }
  if(abs(a)>1e-14){
    typename MPOMat<1,Complex>::iterator o;
    for(o=(*this)[0].begin();o!=(*this)[0].end();++o)
      o->second/=a;
  }

  if(abs(W)>1e-14){
    id[OpKeyType<1>(0,0)]=Complex(1.);id[OpKeyType<1>(1,1)]=Complex(1.);
    const Complex _w=T(W.real()*1.0/N,W.imag()*1.0/N);
    id*=_w;
    id.setSite(0);
    (*this)[0][ukey<2>(0,0)]-=id;
      for(uint site=1;site<N;++site){
	const uint lend=(*this)[site].Dl;
	id.setSite(site);
	(*this)[site][ukey<2>(lend-1,0)]-=id ;
      }
  }
}


template<typename T>
SpinlessFermionsMPO<T>:: SpinlessFermionsMPO(VectorType V,VectorType t,VectorType mu){
  uint N = mu.size();
  MPOMat<1,T> M;
  //  this->resize(N);
  SparseLocalOperator<T> n(2,2),cd(2,2),c(2,2),id(2,2),p(2,2);
  n[OpKeyType<1>(0,0)]=T(0);n[OpKeyType<1>(1,1)]=T(1);
  cd[OpKeyType<1>(1,0)]=T(1.);c[OpKeyType<1>(0,1)]=T(1.);
  id[OpKeyType<1>(0,0)]=T(1.);id[OpKeyType<1>(1,1)]=T(1.);
  p[OpKeyType<1>(0,0)]=T(1.0);p[OpKeyType<1>(1,1)]=T(-1.0);
  n.setSite(0);cd.setSite(0);c.setSite(0);id.setSite(0);p.setSite(0);
  M[ukey<2>(0,0)] = -mu(0)*n;
  M[ukey<2>(0,1)] = t(0)*c*p;
  M[ukey<2>(0,2)] = -t(0)*cd*p;
  M[ukey<2>(0,3)] = V(0)*n;
  M[ukey<2>(0,4)] = id;
  M.Dl=1,M.Dr=5;
  M.squeeze(1e-16);
  (*this)[0]=M;

  M.clear();
  n.setSite(N-1);cd.setSite(N-1);c.setSite(N-1);id.setSite(N-1);p.setSite(N-1);
  M[ukey<2>(0,0)] = id;
  M[ukey<2>(1,0)] = cd;
  M[ukey<2>(2,0)] = c;
  M[ukey<2>(3,0)] = n;
  M[ukey<2>(4,0)] = -mu(N-1)*n;
  M.Dl=5,M.Dr=1;
  M.squeeze(1e-16);
  (*this)[N-1]=M;

  for(uint site = 1;site<(N-1);++site){
    MPOMat<1,T> Mb;
    n.setSite(site);cd.setSite(site);c.setSite(site);id.setSite(site);p.setSite(site);
    Mb[ukey<2>(0,0)] = id;
    Mb[ukey<2>(1,0)] = cd;
    Mb[ukey<2>(2,0)] = c;
    Mb[ukey<2>(3,0)] = n;
    Mb[ukey<2>(4,0)] = -mu(site)*n;

    Mb[ukey<2>(4,1)] = t(site)*(c*p);
    Mb[ukey<2>(4,2)] = -t(site)*(cd*p);
    Mb[ukey<2>(4,3)] = V(site)*n;
    Mb[ukey<2>(4,4)] = id;
    Mb.Dl=5,Mb.Dr=5;
    Mb.squeeze(1e-16);
    (*this)[site]=Mb;
  }
}



template<typename T>
class DMFTHamiltonianMPO:public MPO<T>
{
public:
  DMFTHamiltonianMPO(){};
  DMFTHamiltonianMPO(VectorType t,VectorType mu,uint Norb);
  ~DMFTHamiltonianMPO(){};

};


template<typename T>
DMFTHamiltonianMPO<T>:: DMFTHamiltonianMPO(VectorType t,VectorType mu,uint Norb){
  uint N = mu.size();
  MPOMat<1,T> M;
  //  this->resize(N);
  SparseLocalOperator<T> n(2,2),cd(2,2),c(2,2),id(2,2),p(2,2);
  n[OpKeyType<1>(0,0)]=T(0);n[OpKeyType<1>(1,1)]=T(1);
  cd[OpKeyType<1>(1,0)]=T(1.);c[OpKeyType<1>(0,1)]=T(1.);
  id[OpKeyType<1>(0,0)]=T(1.);id[OpKeyType<1>(1,1)]=T(1.);
  p[OpKeyType<1>(0,0)]=T(1.0);p[OpKeyType<1>(1,1)]=T(-1.0);
  n.setSite(0);cd.setSite(0);c.setSite(0);id.setSite(0);p.setSite(0);


  M[ukey<2>(0,Norb)] = t(0)*(c*p);
  M[ukey<2>(0,2*Norb)] = -t(0)*(cd*p);
  M[ukey<2>(0,2*Norb+1)] = id;
  M.Dl=1;M.Dr=2*Norb+2;
  M.squeeze(1e-16);
  (*this)[0]=M;


  for(uint site = 1;site<(N-1);++site){
    MPOMat<1,T> Mb;
    n.setSite(site);cd.setSite(site);c.setSite(site);id.setSite(site);p.setSite(site);
    Mb[ukey<2>(0,0)] = id;
    Mb[ukey<2>(1,0)] = cd;
    Mb[ukey<2>(1+Norb,0)] = c;
    for(uint orb=0;orb<(Norb-1);++orb){
      Mb[ukey<2>(orb+2,orb+1)] = p;
      Mb[ukey<2>(2+Norb+orb,1+Norb+orb)] = p;
    }
    Mb[ukey<2>(2*Norb+1,Norb)] = t(site)*(c*p);
    Mb[ukey<2>(2*Norb+1,2*Norb)] = -t(site)*(cd*p);
    Mb[ukey<2>(2*Norb+1,2*Norb+1)] = id;

    Mb.Dl=2*Norb+2,Mb.Dr=2*Norb+2;
    Mb.squeeze(1e-14);
    (*this)[site]=Mb;

  }
  M.clear();
  n.setSite(N-1);cd.setSite(N-1);c.setSite(N-1);id.setSite(N-1);p.setSite(N-1);
  M[ukey<2>(0,0)] = id;
  M[ukey<2>(1,0)] = cd;
  M[ukey<2>(1+Norb,0)] = c;
  M.Dl=2*Norb+2,M.Dr=1;
  M.squeeze(1e-16);
  (*this)[N-1]=M;

}


//this class implements time dependent spinless fermions hamiltonian; N and Tmax are linear dimension of the system and the 
//maximum time step, hop V and chpot are files containing the time dependent parameters in matrix form, such that
//number of rows (x-dimension)=Tmax, number of columns (y-dimension)=N (or N-1 in the case of hop).
//to get the mpo at a certain step, call TDSpinlessFermionsMPO.construct(step). This operation materializes the MPO at time step 
//step: mpo.construct(step); now mpo can be used in the common way.
template<typename T>
class TDSpinlessFermionsMPO:public SpinlessFermionsMPO<T>
{
public:
  TDSpinlessFermionsMPO(){};
  TDSpinlessFermionsMPO(const uint N,const uint Tmax,const std::string hop,const std::string V,const std::string chpot);
  // TDSpinlessFermionsMPO(const uint N, const uint Tmax,const VectorType V,const VectorType t,const VectorType mu):SpinlessFermionsMPO<T>(V,t,mu){
  //   _N=N;
  //   _T=Tmax;
  // };
  TDSpinlessFermionsMPO(TDSpinlessFermionsMPO<T>&cpy):SpinlessFermionsMPO<T>(cpy){
    _hop=cpy.hop;
    _V=cpy.V();
    _mu=cpy.mu();
    _N=cpy.sizeN();
    _T=cpy.sizeT();
  
  }
  ~TDSpinlessFermionsMPO(){};
  //virtual T measure(const MPS<AbKeyType,T>&mps)const {return SpinlessFermionsMPO<T>::measure(mps);}
  //virtual T measure(const CanonizedMPS<AbKeyType,T>&mps)const {return SpinlessFermionsMPO<T>::measure(mps);}
  virtual void construct(const uint step);
  MatrixType hop()const {return _hop;}
  MatrixType V()const {return _V;}
  MatrixType sizeN()const {return _N;}
  MatrixType sizeT()const {return _T;}
  
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version){
    ar & boost::serialization::base_object<MPO<T> >(*this);
    ar &_hop;
    ar &_V;
    ar &_mu;
    ar &_N;
    ar &_T;
  }


private:
  MatrixType _hop,_V,_mu;
  uint _N,_T;
};




template <typename T>
void TDSpinlessFermionsMPO<T>::construct(const uint step){
  if(step>_T){
    cerr<<"TDSpinlessFermionsMPO<T>::construct: step>_T"<<endl;
    abort();
  }
  cout<<"her her "<<step<<endl;
  cout<<"_N "<<_N<<endl;
  VectorType h(_N-1),V(_N),mu(_N);

  for(uint n=0;n<_N;++n){
    if(n<_N-1)
      h(n)=_hop(step,n);
    V(n)=_V(step,n);
    mu(n)=_mu(step,n);
  }
  //(*this)=TDSpinlessFermionsMPO<T>(_N,_T,V,h,mu);
  this->clear();
  SpinlessFermionsMPO<T>mpo(V,h,mu);
  for(uint n=0;n<_N;++n){
    (*this)[n]=mpo[n];
  }
}

template <typename T>
TDSpinlessFermionsMPO<T>::TDSpinlessFermionsMPO(const uint N,const uint Tmax,const std::string hop,const std::string V,const std::string mu){
  _N=N,_T=Tmax;
  std::ifstream infile;

  infile.open(hop.c_str());
  cout<<"#TDSpinlessFermionsMPO<T>: reading hopping parameters from file "<<hop<<endl;
  _hop.resize(_T,_N-1);
  boost_setTo(_hop,0.0);
  for(uint t=0;t<_T;++t){
    for(uint n=0;n<_N-1;++n){
      infile>>_hop(t,n);
    }
  }
  infile.close();

  infile.open(V.c_str());
  cout<<"#TDSpinlessFermionsMPO<T>: reading V parameters from file "<<V<<endl;
  _V.resize(_T,_N);
  boost_setTo(_V,0.0);
  for(uint t=0;t<_T;++t){
    for(uint n=0;n<_N;++n){
      infile>>_V(t,n);
    }
  }
  infile.close();

  infile.open(mu.c_str());
  cout<<"#TDSpinlessFermionsMPO<T>: reading chemical potential parameters from file "<<mu<<endl;
  _mu.resize(_T,_N);
  boost_setTo(_mu,0.0);
  for(uint t=0;t<_T;++t){
    for(uint n=0;n<_N;++n){
      infile>>_mu(t,n);
    }
  }
  infile.close();

  cout<<"#done reading in parameters"<<endl;
}



template<typename T>
class FermiHubbardMPO:public MPO<T>
{
public:
  //hop: hopping, U: local interaction; chUp(chDown: chemical potential for up/down-spin electrons;
  FermiHubbardMPO(VectorType hop,VectorType U,VectorType chUp,VectorType chDown);
  //V: next-nearest neighbor interaction
  FermiHubbardMPO(VectorType hop,VectorType U,VectorType chUp,VectorType chDown,VectorType V);
  ~FermiHubbardMPO(){};
  virtual T measure(const MPS<AbKeyType,T>&mps)const {return MPO<T>::measure(mps);}
  virtual T measure(const CanonizedMPS<AbKeyType,T>&mps)const {return MPO<T>::measure(mps);}

  void shiftandrescale(const Real E0,const Real a, const Real W){
    SparseLocalOperator<Real> id(4,4);
    const int N=this->size();
    cout<<"that's N: "<<N<<endl;
    const Real E=E0*1.0/N;
    //cout<<"that's E: "<<E<<endl;
    id[OpKeyType<1>(0,0)]=E;id[OpKeyType<1>(1,1)]=E;id[OpKeyType<1>(2,2)]=E;id[OpKeyType<1>(3,3)]=E;
    id.setSite(0);
    (*this)[0][ukey<2>(0,0)]-=id;
    for(uint site=1;site<N;++site){
      const uint lend=(*this)[site].Dl;
      id.setSite(site);
      (*this)[site][ukey<2>(lend-1,0)]-=id;
    }
    if(abs(W)>1e-14){
      id[OpKeyType<1>(0,0)]=Real(1.);id[OpKeyType<1>(1,1)]=Real(1.);
      const Real _w=W*1.0/N;
      id*=_w;
      id.setSite(0);
      (*this)[0][ukey<2>(0,0)]-=id;
      for(uint site=1;site<N;++site){
	const uint lend=(*this)[site].Dl;
	id.setSite(site);
	(*this)[site][ukey<2>(lend-1,0)]-=id ;
      }
    }
  }
  void shiftandrescale(const Complex E0,const Complex a, const Complex W){
    SparseLocalOperator<Complex> id(4,4);
    const int N=this->size();
    cout<<"that's N: "<<N<<endl;
    const Complex E=Complex(E0.real()*1.0/N,E0.imag()*1.0/N);
    //cout<<"that's E: "<<E<<endl;
    id[OpKeyType<1>(0,0)]=Complex(E);id[OpKeyType<1>(1,1)]=Complex(E);id[OpKeyType<1>(2,2)]=Complex(E);id[OpKeyType<1>(3,3)]=Complex(E);
    id.setSite(0);
    (*this)[0][ukey<2>(0,0)]-=id;
    for(uint site=1;site<N;++site){
      const uint lend=(*this)[site].Dl;
      id.setSite(site);
      (*this)[site][ukey<2>(lend-1,0)]-=id;
    }
    if(abs(W)>1e-14){
      id[OpKeyType<1>(0,0)]=Complex(1.);id[OpKeyType<1>(1,1)]=Complex(1.);
      const Complex _w=Complex(W.real()*1.0/N,W.imag()*1.0/N);
      id*=_w;
      id.setSite(0);
      (*this)[0][ukey<2>(0,0)]-=id;
      for(uint site=1;site<N;++site){
	const uint lend=(*this)[site].Dl;
	id.setSite(site);
	(*this)[site][ukey<2>(lend-1,0)]-=id ;
      }
    }
  }


};
template<typename T>
FermiHubbardMPO<T>::FermiHubbardMPO(VectorType hop,VectorType U,VectorType chUp,VectorType chDown){
  uint N =chUp.size();
  assert(chUp.size()==chDown.size());
  MPOMat<1,T> M;
  //  this->resize(N);
  SparseLocalOperator<T> nu(4,4),nd(4,4),cuD(4,4),cdD(4,4),cu(4,4),cd(4,4),id(4,4),p(4,4);//all statistics are set to bosonic (_statistics=1)
  nu[OpKeyType<1>(2,2)]=T(1.0);nu[OpKeyType<1>(3,3)]=T(1.0);
  nd[OpKeyType<1>(1,1)]=T(1.0);nd[OpKeyType<1>(3,3)]=T(1.0);

  cuD[OpKeyType<1>(2,0)] = T(1.0);cuD[OpKeyType<1>(3,1)] = T(1.0); 
  cu[OpKeyType<1>(0,2)] = T(1.0);cu[OpKeyType<1>(1,3)] = T(1.0);   
								   
  cdD[OpKeyType<1>(1,0)] = T(1.0);cdD[OpKeyType<1>(3,2)] = T(-1.0);
  cd[OpKeyType<1>(0,1)] = T(1.0);cd[OpKeyType<1>(2,3)] = T(-1.0);  

  id[OpKeyType<1>(0,0)]=T(1);id[OpKeyType<1>(1,1)]=T(1);id[OpKeyType<1>(2,2)]=T(1);id[OpKeyType<1>(3,3)]=T(1);
  p[OpKeyType<1>(0,0)]=T(1);p[OpKeyType<1>(1,1)]=T(-1);p[OpKeyType<1>(2,2)]=T(-1);p[OpKeyType<1>(3,3)]=T(1);
  //p[OpKeyType<1>(0,0)]=T(1);p[OpKeyType<1>(1,1)]=T(1);p[OpKeyType<1>(2,2)]=T(1);p[OpKeyType<1>(3,3)]=T(1);
  

  nu.setSite(0);nd.setSite(0);cu.setSite(0);cd.setSite(0);cuD.setSite(0);cdD.setSite(0);id.setSite(0);p.setSite(0);
  //  M[ukey<2>(0,0)] = chUp(0)*nu+chDown(0)*nd+U(0)*((nu-0.5*id)*(nd-0.5*id));
  M[ukey<2>(0,0)] = chUp(0)*nu+chDown(0)*nd+U(0)*(nu*nd);
  M[ukey<2>(0,1)] = -hop(0)*(cuD*p);
  M[ukey<2>(0,2)] = -hop(0)*(cdD*p);
  M[ukey<2>(0,3)] = hop(0)*(cu*p);
  M[ukey<2>(0,4)] = hop(0)*(cd*p);
  M[ukey<2>(0,5)] = id;

  M.Dl=1;M.Dr=6;
  M.squeeze(1e-15);
  (*this)[0]=M;
  M.clear();

  nu.setSite(N-1);nd.setSite(N-1);cu.setSite(N-1);cd.setSite(N-1);cuD.setSite(N-1);cdD.setSite(N-1);id.setSite(N-1);p.setSite(N-1);
  M[ukey<2>(0,0)] = id;
  M[ukey<2>(1,0)] = cu;
  M[ukey<2>(2,0)] = cd;
  M[ukey<2>(3,0)] = cuD;
  M[ukey<2>(4,0)] = cdD;
  //  M[ukey<2>(5,0)] = chUp(N-1)*nu+chDown(N-1)*nd+U(N-1)*((nu-0.5*id)*(nd-0.5*id));
  M[ukey<2>(5,0)] = chUp(N-1)*nu+chDown(N-1)*nd+U(N-1)*(nu*nd);
  M.Dl=6;M.Dr=1;
  M.squeeze(1e-15);
  (*this)[N-1]=M;
  for(uint site = 1;site<(N-1);++site){
    MPOMat<1,T> Mb;
    nu.setSite(site);nd.setSite(site);cu.setSite(site);cd.setSite(site);cuD.setSite(site);cdD.setSite(site);id.setSite(site);p.setSite(site);
    Mb[ukey<2>(0,0)] = id;
    Mb[ukey<2>(1,0)] = cu;
    Mb[ukey<2>(2,0)] = cd;
    Mb[ukey<2>(3,0)] = cuD;
    Mb[ukey<2>(4,0)] = cdD;
    //    Mb[ukey<2>(5,0)] = chUp(site)*nu+chDown(site)*nd+U(site)*((nu-0.5*id)*(nd-0.5*id));
    Mb[ukey<2>(5,0)] = chUp(site)*nu+chDown(site)*nd+U(site)*(nu*nd);
     
    Mb[ukey<2>(5,1)] = -hop(site)*cuD*p;
    Mb[ukey<2>(5,2)] = -hop(site)*cdD*p;
    Mb[ukey<2>(5,3)] = hop(site)*cu*p;
    Mb[ukey<2>(5,4)] = hop(site)*cd*p;
    Mb[ukey<2>(5,5)] = id;
    Mb.Dl=6;Mb.Dr=6;
    Mb.squeeze(1e-15);
    (*this)[site]=Mb;
  }
}

template<typename T>
FermiHubbardMPO<T>::FermiHubbardMPO(VectorType hop,VectorType U,VectorType chUp,VectorType chDown,VectorType V){
  uint N =chUp.size();
  assert(chUp.size()==chDown.size());
  MPOMat<1,T> M;
  //  this->resize(N);
  SparseLocalOperator<T> n(4,4),nu(4,4),nd(4,4),cuD(4,4),cdD(4,4),cu(4,4),cd(4,4),id(4,4),p(4,4);//all statistics are set to bosonic (_statistics=1)
  nu[OpKeyType<1>(2,2)]=T(1.0);nu[OpKeyType<1>(3,3)]=T(1.0);
  nd[OpKeyType<1>(1,1)]=T(1.0);nd[OpKeyType<1>(3,3)]=T(1.0);

  cuD[OpKeyType<1>(2,0)] = T(1.0);cuD[OpKeyType<1>(3,1)] = T(1.0); // cuD.setSign(-1);
  cu[OpKeyType<1>(0,2)] = T(1.0);cu[OpKeyType<1>(1,3)] = T(1.0);   // cu.setSign(-1);
								   // 
  cdD[OpKeyType<1>(1,0)] = T(1.0);cdD[OpKeyType<1>(3,2)] = T(-1.0);// cdD.setSign(-1);
  cd[OpKeyType<1>(0,1)] = T(1.0);cd[OpKeyType<1>(2,3)] = T(-1.0);  // cd.setSign(-1);

  id[OpKeyType<1>(0,0)]=T(1);id[OpKeyType<1>(1,1)]=T(1);id[OpKeyType<1>(2,2)]=T(1);id[OpKeyType<1>(3,3)]=T(1);
  p[OpKeyType<1>(0,0)]=T(1);p[OpKeyType<1>(1,1)]=T(-1);p[OpKeyType<1>(2,2)]=T(-1);p[OpKeyType<1>(3,3)]=T(1);

  nu.setSite(0);nd.setSite(0);cu.setSite(0);cd.setSite(0);cuD.setSite(0);cdD.setSite(0);id.setSite(0);p.setSite(0);
  n=nu*nd;n.setSite(0);
  M[ukey<2>(0,0)] = chUp(0)*nu+chDown(0)*nd+U(0)*(nu*nd);
  M[ukey<2>(0,1)] = -hop(0)*(cuD*p);
  M[ukey<2>(0,2)] = -hop(0)*(cdD*p);
  M[ukey<2>(0,3)] = hop(0)*(cu*p);//trying out ....
  M[ukey<2>(0,4)] = hop(0)*(cd*p);
  M[ukey<2>(0,5)] = V(0)*n;
  M[ukey<2>(0,6)] = id;

  M.Dl=1;M.Dr=7;
  M.squeeze(1e-15);
  (*this)[0]=M;
  M.clear();

  nu.setSite(N-1);nd.setSite(N-1);cu.setSite(N-1);cd.setSite(N-1);cuD.setSite(N-1);cdD.setSite(N-1);id.setSite(N-1);p.setSite(N-1);
  n.setSite(N-1);
  M[ukey<2>(0,0)] = id;
  M[ukey<2>(1,0)] = cu;
  M[ukey<2>(2,0)] = cd;
  M[ukey<2>(3,0)] = cuD;
  M[ukey<2>(4,0)] = cdD;
  M[ukey<2>(5,0)] = n;
  M[ukey<2>(6,0)] = chUp(N-1)*nu+chDown(N-1)*nd+U(N-1)*(nu*nd);
  M.Dl=7;M.Dr=1;
  M.squeeze(1e-15);
  (*this)[N-1]=M;
  for(uint site = 1;site<(N-1);++site){
    MPOMat<1,T> Mb;
    nu.setSite(site);nd.setSite(site);cu.setSite(site);cd.setSite(site);cuD.setSite(site);cdD.setSite(site);id.setSite(site);p.setSite(site);
    n.setSite(site);
    Mb[ukey<2>(0,0)] = id;
    Mb[ukey<2>(1,0)] = cu;
    Mb[ukey<2>(2,0)] = cd;
    Mb[ukey<2>(3,0)] = cuD;
    Mb[ukey<2>(4,0)] = cdD;
    Mb[ukey<2>(5,0)] = n;
    Mb[ukey<2>(6,0)] = chUp(site)*nu+chDown(site)*nd+U(site)*(nu*nd);
     
    Mb[ukey<2>(6,1)] = -hop(site)*cuD*p;
    Mb[ukey<2>(6,2)] = -hop(site)*cdD*p;
    Mb[ukey<2>(6,3)] = hop(site)*cu*p;
    Mb[ukey<2>(6,4)] = hop(site)*cd*p;
    Mb[ukey<2>(6,5)] = V(site)*n;
    Mb[ukey<2>(6,6)] = id;
    Mb.Dl=7;Mb.Dr=7;
    Mb.squeeze(1e-15);
    (*this)[site]=Mb;
  }
}


//this class implements time dependent  fermi hubbard hamiltonian; N and Tmax are linear dimension of the system and the 
//maximum time step, hop V and chpot are files containing the time dependent parameters in matrix form, such that
//number of rows (x-dimension)=Tmax, number of columns (y-dimension)=N (or N-1 in the case of hop).
//to get the mpo at a certain step, call TDSpinlessFermionsMPO.construct(step). This operation materializes the MPO at time step 
//step: mpo.construct(step); now mpo can be used in the common way.
//still under construction
// template<typename T>
// class TDFermiHubbardMPO:public FermiHubbardMPO<T>
// {
// public:
//   TDFermiHubbardMPO(const uint N,const uint Tmax,const std::string hop,const std::string U,const std::string chUp,const std::string chDown);
//   TDFermiHubbardMPO(VectorType hop,VectorType U,VectorType chUp,VectorType chDown);
//   ~TDFermiHubbardMPO(){};
//   virtual T measure(const MPS<AbKeyType,T>&mps)const {return FermiHubbardMPO<T>::measure(mps);}
//   virtual T measure(const CanonizedMPS<AbKeyType,T>&mps)const {return FermiHubbardMPO<T>::measure(mps);}
//   virtual void construct(const uint step);
//   template<class Archive>
//   void serialize(Archive & ar, const unsigned int version){
//     ar & boost::serialization::base_object<MPO<T> >(*this);
//     ar &_hop;
//     ar &_U;
//     ar &_chUp;
//     ar &_chDown;
//     ar &_N;
//     ar &_T;
//   }

// private:

//   MatrixType _hop,_U,_chUp,_chDown;
//   uint _N,_T;
// };


// template <typename T>
// TDFermiHubbardMPO<T>::TDFermiHubbardMPO(const uint N,const uint Tmax,const std::string hop,const std::string U,const std::string chUp,const std::string chDown){
//   _N=N,_T=Tmax;
//   std::ifstream infile;

//   infile.open(hop.c_str());
//   cout<<"#TDFermiHuubardMPO<T>: reading hopping parameters from file "<<hop<<endl;
//   _hop.resize(_T,_N-1);
//   boost_setTo(_hop,0.0);
//   for(uint t=0;t<_T;++t){
//     for(uint n=0;n<_N-1;++n){
//       infile>>_hop(t,n);
//     }
//   }
//   infile.close();

//   infile.open(U.c_str());
//   cout<<"#TDFermiHuubardMPO<T>: reading U parameters from file "<<U<<endl;
//   _U.resize(_T,_N);
//   boost_setTo(_U,0.0);
//   for(uint t=0;t<_T;++t){
//     for(uint n=0;n<_N;++n){
//       infile>>_U(t,n);
//     }
//   }
//   infile.close();

//   infile.open(chUp.c_str());
//   cout<<"#TDFermiHuubardMPO<T>: reading muUp parameters from file "<<chUp<<endl;
//   _muup.resize(_T,_N);
//   boost_setTo(_muup,0.0);
//   for(uint t=0;t<_T;++t){
//     for(uint n=0;n<_N;++n){
//       infile>>_muup(t,n);
//     }
//   }
//   infile.close();

//   infile.open(chDown.c_str());
//   cout<<"#TDFermiHuubardMPO<T>: reading muDown parameters from file "<<chDown<<endl;
//   _mudown.resize(_T,_N);
//   boost_setTo(_mudown,0.0);
//   for(uint t=0;t<_T;++t){
//     for(uint n=0;n<_N;++n){
//       infile>>_mudown(t,n);
//     }
//   }
//   infile.close();
//   cout<<"#done reading in parameters"<<endl;
// }

// template <typename T>
// void TDFermiHubbardMPO<T>::construct(const uint step){
//   if(step>_T){
//     cerr<<"TDFermiHubbardMPO<T>::construct: step>_T"<<endl;
//     abort();
//   }
//   VectorType h(_N-1),U(_N),mu(_N),md(_N);
  
//   for(uint n=0;n<_N;++n){
//     if(n<_N-1)
//       h(n)=_hop(step,n);
//     U(n)=_U(step,n);
//     mu(n)=_muup(step,n);
//     md(n)=_mudown(step,n);
//   }
//   (*this)=TDFermiHubbardMPO<T>(h,U,mu,md);
//}







//p-wave superconductor!
template<typename T>
class FermiHubbardSCMPO:public MPO<T>
{
public:
  FermiHubbardSCMPO(VectorType hop,VectorType U,VectorType chUp,VectorType chDown,VectorType Gamma,VectorType Delta);
  ~FermiHubbardSCMPO(){};
  virtual T measure(const MPS<AbKeyType,T>&mps)const {return MPO<T>::measure(mps);}
  virtual T measure(const CanonizedMPS<AbKeyType,T>&mps)const {return MPO<T>::measure(mps);}
};

template<typename T>
FermiHubbardSCMPO<T>::FermiHubbardSCMPO(VectorType hop,VectorType U,VectorType chUp,VectorType chDown,VectorType Gamma,VectorType Delta){
  uint N =chUp.size();
  assert(chUp.size()==chDown.size());
  MPOMat<1,T> M;

  SparseLocalOperator<T> nu(4,4),nd(4,4),cuD(4,4),cdD(4,4),cu(4,4),cd(4,4),id(4,4),p(4,4);//all statistics are set to bosonic (_statistics=1)
  nu[OpKeyType<1>(2,2)]=T(1.0);nu[OpKeyType<1>(3,3)]=T(1.0);
  nd[OpKeyType<1>(1,1)]=T(1.0);nd[OpKeyType<1>(3,3)]=T(1.0);

  cuD[OpKeyType<1>(2,0)] = T(1.0);cuD[OpKeyType<1>(3,1)] = T(1.0); 
  cu[OpKeyType<1>(0,2)] = T(1.0);cu[OpKeyType<1>(1,3)] = T(1.0);   
								   
  cdD[OpKeyType<1>(1,0)] = T(1.0);cdD[OpKeyType<1>(3,2)] = T(-1.0);
  cd[OpKeyType<1>(0,1)] = T(1.0);cd[OpKeyType<1>(2,3)] = T(-1.0);  

  id[OpKeyType<1>(0,0)]=T(1);id[OpKeyType<1>(1,1)]=T(1);id[OpKeyType<1>(2,2)]=T(1);id[OpKeyType<1>(3,3)]=T(1);
  p[OpKeyType<1>(0,0)]=T(1);p[OpKeyType<1>(1,1)]=T(-1);p[OpKeyType<1>(2,2)]=T(-1);p[OpKeyType<1>(3,3)]=T(1);
  

  nu.setSite(0);nd.setSite(0);cu.setSite(0);cd.setSite(0);cuD.setSite(0);cdD.setSite(0);id.setSite(0);p.setSite(0);
  M[ukey<2>(0,0)] = chUp(0)*nu+chDown(0)*nd+U(0)*(nu*nd)+Gamma(0)*((cuD*cdD)+(cd*cu));
  M[ukey<2>(0,1)] = -hop(0)*(cuD*p);
  M[ukey<2>(0,2)] = -hop(0)*(cdD*p)-Delta(0)*(cu*p);
  M[ukey<2>(0,3)] = hop(0)*(cu*p);
  M[ukey<2>(0,4)] = hop(0)*(cd*p)+Delta(0)*(cuD*p);
  M[ukey<2>(0,5)] = id;

  M.Dl=1;M.Dr=6;
  M.squeeze(1e-15);
  (*this)[0]=M;
  M.clear();

  nu.setSite(N-1);nd.setSite(N-1);cu.setSite(N-1);cd.setSite(N-1);cuD.setSite(N-1);cdD.setSite(N-1);id.setSite(N-1);p.setSite(N-1);
  M[ukey<2>(0,0)] = id;
  M[ukey<2>(1,0)] = cu;
  M[ukey<2>(2,0)] = cd;
  M[ukey<2>(3,0)] = cuD;
  M[ukey<2>(4,0)] = cdD;
  M[ukey<2>(5,0)] = chUp(N-1)*nu+chDown(N-1)*nd+U(N-1)*(nu*nd)+Gamma(N-1)*((cuD*cdD)+(cd*cu));//order of up and down repsects fermionic sign
  M.Dl=6;M.Dr=1;
  M.squeeze(1e-15);
  (*this)[N-1]=M;
  for(uint site = 1;site<(N-1);++site){
    MPOMat<1,T> Mb;
    nu.setSite(site);nd.setSite(site);cu.setSite(site);cd.setSite(site);cuD.setSite(site);cdD.setSite(site);id.setSite(site);p.setSite(site);
    Mb[ukey<2>(0,0)] = id;
    Mb[ukey<2>(1,0)] = cu;
    Mb[ukey<2>(2,0)] = cd;
    Mb[ukey<2>(3,0)] = cuD;
    Mb[ukey<2>(4,0)] = cdD;
    Mb[ukey<2>(5,0)] = chUp(site)*nu+chDown(site)*nd+U(site)*(nu*nd)+Gamma(site)*((cuD*cdD)+(cd*cu));//order of up and down repsects fermionic sign
     
    Mb[ukey<2>(5,1)] = -hop(site)*cuD*p;
    Mb[ukey<2>(5,2)] = -hop(site)*cdD*p-Delta(site)*(cu*p);
    Mb[ukey<2>(5,3)] = hop(site)*cu*p;
    Mb[ukey<2>(5,4)] = hop(site)*cd*p+Delta(site)*(cuD*p);
    Mb[ukey<2>(5,5)] = id;
    Mb.Dl=6;Mb.Dr=6;
    Mb.squeeze(1e-15);
    (*this)[site]=Mb;
  }
}


// template<typename T>
// class UnfoldedSiamMPO:public MPO<T>
// {
// public:
//   UnfoldedSiamMPO(VectorType hop,VectorType U,VectorType chUp,VectorType chDown);
//   UnfoldedSiamMPO(VectorType hop,VectorType U,VectorType chUp,VectorType chDown,VectorType V);
//   ~UnfoldedSiamMPO(){};
//   virtual T measure(const MPS<AbKeyType,T>&mps)const {return MPO<T>::measure(mps);}
//   virtual T measure(const CanonizedMPS<AbKeyType,T>&mps)const {return MPO<T>::measure(mps);}
//   void shiftandrescale(const T E0,const T a, const T W){
//     SparseLocalOperator<T> id(4,4);
//     id[OpKeyType<1>(0,0)]=T(1.);id[OpKeyType<1>(1,1)]=T(1.);id[OpKeyType<1>(2,2)]=T(1.);id[OpKeyType<1>(3,3)]=T(1.);
//     id*=E0;
//     (*this)[0][ukey<2>(0,0)]-=id;
//     typename MPOMat<1,T>::iterator o;
//     for(o=(*this)[0].begin();o!=(*this)[0].end();++o)
//       o->second/=a;
//     id[OpKeyType<1>(0,0)]=T(1.);id[OpKeyType<1>(1,1)]=T(1.);id[OpKeyType<1>(2,2)]=T(1.);id[OpKeyType<1>(3,3)]=T(1.);
//     id*=W;
//     (*this)[0][ukey<2>(0,0)]-=id;
//     //that's it ...
//   }


// };
// template<typename T>
// UnfoldedSiamMPO<T>::UnfoldedSiamMPO(VectorType hop,VectorType U,VectorType chUp,VectorType chDown){
//   uint N =2*chUp.size();
//   assert(chUp.size()==chDown.size());
//   MPOMat<1,T> M;
//   //  this->resize(N);
//   SparseLocalOperator<T> nu(2,2),nd(2,2),cuD(2,2),cdD(2,2),cu(2,2),cd(2,2),id(2,2),p(2,2);//all statistics are set to bosonic (_statistics=1)
//   nu[OpKeyType<1>(0,0)]=T(0.0);nu[OpKeyType<1>(1,1)]=T(1.0);
//   nd[OpKeyType<1>(0,0)]=T(0.0);nd[OpKeyType<1>(1,1)]=T(1.0);

//   cuD[OpKeyType<1>(1,0)] = T(1.0);
//   cu[OpKeyType<1>(0,1)] = T(1.0);
// 								   // 
//   cdD[OpKeyType<1>(1,0)] = T(1.0);
//   cd[OpKeyType<1>(0,1)] = T(1.0);

//   id[OpKeyType<1>(0,0)]=T(1);id[OpKeyType<1>(1,1)]=T(1);id[OpKeyType<1>(2,2)]=T(1);id[OpKeyType<1>(3,3)]=T(1);
//   p[OpKeyType<1>(0,0)]=T(1);p[OpKeyType<1>(1,1)]=T(-1);p[OpKeyType<1>(2,2)]=T(-1);p[OpKeyType<1>(3,3)]=T(1);

//   nu.setSite(0);nd.setSite(0);cu.setSite(0);cd.setSite(0);cuD.setSite(0);cdD.setSite(0);id.setSite(0);p.setSite(0);
//   M[ukey<2>(0,0)] = chUp(0)*nu+chDown(0)*nd+U(0)*(nu*nd);
//   M[ukey<2>(0,1)] = -hop(0)*(cuD*p);
//   M[ukey<2>(0,2)] = -hop(0)*(cdD*p);
//   M[ukey<2>(0,3)] = hop(0)*(cu*p);//trying out ....
//   M[ukey<2>(0,4)] = hop(0)*(cd*p);
//   M[ukey<2>(0,5)] = id;

//   M.Dl=1;M.Dr=6;
//   M.squeeze(1e-15);
//   (*this)[0]=M;
//   M.clear();

//   nu.setSite(N-1);nd.setSite(N-1);cu.setSite(N-1);cd.setSite(N-1);cuD.setSite(N-1);cdD.setSite(N-1);id.setSite(N-1);p.setSite(N-1);
//   M[ukey<2>(0,0)] = id;
//   M[ukey<2>(1,0)] = cu;
//   M[ukey<2>(2,0)] = cd;
//   M[ukey<2>(3,0)] = cuD;
//   M[ukey<2>(4,0)] = cdD;
//   M[ukey<2>(5,0)] = chUp(N-1)*nu+chDown(N-1)*nd+U(N-1)*(nu*nd);
//   M.Dl=6;M.Dr=1;
//   M.squeeze(1e-15);
//   (*this)[N-1]=M;
//   for(uint site = 1;site<(N-1);++site){
//     MPOMat<1,T> Mb;
//     nu.setSite(site);nd.setSite(site);cu.setSite(site);cd.setSite(site);cuD.setSite(site);cdD.setSite(site);id.setSite(site);p.setSite(site);
//     Mb[ukey<2>(0,0)] = id;
//     Mb[ukey<2>(1,0)] = cu;
//     Mb[ukey<2>(2,0)] = cd;
//     Mb[ukey<2>(3,0)] = cuD;
//     Mb[ukey<2>(4,0)] = cdD;
//     Mb[ukey<2>(5,0)] = chUp(site)*nu+chDown(site)*nd+U(site)*(nu*nd);
     
//     Mb[ukey<2>(5,1)] = -hop(site)*cuD*p;
//     Mb[ukey<2>(5,2)] = -hop(site)*cdD*p;
//     Mb[ukey<2>(5,3)] = hop(site)*cu*p;//trying out ...
//     Mb[ukey<2>(5,4)] = hop(site)*cd*p;
//     Mb[ukey<2>(5,5)] = id;
//     Mb.Dl=6;Mb.Dr=6;
//     Mb.squeeze(1e-15);
//     (*this)[site]=Mb;
//   }
// }



//this is a fermi hubbard ladder system, string order operators could be wrong ...
//for completeness,here is the state numbering I use; first one is the lower chain, the secon one the upper
//0  0,0
//1  d,0
//2  u,0
//3  ud,0
//4  0,d
//5  d,d
//6  u,d
//7  ud,d
//8  0,u
//9  d,u
//10 u,u
//11 ud,u
//12 0,ud
//13 d,ud
//14 u,ud
//15 ud,ud
template<typename T>
class FermiHubbardLadderMPO:public MPO<T>
{
public:
  FermiHubbardLadderMPO(VectorType hopL,VectorType hopU,VectorType UL,VectorType UU,VectorType VL,VectorType VU,VectorType chUpL,VectorType chDownL,
			VectorType chUpU,VectorType chDownU,VectorType hopInt,VectorType Vint);
  ~FermiHubbardLadderMPO(){};
  virtual T measure(const MPS<AbKeyType,T>&mps)const {return MPO<T>::measure(mps);}
  virtual T measure(const CanonizedMPS<AbKeyType,T>&mps)const {return MPO<T>::measure(mps);}

};
template<typename T>
FermiHubbardLadderMPO<T>::  FermiHubbardLadderMPO(VectorType hopL,VectorType hopU,VectorType UL,VectorType UU,VectorType VL,VectorType VU,VectorType chUpL,
						  VectorType chDownL,VectorType chUpU,VectorType chDownU,VectorType hopInt,VectorType Vint){
  const uint N =chUpL.size();
  assert((N)==chUpU.size());
  assert((N)==chDownL.size());
  assert((N)==chDownU.size());
  assert(hopInt.size()==(N));
  MPOMat<1,T> M;

  //the Lower chain
  SparseLocalOperator<T> nuU(16,16),ndU(16,16),cuDU(16,16),cdDU(16,16),cuU(16,16),cdU(16,16),id(16,16),p(16,16),pU(16,16),pL(16,16);
  SparseLocalOperator<T> nuL(16,16),ndL(16,16),cuDL(16,16),cdDL(16,16),cuL(16,16),cdL(16,16),nU(16,16),nL(16,16);

//  cuDL[OpKeyType<1>(2,0)] = T(1.0);cuDL[OpKeyType<1>(3,1)] = T(1.0);cuDL[OpKeyType<1>(6,4)] = T(1.0);cuDL[OpKeyType<1>(7,5)] = T(1.0);
//  cuDL[OpKeyType<1>(10,8)] = T(1.0);cuDL[OpKeyType<1>(11,9)] = T(1.0);cuDL[OpKeyType<1>(14,12)] = T(1.0);cuDL[OpKeyType<1>(15,13)] = T(1.0);
// 
//  cuL[OpKeyType<1>(0,2)] = T(1.0);cuL[OpKeyType<1>(1,3)] = T(1.0);cuL[OpKeyType<1>(4,6)] = T(1.0);cuL[OpKeyType<1>(5,7)] = T(1.0);
//  cuL[OpKeyType<1>(8,10)] = T(1.0);cuL[OpKeyType<1>(9,11)] = T(1.0);cuL[OpKeyType<1>(12,14)] = T(1.0);cuL[OpKeyType<1>(13,15)] = T(1.0);
//								    
//  cdL[OpKeyType<1>(0,1)] = T(1.0);cdL[OpKeyType<1>(2,3)] = T(1.0);cdL[OpKeyType<1>(4,5)] = T(1.0);cdL[OpKeyType<1>(6,7)] = T(1.0);
//  cdL[OpKeyType<1>(8,9)] = T(1.0);cdL[OpKeyType<1>(10,11)] = T(1.0);cdL[OpKeyType<1>(12,13)] = T(1.0);cdL[OpKeyType<1>(14,15)] = T(1.0);
//
//  cdDL[OpKeyType<1>(1,0)] = T(1.0);cdDL[OpKeyType<1>(3,2)] = T(1.0);cdDL[OpKeyType<1>(5,4)] = T(1.0);cdDL[OpKeyType<1>(7,6)] = T(1.0);
//  cdDL[OpKeyType<1>(9,8)] = T(1.0);cdDL[OpKeyType<1>(11,10)] = T(1.0);cdDL[OpKeyType<1>(13,12)] = T(1.0);cdDL[OpKeyType<1>(15,14)] = T(1.0);
//
//  //the Upper chain
//  cuDU[OpKeyType<1>(8,0)] = T(1.0);cuDU[OpKeyType<1>(9,1)] = T(1.0);cuDU[OpKeyType<1>(10,2)] = T(1.0);cuDU[OpKeyType<1>(11,3)] = T(1.0);
//  cuDU[OpKeyType<1>(12,4)] = T(1.0);cuDU[OpKeyType<1>(13,5)] = T(1.0);cuDU[OpKeyType<1>(14,6)] = T(1.0);cuDU[OpKeyType<1>(15,7)] = T(1.0);
//     
//  cuU[OpKeyType<1>(0,8)] = T(1.0);cuU[OpKeyType<1>(1,9)] = T(1.0);cuU[OpKeyType<1>(2,10)] = T(1.0);cuU[OpKeyType<1>(3,11)] = T(1.0);
//  cuU[OpKeyType<1>(4,12)] = T(1.0);cuU[OpKeyType<1>(5,13)] = T(1.0);cuU[OpKeyType<1>(6,14)] = T(1.0);cuU[OpKeyType<1>(7,15)] = T(1.0);
//     
//  cdDU[OpKeyType<1>(4,0)] = T(1.0);cdDU[OpKeyType<1>(5,1)] = T(1.0);cdDU[OpKeyType<1>(6,2)] = T(1.0);cdDU[OpKeyType<1>(7,3)] = T(1.0);
//  cdDU[OpKeyType<1>(12,8)] = T(1.0);cdDU[OpKeyType<1>(13,9)] = T(1.0);cdDU[OpKeyType<1>(14,10)] = T(1.0);cdDU[OpKeyType<1>(15,11)] = T(1.0);
//
//  cdU[OpKeyType<1>(0,4)] = T(1.0);cdU[OpKeyType<1>(1,5)] = T(1.0);cdU[OpKeyType<1>(2,6)] = T(1.0);cdU[OpKeyType<1>(3,7)] = T(1.0);
//  cdU[OpKeyType<1>(8,12)] = T(1.0);cdU[OpKeyType<1>(9,13)] = T(1.0);cdU[OpKeyType<1>(10,14)] = T(1.0);cdU[OpKeyType<1>(11,15)] = T(1.0);
//

//

  cuDL[OpKeyType<1>(2,0)] = T(1.0);cuDL[OpKeyType<1>(3,1)] = T(1.0);cuDL[OpKeyType<1>(6,4)] = T(1.0);cuDL[OpKeyType<1>(7,5)] = T(1.0);
  cuDL[OpKeyType<1>(10,8)] = T(1.0);cuDL[OpKeyType<1>(11,9)] = T(1.0);cuDL[OpKeyType<1>(14,12)] = T(1.0);cuDL[OpKeyType<1>(15,13)] = T(1.0);
 
  cuL[OpKeyType<1>(0,2)] = T(1.0);cuL[OpKeyType<1>(1,3)] = T(1.0);cuL[OpKeyType<1>(4,6)] = T(1.0);cuL[OpKeyType<1>(5,7)] = T(1.0);
  cuL[OpKeyType<1>(8,10)] = T(1.0);cuL[OpKeyType<1>(9,11)] = T(1.0);cuL[OpKeyType<1>(12,14)] = T(1.0);cuL[OpKeyType<1>(13,15)] = T(1.0);
								    
  cdL[OpKeyType<1>(0,1)] = T(1.0);cdL[OpKeyType<1>(2,3)] = T(-1.0);cdL[OpKeyType<1>(4,5)] = T(1.0);cdL[OpKeyType<1>(6,7)] = T(-1.0);
  cdL[OpKeyType<1>(8,9)] = T(1.0);cdL[OpKeyType<1>(10,11)] = T(-1.0);cdL[OpKeyType<1>(12,13)] = T(1.0);cdL[OpKeyType<1>(14,15)] = T(-1.0);

  cdDL[OpKeyType<1>(1,0)] = T(1.0);cdDL[OpKeyType<1>(3,2)] = T(-1.0);cdDL[OpKeyType<1>(5,4)] = T(1.0);cdDL[OpKeyType<1>(7,6)] = T(-1.0);
  cdDL[OpKeyType<1>(9,8)] = T(1.0);cdDL[OpKeyType<1>(11,10)] = T(-1.0);cdDL[OpKeyType<1>(13,12)] = T(1.0);cdDL[OpKeyType<1>(15,14)] = T(-1.0);

  //the Upper chain
  cuDU[OpKeyType<1>(8,0)] = T(1.0);cuDU[OpKeyType<1>(9,1)] = T(1.0);cuDU[OpKeyType<1>(10,2)] = T(1.0);cuDU[OpKeyType<1>(11,3)] = T(1.0);
  cuDU[OpKeyType<1>(12,4)] = T(1.0);cuDU[OpKeyType<1>(13,5)] = T(1.0);cuDU[OpKeyType<1>(14,6)] = T(1.0);cuDU[OpKeyType<1>(15,7)] = T(1.0);
     
  cuU[OpKeyType<1>(0,8)] = T(1.0);cuU[OpKeyType<1>(1,9)] = T(1.0);cuU[OpKeyType<1>(2,10)] = T(1.0);cuU[OpKeyType<1>(3,11)] = T(1.0);
  cuU[OpKeyType<1>(4,12)] = T(1.0);cuU[OpKeyType<1>(5,13)] = T(1.0);cuU[OpKeyType<1>(6,14)] = T(1.0);cuU[OpKeyType<1>(7,15)] = T(1.0);
     
  cdDU[OpKeyType<1>(4,0)] = T(1.0);cdDU[OpKeyType<1>(5,1)] = T(1.0);cdDU[OpKeyType<1>(6,2)] = T(1.0);cdDU[OpKeyType<1>(7,3)] = T(1.0);
  cdDU[OpKeyType<1>(12,8)] = T(-1.0);cdDU[OpKeyType<1>(13,9)] = T(-1.0);cdDU[OpKeyType<1>(14,10)] = T(-1.0);cdDU[OpKeyType<1>(15,11)] = T(-1.0);

  cdU[OpKeyType<1>(0,4)] = T(1.0);cdU[OpKeyType<1>(1,5)] = T(1.0);cdU[OpKeyType<1>(2,6)] = T(1.0);cdU[OpKeyType<1>(3,7)] = T(1.0);
  cdU[OpKeyType<1>(8,12)] = T(-1.0);cdU[OpKeyType<1>(9,13)] = T(-1.0);cdU[OpKeyType<1>(10,14)] = T(-1.0);cdU[OpKeyType<1>(11,15)] = T(-1.0);

  SparseLocalOperator<T> bla(16,16);
  
  for(uint ind=0;ind<16;++ind)
    id[OpKeyType<1>(ind,ind)]=T(1);

  //this seems to be right, so DON'T CHANGE ANYTHING HERE
  p[OpKeyType<1>(0,0)]=T(1);p[OpKeyType<1>(1,1)]=T(-1);p[OpKeyType<1>(2,2)]=T(-1);p[OpKeyType<1>(3,3)]=T(1);
  p[OpKeyType<1>(4,4)]=T(-1);p[OpKeyType<1>(5,5)]=T(1);p[OpKeyType<1>(6,6)]=T(1);p[OpKeyType<1>(7,7)]=T(-1);
  p[OpKeyType<1>(8,8)]=T(-1);p[OpKeyType<1>(9,9)]=T(1);p[OpKeyType<1>(10,10)]=T(1);p[OpKeyType<1>(11,11)]=T(-1);
  p[OpKeyType<1>(12,12)]=T(1);p[OpKeyType<1>(13,13)]=T(-1);p[OpKeyType<1>(14,14)]=T(-1);p[OpKeyType<1>(15,15)]=T(1);

  pL[OpKeyType<1>(0,0)]=T(1);pL[OpKeyType<1>(1,1)]=T(-1);pL[OpKeyType<1>(2,2)]=T(-1);pL[OpKeyType<1>(3,3)]=T(1);
  pL[OpKeyType<1>(4,4)]=T(1);pL[OpKeyType<1>(5,5)]=T(-1);pL[OpKeyType<1>(6,6)]=T(-1);pL[OpKeyType<1>(7,7)]=T(1);
  pL[OpKeyType<1>(8,8)]=T(1);pL[OpKeyType<1>(9,9)]=T(-1);pL[OpKeyType<1>(10,10)]=T(-1);pL[OpKeyType<1>(11,11)]=T(1);
  pL[OpKeyType<1>(12,12)]=T(1);pL[OpKeyType<1>(13,13)]=T(-1);pL[OpKeyType<1>(14,14)]=T(-1);pL[OpKeyType<1>(15,15)]=T(1);
 
  pU[OpKeyType<1>(0,0)]=T(1);pU[OpKeyType<1>(1,1)]=T(1);pU[OpKeyType<1>(2,2)]=T(1);pU[OpKeyType<1>(3,3)]=T(1);
  pU[OpKeyType<1>(4,4)]=T(-1);pU[OpKeyType<1>(5,5)]=T(-1);pU[OpKeyType<1>(6,6)]=T(-1);pU[OpKeyType<1>(7,7)]=T(-1);
  pU[OpKeyType<1>(8,8)]=T(1);pU[OpKeyType<1>(9,9)]=T(1);pU[OpKeyType<1>(10,10)]=T(1);pU[OpKeyType<1>(11,11)]=T(1);
  pU[OpKeyType<1>(12,12)]=T(-1);pU[OpKeyType<1>(13,13)]=T(-1);pU[OpKeyType<1>(14,14)]=T(-1);pU[OpKeyType<1>(15,15)]=T(-1);
  
  bla=pL*cuU;cuU=bla;
  bla=pL*cdU;cdU=bla;
  bla=pL*cuDU;cuDU=bla;
  bla=pL*cdDU;cdDU=bla;

    // cout<<"commutators between same spin"<<endl;
    // cout<<"    equal sites"<<endl;
    // cout<<"{cuU,cuDU}"<<endl;
    // bla=cuU*cuDU+cuDU*cuU;
    // cout<<bla<<endl;

    // cout<<"{cdU,cdDU}"<<endl;
    // bla=cdU*cdDU+cdDU*cdU;
    // cout<<bla<<endl;

    // cout<<"{cuL,cuDL}"<<endl;
    // bla=cuL*cuDL+cuDL*cuL;
    // cout<<bla<<endl;

    // cout<<"{cdL,cdDL}"<<endl;
    // bla=cdL*cdDL+cdDL*cdL;
    // cout<<bla<<endl;
    // cout<<"    different sites"<<endl;

    // cout<<"{cuU,cuDL}"<<endl;
    // bla=cuU*cuDL+cuDL*cuU;
    // cout<<bla<<endl;

    // cout<<"{cdU,cdDL}"<<endl;
    // bla=cdU*cdDL+cdDL*cdU;
    // cout<<bla<<endl;
  
    // cout<<"{cuU,cuL}"<<endl;
    // bla=cuU*cuL+cuL*cuU;
    // cout<<bla<<endl;

  
    // cout<<"{cdU,cdL}"<<endl;
    // bla=cdU*cdL+cdL*cdU;
    // cout<<bla<<endl;


    // cout<<"commutators between different spin"<<endl;
    // cout<<"    equal sites"<<endl;
    // cout<<"{cdU,cuDU}"<<endl;
    // bla=cdU*cuDU+cuDU*cdU;
    // cout<<bla<<endl;

    // cout<<"{cuU,cdDU}"<<endl;
    // bla=cuU*cdDU+cdDU*cuU;
    // cout<<bla<<endl;

    // cout<<"{cdL,cuDL}"<<endl;
    // bla=cdL*cuDL+cuDL*cdL;
    // cout<<bla<<endl;

    // cout<<"{cuL,cdDL}"<<endl;
    // bla=cuL*cdDL+cdDL*cuL;
    // cout<<bla<<endl;

    // cout<<"{cdL,cuL}"<<endl;
    // bla=cdL*cuL+cuL*cdL;
    // cout<<bla<<endl;

    // cout<<"{cdU,cuU}"<<endl;
    // bla=cdU*cuU+cuU*cdU;
    // cout<<bla<<endl;

    // cout<<"    different sites"<<endl;

    // cout<<"{cdU,cuDL}"<<endl;
    // bla=cdU*cuDL+cuDL*cdU;
    // cout<<bla<<endl;

    // cout<<"{cuU,cdDL}"<<endl;
    // bla=cuU*cdDL+cdDL*cuU;
    // cout<<bla<<endl;
  
    // cout<<"{cdU,cuL}"<<endl;
    // bla=cdU*cuL+cuL*cdU;
    // cout<<bla<<endl;

  
    // cout<<"{cuU,cdL}"<<endl;
    // bla=cuU*cdL+cdL*cuU;
    // cout<<bla<<endl;

    // my_pause();
  // abort();

  // for(uint ind=0;ind<16;++ind)
  //   p[OpKeyType<1>(ind,ind)]=T(1);

  nuL[OpKeyType<1>(2,2)]=T(1);nuL[OpKeyType<1>(3,3)]=T(1);nuL[OpKeyType<1>(6,6)]=T(1);nuL[OpKeyType<1>(7,7)]=T(1);
  nuL[OpKeyType<1>(10,10)]=T(1);nuL[OpKeyType<1>(11,11)]=T(1);nuL[OpKeyType<1>(14,14)]=T(1);nuL[OpKeyType<1>(15,15)]=T(1);

  ndL[OpKeyType<1>(1,1)]=T(1);ndL[OpKeyType<1>(3,3)]=T(1);ndL[OpKeyType<1>(5,5)]=T(1);ndL[OpKeyType<1>(7,7)]=T(1);
  ndL[OpKeyType<1>(9,9)]=T(1);ndL[OpKeyType<1>(11,11)]=T(1);ndL[OpKeyType<1>(13,13)]=T(1);ndL[OpKeyType<1>(15,15)]=T(1);
  nL=nuL+ndL;
  
  for(uint ind=8;ind<16;++ind)
    nuU[OpKeyType<1>(ind,ind)]=T(1);
  for(uint ind=4;ind<8;++ind)
    ndU[OpKeyType<1>(ind,ind)]=T(1);
  for(uint ind=12;ind<16;++ind)
    ndU[OpKeyType<1>(ind,ind)]=T(1);
  nU=nuU+ndU;

  nuU.setSite(0);ndU.setSite(0);cuU.setSite(0);cdU.setSite(0);cuDU.setSite(0);cdDU.setSite(0);id.setSite(0);p.setSite(0);
  nuL.setSite(0);ndL.setSite(0);cuL.setSite(0);cdL.setSite(0);cuDL.setSite(0);cdDL.setSite(0);pU.setSite(0);pL.setSite(0);
  nL.setSite(0);
  nU.setSite(0);
  // M[ukey<2>(0,0)] = chUpU(0)*nuU+chDownU(0)*ndU+UU(0)*(nuU*ndU)+chUpL(0)*nuL+chDownL(0)*ndL+UL(0)*(nuL*ndL)-
  //   hopInt(0)*((pL*cuDL)*cuU-cuDU*(pL*cuL)+(pL*cdDL)*cdU-cdDU*(pL*cdL))+Vint(0)*((nL-id)*(nU-id)) ;
  M[ukey<2>(0,0)] = chUpU(0)*nuU+chDownU(0)*ndU+UU(0)*(nuU*ndU)+chUpL(0)*nuL+chDownL(0)*ndL+UL(0)*(nuL*ndL)-
    hopInt(0)*(cuDL*cuU+cuDU*cuL+cdDL*cdU+cdDU*cdL)+Vint(0)*((nL-id)*(nU-id)) ;

  M[ukey<2>(0,1)] = -hopL(0)*cuL*p;
  M[ukey<2>(0,2)] = hopL(0)*cuDL*p;//
  M[ukey<2>(0,3)] = -hopU(0)*cuU*p;
  M[ukey<2>(0,4)] = hopU(0)*cuDU*p;//
  M[ukey<2>(0,5)] = -hopL(0)*cdL*p;
  M[ukey<2>(0,6)] = hopL(0)*cdDL*p;//
  M[ukey<2>(0,7)] = -hopU(0)*cdU*p;
  M[ukey<2>(0,8)] = hopU(0)*cdDU*p;//
  M[ukey<2>(0,9)] = VU(0)*(nU-id);
  M[ukey<2>(0,10)] = VL(0)*(nL-id);
  M[ukey<2>(0,11)] = id;
  M.Dl=1,M.Dr=12;
  M.squeeze(1e-15);
  (*this)[0]=M;
  M.clear();

  nuU.setSite(N-1);ndU.setSite(N-1);cuU.setSite(N-1);cdU.setSite(N-1);cuDU.setSite(N-1);cdDU.setSite(N-1);id.setSite(N-1);p.setSite(N-1);
  nuL.setSite(N-1);ndL.setSite(N-1);cuL.setSite(N-1);cdL.setSite(N-1);cuDL.setSite(N-1);cdDL.setSite(N-1);pU.setSite(N-1);pL.setSite(N-1);
  nL.setSite(N-1);
  nU.setSite(N-1);

  M[ukey<2>(0,0)] = id;
  M[ukey<2>(1,0)] = cuDL;
  M[ukey<2>(2,0)] = cuL;
  M[ukey<2>(3,0)] = cuDU;
  M[ukey<2>(4,0)] = cuU;
  M[ukey<2>(5,0)] = cdDL;
  M[ukey<2>(6,0)] = cdL;
  M[ukey<2>(7,0)] = cdDU;
  M[ukey<2>(8,0)] = cdU;
  M[ukey<2>(9,0)] = nU-id;
  M[ukey<2>(10,0)]= nL-id;
  // M[ukey<2>(11,0)] = chUpU(N-1)*nuU+chDownU(N-1)*ndU+UU(N-1)*(nuU*ndU)+chUpL(N-1)*nuL+chDownL(N-1)*ndL+UL(N-1)*(nuL*ndL)-
  //   hopInt(N-1)*((pL*cuDL)*cuU-cuDU*(pL*cuL)+(pL*cdDL)*cdU-cdDU*(pL*cdL))+Vint(N-1)*((nL-id)*(nU-id));
  M[ukey<2>(11,0)] = chUpU(N-1)*nuU+chDownU(N-1)*ndU+UU(N-1)*(nuU*ndU)+chUpL(N-1)*nuL+chDownL(N-1)*ndL+UL(N-1)*(nuL*ndL)-
    hopInt(N-1)*(cuDL*cuU+cuDU*cuL+cdDL*cdU+cdDU*cdL)+Vint(N-1)*((nL-id)*(nU-id));

  M.Dl=12,M.Dr=1;
  M.squeeze(1e-15);
  (*this)[N-1]=M;
  for(uint site = 1;site<(N-1);++site){
    MPOMat<1,T> Mb;
    nuU.setSite(site);ndU.setSite(site);cuU.setSite(site);cdU.setSite(site);cuDU.setSite(site);cdDU.setSite(site);id.setSite(site);p.setSite(site);
    nuL.setSite(site);ndL.setSite(site);cuL.setSite(site);cdL.setSite(site);cuDL.setSite(site);cdDL.setSite(site);pU.setSite(site);pL.setSite(site);
    nL.setSite(site);
    nU.setSite(site);

    Mb[ukey<2>(0,0)] = id;
    Mb[ukey<2>(1,0)] = cuDL;
    Mb[ukey<2>(2,0)] = cuL;
    Mb[ukey<2>(3,0)] = cuDU;
    Mb[ukey<2>(4,0)] = cuU;
    Mb[ukey<2>(5,0)] = cdDL;
    Mb[ukey<2>(6,0)] = cdL;
    Mb[ukey<2>(7,0)] = cdDU;
    Mb[ukey<2>(8,0)] = cdU;
    Mb[ukey<2>(9,0)] = nU-id;
    Mb[ukey<2>(10,0)]= nL-id;
    // Mb[ukey<2>(11,0)] = chUpU(site)*nuU+chDownU(site)*ndU+UU(site)*(nuU*ndU)+chUpL(site)*nuL+chDownL(site)*ndL+UL(site)*(nuL*ndL)
    //   -hopInt(site)*((pL*cuDL)*cuU-cuDU*(pL*cuL)+(pL*cdDL)*cdU-cdDU*(pL*cdL))+Vint(site)*((nL-id)*(nU-id));
    Mb[ukey<2>(11,0)] = chUpU(site)*nuU+chDownU(site)*ndU+UU(site)*(nuU*ndU)+chUpL(site)*nuL+chDownL(site)*ndL+UL(site)*(nuL*ndL)
      -hopInt(site)*(cuDL*cuU+cuDU*cuL+cdDL*cdU+cdDU*cdL)+Vint(site)*((nL-id)*(nU-id));
     
    Mb[ukey<2>(11,1)] = -hopL(site)*cuL*p;
    Mb[ukey<2>(11,2)] = hopL(site)*cuDL*p;//
    Mb[ukey<2>(11,3)] = -hopU(site)*cuU*p;
    Mb[ukey<2>(11,4)] = hopU(site)*cuDU*p;//
    Mb[ukey<2>(11,5)] = -hopL(site)*cdL*p;
    Mb[ukey<2>(11,6)] = hopL(site)*cdDL*p;//
    Mb[ukey<2>(11,7)] = -hopU(site)*cdU*p;
    Mb[ukey<2>(11,8)] = hopU(site)*cdDU*p;//
    Mb[ukey<2>(11,9)] = VU(site)*(nU-id);
    Mb[ukey<2>(11,10)] = VL(site)*(nL-id);
    Mb[ukey<2>(11,11)] = id;
    Mb.Dl=12,Mb.Dr=12;
    Mb.squeeze(1e-15);
    (*this)[site]=Mb;
  }
}




//this is a fermi hubbard ladder system, string order operators could be wrong ...
//for completeness,here is the state numbering I use; first one is the lower chain, the secon one the upper
//0  0,0
//1  d,0
//2  u,0
//3  ud,0
//4  0,d
//5  d,d
//6  u,d
//7  ud,d
//8  0,u
//9  d,u
//10 u,u
//11 ud,u
//12 0,ud
//13 d,ud
//14 u,ud
//15 ud,ud
template<typename T>
class KanamoriMPO:public MPO<T>
{
public:

  KanamoriMPO(VectorType hopL,VectorType hopU,VectorType U1,VectorType U2,VectorType U3,VectorType J1,VectorType J2, VectorType chUpL,VectorType chDownL,VectorType chUpU,VectorType chDownU);
  ~KanamoriMPO(){};
  virtual T measure(const MPS<AbKeyType,T>&mps)const {return MPO<T>::measure(mps);}
  virtual T measure(const CanonizedMPS<AbKeyType,T>&mps)const {return MPO<T>::measure(mps);}

};
template<typename T>

KanamoriMPO<T>::KanamoriMPO(VectorType hopL,VectorType hopU,VectorType U1,VectorType U2,VectorType U3,VectorType J1,VectorType J2, VectorType chUpL,VectorType chDownL,VectorType chUpU,VectorType chDownU){

  const uint N =chUpL.size();
  assert((N)==chUpU.size());
  assert((N)==chDownL.size());
  assert((N)==chDownU.size());
  assert(hopInt.size()==(N));
  MPOMat<1,T> M;

  //the Lower chain
  SparseLocalOperator<T> nuU(16,16),ndU(16,16),cuDU(16,16),cdDU(16,16),cuU(16,16),cdU(16,16),id(16,16),p(16,16),pU(16,16),pL(16,16);
  SparseLocalOperator<T> nuL(16,16),ndL(16,16),cuDL(16,16),cdDL(16,16),cuL(16,16),cdL(16,16),nU(16,16),nL(16,16);


  cuDL[OpKeyType<1>(2,0)] = T(1.0);cuDL[OpKeyType<1>(3,1)] = T(1.0);cuDL[OpKeyType<1>(6,4)] = T(1.0);cuDL[OpKeyType<1>(7,5)] = T(1.0);
  cuDL[OpKeyType<1>(10,8)] = T(1.0);cuDL[OpKeyType<1>(11,9)] = T(1.0);cuDL[OpKeyType<1>(14,12)] = T(1.0);cuDL[OpKeyType<1>(15,13)] = T(1.0);
 
  cuL[OpKeyType<1>(0,2)] = T(1.0);cuL[OpKeyType<1>(1,3)] = T(1.0);cuL[OpKeyType<1>(4,6)] = T(1.0);cuL[OpKeyType<1>(5,7)] = T(1.0);
  cuL[OpKeyType<1>(8,10)] = T(1.0);cuL[OpKeyType<1>(9,11)] = T(1.0);cuL[OpKeyType<1>(12,14)] = T(1.0);cuL[OpKeyType<1>(13,15)] = T(1.0);
								    
  cdL[OpKeyType<1>(0,1)] = T(1.0);cdL[OpKeyType<1>(2,3)] = T(-1.0);cdL[OpKeyType<1>(4,5)] = T(1.0);cdL[OpKeyType<1>(6,7)] = T(-1.0);
  cdL[OpKeyType<1>(8,9)] = T(1.0);cdL[OpKeyType<1>(10,11)] = T(-1.0);cdL[OpKeyType<1>(12,13)] = T(1.0);cdL[OpKeyType<1>(14,15)] = T(-1.0);

  cdDL[OpKeyType<1>(1,0)] = T(1.0);cdDL[OpKeyType<1>(3,2)] = T(-1.0);cdDL[OpKeyType<1>(5,4)] = T(1.0);cdDL[OpKeyType<1>(7,6)] = T(-1.0);
  cdDL[OpKeyType<1>(9,8)] = T(1.0);cdDL[OpKeyType<1>(11,10)] = T(-1.0);cdDL[OpKeyType<1>(13,12)] = T(1.0);cdDL[OpKeyType<1>(15,14)] = T(-1.0);

  //the Upper chain
  cuDU[OpKeyType<1>(8,0)] = T(1.0);cuDU[OpKeyType<1>(9,1)] = T(1.0);cuDU[OpKeyType<1>(10,2)] = T(1.0);cuDU[OpKeyType<1>(11,3)] = T(1.0);
  cuDU[OpKeyType<1>(12,4)] = T(1.0);cuDU[OpKeyType<1>(13,5)] = T(1.0);cuDU[OpKeyType<1>(14,6)] = T(1.0);cuDU[OpKeyType<1>(15,7)] = T(1.0);
     
  cuU[OpKeyType<1>(0,8)] = T(1.0);cuU[OpKeyType<1>(1,9)] = T(1.0);cuU[OpKeyType<1>(2,10)] = T(1.0);cuU[OpKeyType<1>(3,11)] = T(1.0);
  cuU[OpKeyType<1>(4,12)] = T(1.0);cuU[OpKeyType<1>(5,13)] = T(1.0);cuU[OpKeyType<1>(6,14)] = T(1.0);cuU[OpKeyType<1>(7,15)] = T(1.0);
     
  cdDU[OpKeyType<1>(4,0)] = T(1.0);cdDU[OpKeyType<1>(5,1)] = T(1.0);cdDU[OpKeyType<1>(6,2)] = T(1.0);cdDU[OpKeyType<1>(7,3)] = T(1.0);
  cdDU[OpKeyType<1>(12,8)] = T(-1.0);cdDU[OpKeyType<1>(13,9)] = T(-1.0);cdDU[OpKeyType<1>(14,10)] = T(-1.0);cdDU[OpKeyType<1>(15,11)] = T(-1.0);

  cdU[OpKeyType<1>(0,4)] = T(1.0);cdU[OpKeyType<1>(1,5)] = T(1.0);cdU[OpKeyType<1>(2,6)] = T(1.0);cdU[OpKeyType<1>(3,7)] = T(1.0);
  cdU[OpKeyType<1>(8,12)] = T(-1.0);cdU[OpKeyType<1>(9,13)] = T(-1.0);cdU[OpKeyType<1>(10,14)] = T(-1.0);cdU[OpKeyType<1>(11,15)] = T(-1.0);

  SparseLocalOperator<T> bla(16,16);
  
  for(uint ind=0;ind<16;++ind)
    id[OpKeyType<1>(ind,ind)]=T(1);

  //this seems to be right, so DON'T CHANGE ANYTHING HERE
  p[OpKeyType<1>(0,0)]=T(1);p[OpKeyType<1>(1,1)]=T(-1);p[OpKeyType<1>(2,2)]=T(-1);p[OpKeyType<1>(3,3)]=T(1);
  p[OpKeyType<1>(4,4)]=T(-1);p[OpKeyType<1>(5,5)]=T(1);p[OpKeyType<1>(6,6)]=T(1);p[OpKeyType<1>(7,7)]=T(-1);
  p[OpKeyType<1>(8,8)]=T(-1);p[OpKeyType<1>(9,9)]=T(1);p[OpKeyType<1>(10,10)]=T(1);p[OpKeyType<1>(11,11)]=T(-1);
  p[OpKeyType<1>(12,12)]=T(1);p[OpKeyType<1>(13,13)]=T(-1);p[OpKeyType<1>(14,14)]=T(-1);p[OpKeyType<1>(15,15)]=T(1);

  pL[OpKeyType<1>(0,0)]=T(1);pL[OpKeyType<1>(1,1)]=T(-1);pL[OpKeyType<1>(2,2)]=T(-1);pL[OpKeyType<1>(3,3)]=T(1);
  pL[OpKeyType<1>(4,4)]=T(1);pL[OpKeyType<1>(5,5)]=T(-1);pL[OpKeyType<1>(6,6)]=T(-1);pL[OpKeyType<1>(7,7)]=T(1);
  pL[OpKeyType<1>(8,8)]=T(1);pL[OpKeyType<1>(9,9)]=T(-1);pL[OpKeyType<1>(10,10)]=T(-1);pL[OpKeyType<1>(11,11)]=T(1);
  pL[OpKeyType<1>(12,12)]=T(1);pL[OpKeyType<1>(13,13)]=T(-1);pL[OpKeyType<1>(14,14)]=T(-1);pL[OpKeyType<1>(15,15)]=T(1);
 
  pU[OpKeyType<1>(0,0)]=T(1);pU[OpKeyType<1>(1,1)]=T(1);pU[OpKeyType<1>(2,2)]=T(1);pU[OpKeyType<1>(3,3)]=T(1);
  pU[OpKeyType<1>(4,4)]=T(-1);pU[OpKeyType<1>(5,5)]=T(-1);pU[OpKeyType<1>(6,6)]=T(-1);pU[OpKeyType<1>(7,7)]=T(-1);
  pU[OpKeyType<1>(8,8)]=T(1);pU[OpKeyType<1>(9,9)]=T(1);pU[OpKeyType<1>(10,10)]=T(1);pU[OpKeyType<1>(11,11)]=T(1);
  pU[OpKeyType<1>(12,12)]=T(-1);pU[OpKeyType<1>(13,13)]=T(-1);pU[OpKeyType<1>(14,14)]=T(-1);pU[OpKeyType<1>(15,15)]=T(-1);
  
  bla=pL*cuU;cuU=bla;
  bla=pL*cdU;cdU=bla;
  bla=pL*cuDU;cuDU=bla;
  bla=pL*cdDU;cdDU=bla;


  nuL[OpKeyType<1>(2,2)]=T(1);nuL[OpKeyType<1>(3,3)]=T(1);nuL[OpKeyType<1>(6,6)]=T(1);nuL[OpKeyType<1>(7,7)]=T(1);
  nuL[OpKeyType<1>(10,10)]=T(1);nuL[OpKeyType<1>(11,11)]=T(1);nuL[OpKeyType<1>(14,14)]=T(1);nuL[OpKeyType<1>(15,15)]=T(1);

  ndL[OpKeyType<1>(1,1)]=T(1);ndL[OpKeyType<1>(3,3)]=T(1);ndL[OpKeyType<1>(5,5)]=T(1);ndL[OpKeyType<1>(7,7)]=T(1);
  ndL[OpKeyType<1>(9,9)]=T(1);ndL[OpKeyType<1>(11,11)]=T(1);ndL[OpKeyType<1>(13,13)]=T(1);ndL[OpKeyType<1>(15,15)]=T(1);
  nL=nuL+ndL;
  
  for(uint ind=8;ind<16;++ind)
    nuU[OpKeyType<1>(ind,ind)]=T(1);
  for(uint ind=4;ind<8;++ind)
    ndU[OpKeyType<1>(ind,ind)]=T(1);
  for(uint ind=12;ind<16;++ind)
    ndU[OpKeyType<1>(ind,ind)]=T(1);
  nU=nuU+ndU;

  nuU.setSite(0);ndU.setSite(0);cuU.setSite(0);cdU.setSite(0);cuDU.setSite(0);cdDU.setSite(0);id.setSite(0);p.setSite(0);
  nuL.setSite(0);ndL.setSite(0);cuL.setSite(0);cdL.setSite(0);cuDL.setSite(0);cdDL.setSite(0);pU.setSite(0);pL.setSite(0);
  nL.setSite(0);
  nU.setSite(0);
  // M[ukey<2>(0,0)] = chUpU(0)*nuU+chDownU(0)*ndU+UU(0)*(nuU*ndU)+chUpL(0)*nuL+chDownL(0)*ndL+UL(0)*(nuL*ndL)-
  //   hopInt(0)*((pL*cuDL)*cuU-cuDU*(pL*cuL)+(pL*cdDL)*cdU-cdDU*(pL*cdL))+Vint(0)*((nL-id)*(nU-id)) ;
  M[ukey<2>(0,0)] = chUpU(0)*nuU+chDownU(0)*ndU+U1(0)*(nuU*ndU)+chUpL(0)*nuL+chDownL(0)*ndL+U1(0)*(nuL*ndL)+U2(0)*(nuL*ndU+nuU*ndL)+U3(0)*(nuL*nuU+ndL*ndU)+
    J1(0)*(cuDL*cdL*cdDU*cuU+cuDU*cdU*cdDL*cuL)+J2(0)*(cuDL*cdDL*cdU*cuU+cuDU*cdDU*cdL*cuL);

  M[ukey<2>(0,1)] = -hopL(0)*cuL*p;
  M[ukey<2>(0,2)] = hopL(0)*cuDL*p;//
  M[ukey<2>(0,3)] = -hopU(0)*cuU*p;
  M[ukey<2>(0,4)] = hopU(0)*cuDU*p;//
  M[ukey<2>(0,5)] = -hopL(0)*cdL*p;
  M[ukey<2>(0,6)] = hopL(0)*cdDL*p;//
  M[ukey<2>(0,7)] = -hopU(0)*cdU*p;
  M[ukey<2>(0,8)] = hopU(0)*cdDU*p;//
  M[ukey<2>(0,9)] = id;
  M.Dl=1,M.Dr=10;
  M.squeeze(1e-15);
  (*this)[0]=M;
  M.clear();

  nuU.setSite(N-1);ndU.setSite(N-1);cuU.setSite(N-1);cdU.setSite(N-1);cuDU.setSite(N-1);cdDU.setSite(N-1);id.setSite(N-1);p.setSite(N-1);
  nuL.setSite(N-1);ndL.setSite(N-1);cuL.setSite(N-1);cdL.setSite(N-1);cuDL.setSite(N-1);cdDL.setSite(N-1);pU.setSite(N-1);pL.setSite(N-1);
  nL.setSite(N-1);
  nU.setSite(N-1);

  M[ukey<2>(0,0)] = id;
  M[ukey<2>(1,0)] = cuDL;
  M[ukey<2>(2,0)] = cuL;
  M[ukey<2>(3,0)] = cuDU;
  M[ukey<2>(4,0)] = cuU;
  M[ukey<2>(5,0)] = cdDL;
  M[ukey<2>(6,0)] = cdL;
  M[ukey<2>(7,0)] = cdDU;
  M[ukey<2>(8,0)] = cdU;
  M[ukey<2>(9,0)] = chUpU(N-1)*nuU+chDownU(N-1)*ndU+U1(N-1)*(nuU*ndU)+chUpL(N-1)*nuL+chDownL(N-1)*ndL+U1(N-1)*(nuL*ndL)+U2(N-1)*(nuL*ndU+nuU*ndL)+U3(N-1)*(nuL*nuU+ndL*ndU)+
    J1(N-1)*(cuDL*cdL*cdDU*cuU+cuDU*cdU*cdDL*cuL)+J2(N-1)*(cuDL*cdDL*cdU*cuU+cuDU*cdDU*cdL*cuL);


  M.Dl=10,M.Dr=1;
  M.squeeze(1e-15);
  (*this)[N-1]=M;

  for(uint site = 1;site<(N-1);++site){
    MPOMat<1,T> Mb;
    nuU.setSite(site);ndU.setSite(site);cuU.setSite(site);cdU.setSite(site);cuDU.setSite(site);cdDU.setSite(site);id.setSite(site);p.setSite(site);
    nuL.setSite(site);ndL.setSite(site);cuL.setSite(site);cdL.setSite(site);cuDL.setSite(site);cdDL.setSite(site);pU.setSite(site);pL.setSite(site);
    nL.setSite(site);
    nU.setSite(site);

    Mb[ukey<2>(0,0)] = id;
    Mb[ukey<2>(1,0)] = cuDL;
    Mb[ukey<2>(2,0)] = cuL;
    Mb[ukey<2>(3,0)] = cuDU;
    Mb[ukey<2>(4,0)] = cuU;
    Mb[ukey<2>(5,0)] = cdDL;
    Mb[ukey<2>(6,0)] = cdL;
    Mb[ukey<2>(7,0)] = cdDU;
    Mb[ukey<2>(8,0)] = cdU;
    Mb[ukey<2>(9,0)] = chUpU(site)*nuU+chDownU(site)*ndU+U1(site)*(nuU*ndU)+chUpL(site)*nuL+chDownL(site)*ndL+U1(site)*(nuL*ndL)+U2(site)*(nuL*ndU+nuU*ndL)+U3(site)*(nuL*nuU+ndL*ndU)+
      J1(site)*(cuDL*cdL*cdDU*cuU+cuDU*cdU*cdDL*cuL)+J2(site)*(cuDL*cdDL*cdU*cuU+cuDU*cdDU*cdL*cuL);
     
    Mb[ukey<2>(9,1)] = -hopL(site)*cuL*p;
    Mb[ukey<2>(9,2)] = hopL(site)*cuDL*p;//
    Mb[ukey<2>(9,3)] = -hopU(site)*cuU*p;
    Mb[ukey<2>(9,4)] = hopU(site)*cuDU*p;//
    Mb[ukey<2>(9,5)] = -hopL(site)*cdL*p;
    Mb[ukey<2>(9,6)] = hopL(site)*cdDL*p;//
    Mb[ukey<2>(9,7)] = -hopU(site)*cdU*p;
    Mb[ukey<2>(9,8)] = hopU(site)*cdDU*p;//
    Mb[ukey<2>(9,9)] = id;
    Mb.Dl=10,Mb.Dr=10;
    Mb.squeeze(1e-15);
    (*this)[site]=Mb;
  }
};


template<typename T>
class SpinlessFermNN:public MPO<T> {
public:
	SpinlessFermNN(VectorType t, VectorType eps, double U, double W){

		uint N=eps.size();

		MPOMat<1,T> M;
		//this->resize(N);
		SparseLocalOperator<T> c(2,2),cD(2,2),n(2,2),id(2,2), p(2,2);
		c[OpKeyType<1>(0,1)]=T(1.);
		cD[OpKeyType<1>(1,0)]=T(1.);
		n[OpKeyType<1>(1,1)]=T(1.);
		id[OpKeyType<1>(0,0)]=T(1.);id[OpKeyType<1>(1,1)]=T(1.);
		p[OpKeyType<1>(0,0)]=T(1.);p[OpKeyType<1>(1,1)]=T(-1.);

		c.setSite(0); cD.setSite(0); n.setSite(0); id.setSite(0); p.setSite(0);
		M[ukey<2>(0,0)] = eps[0]*n;
		M[ukey<2>(0,1)] = -t[0]*cD*p;
		M[ukey<2>(0,2)] = t[0]*c*p;//-
		M[ukey<2>(0,3)] = U*n;
		M[ukey<2>(0,4)] = W*n;
		M[ukey<2>(0,5)] = id;
		M.Dl=1,M.Dr=6;
		M.squeeze(1e-16);
		(*this)[0]=M;
		M.clear();

		c.setSite(N-1); cD.setSite(N-1); n.setSite(N-1); id.setSite(N-1);
		M[ukey<2>(0,0)] = id;
		M[ukey<2>(1,0)] = c;
		M[ukey<2>(2,0)] = cD;
		M[ukey<2>(3,0)] = n;
		M[ukey<2>(5,0)] = eps[N-1]*n;
		M.Dl=6, M.Dr=1;
		M.squeeze(1e-16);
		(*this)[N-1]=M;
		M.clear();

		// intermediate sites
		for (uint site = 1; site < N-1; site++) {
			MPOMat<1,T> Mb;
			c.setSite(site); cD.setSite(site); n.setSite(site); id.setSite(site);p.setSite(site);
			Mb[ukey<2>(0,0)] = id;
			Mb[ukey<2>(1,0)] = c;
			Mb[ukey<2>(2,0)] = cD;
			Mb[ukey<2>(3,0)] = n;
			Mb[ukey<2>(4,3)] = id;
			Mb[ukey<2>(5,0)] = eps[site]*n;
			Mb[ukey<2>(5,1)] = -t[site]*cD*p;
			Mb[ukey<2>(5,2)] = t[site]*c*p;//-
			Mb[ukey<2>(5,3)] = U*n;
			Mb[ukey<2>(5,4)] = W*n;
			Mb[ukey<2>(5,5)] = id;

			Mb.Dl=6,Mb.Dr=6;

		    Mb.squeeze(1e-16);
		    (*this)[site]=Mb;
		}

	}

	~SpinlessFermNN(){}

};


template<typename T>
class SpinlessFermNNImp:public MPO<T> {
public:
  SpinlessFermNNImp(VectorType t, VectorType eps, double V, double U, double W, double delta, double omega, uint impSite){

    std::cout << "IMPURITY SITE="<<impSite << std::endl;

    uint N=eps.size();

    MPOMat<1,T> M;
    //this->resize(N);
    SparseLocalOperator<T> c(2,2),cD(2,2),n(2,2),id(2,2), p(2,2);
    SparseLocalOperator<T> cI(4,4),cDI(4,4),nI(4,4),idI(4,4),szI(4,4),sxI(4,4), pI(4,4);

    c[OpKeyType<1>(0,1)]=T(1.);
    cD[OpKeyType<1>(1,0)]=T(1.);
    n[OpKeyType<1>(1,1)]=T(1.);
    id[OpKeyType<1>(0,0)]=T(1.);id[OpKeyType<1>(1,1)]=T(1.);
    p[OpKeyType<1>(0,0)]=T(1.);p[OpKeyType<1>(1,1)]=T(-1.);

    cI[OpKeyType<1>(0,2)]=T(1.); cI[OpKeyType<1>(1,3)]=T(1.);
    cDI[OpKeyType<1>(2,0)]=T(1.); cDI[OpKeyType<1>(3,1)]=T(1.);
    nI[OpKeyType<1>(2,2)]=T(1.); nI[OpKeyType<1>(3,3)]=T(1.);
    idI[OpKeyType<1>(0,0)]=T(1.); idI[OpKeyType<1>(1,1)]=T(1.); idI[OpKeyType<1>(2,2)]=T(1.); idI[OpKeyType<1>(3,3)]=T(1.);
    szI[OpKeyType<1>(0,0)]=T(-0.5); szI[OpKeyType<1>(1,1)]=T(0.5); szI[OpKeyType<1>(2,2)]=T(-0.5); szI[OpKeyType<1>(3,3)]=T(0.5);
    sxI[OpKeyType<1>(0,1)]=T(0.5); sxI[OpKeyType<1>(1,0)]=T(0.5); sxI[OpKeyType<1>(2,3)]=T(0.5); sxI[OpKeyType<1>(3,2)]=T(0.5);
    pI[OpKeyType<1>(0,0)]=T(1.);pI[OpKeyType<1>(1,1)]=T(1.); pI[OpKeyType<1>(2,2)]=T(-1.);pI[OpKeyType<1>(3,3)]=T(-1.);

    if (impSite!=0){
      c.setSite(0); cD.setSite(0); n.setSite(0); id.setSite(0); p.setSite(0);
      M[ukey<2>(0,0)] = eps[0]*n;
      M[ukey<2>(0,1)] = -t[0]*cD*p;
      M[ukey<2>(0,2)] = t[0]*c*p;//-
      M[ukey<2>(0,3)] = U*n;
      M[ukey<2>(0,4)] = W*n;
      M[ukey<2>(0,5)] = id;
    }else{
      cI.setSite(0); cDI.setSite(0); nI.setSite(0); idI.setSite(0); szI.setSite(0); sxI.setSite(0); pI.setSite(0);
      M[ukey<2>(0,0)] = (eps[0]*nI)+(omega*sxI)+(delta*szI)+(V*szI*nI);
      M[ukey<2>(0,1)] = -t[0]*cDI*pI;
      M[ukey<2>(0,2)] = t[0]*cI*pI;//-
      M[ukey<2>(0,3)] = U*nI;
      M[ukey<2>(0,4)] = W*nI;
      M[ukey<2>(0,5)] = idI;
    }
    M.Dl=1,M.Dr=6;
    M.squeeze(1e-16);
    (*this)[0]=M;
    M.clear();

    if (impSite!=N-1){
      c.setSite(N-1); cD.setSite(N-1); n.setSite(N-1); id.setSite(N-1);
      M[ukey<2>(0,0)] = id;
      M[ukey<2>(1,0)] = c;
      M[ukey<2>(2,0)] = cD;
      M[ukey<2>(3,0)] = n;
      M[ukey<2>(5,0)] = eps[N-1]*n;
    }else{
      cI.setSite(N-1); cDI.setSite(N-1); nI.setSite(N-1); idI.setSite(N-1); szI.setSite(N-1); sxI.setSite(N-1);
      M[ukey<2>(0,0)] = idI;
      M[ukey<2>(1,0)] = cI;
      M[ukey<2>(2,0)] = cDI;
      M[ukey<2>(3,0)] = nI;
      M[ukey<2>(5,0)] = (eps[N-1]*nI)+(omega*sxI)+(delta*szI)+(V*szI*nI);
    }
    M.Dl=6, M.Dr=1;
    M.squeeze(1e-16);
    (*this)[N-1]=M;
    M.clear();

    // intermediate sites
    for (uint site = 1; site < N-1; site++) {
      MPOMat<1,T> Mb;
      if(site!=impSite){ // normal site, local Hilbert space=2
	c.setSite(site); cD.setSite(site); n.setSite(site); id.setSite(site);p.setSite(site);
	Mb[ukey<2>(0,0)] = id;
	Mb[ukey<2>(1,0)] = c;
	Mb[ukey<2>(2,0)] = cD;
	Mb[ukey<2>(3,0)] = n;
	Mb[ukey<2>(4,3)] = id;
	Mb[ukey<2>(5,0)] = eps[site]*n;
	Mb[ukey<2>(5,1)] = -t[site]*cD*p;
	Mb[ukey<2>(5,2)] = t[site]*c*p;//-
	Mb[ukey<2>(5,3)] = U*n;
	Mb[ukey<2>(5,4)] = W*n;
	Mb[ukey<2>(5,5)] = id;
      }else{ // impurity site, local Hilbert space =4
	cI.setSite(site); cDI.setSite(site); nI.setSite(site); idI.setSite(site); szI.setSite(site); sxI.setSite(site);pI.setSite(site);
	Mb[ukey<2>(0,0)] = idI;
	Mb[ukey<2>(1,0)] = cI;
	Mb[ukey<2>(2,0)] = cDI;
	Mb[ukey<2>(3,0)] = nI;
	Mb[ukey<2>(4,3)] = idI;
	Mb[ukey<2>(5,0)] = (eps[site]*nI)+(omega*sxI)+(delta*szI)+V*(szI*nI);
	Mb[ukey<2>(5,1)] = -t[site]*cDI*pI;
	Mb[ukey<2>(5,2)] = t[site]*cI*pI;//-
	Mb[ukey<2>(5,3)] = U*nI;
	Mb[ukey<2>(5,4)] = W*nI;
	Mb[ukey<2>(5,5)] = idI;
      }

      Mb.Dl=6,Mb.Dr=6;

      Mb.squeeze(1e-16);
      (*this)[site]=Mb;
    }

  }

  ~SpinlessFermNNImp(){}

};


template<typename T>
class CondImp:public MPO<T> {
public:
	CondImp(uint N, VectorType t, double delta, double omega, double V, uint impSite, VectorType eps){

		std::cout << "IMPURITY SITE="<<impSite << std::endl;

		MPOMat<1,T> M;
		//this->resize(N);
		SparseLocalOperator<T> c(2,2),cD(2,2),n(2,2),id(2,2), p(2,2);//all statistics are set to bosonic (_statistics=1)
		SparseLocalOperator<T> cI(4,4),cDI(4,4),nI(4,4),idI(4,4),szI(4,4),sxI(4,4), pI(4,4);
		c[OpKeyType<1>(0,1)]=T(1.);
		cD[OpKeyType<1>(1,0)]=T(1.);
		n[OpKeyType<1>(1,1)]=T(1.);
		id[OpKeyType<1>(0,0)]=T(1.);id[OpKeyType<1>(1,1)]=T(1.);
		p[OpKeyType<1>(0,0)]=T(1.);p[OpKeyType<1>(1,1)]=T(-1.);

		cI[OpKeyType<1>(0,2)]=T(1.); cI[OpKeyType<1>(1,3)]=T(1.);
		cDI[OpKeyType<1>(2,0)]=T(1.); cDI[OpKeyType<1>(3,1)]=T(1.);
		nI[OpKeyType<1>(2,2)]=T(1.); nI[OpKeyType<1>(3,3)]=T(1.);
		idI[OpKeyType<1>(0,0)]=T(1.); idI[OpKeyType<1>(1,1)]=T(1.); idI[OpKeyType<1>(2,2)]=T(1.); idI[OpKeyType<1>(3,3)]=T(1.);
		szI[OpKeyType<1>(0,0)]=T(-0.5); szI[OpKeyType<1>(1,1)]=T(0.5); szI[OpKeyType<1>(2,2)]=T(-0.5); szI[OpKeyType<1>(3,3)]=T(0.5);
		sxI[OpKeyType<1>(0,1)]=T(0.5); sxI[OpKeyType<1>(1,0)]=T(0.5); sxI[OpKeyType<1>(2,3)]=T(0.5); sxI[OpKeyType<1>(3,2)]=T(0.5);
		pI[OpKeyType<1>(0,0)]=T(1.);pI[OpKeyType<1>(1,1)]=T(1.); pI[OpKeyType<1>(2,2)]=T(-1.);pI[OpKeyType<1>(3,3)]=T(-1.);

		if (impSite!=0){
			c.setSite(0); cD.setSite(0); n.setSite(0); id.setSite(0); p.setSite(0);
			M[ukey<2>(0,0)] = eps[0]*n;
			M[ukey<2>(0,1)] = -t[0]*cD*p;
			M[ukey<2>(0,2)] = t[0]*c*p;//-
			M[ukey<2>(0,3)] = id;
			M.Dl=1,M.Dr=4;
		}else{
			cI.setSite(0); cDI.setSite(0); nI.setSite(0); idI.setSite(0); szI.setSite(0); sxI.setSite(0); pI.setSite(0);
			M[ukey<2>(0,0)] = (eps[0]*nI)+(omega*sxI)+(delta*szI)+(V*szI*nI);
			M[ukey<2>(0,1)] = -t[0]*cDI*pI;
			M[ukey<2>(0,2)] = t[0]*cI*pI;//-
			M[ukey<2>(0,3)] = idI;
			M.Dl=1,M.Dr=4;
		}
		M.squeeze(1e-16);
		(*this)[0]=M;
		M.clear();
		if (impSite!=N-1){
			c.setSite(N-1); cD.setSite(N-1); n.setSite(N-1); id.setSite(N-1);
			M[ukey<2>(0,0)] = id;
			M[ukey<2>(1,0)] = c;
			M[ukey<2>(2,0)] = cD;
			M[ukey<2>(3,0)] = eps[N-1]*n;
			M.Dl=4,M.Dr=1;
		}else{
			cI.setSite(N-1); cDI.setSite(N-1); nI.setSite(N-1); idI.setSite(N-1); szI.setSite(N-1); sxI.setSite(N-1);
			M[ukey<2>(0,0)] = idI;
			M[ukey<2>(1,0)] = cI;
			M[ukey<2>(2,0)] = cDI;
			M[ukey<2>(3,0)] = (eps[N-1]*nI)+(omega*sxI)+(delta*szI)+(V*szI*nI);
			M.Dl=4,M.Dr=1;
		}
		M.squeeze(1e-16);
		(*this)[N-1]=M;
		M.clear();

		// intermediate sites
		for (uint site = 1; site < N-1; site++) {
			MPOMat<1,T> Mb;
	    	if(site!=impSite){ // normal site, local Hilbert space=2
		  c.setSite(site); cD.setSite(site); n.setSite(site); id.setSite(site);p.setSite(site);
				Mb[ukey<2>(0,0)] = id;
				Mb[ukey<2>(1,0)] = c;
				Mb[ukey<2>(2,0)] = cD;
				Mb[ukey<2>(3,0)] = eps[site]*n;
				Mb[ukey<2>(3,1)] = -t[site]*cD*p;
				Mb[ukey<2>(3,2)] = t[site]*c*p;//-
				Mb[ukey<2>(3,3)] = id;
			}else{ // impurity site, local Hilbert space =4
		  cI.setSite(site); cDI.setSite(site); nI.setSite(site); idI.setSite(site); szI.setSite(site); sxI.setSite(site);pI.setSite(site);
				Mb[ukey<2>(0,0)] = idI;
				Mb[ukey<2>(1,0)] = cI;
				Mb[ukey<2>(2,0)] = cDI;
				Mb[ukey<2>(3,0)] = (eps[site]*nI)+(omega*sxI)+(delta*szI)+V*(szI*nI);
				Mb[ukey<2>(3,1)] = -t[site]*cDI*pI;
				Mb[ukey<2>(3,2)] = t[site]*cI*pI;//-
				Mb[ukey<2>(3,3)] = idI;
			}
			Mb.Dl=4,Mb.Dr=4;

		    Mb.squeeze(1e-15);
		    (*this)[site]=Mb;
		}

	}

	~CondImp(){}

};



template<typename T>
class CondImpInt:public MPO<T> {
public:
	CondImpInt(uint N, VectorType t, double delta, double omega, double V, double U, uint impSite, VectorType eps){

		std::cout << "IMPURITY SITE="<<impSite << std::endl;

		MPOMat<1,T> M;
		//this->resize(N);
		SparseLocalOperator<T> c(2,2),cD(2,2),n(2,2),id(2,2), p(2,2);//all statistics are set to bosonic (_statistics=1)
		SparseLocalOperator<T> cI(4,4),cDI(4,4),nI(4,4),idI(4,4),szI(4,4),sxI(4,4), pI(4,4);
		c[OpKeyType<1>(0,1)]=T(1.);
		cD[OpKeyType<1>(1,0)]=T(1.);
		n[OpKeyType<1>(1,1)]=T(1.);
		id[OpKeyType<1>(0,0)]=T(1.);id[OpKeyType<1>(1,1)]=T(1.);
		p[OpKeyType<1>(0,0)]=T(1.);p[OpKeyType<1>(1,1)]=T(-1.);

		cI[OpKeyType<1>(0,2)]=T(1.); cI[OpKeyType<1>(1,3)]=T(1.);
		cDI[OpKeyType<1>(2,0)]=T(1.); cDI[OpKeyType<1>(3,1)]=T(1.);
		nI[OpKeyType<1>(2,2)]=T(1.); nI[OpKeyType<1>(3,3)]=T(1.);
		idI[OpKeyType<1>(0,0)]=T(1.); idI[OpKeyType<1>(1,1)]=T(1.); idI[OpKeyType<1>(2,2)]=T(1.); idI[OpKeyType<1>(3,3)]=T(1.);
		szI[OpKeyType<1>(0,0)]=T(-0.5); szI[OpKeyType<1>(1,1)]=T(0.5); szI[OpKeyType<1>(2,2)]=T(-0.5); szI[OpKeyType<1>(3,3)]=T(0.5);
		sxI[OpKeyType<1>(0,1)]=T(0.5); sxI[OpKeyType<1>(1,0)]=T(0.5); sxI[OpKeyType<1>(2,3)]=T(0.5); sxI[OpKeyType<1>(3,2)]=T(0.5);
		pI[OpKeyType<1>(0,0)]=T(1.);pI[OpKeyType<1>(1,1)]=T(1.); pI[OpKeyType<1>(2,2)]=T(-1.);pI[OpKeyType<1>(3,3)]=T(-1.);

		if (impSite!=0){
			c.setSite(0); cD.setSite(0); n.setSite(0); id.setSite(0); p.setSite(0);
			M[ukey<2>(0,0)] = eps[0]*n;
			M[ukey<2>(0,1)] = -t[0]*cD*p;
			M[ukey<2>(0,2)] = t[0]*c*p;//-
			M[ukey<2>(0,3)] = U*n;
			M[ukey<2>(0,4)] = id;
			M.Dl=1,M.Dr=5;
		}else{
			cI.setSite(0); cDI.setSite(0); nI.setSite(0); idI.setSite(0); szI.setSite(0); sxI.setSite(0); pI.setSite(0);
			M[ukey<2>(0,0)] = (eps[0]*nI)+(omega*sxI)+(delta*szI)+(V*szI*nI);
			M[ukey<2>(0,1)] = -t[0]*cDI*pI;
			M[ukey<2>(0,2)] = t[0]*cI*pI;//-
			M[ukey<2>(0,3)] = U*nI;
			M[ukey<2>(0,4)] = idI;
			M.Dl=1,M.Dr=5;
		}
		M.squeeze(1e-16);
		(*this)[0]=M;
		M.clear();
		if (impSite!=N-1){
			c.setSite(N-1); cD.setSite(N-1); n.setSite(N-1); id.setSite(N-1);
			M[ukey<2>(0,0)] = id;
			M[ukey<2>(1,0)] = c;
			M[ukey<2>(2,0)] = cD;
			M[ukey<2>(3,0)] = n;
			M[ukey<2>(4,0)] = eps[N-1]*n;
			M.Dl=5,M.Dr=1;
		}else{
			cI.setSite(N-1); cDI.setSite(N-1); nI.setSite(N-1); idI.setSite(N-1); szI.setSite(N-1); sxI.setSite(N-1);
			M[ukey<2>(0,0)] = idI;
			M[ukey<2>(1,0)] = cI;
			M[ukey<2>(2,0)] = cDI;
			M[ukey<2>(3,0)] = nI;
			M[ukey<2>(4,0)] = (eps[N-1]*nI)+(omega*sxI)+(delta*szI)+(V*szI*nI);
			M.Dl=5,M.Dr=1;
		}
		M.squeeze(1e-16);
		(*this)[N-1]=M;
		M.clear();

		// intermediate sites
		for (uint site = 1; site < N-1; site++) {
			MPOMat<1,T> Mb;
	    	if(site!=impSite){ // normal site, local Hilbert space=2
	    		c.setSite(site); cD.setSite(site); n.setSite(site); id.setSite(site);p.setSite(site);
				Mb[ukey<2>(0,0)] = id;
				Mb[ukey<2>(1,0)] = c;
				Mb[ukey<2>(2,0)] = cD;
				Mb[ukey<2>(3,0)] = n;
				Mb[ukey<2>(4,0)] = eps[site]*n;
				Mb[ukey<2>(4,1)] = -t[site]*cD*p;
				Mb[ukey<2>(4,2)] = t[site]*c*p;//-
				Mb[ukey<2>(4,3)] = U*n;
				Mb[ukey<2>(4,4)] = id;
			}else{ // impurity site, local Hilbert space =4
				cI.setSite(site); cDI.setSite(site); nI.setSite(site); idI.setSite(site); szI.setSite(site); sxI.setSite(site);pI.setSite(site);
				Mb[ukey<2>(0,0)] = idI;
				Mb[ukey<2>(1,0)] = cI;
				Mb[ukey<2>(2,0)] = cDI;
				Mb[ukey<2>(3,0)] = nI;
				Mb[ukey<2>(4,0)] = (eps[site]*nI)+(omega*sxI)+(delta*szI)+V*(szI*nI);
				Mb[ukey<2>(4,1)] = -t[site]*cDI*pI;
				Mb[ukey<2>(4,2)] = t[site]*cI*pI;//-
				Mb[ukey<2>(4,3)] = U*nI;
				Mb[ukey<2>(4,4)] = idI;
			}
			Mb.Dl=5,Mb.Dr=5;

		    Mb.squeeze(1e-15);
		    (*this)[site]=Mb;
		}
	}

	~CondImpInt(){}

};


template<typename T>
class CondMvImp:public MPO<T> {
public:
	CondMvImp(uint N, VectorType J, VectorType eps, VectorType epsImp, double K, double omega, double delta, double V){

		assert(eps.size() == N);
		assert(J.size() == N-1);

		MPOMat<1,T> M;
		//this->resize(N);
		SparseLocalOperator<T> c(6,6),cD(6,6),fu(6,6),fuD(6,6),fd(6,6),fdD(6,6),n(6,6),id(6,6),sz(6,6),sx(6,6),nI(6,6),p(6,6);

		c[OpKeyType<1>(0,3)]=T(1.); c[OpKeyType<1>(1,4)]=T(1.); c[OpKeyType<1>(2,5)]=T(1.);
		cD[OpKeyType<1>(3,0)]=T(1.); cD[OpKeyType<1>(4,1)]=T(1.); cD[OpKeyType<1>(5,2)]=T(1.);
		n[OpKeyType<1>(3,3)]=T(1.); n[OpKeyType<1>(4,4)]=T(1.); n[OpKeyType<1>(5,5)]=T(1.);

		fd[OpKeyType<1>(0,1)]=T(1.); fd[OpKeyType<1>(3,4)]=T(1.);
		fu[OpKeyType<1>(0,2)]=T(1.); fu[OpKeyType<1>(3,5)]=T(1.);
		fdD[OpKeyType<1>(1,0)]=T(1.); fdD[OpKeyType<1>(4,3)]=T(1.);
		fuD[OpKeyType<1>(2,0)]=T(1.); fuD[OpKeyType<1>(5,3)]=T(1.);

		id[OpKeyType<1>(0,0)]=T(1.); id[OpKeyType<1>(1,1)]=T(1.); id[OpKeyType<1>(2,2)]=T(1.);
		id[OpKeyType<1>(3,3)]=T(1.); id[OpKeyType<1>(4,4)]=T(1.); id[OpKeyType<1>(5,5)]=T(1.);
		sz[OpKeyType<1>(1,1)]=T(-0.5); sz[OpKeyType<1>(2,2)]=T(0.5); sz[OpKeyType<1>(4,4)]=T(-0.5); sz[OpKeyType<1>(5,5)]=T(0.5);
		sx[OpKeyType<1>(1,2)]=T(0.5); sx[OpKeyType<1>(2,1)]=T(0.5); sx[OpKeyType<1>(4,5)]=T(0.5); sx[OpKeyType<1>(5,4)]=T(0.5);
		nI[OpKeyType<1>(1,1)]=T(1); nI[OpKeyType<1>(2,2)]=T(1); nI[OpKeyType<1>(4,4)]=T(1); nI[OpKeyType<1>(5,5)]=T(1);

		p[OpKeyType<1>(0,0)]=T(1.); p[OpKeyType<1>(1,1)]=T(1.); p[OpKeyType<1>(2,2)]=T(1.);
		p[OpKeyType<1>(3,3)]=T(-1.); p[OpKeyType<1>(4,4)]=T(-1.); p[OpKeyType<1>(5,5)]=T(-1.);

		c.setSite(0); cD.setSite(0); fu.setSite(0); fuD.setSite(0); fd.setSite(0); fdD.setSite(0);
		n.setSite(0); id.setSite(0); sz.setSite(0); sx.setSite(0); p.setSite(0); nI.setSite(0);
		M[ukey<2>(0,0)] = (eps[0]*n)+(epsImp[0]*nI)+(omega*sx)+(delta*sz)+(V*sz*n);
		M[ukey<2>(0,1)] = -J[0]*cD*p;
		M[ukey<2>(0,2)] = J[0]*c*p;
		M[ukey<2>(0,3)] = -K*fuD;
		M[ukey<2>(0,4)] = -K*fu;
		M[ukey<2>(0,5)] = -K*fdD;
		M[ukey<2>(0,6)] = -K*fd;
		M[ukey<2>(0,7)] = id;
		M.Dl=1,M.Dr=8;


		M.squeeze(1e-16);
		(*this)[0]=M;
		M.clear();

		c.setSite(N-1); cD.setSite(N-1); fu.setSite(N-1); fuD.setSite(N-1); fd.setSite(N-1); fdD.setSite(N-1);
		n.setSite(N-1); id.setSite(N-1); sz.setSite(N-1); sx.setSite(N-1); nI.setSite(N-1);
		M[ukey<2>(0,0)] = id;
		M[ukey<2>(1,0)] = c;
		M[ukey<2>(2,0)] = cD;
		M[ukey<2>(3,0)] = fu;
		M[ukey<2>(4,0)] = fuD;
		M[ukey<2>(5,0)] = fd;
		M[ukey<2>(6,0)] = fdD;
		M[ukey<2>(7,0)] = (eps[N-1]*n)+(epsImp[N-1]*nI)+(omega*sx)+(delta*sz)+(V*sz*n);
		M.Dl=8,M.Dr=1;

		M.squeeze(1e-16);
		(*this)[N-1]=M;
		M.clear();

		// intermediate sites
		for (uint site = 1; site < N-1; site++) {
			MPOMat<1,T> Mb;
			c.setSite(site); cD.setSite(site); fu.setSite(site); fuD.setSite(site); fd.setSite(site); fdD.setSite(site);
			n.setSite(site); id.setSite(site); sz.setSite(site); sx.setSite(site); p.setSite(site); nI.setSite(site);
			Mb[ukey<2>(0,0)] = id;
			Mb[ukey<2>(1,0)] = c;
			Mb[ukey<2>(2,0)] = cD;
			Mb[ukey<2>(3,0)] = fu;
			Mb[ukey<2>(4,0)] = fuD;
			Mb[ukey<2>(5,0)] = fd;
			Mb[ukey<2>(6,0)] = fdD;
			Mb[ukey<2>(7,0)] = (eps[site]*n)+(epsImp[site]*nI)+(omega*sx)+(delta*sz)+(V*sz*n);
			Mb[ukey<2>(7,1)] = -J[site]*cD*p;
			Mb[ukey<2>(7,2)] = J[site]*c*p;
			Mb[ukey<2>(7,3)] = -K*fuD;
			Mb[ukey<2>(7,4)] = -K*fu;
			Mb[ukey<2>(7,5)] = -K*fdD;
			Mb[ukey<2>(7,6)] = -K*fd;
			Mb[ukey<2>(7,7)] = id;

			Mb.Dl=8,Mb.Dr=8;
		    Mb.squeeze(1e-15);
		    (*this)[site]=Mb;
		}
	}

	~CondMvImp(){}

};



template<typename T>
class HubbImp:public MPO<T> {
public:
	HubbImp(uint N, VectorType t, VectorType eps, double delta, double omega, double V, double U, uint impSite){

		std::cout << "IMPURITY SITE="<<impSite << std::endl;

		assert(eps.size() == this->num_sites_);
		assert(t.size() == this->num_sites_ - 1);
		assert(impSite>=0&&impSite<this->num_sites_-1);

		MPOMat<1,T> M;
		//this->resize(N);
		SparseLocalOperator<T> c1(4,4),c1D(4,4),c2(4,4),c2D(4,4),n1(4,4),n2(4,4),id(4,4), p(4,4);//all statistics are set to bosonic (_statistics=1)
		SparseLocalOperator<T> c1I(8,8),c1DI(8,8),c2I(8,8),c2DI(8,8),n1I(8,8),n2I(8,8),idI(8,8),szI(8,8),sxI(8,8), pI(8,8);
		c1[OpKeyType<1>(0,1)]=T(1.); c1[OpKeyType<1>(2,3)]=T(1.);
		c2[OpKeyType<1>(0,2)]=T(1.); c2[OpKeyType<1>(1,3)]=T(-1.);
		c1D[OpKeyType<1>(1,0)]=T(1.); c1D[OpKeyType<1>(3,2)]=T(1.);
		c2D[OpKeyType<1>(2,0)]=T(1.); c2D[OpKeyType<1>(3,1)]=T(-1.);
		n1[OpKeyType<1>(1,1)]=T(1.); n1[OpKeyType<1>(3,3)]=T(1.);
		n2[OpKeyType<1>(2,2)]=T(1.); n2[OpKeyType<1>(3,3)]=T(1.);
		id[OpKeyType<1>(0,0)]=T(1.);id[OpKeyType<1>(1,1)]=T(1.); id[OpKeyType<1>(2,2)]=T(1.);id[OpKeyType<1>(3,3)]=T(1.);
		p[OpKeyType<1>(0,0)]=T(1.); p[OpKeyType<1>(1,1)]=T(-1.); p[OpKeyType<1>(2,2)]=T(-1.); p[OpKeyType<1>(3,3)]=T(1.);

		c1I[OpKeyType<1>(0,1)]=T(1.); c1I[OpKeyType<1>(2,3)]=T(1.); c1I[OpKeyType<1>(4,5)]=T(1.); c1I[OpKeyType<1>(6,7)]=T(1.);
		c2I[OpKeyType<1>(0,2)]=T(1.); c2I[OpKeyType<1>(1,3)]=T(-1.); c2I[OpKeyType<1>(4,6)]=T(1.); c2I[OpKeyType<1>(5,7)]=T(-1.);
		c1DI[OpKeyType<1>(1,0)]=T(1.); c1DI[OpKeyType<1>(3,2)]=T(1.); c1DI[OpKeyType<1>(5,4)]=T(1.); c1DI[OpKeyType<1>(7,6)]=T(1.);
		c2DI[OpKeyType<1>(2,0)]=T(1.); c2DI[OpKeyType<1>(3,1)]=T(-1.); c2DI[OpKeyType<1>(6,4)]=T(1.); c2DI[OpKeyType<1>(7,5)]=T(-1.);
		n1I[OpKeyType<1>(1,1)]=T(1.); n1I[OpKeyType<1>(3,3)]=T(1.); n1I[OpKeyType<1>(5,5)]=T(1.); n1I[OpKeyType<1>(7,7)]=T(1.);
		n2I[OpKeyType<1>(2,2)]=T(1.); n2I[OpKeyType<1>(3,3)]=T(1.); n2I[OpKeyType<1>(6,6)]=T(1.); n2I[OpKeyType<1>(7,7)]=T(1.);
		idI[OpKeyType<1>(0,0)]=T(1.);idI[OpKeyType<1>(1,1)]=T(1.); idI[OpKeyType<1>(2,2)]=T(1.);idI[OpKeyType<1>(3,3)]=T(1.);
		idI[OpKeyType<1>(4,4)]=T(1.);idI[OpKeyType<1>(5,5)]=T(1.); idI[OpKeyType<1>(6,6)]=T(1.);idI[OpKeyType<1>(7,7)]=T(1.);

		szI[OpKeyType<1>(0,0)]=T(-0.5); szI[OpKeyType<1>(1,1)]=T(-0.5); szI[OpKeyType<1>(2,2)]=T(-0.5); szI[OpKeyType<1>(3,3)]=T(-0.5);
		szI[OpKeyType<1>(4,4)]=T(0.5); szI[OpKeyType<1>(5,5)]=T(0.5); szI[OpKeyType<1>(6,6)]=T(0.5); szI[OpKeyType<1>(7,7)]=T(0.5);
		sxI[OpKeyType<1>(0,4)]=T(0.5); sxI[OpKeyType<1>(1,5)]=T(0.5); sxI[OpKeyType<1>(2,6)]=T(0.5); sxI[OpKeyType<1>(3,7)]=T(0.5);
		sxI[OpKeyType<1>(4,0)]=T(0.5); sxI[OpKeyType<1>(5,1)]=T(0.5); sxI[OpKeyType<1>(6,2)]=T(0.5); sxI[OpKeyType<1>(7,3)]=T(0.5);

		pI[OpKeyType<1>(0,0)]=T(1.); pI[OpKeyType<1>(1,1)]=T(-1.); pI[OpKeyType<1>(2,2)]=T(-1.); pI[OpKeyType<1>(3,3)]=T(1.);
		pI[OpKeyType<1>(4,4)]=T(1.); pI[OpKeyType<1>(5,5)]=T(-1.); pI[OpKeyType<1>(6,6)]=T(-1.); pI[OpKeyType<1>(7,7)]=T(1.);

		std::cout << t << std::endl << eps << std::endl;

		if (impSite!=0){
			c1.setSite(0); c1D.setSite(0); n1.setSite(0);
			c2.setSite(0); c2D.setSite(0); n2.setSite(0); id.setSite(0); p.setSite(0);
			M[ukey<2>(0,0)] = eps[0]*n1+eps[0]*n2+U*(n1*n2);
			M[ukey<2>(0,1)] = -t[0]*c1D*p;
			M[ukey<2>(0,2)] = t[0]*c1*p;
			M[ukey<2>(0,3)] = -t[0]*c2D*p;
			M[ukey<2>(0,4)] = t[0]*c2*p;
			M[ukey<2>(0,5)] = id;
		}else{
			c1I.setSite(0); c1DI.setSite(0); n1I.setSite(0);
			c2I.setSite(0); c2DI.setSite(0); n2I.setSite(0);
			idI.setSite(0); szI.setSite(0); sxI.setSite(0); pI.setSite(0);
			M[ukey<2>(0,0)] = eps[0]*n1I+eps[0]*n2I+U*(n1I*n2I)+omega*sxI+delta*szI+V*szI*(n1I+n2I);
			M[ukey<2>(0,1)] = -t[0]*c1DI*pI;
			M[ukey<2>(0,2)] = t[0]*c1I*pI;
			M[ukey<2>(0,3)] = -t[0]*c2DI*pI;
			M[ukey<2>(0,4)] = t[0]*c2I*pI;
			M[ukey<2>(0,5)] = idI;
		}
		M.Dl=1,M.Dr=6;
		M.squeeze(1e-16);
		(*this)[0]=M;
		M.clear();

		if (impSite!=N-1){
			c1.setSite(N-1); c1D.setSite(N-1); n1.setSite(N-1);
			c2.setSite(N-1); c2D.setSite(N-1); n2.setSite(N-1);
			id.setSite(N-1);
			M[ukey<2>(0,0)] = id;
			M[ukey<2>(1,0)] = c1;
			M[ukey<2>(2,0)] = c1D;
			M[ukey<2>(3,0)] = c2;
			M[ukey<2>(4,0)] = c2D;
			M[ukey<2>(5,0)] = eps[N-1]*n1+eps[N-1]*n2+U*(n1*n2);
		}else{
			c1I.setSite(N-1); c1DI.setSite(N-1); n1I.setSite(N-1);
			c2I.setSite(N-1); c2DI.setSite(N-1); n2I.setSite(N-1);
			idI.setSite(N-1); szI.setSite(N-1); sxI.setSite(N-1);
			M[ukey<2>(0,0)] = idI;
			M[ukey<2>(1,0)] = c1I;
			M[ukey<2>(2,0)] = c1DI;
			M[ukey<2>(3,0)] = c2I;
			M[ukey<2>(4,0)] = c2DI;
			M[ukey<2>(5,0)] = eps[N-1]*n1I+eps[N-1]*n2I+U*(n1I*n2I)+omega*sxI+delta*szI+V*szI*(n1I+n2I);
		}
		M.Dl=6,M.Dr=1;
		M.squeeze(1e-16);
		(*this)[N-1]=M;
		M.clear();

		// intermediate sites
		for (uint site = 1; site < N-1; site++) {
			MPOMat<1,T> Mb;
	    	if(site!=impSite){ // normal site, local Hilbert space=4
				c1.setSite(site); c1D.setSite(site); n1.setSite(site);
				c2.setSite(site); c2D.setSite(site); n2.setSite(site);
				id.setSite(site); p.setSite(site);
				Mb[ukey<2>(0,0)] = id;
				Mb[ukey<2>(1,0)] = c1;
				Mb[ukey<2>(2,0)] = c1D;
				Mb[ukey<2>(3,0)] = c2;
				Mb[ukey<2>(4,0)] = c2D;
				Mb[ukey<2>(5,0)] = eps[site]*n1+eps[site]*n2+U*(n1*n2);
				Mb[ukey<2>(5,1)] = -t[site]*c1D*p;
				Mb[ukey<2>(5,2)] = t[site]*c1*p;
				Mb[ukey<2>(5,3)] = -t[site]*c2D*p;
				Mb[ukey<2>(5,4)] = t[site]*c2*p;
				Mb[ukey<2>(5,5)] = id;
			}else{ // impurity site, local Hilbert space=8
				c1I.setSite(site); c1DI.setSite(site); n1I.setSite(site);
				c2I.setSite(site); c2DI.setSite(site); n2I.setSite(site);
				idI.setSite(site); szI.setSite(site); sxI.setSite(site);  pI.setSite(site);
				Mb[ukey<2>(0,0)] = idI;
				Mb[ukey<2>(1,0)] = c1I;
				Mb[ukey<2>(2,0)] = c1DI;
				Mb[ukey<2>(3,0)] = c2I;
				Mb[ukey<2>(4,0)] = c2DI;
				Mb[ukey<2>(5,0)] = eps[site]*n1I+eps[site]*n2I+U*(n1I*n2I)+omega*sxI+delta*szI+V*szI*(n1I+n2I);
				Mb[ukey<2>(5,1)] = -t[site]*c1DI*pI;
				Mb[ukey<2>(5,2)] = t[site]*c1I*pI;
				Mb[ukey<2>(5,3)] = -t[site]*c2DI*pI;
				Mb[ukey<2>(5,4)] = t[site]*c2I*pI;
				Mb[ukey<2>(5,5)] = idI;
			}
			Mb.Dl=6,Mb.Dr=6;

		    Mb.squeeze(1e-15);
		    (*this)[site]=Mb;
		}
	}

	~HubbImp(){}

};

template<typename T>
class CondImpChannels:public MPO<T> {
public:
	CondImpChannels(uint N, uint channels, uint impSite, double delta, double omega, double V, VectorType t, VectorType eps){
		std::cout << "IMPURITY SITE="<<impSite << std::endl;

		if (impSite==0||impSite==N-1){
			throw std::runtime_error("CondImpChannels impurity is not allowed to be on the first or on the last site!");
		}

		MPOMat<1,T> M;
		SparseLocalOperator<T> c(2,2),cD(2,2),n(2,2),id(2,2), p(2,2);
		SparseLocalOperator<T> cI(4,4),cDI(4,4),nI(4,4),idI(4,4),szI(4,4),sxI(4,4), pI(4,4);
		c[OpKeyType<1>(0,1)]=T(1.);
		cD[OpKeyType<1>(1,0)]=T(1.);
		n[OpKeyType<1>(1,1)]=T(1.);
		id[OpKeyType<1>(0,0)]=T(1.);id[OpKeyType<1>(1,1)]=T(1.);
		p[OpKeyType<1>(0,0)]=T(1.);p[OpKeyType<1>(1,1)]=T(-1.);

		cI[OpKeyType<1>(0,2)]=T(1.); cI[OpKeyType<1>(1,3)]=T(1.);
		cDI[OpKeyType<1>(2,0)]=T(1.); cDI[OpKeyType<1>(3,1)]=T(1.);
		nI[OpKeyType<1>(2,2)]=T(1.); nI[OpKeyType<1>(3,3)]=T(1.);
		idI[OpKeyType<1>(0,0)]=T(1.); idI[OpKeyType<1>(1,1)]=T(1.); idI[OpKeyType<1>(2,2)]=T(1.); idI[OpKeyType<1>(3,3)]=T(1.);
		szI[OpKeyType<1>(0,0)]=T(-0.5); szI[OpKeyType<1>(1,1)]=T(0.5); szI[OpKeyType<1>(2,2)]=T(-0.5); szI[OpKeyType<1>(3,3)]=T(0.5);
		sxI[OpKeyType<1>(0,1)]=T(0.5); sxI[OpKeyType<1>(1,0)]=T(0.5); sxI[OpKeyType<1>(2,3)]=T(0.5); sxI[OpKeyType<1>(3,2)]=T(0.5);
		pI[OpKeyType<1>(0,0)]=T(1.);pI[OpKeyType<1>(1,1)]=T(1.); pI[OpKeyType<1>(2,2)]=T(-1.);pI[OpKeyType<1>(3,3)]=T(-1.);

		VectorType tTmp=t, epsTmp=eps;
		t.resize(N*channels-1); t.clear();
		eps.resize(N*channels); eps.clear();
		for (uint i=0; i<channels; i++){
			for (uint j=0; j<N; j++){
				eps(j+i*N)=epsTmp(j);
				if (j!=N-1){
					t(j+i*N)=tTmp(j);
				}
			}
		}
		std::cout << "e="<<eps <<std::endl;
		std::cout << "t="<<t<<std::endl;

		c.setSite(0); cD.setSite(0); n.setSite(0); id.setSite(0); p.setSite(0);
		M[ukey<2>(0,0)] = eps[0]*n;
		M[ukey<2>(0,1)] = -t[0]*cD*p;
		M[ukey<2>(0,2)] = t[0]*c*p;
		M[ukey<2>(0,3)] = id;
		M.Dl=1,M.Dr=4;
		M.squeeze(1e-16);
		(*this)[0]=M;
		M.clear();

		c.setSite(channels*N-1); cD.setSite(channels*N-1); n.setSite(channels*N-1); id.setSite(channels*N-1);
		M[ukey<2>(0,0)] = id;
		M[ukey<2>(1,0)] = c;
		M[ukey<2>(2,0)] = cD;
		M[ukey<2>(3,0)] = eps[channels*N-1]*n;
		M.Dl=4,M.Dr=1;
		M.squeeze(1e-16);
		(*this)[channels*N-1]=M;
		M.clear();

		// intermediate sites
		uint nImp=0;
		for (uint i=0; i<channels; i++){
			for (uint j=0; j<N; j++){
				uint site=j+i*N;
				if (site!=0 && site!=channels*N-1){
					MPOMat<1,T> Mb;
					if (j!=impSite){
						c.setSite(site); cD.setSite(site); n.setSite(site); id.setSite(site); p.setSite(site);
						Mb[ukey<2>(0,0)] = id;
						Mb[ukey<2>(1,0)] = c;
						Mb[ukey<2>(2,0)] = cD;
						for (uint k=0; k<nImp; k++){
							Mb[ukey<2>(3+k,3+k)] = id;
						}
						Mb[ukey<2>(3+nImp,0)] = eps[site]*n;
						Mb[ukey<2>(3+nImp,1)] = -t[site]*cD*p;
						Mb[ukey<2>(3+nImp,2)] = t[site]*c*p;//-
						Mb[ukey<2>(3+nImp,3+nImp)] = id;
						Mb.Dl=4+nImp,Mb.Dr=4+nImp;
					}else{
						if (i==0){
							cI.setSite(site); cDI.setSite(site); nI.setSite(site); idI.setSite(site); szI.setSite(site); sxI.setSite(site);pI.setSite(site);
							Mb[ukey<2>(0,0)] = idI;
							Mb[ukey<2>(1,0)] = cI;
							Mb[ukey<2>(2,0)] = cDI;
							Mb[ukey<2>(3,0)] = (nI*eps[site])+(sxI*omega)+(szI*delta)+(szI*nI)*V;
							Mb[ukey<2>(3,1)] = -t[site]*cDI*pI;
							Mb[ukey<2>(3,2)] = t[site]*cI*pI;//-
							for (uint k=0; k<channels-1; k++){
								Mb[ukey<2>(3,3+k)] = V*szI;
							}
							Mb[ukey<2>(3,2+channels)] = idI;

							nImp=channels-1;
							Mb.Dl=4,Mb.Dr=4+nImp;
						}else{
							cI.setSite(site); cDI.setSite(site); nI.setSite(site); idI.setSite(site); pI.setSite(site);
							Mb[ukey<2>(0,0)] = idI;
							Mb[ukey<2>(1,0)] = cI;
							Mb[ukey<2>(2,0)] = cDI;
							Mb[ukey<2>(3,0)] = nI;
							nImp--;
							for (uint k=0; k<nImp; k++){
								Mb[ukey<2>(4+k,3+k)] = idI;
							}

							Mb[ukey<2>(4+nImp,0)] = eps[site]*nI;
							Mb[ukey<2>(4+nImp,1)] = -t[site]*cDI*pI;
							Mb[ukey<2>(4+nImp,2)] = t[site]*cI*pI;//-
							Mb[ukey<2>(4+nImp,3+nImp)] = idI;
							Mb.Dl=5+nImp,Mb.Dr=4+nImp;
						}
					}
					Mb.squeeze(1e-15);
					(*this)[site]=Mb;
				}
			}
		}

		for (uint i=0; i<channels*N; i++){
			std::cout << (*this)[i] << std::endl;
		}
	}

	~CondImpChannels(){}

};

template<typename T>
class BHMvImp:public MPO<T>{
public:
	BHMvImp(VectorType t1, VectorType t2, VectorType eps1, VectorType eps2, double V, double U, uint maxdim){
		uint N =eps1.size();
		assert(t1.size()==t2.size());
		assert(t1.size()==eps1.size()-1);
		assert(t1.size()==eps2.size()-1);

		MPOMat<1,T> M;
		SparseLocalOperator<T> b(2*maxdim,2*maxdim),bd(2*maxdim,2*maxdim),id(2*maxdim,2*maxdim),
				n1(2*maxdim,2*maxdim),n2(2*maxdim,2*maxdim),c(2*maxdim,2*maxdim),cd(2*maxdim,2*maxdim);
		for(uint ind = 0;ind<maxdim;++ind){
			id[OpKeyType<1>(ind,ind)]=T(1.0);
			id[OpKeyType<1>(ind+maxdim,ind+maxdim)]=T(1.0);
			if(ind!=0){
				n1[OpKeyType<1>(ind,ind)]=T(ind*1.0);
				n1[OpKeyType<1>(ind+maxdim,ind+maxdim)]=T(ind*1.0);
				b[OpKeyType<1>(ind-1,ind)] = T(sqrt(ind));
				b[OpKeyType<1>(ind-1+maxdim,ind+maxdim)] = T(sqrt(ind));
			}
			if(ind<(maxdim-1)){
				bd[OpKeyType<1>(ind+1,ind)] = T(sqrt(ind+1.));
				bd[OpKeyType<1>(ind+1+maxdim,ind+maxdim)] = T(sqrt(ind+1.));
			}
		}
		for (uint ind=0; ind<maxdim;++ind){
			n2[OpKeyType<1>(ind+maxdim,ind+maxdim)]=T(1.0);
			c[OpKeyType<1>(ind,ind+maxdim)]=T(1.);
			cd[OpKeyType<1>(ind+maxdim,ind)]=T(1.);
		}
		b.setSite(0);bd.setSite(0);n1.setSite(0);
		c.setSite(0);cd.setSite(0);n2.setSite(0);
		id.setSite(0);

		M[ukey<2>(0,0)] = eps1[0]*n1+eps2[0]*n2+V*(n1*n2)+U/2.*(n1*(n1-id));
		M[ukey<2>(0,1)] = -t1[0]*bd;//p
		M[ukey<2>(0,2)] = -t1[0]*b;
		M[ukey<2>(0,3)] = -t2[0]*cd;
		M[ukey<2>(0,4)] = -t2[0]*c;
		M[ukey<2>(0,5)] = id;

		M.Dl=1,M.Dr=6;
		M.squeeze(1e-16);
		(*this)[0]=M;
		M.clear();

		b.setSite(N-1); bd.setSite(N-1); n1.setSite(N-1);
		c.setSite(N-1); cd.setSite(N-1); n2.setSite(N-1);
		id.setSite(N-1);

		M[ukey<2>(0,0)] = id;
		M[ukey<2>(1,0)] = b;
		M[ukey<2>(2,0)] = bd;
		M[ukey<2>(3,0)] = c;
		M[ukey<2>(4,0)] = cd;
		M[ukey<2>(5,0)] = eps1[N-1]*n1+eps2[N-1]*n2+V*(n1*n2)+U/2.*(n1*(n1-id));
		M.Dl=6,M.Dr=1;
		M.squeeze(1e-16);
		(*this)[N-1]=M;
		M.clear();

		// intermediate sites
		for (uint site = 1; site < N-1; site++) {
			MPOMat<1,T> Mb;
			b.setSite(site); bd.setSite(site); n1.setSite(site);
			c.setSite(site); cd.setSite(site); n2.setSite(site);
			id.setSite(site);
			Mb[ukey<2>(0,0)] = id;
			Mb[ukey<2>(1,0)] = b;
			Mb[ukey<2>(2,0)] = bd;
			Mb[ukey<2>(3,0)] = c;
			Mb[ukey<2>(4,0)] = cd;
			Mb[ukey<2>(5,0)] = eps1[site]*n1+eps2[site]*n2+V*(n1*n2)+U/2.*(n1*(n1-id));
			Mb[ukey<2>(5,1)] = -t1[site]*bd;//p
			Mb[ukey<2>(5,2)] = -t1[site]*b;
			Mb[ukey<2>(5,3)] = -t2[site]*cd;
			Mb[ukey<2>(5,4)] = -t2[site]*c;
			Mb[ukey<2>(5,5)] = id;
			Mb.Dl=6,Mb.Dr=6;

		    Mb.squeeze(1e-15);
		    (*this)[site]=Mb;
		}
	}
	~BHMvImp(){};
};


template<typename T>
class IntMvImp:public MPO<T> {
public:
	IntMvImp(uint N, VectorType t1, VectorType t2, VectorType eps1, VectorType eps2, double V){

		assert(eps1.size() == N);
		assert(t1.size() == N-1);
		assert(eps2.size() == N);
		assert(t2.size() == N-1);
		double sign=-1.;

		MPOMat<1,T> M;
		//this->resize(N);
		SparseLocalOperator<T> c1(4,4),c1D(4,4),c2(4,4),c2D(4,4),n1(4,4),n2(4,4),id(4,4), p(4,4);
		c1[OpKeyType<1>(0,1)]=T(1.); c1[OpKeyType<1>(2,3)]=T(1.);
		c2[OpKeyType<1>(0,2)]=T(1.); c2[OpKeyType<1>(1,3)]=T(sign);
		c1D[OpKeyType<1>(1,0)]=T(1.); c1D[OpKeyType<1>(3,2)]=T(1.);
		c2D[OpKeyType<1>(2,0)]=T(1.); c2D[OpKeyType<1>(3,1)]=T(sign);
		n1[OpKeyType<1>(1,1)]=T(1.); n1[OpKeyType<1>(3,3)]=T(1.);
		n2[OpKeyType<1>(2,2)]=T(1.); n2[OpKeyType<1>(3,3)]=T(1.);
		id[OpKeyType<1>(0,0)]=T(1.);id[OpKeyType<1>(1,1)]=T(1.); id[OpKeyType<1>(2,2)]=T(1.);id[OpKeyType<1>(3,3)]=T(1.);
		p[OpKeyType<1>(0,0)]=T(1.); p[OpKeyType<1>(1,1)]=T(sign); p[OpKeyType<1>(2,2)]=T(sign); p[OpKeyType<1>(3,3)]=T(1.);

		std::cout << t1 << std::endl << t2 << std::endl<< eps1 << std::endl << eps2 << std::endl;

		c1.setSite(0); c1D.setSite(0); n1.setSite(0);
		c2.setSite(0); c2D.setSite(0); n2.setSite(0);
		id.setSite(0); p.setSite(0);
		M[ukey<2>(0,0)] = eps1[0]*n1+eps2[0]*n2+V*(n1*n2);
		M[ukey<2>(0,1)] = sign*t1[0]*c1D;//p
		M[ukey<2>(0,2)] = t1[0]*c1;
		M[ukey<2>(0,3)] = sign*t2[0]*c2D;
		M[ukey<2>(0,4)] = t2[0]*c2;
		M[ukey<2>(0,5)] = id;

		M.Dl=1,M.Dr=6;
		M.squeeze(1e-16);
		(*this)[0]=M;
		M.clear();

		c1.setSite(N-1); c1D.setSite(N-1); n1.setSite(N-1);
		c2.setSite(N-1); c2D.setSite(N-1); n2.setSite(N-1);
		id.setSite(N-1);
		M[ukey<2>(0,0)] = id;
		M[ukey<2>(1,0)] = p*c1;
		M[ukey<2>(2,0)] = p*c1D;
		M[ukey<2>(3,0)] = p*c2;
		M[ukey<2>(4,0)] = p*c2D;
		M[ukey<2>(5,0)] = eps1[N-1]*n1+eps2[N-1]*n2+V*(n1*n2);
		M.Dl=6,M.Dr=1;
		M.squeeze(1e-16);
		(*this)[N-1]=M;
		M.clear();

		// intermediate sites
		for (uint site = 1; site < N-1; site++) {
			MPOMat<1,T> Mb;
			c1.setSite(site); c1D.setSite(site); n1.setSite(site);
			c2.setSite(site); c2D.setSite(site); n2.setSite(site);
			id.setSite(site); p.setSite(site);
			Mb[ukey<2>(0,0)] = id;
			Mb[ukey<2>(1,0)] = p*c1;
			Mb[ukey<2>(2,0)] = p*c1D;
			Mb[ukey<2>(3,0)] = p*c2;
			Mb[ukey<2>(4,0)] = p*c2D;
			Mb[ukey<2>(5,0)] = eps1[site]*n1+eps2[site]*n2+V*(n1*n2);
			Mb[ukey<2>(5,1)] = sign*t1[site]*c1D;//p
			Mb[ukey<2>(5,2)] = t1[site]*c1;
			Mb[ukey<2>(5,3)] = sign*t2[site]*c2D;
			Mb[ukey<2>(5,4)] = t2[site]*c2;
			Mb[ukey<2>(5,5)] = id;
			Mb.Dl=6,Mb.Dr=6;

		    Mb.squeeze(1e-15);
		    (*this)[site]=Mb;
		}
	}

	~IntMvImp(){}

};


template<typename T>
class IntMvImpInt:public MPO<T> {
public:
	IntMvImpInt(uint N, VectorType t1, VectorType t2, VectorType eps1, VectorType eps2, double V, double U){
		double sign=-1.;

		MPOMat<1,T> M;
		//this->resize(N);
		SparseLocalOperator<T> c1(4,4),c1D(4,4),c2(4,4),c2D(4,4),n1(4,4),n2(4,4),id(4,4), p(4,4);
		c1[OpKeyType<1>(0,1)]=T(1.); c1[OpKeyType<1>(2,3)]=T(1.);
		c2[OpKeyType<1>(0,2)]=T(1.); c2[OpKeyType<1>(1,3)]=T(sign);
		c1D[OpKeyType<1>(1,0)]=T(1.); c1D[OpKeyType<1>(3,2)]=T(1.);
		c2D[OpKeyType<1>(2,0)]=T(1.); c2D[OpKeyType<1>(3,1)]=T(sign);
		n1[OpKeyType<1>(1,1)]=T(1.); n1[OpKeyType<1>(3,3)]=T(1.);
		n2[OpKeyType<1>(2,2)]=T(1.); n2[OpKeyType<1>(3,3)]=T(1.);
		id[OpKeyType<1>(0,0)]=T(1.);id[OpKeyType<1>(1,1)]=T(1.); id[OpKeyType<1>(2,2)]=T(1.);id[OpKeyType<1>(3,3)]=T(1.);
		p[OpKeyType<1>(0,0)]=T(1.); p[OpKeyType<1>(1,1)]=T(sign); p[OpKeyType<1>(2,2)]=T(sign); p[OpKeyType<1>(3,3)]=T(1.);

		std::cout << t1 << std::endl << t2 << std::endl<< eps1 << std::endl << eps2 << std::endl;

		c1.setSite(0); c1D.setSite(0); n1.setSite(0);
		c2.setSite(0); c2D.setSite(0); n2.setSite(0);
		id.setSite(0); p.setSite(0);
		M[ukey<2>(0,0)] = eps1[0]*n1+eps2[0]*n2+V*(n1*n2);
		M[ukey<2>(0,1)] = sign*t1[0]*c1D;//p
		M[ukey<2>(0,2)] = t1[0]*c1;
		M[ukey<2>(0,3)] = sign*t2[0]*c2D;
		M[ukey<2>(0,4)] = t2[0]*c2;
		M[ukey<2>(0,5)] = U*n1;
		M[ukey<2>(0,6)] = id;

		M.Dl=1,M.Dr=7;
		M.squeeze(1e-16);
		(*this)[0]=M;
		M.clear();

		c1.setSite(N-1); c1D.setSite(N-1); n1.setSite(N-1);
		c2.setSite(N-1); c2D.setSite(N-1); n2.setSite(N-1);
		id.setSite(N-1);
		M[ukey<2>(0,0)] = id;
		M[ukey<2>(1,0)] = p*c1;
		M[ukey<2>(2,0)] = p*c1D;
		M[ukey<2>(3,0)] = p*c2;
		M[ukey<2>(4,0)] = p*c2D;
		M[ukey<2>(5,0)] = n1;
		M[ukey<2>(6,0)] = eps1[N-1]*n1+eps2[N-1]*n2+V*(n1*n2);
		M.Dl=7,M.Dr=1;
		M.squeeze(1e-16);
		(*this)[N-1]=M;
		M.clear();

		// intermediate sites
		for (uint site = 1; site < N-1; site++) {
			MPOMat<1,T> Mb;
			c1.setSite(site); c1D.setSite(site); n1.setSite(site);
			c2.setSite(site); c2D.setSite(site); n2.setSite(site);
			id.setSite(site); p.setSite(site);
			Mb[ukey<2>(0,0)] = id;
			Mb[ukey<2>(1,0)] = p*c1;
			Mb[ukey<2>(2,0)] = p*c1D;
			Mb[ukey<2>(3,0)] = p*c2;
			Mb[ukey<2>(4,0)] = p*c2D;
			Mb[ukey<2>(5,0)] = n1;
			Mb[ukey<2>(6,0)] = eps1[site]*n1+eps2[site]*n2+V*(n1*n2);
			Mb[ukey<2>(6,1)] = sign*t1[site]*c1D;//p
			Mb[ukey<2>(6,2)] = t1[site]*c1;
			Mb[ukey<2>(6,3)] = sign*t2[site]*c2D;
			Mb[ukey<2>(6,4)] = t2[site]*c2;
			Mb[ukey<2>(6,5)] = U*n1;
			Mb[ukey<2>(6,6)] = id;
			Mb.Dl=7,Mb.Dr=7;

		    Mb.squeeze(1e-15);
		    (*this)[site]=Mb;
		}
	}

	~IntMvImpInt(){}
};

template<typename T>
class LRIsingExpMPO:public MPO<T>
{
public:
  LRIsingExpMPO(uint N, VectorType p, double h);
  ~LRIsingExpMPO(){};
};

template<typename T>
LRIsingExpMPO<T>:: LRIsingExpMPO(uint N, VectorType p, double h){

  MPOMat<1,T> M;
  //  this->resize(N);
  SparseLocalOperator<T> sz(2,2),sx(2,2),id(2,2);
  sz[OpKeyType<1>(0,0)]=T(-1.); sz[OpKeyType<1>(1,1)]=T(1.);
  sx[OpKeyType<1>(1,0)]=T(1.);  sx[OpKeyType<1>(0,1)]=T(1.);
  id[OpKeyType<1>(0,0)]=T(1.);  id[OpKeyType<1>(1,1)]=T(1.);

  int l=p.size()/2;

  sz.setSite(0);sx.setSite(0);id.setSite(0);
  M[ukey<2>(0,0)] = -h*sx;
  for (uint i=0; i<l; ++i){
	  M[ukey<2>(0,i+1)]=-p(i)*sz;
  }
  M[ukey<2>(0,l+1)] = id;
  M.Dl=1,M.Dr=l+2;
  M.squeeze(1e-16);
  (*this)[0]=M;

  M.clear();
  sz.setSite(N-1);sx.setSite(N-1);id.setSite(N-1);
  M[ukey<2>(0,0)] = id;
  for (uint i=0; i<l; ++i){
	  M[ukey<2>(i+1,0)]=sz;
  }
  M[ukey<2>(l+1,0)] = -h*sx;
  M.Dl=l+2,M.Dr=1;

  M.squeeze(1e-16);
  (*this)[N-1]=M;

  for(uint site = 1; site<(N-1); ++site){
    MPOMat<1,T> Mb;
    sz.setSite(site);sx.setSite(site);id.setSite(site);
    Mb[ukey<2>(0,0)] = id;
    Mb[ukey<2>(l+1,0)] = -h*sx;
    for (uint i=0; i<l; ++i){
  	  Mb[ukey<2>(i+1,0)]=sz;
  	  Mb[ukey<2>(i+1,i+1)]=p(l+i)*id;
	  Mb[ukey<2>(l+1,i+1)]=-p(i)*sz;
    }
    Mb[ukey<2>(l+1,l+1)] = id;

    Mb.Dl=l+2,Mb.Dr=l+2;
    Mb.squeeze(1e-16);
    (*this)[site]=Mb;
  }
}



template<typename T>
class LRIsingMPO:public MPO<T>
{
public:
  LRIsingMPO(uint N, double J, double ia, double h);
  ~LRIsingMPO(){};
};

template<typename T>
LRIsingMPO<T>:: LRIsingMPO(uint N, double J, double ia, double h){
	VectorType Jv(N-1);
	for (uint i=1; i<N; ++i){
		Jv(i-1)=J/(std::pow((double)i,1./ia));
	}
	std::cout << "J="<<Jv << std::endl;
  MPOMat<1,T> M;
  //  this->resize(N);
  SparseLocalOperator<T> sz(2,2),sx(2,2),id(2,2);
  sz[OpKeyType<1>(0,0)]=T(-1.); sz[OpKeyType<1>(1,1)]=T(1.);
  sx[OpKeyType<1>(1,0)]=T(1.);  sx[OpKeyType<1>(0,1)]=T(1.);
  id[OpKeyType<1>(0,0)]=T(1.);  id[OpKeyType<1>(1,1)]=T(1.);

  sz.setSite(0);sx.setSite(0);id.setSite(0);
  M[ukey<2>(0,0)] = -h*sx;
  for (uint i=0; i<Jv.size(); ++i){
	  M[ukey<2>(0,i+1)]=-Jv(i)*sz;
  }
  M[ukey<2>(0,N)] = id;
  M.Dl=1,M.Dr=N+1;
  M.squeeze(1e-16);
  (*this)[0]=M;

  M.clear();
  sz.setSite(N-1);sx.setSite(N-1);id.setSite(N-1);
  M[ukey<2>(0,0)] = id;
  M[ukey<2>(1,0)] = sz;
  M[ukey<2>(N,0)] = -h*sx;
  M.Dl=N+1,M.Dr=1;
  M.squeeze(1e-16);
  (*this)[N-1]=M;

  for(uint site = 1; site<(N-1); ++site){
    MPOMat<1,T> Mb;
    sz.setSite(site);sx.setSite(site);id.setSite(site);
    Mb[ukey<2>(0,0)] = id;
    Mb[ukey<2>(1,0)] = sz;
    Mb[ukey<2>(N,0)] = -h*sx;
    for (uint i=0; i<Jv.size(); ++i){
      Mb[ukey<2>(i+2,i+1)] = id;
  	  Mb[ukey<2>(N,i+1)]=-Jv(i)*sz;
    }
    Mb[ukey<2>(N,N)] = id;

    Mb.Dl=N+1,Mb.Dr=N+1;
    Mb.squeeze(1e-16);
    (*this)[site]=Mb;
  }
}


template<typename T>
class TransverseIsingMPO:public MPO<T>
{
public:
  TransverseIsingMPO(VectorType J,VectorType h);
  ~TransverseIsingMPO(){};
};

template<typename T>
TransverseIsingMPO<T>:: TransverseIsingMPO(VectorType J,VectorType h){
  uint N = h.size();
  MPOMat<1,T> M;
  //  this->resize(N);
  SparseLocalOperator<T> sz(2,2),sx(2,2),id(2,2);
  sz[OpKeyType<1>(0,0)]=T(-1.); sz[OpKeyType<1>(1,1)]=T(1.);
  sx[OpKeyType<1>(1,0)]=T(1.);  sx[OpKeyType<1>(0,1)]=T(1.);
  id[OpKeyType<1>(0,0)]=T(1.);  id[OpKeyType<1>(1,1)]=T(1.);

  sz.setSite(0);sx.setSite(0);id.setSite(0);
  M[ukey<2>(0,0)] = -h(0)*sx;
  M[ukey<2>(0,1)] = -J(0)*sz;
  M[ukey<2>(0,2)] = id;
  M.Dl=1,M.Dr=3;
  M.squeeze(1e-16);
  (*this)[0]=M;

  M.clear();
  sz.setSite(N-1);sx.setSite(N-1);id.setSite(N-1);
  M[ukey<2>(0,0)] = id;
  M[ukey<2>(1,0)] = sz;
  M[ukey<2>(2,0)] = -h(N-1)*sx;
  M.Dl=3,M.Dr=1;
  M.squeeze(1e-16);
  (*this)[N-1]=M;

  for(uint site = 1;site<(N-1);++site){
    MPOMat<1,T> Mb;
    sz.setSite(site);sx.setSite(site);id.setSite(site);
    Mb[ukey<2>(0,0)] = id;
    Mb[ukey<2>(1,0)] = sz;
    Mb[ukey<2>(2,0)] = -h(site)*sx;
    Mb[ukey<2>(2,1)] = -J(site)*sz;
    Mb[ukey<2>(2,2)] = id;

    Mb.Dl=3,Mb.Dr=3;
    Mb.squeeze(1e-16);
    (*this)[site]=Mb;
  }
}


//definition of bose hubbard model:
//-\sum_i hop(i)b_i^{\dagger}b_{i+1}+h.c. + \sum_i U(i)/2 n_i(n_i-1)
template<typename T>
class BoseHubbardMPO:public MPO<T>
{
public:
  BoseHubbardMPO(const VectorType hop,const VectorType U,const VectorType mu,const UintVectorType maxdim);
  ~BoseHubbardMPO(){};
  virtual T measure(const MPS<AbKeyType,T>&mps)const {return MPO<T>::measure(mps);}
  virtual T measure(const CanonizedMPS<AbKeyType,T>&mps)const {return MPO<T>::measure(mps);}

};

template<typename T>
BoseHubbardMPO<T>::BoseHubbardMPO(const VectorType hop,const VectorType U,const VectorType mu,const UintVectorType maxdim){
  uint N =mu.size();
  assert(mu.size()==U.size());
  assert(maxdim.size()==mu.size());
  MPOMat<1,T> M;
  //  this->resize(N);
  //all statistics are set to bosonic (_statistics=1)
  SparseLocalOperator<T> b(maxdim(0),maxdim(0)),bd(maxdim(0),maxdim(0)),id(maxdim(0),maxdim(0)),n(maxdim(0),maxdim(0));
  for(int ind = 0;ind<maxdim(0);++ind){
    n[OpKeyType<1>(ind,ind)]=T(ind*1.0);
    id[OpKeyType<1>(ind,ind)]=T(1.0);
    if(ind!=0)
      b[OpKeyType<1>(ind-1,ind)] = T(sqrt(ind));
    if(ind<(maxdim(0)-1))
      bd[OpKeyType<1>(ind+1,ind)] = T(sqrt(ind+1));
  }
  b.setSite(0);bd.setSite(0);id.setSite(0);n.setSite(0);
  

  M[ukey<2>(0,0)] = -mu(0)*n+(U(0)/2.*(n*(n-id)));
  M[ukey<2>(0,1)] = -hop(0)*b;
  M[ukey<2>(0,2)] = -hop(0)*bd;
  M[ukey<2>(0,3)] = id;
  M.Dl=1,M.Dr=4;
  M.squeeze(1e-15);
  (*this)[0]=M;
  M.clear();


  b.clear();bd.clear();n.clear();id.clear();
  b.resize(maxdim(N-1),maxdim(N-1)),bd.resize(maxdim(N-1),maxdim(N-1)),id.resize(maxdim(N-1),maxdim(N-1)),n.resize(maxdim(N-1),maxdim(N-1));
  for(int ind = 0;ind<maxdim(0);++ind){
    n[OpKeyType<1>(ind,ind)]=T(ind*1.0);
    id[OpKeyType<1>(ind,ind)]=T(1.0);
    if(ind!=0)
      b[OpKeyType<1>(ind-1,ind)] = T(sqrt(ind));
    if(ind<(maxdim(N-1)-1))
      bd[OpKeyType<1>(ind+1,ind)] = T(sqrt(ind+1));
  }
  b.setSite(N-1);bd.setSite(N-1);id.setSite(N-1);n.setSite(N-1);

  M[ukey<2>(0,0)] = id;
  M[ukey<2>(1,0)] = bd;
  M[ukey<2>(2,0)] = b;
  M[ukey<2>(3,0)] = -mu(N-1)*n+(U(N-1)/2.*(n*(n-id)));
  M[ukey<2>(3,0)].resize(maxdim(N-1),maxdim(N-1));

  M.Dl=4,M.Dr=1;
  M.squeeze(1e-15);
  (*this)[N-1]=M;
  for(uint site = 1;site<(N-1);++site){
    MPOMat<1,T> Mb;

    b.clear();bd.clear();n.clear();id.clear();
    b.resize(maxdim(site),maxdim(site)),bd.resize(maxdim(site),maxdim(site)),id.resize(maxdim(site),maxdim(site)),n.resize(maxdim(site),maxdim(site));
    for(int ind = 0;ind<maxdim(0);++ind){
      n[OpKeyType<1>(ind,ind)]=T(ind*1.0);
      id[OpKeyType<1>(ind,ind)]=T(1.0);
      if(ind!=0)
	b[OpKeyType<1>(ind-1,ind)] = T(sqrt(ind));
      if(ind<(maxdim(site)-1))
	bd[OpKeyType<1>(ind+1,ind)] = T(sqrt(ind+1));
    }
    b.setSite(site);bd.setSite(site);id.setSite(site);n.setSite(site);

    Mb[ukey<2>(0,0)]= id;
    Mb[ukey<2>(1,0)]= bd;
    Mb[ukey<2>(2,0)] = b;
    Mb[ukey<2>(3,0)] = -mu(site)*n+U(site)/2.*(n*(n-id));
    Mb[ukey<2>(3,0)].resize(maxdim(site),maxdim(site));//this resize is necessary if you run the program in debug mode; the mat sizes are checked for 
                                                       //consisteny, though actually it should't matter cause I never do explicit indexing
     
    Mb[ukey<2>(3,1)] = -hop(site)*b;
    Mb[ukey<2>(3,2)] = -hop(site)*bd;
    Mb[ukey<2>(3,3)] = id;
    Mb.Dl=4,Mb.Dr=4;
    Mb.squeeze(1e-15);
    (*this)[site]=Mb;
  }
}


template<typename T>
class BoseHubbard2ChMPO:public MPO<T>
{
public:
  BoseHubbard2ChMPO(const VectorType J1,const VectorType J2,const VectorType U1,const VectorType U2,const VectorType U12,
		    const VectorType mu1,const VectorType mu2,const UintVectorType maxdim1,
		    const UintVectorType maxdim2);
  ~BoseHubbard2ChMPO(){};
  virtual T measure(const MPS<AbKeyType,T>&mps)const {return MPO<T>::measure(mps);}
  virtual T measure(const CanonizedMPS<AbKeyType,T>&mps)const {return MPO<T>::measure(mps);}

};

template<typename T>
BoseHubbard2ChMPO<T>::BoseHubbard2ChMPO(const VectorType J1,const VectorType J2,const VectorType U1,const VectorType U2,const VectorType U12,
					const VectorType mu1,const VectorType mu2,const UintVectorType maxdim1,
					const UintVectorType maxdim2){
  uint N =mu1.size();
  assert(N==mu2.size());
  assert(N==U1.size());
  assert(N==U2.size());
  assert(N==U12.size());
  assert((N-1)==J1.size());
  assert((N-1)==J2.size());
  assert(maxdim1.size()==N);
  assert(maxdim2.size()==N);

  MPOMat<1,T> M;
  //  this->resize(N);
  //all statistics are set to bosonic (_statistics=1)
  SparseLocalOperator<T> b1_(maxdim1(0)*maxdim2(0),maxdim1(0)*maxdim2(0)),bd1_(maxdim1(0)*maxdim2(0),maxdim1(0)*maxdim2(0)),
    id_(maxdim1(0)*maxdim2(0),maxdim1(0)*maxdim2(0)),n1_(maxdim1(0)*maxdim2(0),maxdim1(0)*maxdim2(0)),
    b2_(maxdim1(0)*maxdim2(0),maxdim1(0)*maxdim2(0)),bd2_(maxdim1(0)*maxdim2(0),maxdim1(0)*maxdim2(0)),n2_(maxdim1(0)*maxdim2(0),maxdim1(0)*maxdim2(0));

  SparseLocalOperator<T> b1(maxdim1(0),maxdim1(0)),bd1(maxdim1(0),maxdim1(0)),
    id1(maxdim1(0),maxdim1(0)),n1(maxdim1(0),maxdim1(0)),b2(maxdim2(0),maxdim2(0)),bd2(maxdim2(0),maxdim2(0)),
    id2(maxdim2(0),maxdim2(0)),n2(maxdim2(0),maxdim2(0));


  for(int ind = 0;ind<maxdim1(0);++ind){
    n1[OpKeyType<1>(ind,ind)]=T(ind*1.0);
    id1[OpKeyType<1>(ind,ind)]=T(1.0);
    if(ind!=0)
      b1[OpKeyType<1>(ind-1,ind)] = T(sqrt(ind));
    if(ind<(maxdim1(0)-1))
      bd1[OpKeyType<1>(ind+1,ind)] = T(sqrt(ind+1));
  }
  b1.setSite(0);bd1.setSite(0);id1.setSite(0);n1.setSite(0);
  
  for(int ind = 0;ind<maxdim2(0);++ind){
    n2[OpKeyType<1>(ind,ind)]=T(ind*1.0);
    id2[OpKeyType<1>(ind,ind)]=T(1.0);
    if(ind!=0)
      b2[OpKeyType<1>(ind-1,ind)] = T(sqrt(ind));
    if(ind<(maxdim2(0)-1))
      bd2[OpKeyType<1>(ind+1,ind)] = T(sqrt(ind+1));
  }
  b2.setSite(0);bd2.setSite(0);id2.setSite(0);n2.setSite(0);

//  cout<<"#######################this is the first test: ###########################"<<endl;
//  cout<<endl;
//
//  cout<<"[b1,b1']:"<<endl;
//  cout<<((b1*bd1)-(bd1*b1))<<endl;
//  cout<<endl;
//
//  cout<<"[b2,b2']:"<<endl;
//  cout<<((b2*bd2)-(bd2*b2))<<endl;
//  cout<<endl;
//
//  cin.get();

  
  b1_=kronecker(b1,id2);
  bd1_=kronecker(bd1,id2);
  n1_=kronecker(n1,id2);


  b2_=kronecker(id1,b2);
  bd2_=kronecker(id1,bd2);
  n2_=kronecker(id1,n2);

  id_=kronecker(id1,id2);  


//  cout<<"#######################this is a test: ###########################"<<endl;
//  cout<<endl;
//  cout<<"[n1_,n2_]:"<<endl;
//  cout<<(n1_*n2_-n2_*n1_)<<endl;
//  cout<<endl;
//  cin.get();
//
//
//  cout<<"[b1,b1']:"<<endl;
//  cout<<(b1_*bd1_-bd1_*b1_)<<endl;
//  cout<<endl;
//
//  cout<<"[b1,b1]:"<<endl;
//  cout<<(b1_*b1_-b1_*b1_)<<endl;
//  cout<<endl;
//
//
//  cout<<"[b2,b2']:"<<endl;
//  cout<<(b2_*bd2_-bd2_*b2_)<<endl;
//  cout<<endl;
//
//  cout<<"[b2,b2]:"<<endl;
//  cout<<(b2_*b2_-b2_*b2_)<<endl;
//  cout<<endl;
//
//
//  cout<<"[b1,b2]:"<<endl;
//  cout<<(b1_*b2_-b2_*b1_)<<endl;
//  cout<<endl;
//
//  cout<<"[b1,b2']:"<<endl;
//  cout<<(b1_*bd2_-bd2_*b1_)<<endl;
//  cout<<endl;
//
//
//  cout<<"[b1',b2]:"<<endl;
//  cout<<(bd1_*b2_-b2_*bd1_)<<endl;
//  cout<<endl;
//
//  cin.get();


  b1_.setSite(0);bd1_.setSite(0);n1_.setSite(0);b2_.setSite(0);bd2_.setSite(0);id_.setSite(0);n2_.setSite(0);

  //  M[ukey<2>(0,0)] = U1(0)/2.*(n1_*(n1_-id_))+U2(0)/2.*(n2_*(n2_-id_))+U12(0)*((n1_-(0.5*id_))*(n2_-(0.5*id_)))-mu1(0)*n1_-mu2(0)*n2_;
  M[ukey<2>(0,0)] = U1(0)/2.*(n1_*(n1_-id_))+U2(0)/2.*(n2_*(n2_-id_))+U12(0)*(n1_*n2_)+mu1(0)*n1_+mu2(0)*n2_;
  M[ukey<2>(0,1)] = -J1(0)*bd1_;
  M[ukey<2>(0,2)] = -J1(0)*b1_;
  M[ukey<2>(0,3)] = -J2(0)*bd2_;
  M[ukey<2>(0,4)] = -J2(0)*b2_;
  M[ukey<2>(0,5)] = id_;
  M.Dl=1,M.Dr=6;
  M.squeeze(1e-15);
  (*this)[0]=M;
  M.clear();

//  cout<<n1_<<endl;
//  cout<<n2_<<endl;
//  cout<<n1_*n2_<<endl;
//  cout<<n2_*n1_<<endl;
//  cin.get();
//
  b1_.clear();bd1_.clear();n1_.clear();id_.clear();b2_.clear();bd2_.clear();n2_.clear();
  b1.clear();bd1.clear();n1.clear();id1.clear();b2.clear();bd2.clear();n2.clear();id2.clear();

  b1.resize(maxdim1(N-1),maxdim1(N-1)),bd1.resize(maxdim1(N-1),maxdim1(N-1)),id1.resize(maxdim1(N-1),maxdim1(N-1)),n1.resize(maxdim1(N-1),maxdim1(N-1));
  b2.resize(maxdim2(N-1),maxdim2(N-1)),bd2.resize(maxdim2(N-1),maxdim2(N-1)),id2.resize(maxdim2(N-1),maxdim2(N-1)),n2.resize(maxdim2(N-1),maxdim2(N-1));

  b1_.resize(maxdim1(N-1)*maxdim2(N-1),maxdim1(N-1)*maxdim2(N-1)),bd1_.resize(maxdim1(N-1)*maxdim2(N-1),maxdim1(N-1)*maxdim2(N-1)),
    id_.resize(maxdim1(N-1)*maxdim2(N-1),maxdim1(N-1)*maxdim2(N-1)),n1_.resize(maxdim1(N-1)*maxdim2(N-1),maxdim1(N-1)*maxdim2(N-1)),
    b2_.resize(maxdim1(N-1)*maxdim2(N-1),maxdim1(N-1)*maxdim2(N-1)),bd2_.resize(maxdim1(N-1)*maxdim2(N-1),maxdim1(N-1)*maxdim2(N-1)),
    n2_.resize(maxdim1(N-1)*maxdim2(N-1),maxdim1(N-1)*maxdim2(N-1));

  for(int ind = 0;ind<maxdim1(N-1);++ind){
    n1[OpKeyType<1>(ind,ind)]=T(ind*1.0);
    id1[OpKeyType<1>(ind,ind)]=T(1.0);
    if(ind!=0)
      b1[OpKeyType<1>(ind-1,ind)] = T(sqrt(ind));
    if(ind<(maxdim1(N-1)-1))
      bd1[OpKeyType<1>(ind+1,ind)] = T(sqrt(ind+1));
  }
  b1.setSite(N-1);bd1.setSite(N-1);id1.setSite(N-1);n1.setSite(N-1);
  
  for(int ind = 0;ind<maxdim2(N-1);++ind){
    n2[OpKeyType<1>(ind,ind)]=T(ind*1.0);
    id2[OpKeyType<1>(ind,ind)]=T(1.0);
    if(ind!=0)
      b2[OpKeyType<1>(ind-1,ind)] = T(sqrt(ind));
    if(ind<(maxdim2(N-1)-1))
      bd2[OpKeyType<1>(ind+1,ind)] = T(sqrt(ind+1));
  }
  b2.setSite(N-1);bd2.setSite(N-1);id2.setSite(N-1);n2.setSite(N-1);

  b1_=kronecker(b1,id2);
  bd1_=kronecker(bd1,id2);
  n1_=kronecker(n1,id2);
  b2_=kronecker(id1,b2);
  bd2_=kronecker(id1,bd2);
  n2_=kronecker(id1,n2);

  id_=kronecker(id1,id2);  

  b1_.setSite(N-1);bd1_.setSite(N-1);n1_.setSite(N-1);b2_.setSite(N-1);bd2_.setSite(N-1);id_.setSite(N-1);n2_.setSite(N-1);


  M[ukey<2>(0,0)] = id_;
  M[ukey<2>(1,0)] = b1_;
  M[ukey<2>(2,0)] = bd1_;
  M[ukey<2>(3,0)] = b2_;
  M[ukey<2>(4,0)] = bd2_;
  //  M[ukey<2>(5,0)] = U1(N-1)/2.*(n1_*(n1_-id_))+U2(N-1)/2.*(n2_*(n2_-id_))+U12(N-1)*((n1_-(0.5*id_))*(n2_-(0.5*id_)))-mu1(N-1)*n1_-mu2(N-1)*n2_;
  M[ukey<2>(5,0)] = U1(N-1)/2.*(n1_*(n1_-id_))+U2(N-1)/2.*(n2_*(n2_-id_))+U12(N-1)*(n1_*n2_)+mu1(N-1)*n1_+mu2(N-1)*n2_;

  M.Dl=6,M.Dr=1;
  M.squeeze(1e-15);
  (*this)[N-1]=M;
  for(uint site = 1;site<(N-1);++site){
    MPOMat<1,T> Mb;
    
    //clear the operators
    b1_.clear();bd1_.clear();n1_.clear();b2_.clear();bd2_.clear();n2_.clear();id_.clear();
    b1.clear();bd1.clear();n1.clear();id1.clear();b2.clear();bd2.clear();n2.clear();id2.clear();

    //resize the operators
    b1.resize(maxdim1(site),maxdim1(site)),bd1.resize(maxdim1(site),maxdim1(site)),id1.resize(maxdim1(site),maxdim1(site)),n1.resize(maxdim1(site),maxdim1(site));
    b2.resize(maxdim2(site),maxdim2(site)),bd2.resize(maxdim2(site),maxdim2(site)),id2.resize(maxdim2(site),maxdim2(site)),n2.resize(maxdim2(site),maxdim2(site));

    b1_.resize(maxdim1(site)*maxdim2(site),maxdim1(site)*maxdim2(site)),bd1_.resize(maxdim1(site)*maxdim2(site),maxdim1(site)*maxdim2(site)),
      id_.resize(maxdim1(site)*maxdim2(site),maxdim1(site)*maxdim2(site)),n1_.resize(maxdim1(site)*maxdim2(site),maxdim1(site)*maxdim2(site)),
      b2_.resize(maxdim1(site)*maxdim2(site),maxdim1(site)*maxdim2(site)),bd2_.resize(maxdim1(site)*maxdim2(site),maxdim1(site)*maxdim2(site)),
      n2_.resize(maxdim1(site)*maxdim2(site),maxdim1(site)*maxdim2(site));

    for(int ind = 0;ind<maxdim1(site);++ind){
      n1[OpKeyType<1>(ind,ind)]=T(ind*1.0);
      id1[OpKeyType<1>(ind,ind)]=T(1.0);
      if(ind!=0)
	b1[OpKeyType<1>(ind-1,ind)] = T(sqrt(ind));
      if(ind<(maxdim1(site)-1))
	bd1[OpKeyType<1>(ind+1,ind)] = T(sqrt(ind+1));
    }
    b1.setSite(site);bd1.setSite(site);id1.setSite(site);n1.setSite(site);
  
    for(int ind = 0;ind<maxdim2(site);++ind){
      n2[OpKeyType<1>(ind,ind)]=T(ind*1.0);
      id2[OpKeyType<1>(ind,ind)]=T(1.0);
      if(ind!=0)
	b2[OpKeyType<1>(ind-1,ind)] = T(sqrt(ind));
      if(ind<(maxdim2(site)-1))
	bd2[OpKeyType<1>(ind+1,ind)] = T(sqrt(ind+1));
    }
    b2.setSite(site);bd2.setSite(site);id2.setSite(site);n2.setSite(site);

    b1_=kronecker(b1,id2);
    bd1_=kronecker(bd1,id2);
    n1_=kronecker(n1,id2);
    b2_=kronecker(id1,b2);
    bd2_=kronecker(id1,bd2);
    n2_=kronecker(id1,n2);

    id_=kronecker(id1,id2);  
    b1_.setSite(site);bd1_.setSite(site);n1_.setSite(site);b2_.setSite(site);bd2_.setSite(site);id_.setSite(site);n2_.setSite(site);


    Mb[ukey<2>(0,0)] = id_;
    Mb[ukey<2>(1,0)] = b1_;
    Mb[ukey<2>(2,0)] = bd1_;
    Mb[ukey<2>(3,0)] = b2_;
    Mb[ukey<2>(4,0)] = bd2_;
    Mb[ukey<2>(5,0)] = U1(site)/2.0*(n1_*(n1_-id_))+U2(site)/2.*(n2_*(n2_-id_))+U12(site)*(n1_*n2_)+mu1(site)*n1_+mu2(site)*n2_;

    Mb[ukey<2>(5,1)] = -J1(site)*bd1_;
    Mb[ukey<2>(5,2)] = -J1(site)*b1_;
    Mb[ukey<2>(5,3)] = -J2(site)*bd2_;
    Mb[ukey<2>(5,4)] = -J2(site)*b2_;
    Mb[ukey<2>(5,5)] = id_;

    Mb.Dl=6,Mb.Dr=6;
    Mb.squeeze(1e-15);
    (*this)[site]=Mb;
  }
}






//this is the siam mpo from iztok's nrg;
double t_n(const double lambda, const int n, bool rescale)
{
  double x = ( (1+1./lambda)*(1-pow(lambda, -n-1)) )/( 2*sqrt(1- pow(lambda,-2*n-1) )*sqrt(1-pow(lambda,-2*n-3) ) );
  if (!rescale)
    x *= pow(lambda, -n/2.);
  return x;
}


template<typename T>
class siammpo: public MPO<T>
{
public:
  siammpo(const unsigned int n, const double u, const double xi, const double ef, const double lambda, const bool rescale,const bool minus=false)
  {
    // first generate the single basis
    SparseOperator<1,T> O_11(4,4),O_mZ(4,4), O_pZ(4,4), O_Zm(4,4), O_Zp(4,4), O_1m(4,4), O_m1(4,4), O_1p(4,4), O_p1(4,4);
    //O_11 = { ((0,0),(0,0)):1.0, ((1,1),(1,1)):1.0, ((0,1),(0,1)):1.0, ((1,0),(1,0)):1.0}
    O_11[ OpKeyType<1>(0,0)] = T(1.0); //(0,0)=state 0
    O_11[ OpKeyType<1>(1,1)] = T(1.0); //(0,1)=state 1
    O_11[ OpKeyType<1>(2,2)] = T(1.0); //(1,0)=state 2
    O_11[ OpKeyType<1>(3,3)] = T(1.0); //(1,1)= state 3

    //c-up dagger
    //O_mZ = { ((1,0),(0,0)):1.0, ((1,1),(0,1)):-1.0}
    O_mZ[ OpKeyType<1>(2,0)] = T(1.0);
    O_mZ[ OpKeyType<1>(3,1)] = T(-1.0);

    //c-up
    //O_pZ = { ((0,0),(1,0)):1.0, ((0,1),(1,1)):-1.0}
    O_pZ[ OpKeyType<1>(0,2)] = T(1.0);
    O_pZ[ OpKeyType<1>(1,3)] = T(-1.0);

    //c-down dagger
    //O_Zm = { ((0,1),(0,0)):1.0, ((1,1),(1,0)):-1.0}
    O_Zm[ OpKeyType<1>(1,0)] = T(1.0);
    O_Zm[ OpKeyType<1>(3,2)] = T(-1.0);

    //c-down
    //O_Zp = { ((0,0),(0,1)):1.0, ((1,0),(1,1)):-1.0}
    O_Zp[ OpKeyType<1>(0,1)] = T(1.0);   
    O_Zp[ OpKeyType<1>(2,3)] = T(-1.0);
    
    //c-down dagger
    //O_1m = { ((0,1),(0,0)):1.0, ((1,1),(1,0)):1.0}
    O_1m[ OpKeyType<1>(1,0)] =  T(1.0);
    O_1m[ OpKeyType<1>(3,2)] =  T(1.0);

    //c-up dagger
    //O_m1 = { ((1,0),(0,0)):1.0, ((1,1),(0,1)):1.0}
    O_m1[ OpKeyType<1>(2,0)] =  T(1.0);
    O_m1[ OpKeyType<1>(3,1)] =  T(1.0);

    //c-down
    //O_1p = { ((0,0),(0,1)):1.0, ((1,0),(1,1)):1.0}
    O_1p[ OpKeyType<1>(0,1)] =  T(1.0);
    O_1p[ OpKeyType<1>(2,3)] =  T(1.0);

    //c-up
    //O_p1 = { ((0,0),(1,0)):1.0, ((0,1),(1,1)):1.0}
    O_p1[ OpKeyType<1>(0,2)] =  T(1.0);
    O_p1[ OpKeyType<1>(1,3)] =  T(1.0);

    SparseOperator<1,T> H_imp(4,4);
    H_imp[ OpKeyType<1>(1,1)] = ef;
    H_imp[ OpKeyType<1>(3,3)] = (2*ef + u);
    H_imp[ OpKeyType<1>(2,2)] = ef;
        
    if (rescale)
      H_imp *= (1.0/lambda);
    if(minus)
      H_imp*=(-1.0);
    
    if (n==0) return;
    if (n==1) {
      MPOMat<1,T> A;
      A[ukey<2>(0,0)] = H_imp;
      A.Dl=1;A.Dr=1;
      (*this)[0] = A;
      return;
    }
    const double sql = (rescale ? sqrt(lambda) : 1.0);
    //std::cout  << "sql: " << sql << std::endl;
    {
      MPOMat<1,T> A;
      A[ukey<2>(0,0)] = O_11;
      A[ukey<2>(0,1)] = H_imp;
      double xfac = sqrt(xi/M_PI);
      //std::cout << "xfac: " << xfac << std::endl;
      if (rescale)
	xfac *= (1./sql);
      //std::cout << "xfac--> " << xfac << std::endl;
      if(minus)
	xfac*=(-1.0);
      A[ukey<2>(0,2)] = xfac * O_mZ;
      A[ukey<2>(0,3)] = xfac * O_pZ;
      A[ukey<2>(0,4)] = xfac * O_1m;
      A[ukey<2>(0,5)] = xfac * O_1p;
      A.Dl=1;A.Dr=6;
      (*this)[0] = A;
    }
    for(int j=1;j<n-1;++j) {
      MPOMat<1,T> A;
      //const double bx = pow(lambda, (j-1.)/2.0);
      double xfac = t_n(lambda, j-1, rescale);
      if(minus){
      	xfac*=(-1.0);
      }
      A[ukey<2>(0,0)] = O_11;
      A[ukey<2>(1,1)] = sql*O_11;
      A[ukey<2>(0,2)] = xfac*O_mZ;
      A[ukey<2>(0,3)] = xfac*O_pZ;
      A[ukey<2>(0,4)] = xfac*O_1m;
      A[ukey<2>(0,5)] = xfac*O_1p;
      //
      A[ukey<2>(2,1)] = O_p1;
      A[ukey<2>(3,1)] = O_m1;
      A[ukey<2>(4,1)] = O_Zp;
      A[ukey<2>(5,1)] = O_Zm;
      A.Dl=6;A.Dr=6;
      (*this)[j] = A;
    }
    {
      MPOMat<1,T> A;     
      A[ukey<2>(1,0)] = sql*O_11;
      A[ukey<2>(2,0)] = O_p1;
      A[ukey<2>(3,0)] = O_m1;
      A[ukey<2>(4,0)] = O_Zp;
      A[ukey<2>(5,0)] = O_Zm;
      A.Dl=6;A.Dr=1;
      (*this)[n-1] = A;
    }
        
  }
  //  template<typename T>
  void shiftandrescale(const T E0,const T a, const T W){
    SparseLocalOperator<T> id(4,4);
    id[OpKeyType<1>(0,0)]=T(1.);id[OpKeyType<1>(1,1)]=T(1.);id[OpKeyType<1>(2,2)]=T(1.);id[OpKeyType<1>(3,3)]=T(1.);
    id*=E0;
    (*this)[0][ukey<2>(0,0)]=(-1.0)*id;
    typename MPOMat<1,T>::iterator o;
    for(o=(*this)[0].begin();o!=(*this)[0].end();++o)
      o->second/=a;
    id[OpKeyType<1>(0,0)]=T(1.);id[OpKeyType<1>(1,1)]=T(1.);id[OpKeyType<1>(2,2)]=T(1.);id[OpKeyType<1>(3,3)]=T(1.);
    id*=W;
    (*this)[0][ukey<2>(0,0)]-=id;
    //that's it ...
  }

};

class Sk_dagger:public MPO<Complex>{
public:
  Sk_dagger(){};
  //k is the quantum number, it has integer values
  Sk_dagger(uint N,Real k,bool obc=true){
    std::vector<Complex> coeff(N);
    if(obc==true){
      for(uint s=0;s<N;++s){
	coeff[s]=1.0/sqrt(2*(N+1))*(exp(Complex(0,-k*1.0*M_PI/(1.0*(N+1))*(s+1)))-exp(Complex(0,k*1.0*M_PI/(1.0*(N+1))*(s+1))));
      }
    }else{
      for(uint s=0;s<N;++s)
	coeff[s]=1.0/sqrt(2*N)*exp(Complex(0,-k*1.0*M_PI/(1.0*(N+1))*s))+exp(Complex(0,k*1.0*M_PI/(1.0*(N+1))*s));
    }
    MPOMat<1,Complex> mat;
    SparseOperator<1,Complex> O_11,Sp;
    Sp[OpKeyType<1>(1,0)]=Complex(1.);
    O_11[OpKeyType<1>(0,0)]=Complex(1.);
    O_11[OpKeyType<1>(1,1)]=Complex(1.);
    
    mat[ukey<2>(0,0)] = coeff[0]*Sp;
    mat[ukey<2>(0,1)] = O_11;
    mat.Dl=1;mat.Dr=2;
    (*this)[0]=mat;
    mat.clear();

    for(uint s=1;s<(N-1);++s){
      mat[ukey<2>(0,0)] = O_11;
      mat[ukey<2>(1,0)] = coeff[s]*Sp;
      mat[ukey<2>(1,1)] = O_11;
      mat.Dl=2;mat.Dr=2;
      (*this)[s]=mat;
      mat.clear();
    }

    mat[ukey<2>(0,0)] = O_11;
    mat[ukey<2>(1,0)] = coeff[N-1]*Sp;
    mat.Dl=2;mat.Dr=1;
    (*this)[N-1]=mat;
  }
private:
  
};

class ck_dagger:public MPO<Complex>{
public:
  ck_dagger(){};
  //k is the quantum number, it has integer values
  ck_dagger(uint N,Real k,bool obc=true){
    std::vector<Complex> coeff(N);
    if(obc==true){
      for(uint s=0;s<N;++s){
	coeff[s]=1.0/sqrt(2*(N+1))*(exp(Complex(0,-k*1.0*M_PI/(1.0*(N+1))*(s+1)))-exp(Complex(0,k*1.0*M_PI/(1.0*(N+1))*(s+1))));
      }
    }else{
      for(uint s=0;s<N;++s)
	coeff[s]=1.0/sqrt(2*N)*exp(Complex(0,-k*1.0*M_PI/(1.0*(N+1))*s))+exp(Complex(0,k*1.0*M_PI/(1.0*(N+1))*s));
    }
    MPOMat<1,Complex> mat;
    SparseOperator<1,Complex> p(2,2),cd(2,2),O_11(2,2);
    cd[OpKeyType<1>(1,0)]=Complex(1.);
    p[OpKeyType<1>(0,0)]=Complex(1.);
    p[OpKeyType<1>(1,1)]=Complex(-1.);
    O_11[OpKeyType<1>(0,0)]=Complex(1.);
    O_11[OpKeyType<1>(1,1)]=Complex(1.);
    

    mat[ukey<2>(0,0)] = coeff[0]*cd;
    mat[ukey<2>(0,1)] = p;
    mat.Dl=1;mat.Dr=2;
    (*this)[0]=mat;
    mat.clear();

    for(uint s=1;s<(N-1);++s){
      mat[ukey<2>(0,0)] = O_11;
      mat[ukey<2>(1,0)] = coeff[s]*cd;
      mat[ukey<2>(1,1)] = p;
      mat.Dl=2;mat.Dr=2;
      (*this)[s]=mat;
      mat.clear();
    }

    mat[ukey<2>(0,0)] = O_11;
    mat[ukey<2>(1,0)] = coeff[N-1]*cd;
    mat.Dl=2;mat.Dr=1;
    (*this)[N-1]=mat;
  }
private:
  
};

class ck_daggerReal:public MPO<Real>{
public:
  ck_daggerReal(){};
  //k is the quantum number, it has integer values
  ck_daggerReal(uint N,Real k,bool obc=true){
    std::vector<Real> coeff(N);
    if(obc==true){
      for(uint s=0;s<N;++s){
	coeff[s]=1.0/sqrt(2*N)*2*(sin(k*1.0*M_PI/(1.0*(N+1))*(s+1)));
      }
    }else{
      for(uint s=0;s<N;++s)
	coeff[s]=1.0/sqrt(2*N)*2*(cos(k*1.0*M_PI/(1.0*(N+1))*(s+1)));
    }
    MPOMat<1,Real> mat;
    SparseOperator<1,Real> p(2,2),cd(2,2),O_11(2,2);
    cd[OpKeyType<1>(1,0)]=1.;
    p[OpKeyType<1>(0,0)]=1.;
    p[OpKeyType<1>(1,1)]=-1.;
    O_11[OpKeyType<1>(0,0)]=1.;
    O_11[OpKeyType<1>(1,1)]=1.;
    

    mat[ukey<2>(0,0)] = coeff[0]*cd;
    mat[ukey<2>(0,1)] = p;
    mat.Dl=1;mat.Dr=2;
    (*this)[0]=mat;
    mat.clear();

    for(uint s=1;s<(N-1);++s){
      mat[ukey<2>(0,0)] = O_11;
      mat[ukey<2>(1,0)] = coeff[s]*cd;
      mat[ukey<2>(1,1)] = p;
      mat.Dl=2;mat.Dr=2;
      (*this)[s]=mat;
      mat.clear();
    }

    mat[ukey<2>(0,0)] = O_11;
    mat[ukey<2>(1,0)] = coeff[N-1]*cd;
    mat.Dl=2;mat.Dr=1;
    (*this)[N-1]=mat;
  }
private:
  
};


class ck:public MPO<Complex>{
public:
  ck(){};
  //k is the quantum number, it has integer values
  ck(uint N,Real k,bool obc=true){
    std::vector<Complex> coeff(N);
    if(obc==true){
      for(uint s=0;s<N;++s){
	coeff[s]=1.0/sqrt(2*(N+1))*(exp(Complex(0,k*1.0*M_PI/(1.0*(N+1))*(s+1)))-exp(Complex(0,-k*1.0*M_PI/(1.0*(N+1))*(s+1))));
      }
    }else{
      for(uint s=0;s<N;++s)
	coeff[s]=1.0/sqrt(2*N)*exp(Complex(0,+k*1.0*M_PI/(1.0*(N+1))*s))+exp(Complex(0,-k*1.0*M_PI/(1.0*(N+1))*s));
    }
    MPOMat<1,Complex> mat;
    SparseOperator<1,Complex> p(2,2),c(2,2),O_11(2,2);
    c[OpKeyType<1>(0,1)]=Complex(1.);
    p[OpKeyType<1>(0,0)]=Complex(1.);
    p[OpKeyType<1>(1,1)]=Complex(-1.);
    O_11[OpKeyType<1>(0,0)]=Complex(1.);
    O_11[OpKeyType<1>(1,1)]=Complex(1.);
    

    mat[ukey<2>(0,0)] = coeff[0]*c;
    mat[ukey<2>(0,1)] = p;
    mat.Dl=1;mat.Dr=2;
    (*this)[0]=mat;
    mat.clear();

    for(uint s=1;s<(N-1);++s){
      mat[ukey<2>(0,0)] = O_11;
      mat[ukey<2>(1,0)] = coeff[s]*c;
      mat[ukey<2>(1,1)] = p;
      mat.Dl=2;mat.Dr=2;
      (*this)[s]=mat;
      mat.clear();
    }

    mat[ukey<2>(0,0)] = O_11;
    mat[ukey<2>(1,0)] = coeff[N-1]*c;
    mat.Dl=2;mat.Dr=1;
    (*this)[N-1]=mat;
  }
private:
  
};


class ckDagIntMvImp:public MPO<Complex>{
public:
  ckDagIntMvImp(){};
  //k is the quantum number, it has integer values
  ckDagIntMvImp(uint N,Real k, uint type){
    std::vector<Complex> coeff(N);
    double sign=-1.;
    /*if(obc==true){
      for(uint s=0;s<N;++s){
      coeff[s]=1.0/sqrt(2*N)*(exp(Complex(0,-k*1.0*M_PI/(1.0*(N+1))*(s+1)))-exp(Complex(0,k*1.0*M_PI/(1.0*(N+1))*(s+1))));
      }
      }else{
      for(uint s=0;s<N;++s)
      coeff[s]=1.0/sqrt(2*N)*exp(Complex(0,-k*1.0*M_PI/(1.0*(N+1))*s))+exp(Complex(0,k*1.0*M_PI/(1.0*(N+1))*s));
      }*/
    for(uint s=0;s<N;++s){
    	coeff[s]=1.0/sqrt(N)*(exp(Complex(0,k*M_PI/(N+1.)*(s+1.))));
    	//coeff[s]=1.0/sqrt(N)*(exp(Complex(0,k*M_PI/(N)*(s))));
    	//coeff[s]=1.0/sqrt(2.*N)*(exp(Complex(0,k*M_PI/(N+1.)*(s+1.)))-exp(Complex(0,-k*M_PI/(N+1.)*(s+1.))));
    }
    MPOMat<1,Complex> mat;

    SparseLocalOperator<Complex> cD(4,4),id(4,4),p(4,4);

    if (type==1){
		cD[OpKeyType<1>(1,0)]=1.; cD[OpKeyType<1>(3,2)]=1.;
    }else{
    	if (type==2){
    		cD[OpKeyType<1>(2,0)]=1.; cD[OpKeyType<1>(3,1)]=sign;
    	}else{
    		throw std::runtime_error("ckIntMvImp error");
    	}
    }

	id[OpKeyType<1>(0,0)]=1.;id[OpKeyType<1>(1,1)]=1.; id[OpKeyType<1>(2,2)]=1.;id[OpKeyType<1>(3,3)]=1.;
	p[OpKeyType<1>(0,0)]=1.; p[OpKeyType<1>(1,1)]=sign; p[OpKeyType<1>(2,2)]=sign; p[OpKeyType<1>(3,3)]=1.;

    cD.setSite(0); p.setSite(0); id.setSite(0);
    mat[ukey<2>(0,0)] = coeff[0]*cD; //p
    mat[ukey<2>(0,1)] = id;//p;
    mat.Dl=1;mat.Dr=2;
    (*this)[0]=mat;
    mat.clear();

    for(uint s=1;s<(N-1);++s){
      cD.setSite(s); p.setSite(s); id.setSite(s);
      mat[ukey<2>(0,0)] = p;//id;
      mat[ukey<2>(1,0)] = coeff[s]*cD;
      mat[ukey<2>(1,1)] = id;//p;
      mat.Dl=2;mat.Dr=2;
      (*this)[s]=mat;
      mat.clear();
    }

    cD.setSite(N-1); p.setSite(N-1); id.setSite(N-1);
    mat[ukey<2>(0,0)] = p;//id;
    mat[ukey<2>(1,0)] = coeff[N-1]*cD;
    mat.Dl=2;mat.Dr=1;
    (*this)[N-1]=mat;
  }
private:

};


class ckIntMvImp:public MPO<Complex>{
public:
  ckIntMvImp(){};
  //k is the quantum number, it has integer values
  //type 1: majority, 2 minority
  ckIntMvImp(uint N, Real k, uint type){
    std::vector<Complex> coeff(N);
    double sign=-1.;

    for(uint s=0; s<N; ++s){
    	coeff[s]=1.0/sqrt(N)*(exp(Complex(0,-k*M_PI/(N+1.)*(s+1.))));
    }
    MPOMat<1,Complex> mat;

    SparseLocalOperator<Complex> c(4,4),id(4,4),p(4,4);
    if (type==1){
        c[OpKeyType<1>(0,1)]=(1.); c[OpKeyType<1>(2,3)]=(1.);
    }else{
    	if (type==2){
    		c[OpKeyType<1>(0,2)]=(1.); c[OpKeyType<1>(1,3)]=(sign);
    	}else{
    		throw std::runtime_error("ckIntMvImp error");
    	}
    }

    id[OpKeyType<1>(0,0)]=1.;id[OpKeyType<1>(1,1)]=1.; id[OpKeyType<1>(2,2)]=1.;id[OpKeyType<1>(3,3)]=1.;
	p[OpKeyType<1>(0,0)]=1.; p[OpKeyType<1>(1,1)]=sign; p[OpKeyType<1>(2,2)]=sign; p[OpKeyType<1>(3,3)]=1.;

    c.setSite(0); p.setSite(0); id.setSite(0);
    mat[ukey<2>(0,0)] = coeff[0]*c; //p
    mat[ukey<2>(0,1)] = id;//p;
    mat.Dl=1;mat.Dr=2;
    (*this)[0]=mat;
    mat.clear();

    for(uint s=1;s<(N-1);++s){
      c.setSite(s); p.setSite(s); id.setSite(s);
      mat[ukey<2>(0,0)] = p;//id;
      mat[ukey<2>(1,0)] = coeff[s]*c;
      mat[ukey<2>(1,1)] = id;//p;
      mat.Dl=2;mat.Dr=2;
      (*this)[s]=mat;
      mat.clear();
    }

    c.setSite(N-1); p.setSite(N-1); id.setSite(N-1);
    mat[ukey<2>(0,0)] = p;//id;
    mat[ukey<2>(1,0)] = coeff[N-1]*c;
    mat.Dl=2;mat.Dr=1;
    (*this)[N-1]=mat;
  }
private:

};
class ckDagBHMvImp:public MPO<Complex>{
public:
	ckDagBHMvImp(){};
  //k is the quantum number, it has integer values
	ckDagBHMvImp(uint N,Real k, uint maxdim){
    std::vector<Complex> coeff(N);

    for(uint s=0;s<N;++s){
    	coeff[s]=1.0/sqrt(N)*(exp(Complex(0,k*M_PI/(N+1.)*(s+1.))));
    }
    MPOMat<1,Complex> mat;

    SparseLocalOperator<Complex> cd(2*maxdim,2*maxdim),id(2*maxdim,2*maxdim);
	for (uint ind=0; ind<maxdim;++ind){
		id[OpKeyType<1>(ind,ind)]=(1.0);
		id[OpKeyType<1>(ind+maxdim,ind+maxdim)]=(1.0);
		cd[OpKeyType<1>(ind+maxdim,ind)]=(1.);
	}

    cd.setSite(0); id.setSite(0);
    mat[ukey<2>(0,0)] = coeff[0]*cd;
    mat[ukey<2>(0,1)] = id;
    mat.Dl=1;mat.Dr=2;
    (*this)[0]=mat;
    mat.clear();

    for(uint s=1;s<(N-1);++s){
      cd.setSite(s); id.setSite(s);
      mat[ukey<2>(0,0)] = id;
      mat[ukey<2>(1,0)] = coeff[s]*cd;
      mat[ukey<2>(1,1)] = id;
      mat.Dl=2;mat.Dr=2;
      (*this)[s]=mat;
      mat.clear();
    }

    cd.setSite(N-1); id.setSite(N-1);
    mat[ukey<2>(0,0)] = id;
    mat[ukey<2>(1,0)] = coeff[N-1]*cd;
    mat.Dl=2;mat.Dr=1;
    (*this)[N-1]=mat;
  }
private:

};

class ckNDagBHMvImp:public MPO<Complex>{
public:
	ckNDagBHMvImp(){};
  //k is the quantum number, it has integer values
	ckNDagBHMvImp(uint N,Real k, uint maxdim,bool create,bool majority){
    std::vector<Complex> coeff(N);

    double s=create?1.:-1.;
    for(uint s=0;s<N;++s){
    	coeff[s]=1.0/sqrt(N)*(exp(Complex(0,s*k*M_PI/(N+1.)*(s+1.))));
    }
    MPOMat<1,Complex> mat;

    SparseLocalOperator<Complex> cnd(2*maxdim,2*maxdim),id(2*maxdim,2*maxdim);
	for (uint ind=0; ind<maxdim;++ind){
		id[OpKeyType<1>(ind,ind)]=(1.0);
		id[OpKeyType<1>(ind+maxdim,ind+maxdim)]=(1.0);
	}
	if (create){
		if (majority){
			for(int ind = 0;ind<maxdim-1;++ind){
				cnd[OpKeyType<1>(ind+1,ind)] = (sqrt(ind+1));
				cnd[OpKeyType<1>(ind+1+maxdim,ind+maxdim)] = (sqrt(ind+1));
			}
		}else{
			for (uint ind=0; ind<maxdim;++ind){
				cnd[OpKeyType<1>(ind+maxdim,ind)]=(1.);
			}
		}
	}else{
		if (majority){
			for(int ind = 1;ind<maxdim-1;++ind){
				cnd[OpKeyType<1>(ind-1,ind)] = (sqrt(ind));
				cnd[OpKeyType<1>(ind-1+maxdim,ind+maxdim)] = (sqrt(ind));
			}
		}else{
			for (uint ind=0; ind<maxdim;++ind){
				cnd[OpKeyType<1>(ind,ind+maxdim)]=(1.);
			}
		}
	}

    cnd.setSite(0); id.setSite(0);
    mat[ukey<2>(0,0)] = coeff[0]*cnd;
    mat[ukey<2>(0,1)] = id;
    mat.Dl=1;mat.Dr=2;
    (*this)[0]=mat;
    mat.clear();

    for(uint s=1;s<(N-1);++s){
      cnd.setSite(s); id.setSite(s);
      mat[ukey<2>(0,0)] = id;
      mat[ukey<2>(1,0)] = coeff[s]*cnd;
      mat[ukey<2>(1,1)] = id;
      mat.Dl=2;mat.Dr=2;
      (*this)[s]=mat;
      mat.clear();
    }

    cnd.setSite(N-1); id.setSite(N-1);
    mat[ukey<2>(0,0)] = id;
    mat[ukey<2>(1,0)] = coeff[N-1]*cnd;
    mat.Dl=2;mat.Dr=1;
    (*this)[N-1]=mat;
  }
private:

};


class nk:public MPO<Complex>{
public:
  nk(){};
  //k is the quantum number, it has integer values
  nk(uint N,Real k){
    std::vector<Complex> coeff(N);
    //that's the ck^dagger coefficient ...
    for(uint s=0;s<N;++s){
      coeff[s]=1.0/sqrt(2*(N+1))*(exp(Complex(0,-k*1.0*M_PI/(1.0*(N+1))*(s+1)))-exp(Complex(0,k*1.0*M_PI/(1.0*(N+1))*(s+1))));
    }
    
    MPOMat<1,Complex> mat;
    SparseLocalOperator<Complex> p(2,2),cd(2,2),c(2,2),O_11(2,2),n(2,2);
    cd[OpKeyType<1>(1,0)]=Complex(1.);
    c[OpKeyType<1>(0,1)]=Complex(1.);
    p[OpKeyType<1>(0,0)]=Complex(1.);
    p[OpKeyType<1>(1,1)]=Complex(-1.);
    O_11[OpKeyType<1>(0,0)]=Complex(1.);
    O_11[OpKeyType<1>(1,1)]=Complex(1.);
    n[OpKeyType<1>(1,1)]=Complex(1.);    
    
    n.setSite(0);
    p.setSite(0);
    O_11.setSite(0);
    cd.setSite(0);
    c.setSite(0);
    
    mat[ukey<2>(0,0)] = pow(abs(coeff[0]),2)*n;
    mat[ukey<2>(0,1)] = coeff[0]*cd;
    mat[ukey<2>(0,2)] = conj(coeff[0])*c;
    mat[ukey<2>(0,3)] = O_11;
    mat.Dl=1;mat.Dr=4;
    (*this)[0]=mat;
    mat.clear();

    for(uint s=1;s<(N-1);++s){
      n.setSite(s);
      p.setSite(s);
      O_11.setSite(s);
      cd.setSite(s);
      c.setSite(s);

      mat[ukey<2>(0,0)] = O_11;
      mat[ukey<2>(1,1)] = p;
      mat[ukey<2>(2,2)] = p;
      mat[ukey<2>(3,3)] = O_11;

      mat[ukey<2>(1,0)] = conj(coeff[s])*c;
      mat[ukey<2>(2,0)] = coeff[s]*cd;
      mat[ukey<2>(3,0)] = pow(abs(coeff[s]),2)*n;

      mat[ukey<2>(3,1)] = coeff[s]*cd;
      mat[ukey<2>(3,2)] = conj(coeff[s])*c;

      mat.Dl=4;mat.Dr=4;
      (*this)[s]=mat;
      mat.clear();
    }
    uint s=N-1;
    n.setSite(s);
    p.setSite(s);
    O_11.setSite(s);
    cd.setSite(s);
    c.setSite(s);

    mat[ukey<2>(0,0)] = O_11;
    mat[ukey<2>(1,0)] = conj(coeff[s])*c;
    mat[ukey<2>(2,0)] = coeff[s]*cd;
    mat[ukey<2>(3,0)] = pow(abs(coeff[s]),2)*n;
    mat.Dl=4;mat.Dr=1;
    (*this)[N-1]=mat;
  }
private:
  
};


//is spin > 0, an up spin gaussian is created;
//else, its a downspin gaussian
//k0 is an integer number quantising the momentum of the operator
//N is the chain length, x0 the mean value in real space, sigma its width in k-space
class gaussMPOfermihub: public MPO<Complex> {
public:
  gaussMPOfermihub():MPO<Complex>(){};
  gaussMPOfermihub(uint N,Real x0,Real k0,Real sigma,int spin){
    CVectorType coeff(N);
    for(uint n=0;n<N;++n){
      coeff(n)=abs(exp(Complex(-2*(pow(sigma,2)*pow(M_PI,2)*pow(n- x0,2)),M_PI*k0/N*(n - x0))) )>1e-14?
	exp(Complex(-2*(pow(sigma,2)*pow(M_PI,2)*pow(n- x0,2)),M_PI*k0/N*(n - x0)) ):Complex(0.0,0.0);
    }
    MPOMat<1,Complex> mat;

    SparseLocalOperator<Complex> cr(4,4),cuD(4,4),cdD(4,4),O_11(4,4),p(4,4);//all statistics are set to bosonic (_statistics=1)

    cuD[OpKeyType<1>(2,0)] = Complex(1.0);cuD[OpKeyType<1>(3,1)] = Complex(1.0);

    cdD[OpKeyType<1>(1,0)] = Complex(1.0);cdD[OpKeyType<1>(3,2)] = Complex(-1.0);

    O_11[OpKeyType<1>(0,0)]=Complex(1.0);O_11[OpKeyType<1>(1,1)]=Complex(1.0);O_11[OpKeyType<1>(2,2)]=Complex(1.0);O_11[OpKeyType<1>(3,3)]=Complex(1.0);

    p[OpKeyType<1>(0,0)]=Complex(1.0);p[OpKeyType<1>(1,1)]=Complex(-1.0);p[OpKeyType<1>(2,2)]=Complex(-1.0);p[OpKeyType<1>(3,3)]=Complex(1.0);
    if(spin>0)
      cr=cuD;
    else if(spin<=0)
      cr=cdD;
    p.setSite(0);cr.setSite(0);O_11.setSite(0);
    //you may wonder about this, but if the coefficients get too small (much smaller than machine precision), then, after application of this operator to
    //the GS at halffilling (which btw works fine), doing an orthonormalization gives pretty strange results. I THINK that this is an degeneracy issue.
    //this works now: if the non diagonal entries are very very small, the mpo matrix is 1x1. if the coeff is larger than a threshold, the usual
    //representation is used ... i think that before, i had multiple superpositions of the state with itself, which led to problems with the svd in orthonormalize
    bool first=true,last=true;
    for(uint s=0;s<(N);++s){
      p.setSite(s);cr.setSite(s);O_11.setSite(s);
      if(abs(coeff(s))>1e-14){
	if(first){
	  mat[ukey<2>(0,0)] = coeff(s)*cr;
	  mat[ukey<2>(0,1)] = p;
	  mat.Dl=1;mat.Dr=2;
	  (*this)[s]=mat;
	  mat.clear();
	  first=false;
	}else{
	  mat[ukey<2>(0,0)] = O_11;
	  mat[ukey<2>(1,0)] = coeff(s)*cr;
	  mat[ukey<2>(1,1)] = p;
	  mat.Dl=2;mat.Dr=2;
	  (*this)[s]=mat;
	  mat.clear();
	}
      }else{
	if(first){
	  mat[ukey<2>(0,0)] = p;
	  mat.Dl=1;mat.Dr=1;
	  (*this)[s]=mat;
	  mat.clear();
	}
	if(!first&&last){
	  mat[ukey<2>(0,0)] = O_11;
	  mat[ukey<2>(1,0)] = coeff(s-1)*cr;
	  mat.Dl=2;mat.Dr=1;
	  (*this)[s-1]=mat;
	  mat.clear();
	  last=false;
	}
	if(!first&&!last){
	  mat[ukey<2>(0,0)] = O_11;
	  mat.Dl=1;mat.Dr=1;
	  (*this)[s]=mat;
	  mat.clear();
	}
      }
    }
  }
  virtual ~gaussMPOfermihub(){};
};

//that's a gauss operator for the fermi hubbard ladder
//with "spin" you control the type of operator which is used in the gaussian:
// 2: create up-spin on upper chain
//-2: create down-spin in upper chain
// 1: create up-spin on lower chain
//-1: create down-spin on lower chain
//k0 is an integer number quantising the momentum of the operator
//N is the chain length, x0 the mean value in real space, sigma its width in k-space
class gaussMPOfermihubladder: public MPO<Complex> {
public:
  gaussMPOfermihubladder():MPO<Complex>(){};
  gaussMPOfermihubladder(uint N,Real x0,Real k0,Real sigma,int spin){
    CVectorType coeff(N);
    for(uint n=0;n<N;++n){
      coeff(n)=abs( exp(Complex(-2*(pow(sigma,2)*pow(M_PI,2)*pow(n- x0,2)),M_PI*k0/N*(n - x0))) )>1e-14?
	exp(Complex(-2*(pow(sigma,2)*pow(M_PI,2)*pow(n- x0,2)),M_PI*k0/N*(n - x0)) ):Complex(0.0,0.0);
    }

    MPOMat<1,Complex> mat;
    SparseLocalOperator<Complex> cuDU(16,16),cdDU(16,16),cuU(16,16),cdU(16,16),O_11(16,16),p(16,16),pU(16,16),pL(16,16);
    SparseLocalOperator<Complex> cuDL(16,16),cdDL(16,16),cuL(16,16),cdL(16,16);

    cuDL[OpKeyType<1>(2,0)] = 1.0;cuDL[OpKeyType<1>(3,1)] = 1.0;cuDL[OpKeyType<1>(6,4)] = 1.0;cuDL[OpKeyType<1>(7,5)] = 1.0;
    cuDL[OpKeyType<1>(10,8)] = 1.0;cuDL[OpKeyType<1>(11,9)] = 1.0;cuDL[OpKeyType<1>(14,12)] = 1.0;cuDL[OpKeyType<1>(15,13)] = 1.0;
 
    cuL[OpKeyType<1>(0,2)] = 1.0;cuL[OpKeyType<1>(1,3)] = 1.0;cuL[OpKeyType<1>(4,6)] = 1.0;cuL[OpKeyType<1>(5,7)] = 1.0;
    cuL[OpKeyType<1>(8,10)] = 1.0;cuL[OpKeyType<1>(9,11)] = 1.0;cuL[OpKeyType<1>(12,14)] = 1.0;cuL[OpKeyType<1>(13,15)] = 1.0;
								    
    cdL[OpKeyType<1>(0,1)] = 1.0;cdL[OpKeyType<1>(2,3)] = -1.0;cdL[OpKeyType<1>(4,5)] = 1.0;cdL[OpKeyType<1>(6,7)] = -1.0;
    cdL[OpKeyType<1>(8,9)] = 1.0;cdL[OpKeyType<1>(10,11)] = -1.0;cdL[OpKeyType<1>(12,13)] = 1.0;cdL[OpKeyType<1>(14,15)] = -1.0;

    cdDL[OpKeyType<1>(1,0)] = 1.0;cdDL[OpKeyType<1>(3,2)] = -1.0;cdDL[OpKeyType<1>(5,4)] = 1.0;cdDL[OpKeyType<1>(7,6)] = -1.0;
    cdDL[OpKeyType<1>(9,8)] = 1.0;cdDL[OpKeyType<1>(11,10)] = -1.0;cdDL[OpKeyType<1>(13,12)] = 1.0;cdDL[OpKeyType<1>(15,14)] = -1.0;

    //the Upper chain
    cuDU[OpKeyType<1>(8,0)] = 1.0;cuDU[OpKeyType<1>(9,1)] = 1.0;cuDU[OpKeyType<1>(10,2)] = 1.0;cuDU[OpKeyType<1>(11,3)] = 1.0;
    cuDU[OpKeyType<1>(12,4)] = 1.0;cuDU[OpKeyType<1>(13,5)] = 1.0;cuDU[OpKeyType<1>(14,6)] = 1.0;cuDU[OpKeyType<1>(15,7)] = 1.0;
     
    cuU[OpKeyType<1>(0,8)] = 1.0;cuU[OpKeyType<1>(1,9)] = 1.0;cuU[OpKeyType<1>(2,10)] = 1.0;cuU[OpKeyType<1>(3,11)] = 1.0;
    cuU[OpKeyType<1>(4,12)] = 1.0;cuU[OpKeyType<1>(5,13)] = 1.0;cuU[OpKeyType<1>(6,14)] = 1.0;cuU[OpKeyType<1>(7,15)] = 1.0;
     
    cdDU[OpKeyType<1>(4,0)] = 1.0;cdDU[OpKeyType<1>(5,1)] = 1.0;cdDU[OpKeyType<1>(6,2)] = 1.0;cdDU[OpKeyType<1>(7,3)] = 1.0;
    cdDU[OpKeyType<1>(12,8)] = -1.0;cdDU[OpKeyType<1>(13,9)] = -1.0;cdDU[OpKeyType<1>(14,10)] = -1.0;cdDU[OpKeyType<1>(15,11)] = -1.0;

    cdU[OpKeyType<1>(0,4)] = 1.0;cdU[OpKeyType<1>(1,5)] = 1.0;cdU[OpKeyType<1>(2,6)] = 1.0;cdU[OpKeyType<1>(3,7)] = 1.0;
    cdU[OpKeyType<1>(8,12)] = -1.0;cdU[OpKeyType<1>(9,13)] = -1.0;cdU[OpKeyType<1>(10,14)] = -1.0;cdU[OpKeyType<1>(11,15)] = -1.0;

    SparseLocalOperator<Complex> bla(16,16),cr(16,16);

  
    for(uint ind=0;ind<16;++ind)
      O_11[OpKeyType<1>(ind,ind)]=1;

    p[OpKeyType<1>(0,0)]=1;p[OpKeyType<1>(1,1)]=-1;p[OpKeyType<1>(2,2)]=-1;p[OpKeyType<1>(3,3)]=1;
    p[OpKeyType<1>(4,4)]=-1;p[OpKeyType<1>(5,5)]=1;p[OpKeyType<1>(6,6)]=1;p[OpKeyType<1>(7,7)]=-1;
    p[OpKeyType<1>(8,8)]=-1;p[OpKeyType<1>(9,9)]=1;p[OpKeyType<1>(10,10)]=1;p[OpKeyType<1>(11,11)]=-1;
    p[OpKeyType<1>(12,12)]=1;p[OpKeyType<1>(13,13)]=-1;p[OpKeyType<1>(14,14)]=-1;p[OpKeyType<1>(15,15)]=1;

    pL[OpKeyType<1>(0,0)]=1;pL[OpKeyType<1>(1,1)]=-1;pL[OpKeyType<1>(2,2)]=-1;pL[OpKeyType<1>(3,3)]=1;
    pL[OpKeyType<1>(4,4)]=1;pL[OpKeyType<1>(5,5)]=-1;pL[OpKeyType<1>(6,6)]=-1;pL[OpKeyType<1>(7,7)]=1;
    pL[OpKeyType<1>(8,8)]=1;pL[OpKeyType<1>(9,9)]=-1;pL[OpKeyType<1>(10,10)]=-1;pL[OpKeyType<1>(11,11)]=1;
    pL[OpKeyType<1>(12,12)]=1;pL[OpKeyType<1>(13,13)]=-1;pL[OpKeyType<1>(14,14)]=-1;pL[OpKeyType<1>(15,15)]=1;
 
    pU[OpKeyType<1>(0,0)]=1;pU[OpKeyType<1>(1,1)]=1;pU[OpKeyType<1>(2,2)]=1;pU[OpKeyType<1>(3,3)]=1;
    pU[OpKeyType<1>(4,4)]=-1;pU[OpKeyType<1>(5,5)]=-1;pU[OpKeyType<1>(6,6)]=-1;pU[OpKeyType<1>(7,7)]=-1;
    pU[OpKeyType<1>(8,8)]=1;pU[OpKeyType<1>(9,9)]=1;pU[OpKeyType<1>(10,10)]=1;pU[OpKeyType<1>(11,11)]=1;
    pU[OpKeyType<1>(12,12)]=-1;pU[OpKeyType<1>(13,13)]=-1;pU[OpKeyType<1>(14,14)]=-1;pU[OpKeyType<1>(15,15)]=-1;
  
    bla=pL*cuU;cuU=bla;
    bla=pL*cdU;cdU=bla;
    bla=pL*cuDU;cuDU=bla;
    bla=pL*cdDU;cdDU=bla;
  
    if(spin==2)
      cr=cuDU;
    else if(spin==-2)
      cr=cdDU;
    else if(spin==1)
      cr=cuDL;
    else if(spin==-1)
      cr=cdDL;

    p.setSite(0);cr.setSite(0);O_11.setSite(0);
    //you may wonder about this, but if the coefficients get too small (much smaller than machine precision), then, after application of this operator to
    //the GS at halffilling (which btw works fine), doing an orthonormalization gives pretty strange results. I THINK that this is an degeneracy issue.
    //this works now: if the non diagonal entries are very very small, the mpo matrix is 1x1. if the coeff is larger than a threshold, the usual
    //representation is used ... i think that before, i had multiple superpositions of the state with itself, which led to problems with the svd in orthonormalize
    bool first=true,last=true;
    for(uint s=0;s<(N);++s){
      p.setSite(s);cr.setSite(s);O_11.setSite(s);
      if(abs(coeff(s))>1e-14){
	if(first){
	  mat[ukey<2>(0,0)] = coeff(s)*cr;
	  mat[ukey<2>(0,1)] = p;
	  mat.Dl=1;mat.Dr=2;
	  (*this)[s]=mat;
	  mat.clear();
	  first=false;
	}else{
	  mat[ukey<2>(0,0)] = O_11;
	  mat[ukey<2>(1,0)] = coeff(s)*cr;
	  mat[ukey<2>(1,1)] = p;
	  mat.Dl=2;mat.Dr=2;
	  (*this)[s]=mat;
	  mat.clear();
	}
      }else{
	if(first){
	  mat[ukey<2>(0,0)] = p;
	  mat.Dl=1;mat.Dr=1;
	  (*this)[s]=mat;
	  mat.clear();
	}
	if(!first&&last){
	  mat[ukey<2>(0,0)] = O_11;
	  mat[ukey<2>(1,0)] = coeff(s-1)*cr;
	  mat.Dl=2;mat.Dr=1;
	  (*this)[s-1]=mat;
	  mat.clear();
	  last=false;
	}
	if(!first&&!last){
	  mat[ukey<2>(0,0)] = O_11;
	  mat.Dl=1;mat.Dr=1;
	  (*this)[s]=mat;
	  mat.clear();
	}
      }
    }
  }
  virtual ~gaussMPOfermihubladder(){};
};


//that's for spinless fermions ..
class gaussMPOfermi: public MPO<Complex> {
public:
  gaussMPOfermi():MPO<Complex>(){};
  gaussMPOfermi(uint N,Real x0,Real k0,Real sigma){
    CVectorType coeff(N);
    for(uint n=0;n<N;++n){
      coeff(n)=abs( exp(Complex(-2.*(pow(sigma,2.0)*pow(M_PI,2.0)*pow(n- x0,2.0)),M_PI*k0/N*(n - x0))) )>1e-12?//should be i*2*pi, and k0=N/4
	exp(Complex(-2*(pow(sigma,2)*pow(M_PI,2)*pow(n- x0,2)),M_PI*k0/N*(n - x0)) ):Complex(0.0,0.0);
    }
    MPOMat<1,Complex> mat;
    SparseLocalOperator<Complex> p(2,2),cd(2,2),O_11(2,2);
    p.setSite(0);cd.setSite(0);O_11.setSite(0);
    cd[OpKeyType<1>(1,0)]=Complex(1.);
    p[OpKeyType<1>(0,0)]=Complex(1.);
    p[OpKeyType<1>(1,1)]=Complex(-1.);
    O_11[OpKeyType<1>(0,0)]=Complex(1.);
    O_11[OpKeyType<1>(1,1)]=Complex(1.);

    //you may wonder about this, but if the coefficients get too small (much smaller than machine precision), then, after application of this operator to
    //the GS at halffilling (which btw works fine), doing an orthonormalization gives pretty strange results. I THINK that this is an degeneracy issue.
    //this works now: if the non diagonal entries are very very small, the mpo matrix is 1x1. if the coeff is larger than a threshold, the usual
    //representation is used ... i think that before, i had multiple superpositions of the state with itself, which led to problems with the svd in orthonormalize
    bool first=true,last=true;
    for(uint s=0;s<(N);++s){
      p.setSite(s);cd.setSite(s);O_11.setSite(s);
      if(abs(coeff(s))>1e-12){
	if(first){
	  mat[ukey<2>(0,0)] = coeff(s)*cd;
	  mat[ukey<2>(0,1)] = p;
	  mat.Dl=1;mat.Dr=2;
	  (*this)[s]=mat;
	  mat.clear();
	  first=false;
	}else{
	  mat[ukey<2>(0,0)] = O_11;
	  mat[ukey<2>(1,0)] = coeff(s)*cd;
	  mat[ukey<2>(1,1)] = p;
	  mat.Dl=2;mat.Dr=2;
	  (*this)[s]=mat;
	  mat.clear();
	}
      }else{
	if(first){
	  mat[ukey<2>(0,0)] = p;
	  mat.Dl=1;mat.Dr=1;
	  (*this)[s]=mat;
	  mat.clear();
	}
	if(!first&&last){
	  mat[ukey<2>(0,0)] = O_11;
	  mat[ukey<2>(1,0)] = coeff(s-1)*cd;
	  mat.Dl=2;mat.Dr=1;
	  (*this)[s-1]=mat;
	  mat.clear();
	  last=false;
	}
	if(!first&&!last){
	  mat[ukey<2>(0,0)] = O_11;
	  mat.Dl=1;mat.Dr=1;
	  (*this)[s]=mat;
	  mat.clear();
	}
      }
    }
  }
  virtual ~gaussMPOfermi(){};
};

class gaussMPOXXZ: public MPO<Complex> {
public:
  gaussMPOXXZ():MPO<Complex>(){};
  gaussMPOXXZ(uint N,Real x0,Real k0,Real sigma){
    CVectorType coeff(N);
    for(uint n=0;n<N;++n){
      coeff(n)=abs( exp(Complex(-2*(pow(sigma,2)*pow(M_PI,2)*pow(n- x0,2)),M_PI*k0/(N+1)*(n - x0))) )>1e-15?
	exp(Complex(-2*(pow(sigma,2)*pow(M_PI,2)*pow(n- x0,2)),M_PI*k0/(N+1)*(n - x0)) ):1e-15;
    }

    MPOMat<1,Complex> mat;
    SparseLocalOperator<Complex> Sp(2,2),O_11(2,2);
    Sp.setSite(0);O_11.setSite(0);
    Sp[OpKeyType<1>(1,0)]=Complex(1.);
    O_11[OpKeyType<1>(0,0)]=Complex(1.);
    O_11[OpKeyType<1>(1,1)]=Complex(1.);

    //you may wonder about this, but if the coefficients get too small (much smaller than machine precision), then, after application of this operator to
    //the GS at halffilling (which btw works fine), doing an orthonormalization gives pretty strange results. I THINK that this is an degeneracy issue.
    //this works now: if the non diagonal entries are very very small, the mpo matrix is 1x1. if the coeff is larger than a threshold, the usual
    //representation is used ... i think that before, i had multiple superpositions of the state with itself, which led to problems with the svd in orthonormalize
    bool first=true,last=true;
    for(uint s=0;s<(N);++s){
      Sp.setSite(s);O_11.setSite(s);
      if(abs(coeff(s))>1e-8){
	if(first){
	  mat[ukey<2>(0,0)] = coeff(s)*Sp;
	  mat[ukey<2>(0,1)] = O_11;
	  mat.Dl=1;mat.Dr=2;
	  (*this)[s]=mat;
	  mat.clear();
	  first=false;
	}else{
	  mat[ukey<2>(0,0)] = O_11;
	  mat[ukey<2>(1,0)] = coeff(s)*Sp;
	  mat[ukey<2>(1,1)] = O_11;
	  mat.Dl=2;mat.Dr=2;
	  (*this)[s]=mat;
	  mat.clear();
	}
      }else{
	if(first){
	  mat[ukey<2>(0,0)]=O_11;
	  mat.Dl=1;mat.Dr=1;
	  (*this)[s]=mat;
	  mat.clear();
	}
	if(!first&&last){
	  mat[ukey<2>(0,0)] = O_11;
	  mat[ukey<2>(1,0)] = coeff(s-1)*Sp;
	  mat.Dl=2;mat.Dr=1;
	  (*this)[s-1]=mat;
	  mat.clear();
	  last=false;
	}
	if(!first&&!last){
	  mat[ukey<2>(0,0)] = O_11;
	  mat.Dl=1;mat.Dr=1;
	  (*this)[s]=mat;
	  mat.clear();
	}
      }
    }
  }
  virtual ~gaussMPOXXZ(){};
};

//that's giving not the desired results ... . It doesn't produce a travelling two-string


class gaussTwoStringMPOXXZ: public MPO<Complex> {
public:
  gaussTwoStringMPOXXZ():MPO<Complex>(){};
  gaussTwoStringMPOXXZ(uint N,Real x0,Real k0,Real sigma){
    CVectorType coeff(N);
    for(uint n=0;n<N;++n){
      coeff(n)=abs( exp(Complex(-2*(pow(sigma,2)*pow(M_PI,2)*pow(n- x0,2)),M_PI*k0/(N+1)*(n - x0))) )>1e-16?
	exp(Complex(-2*(pow(sigma,2)*pow(M_PI,2)*pow(n- x0,2)),M_PI*k0/(N+1)*(n - x0)) ):1e-16;
    }

    MPOMat<1,Complex> mat;
    SparseLocalOperator<Complex> Sp(2,2),O_11(2,2);
    Sp.setSite(0);O_11.setSite(0);
    Sp[OpKeyType<1>(1,0)]=Complex(1.);
    O_11[OpKeyType<1>(0,0)]=Complex(1.);
    O_11[OpKeyType<1>(1,1)]=Complex(1.);

    //you may wonder about this, but if the coefficients get too small (much smaller than machine precision), then, after application of this operator to
    //the GS at halffilling (which btw works fine), doing an orthonormalization gives pretty strange results. I THINK that this is an degeneracy issue.
    //this works now: if the non diagonal entries are very very small, the mpo matrix is 1x1. if the coeff is larger than a threshold, the usual
    //representation is used ... i think that before, i had multiple superpositions of the state with itself, which led to problems with the svd in orthonormalize
    bool first=true,last=true;
    for(uint s=0;s<(N);++s){
      Sp.setSite(s);O_11.setSite(s);
      if(abs(coeff(s))>1e-15){
	if(first){
	  mat[ukey<2>(0,1)] = coeff(s)*Sp;
	  mat[ukey<2>(0,2)] = O_11;
	  mat.Dl=1;mat.Dr=3;
	  (*this)[s]=mat;
	  mat.clear();
	  first=false;
	}else{
	  mat[ukey<2>(0,0)] = O_11;
	  mat[ukey<2>(1,0)] = Sp;
	  mat[ukey<2>(2,1)] = coeff(s)*Sp;
	  mat[ukey<2>(2,2)] = O_11;
	  mat.Dl=3;mat.Dr=3;
	  (*this)[s]=mat;
	  mat.clear();
	}
      }else{
	if(first){
	  mat[ukey<2>(0,0)] = O_11;
	  mat.Dl=1;mat.Dr=1;
	  (*this)[s]=mat;
	  mat.clear();
	}
	if(!first&&last){
	  mat[ukey<2>(0,0)] = O_11;
	  mat[ukey<2>(1,0)] = coeff(s-1)*Sp;
	  mat.Dl=3;mat.Dr=1;
	  (*this)[s-1]=mat;
	  mat.clear();
	  last=false;
	}
	if(!first&&!last){
	  mat[ukey<2>(0,0)] = O_11;
	  mat.Dl=1;mat.Dr=1;
	  (*this)[s]=mat;
	  mat.clear();
	}
      }
    }
  }
  virtual ~gaussTwoStringMPOXXZ(){};
};

class gaussMPOBoseHub: public MPO<Complex> {
public:
  gaussMPOBoseHub():MPO<Complex>(){};
  gaussMPOBoseHub(uint N,Real x0,Real k0,Real sigma,uint maxdim){
    CVectorType coeff(N);
    for(uint n=0;n<N;++n){
      coeff(n)=abs( exp(Complex(-2*(pow(sigma,2)*pow(M_PI,2)*pow(n- x0,2)),M_PI*k0/N*(n - x0))) )>1e-15?
	exp(Complex(-2*(pow(sigma,2)*pow(M_PI,2)*pow(n- x0,2)),M_PI*k0/N*(n - x0)) ):Complex(0.0,0.0);
    }
    MPOMat<1,Complex> mat;
    //  this->resize(N);
    //all statistics are set to bosonic (_statistics=1)
    SparseLocalOperator<Complex> bd(maxdim,maxdim),O_11(maxdim,maxdim);
    for(int ind = 0;ind<maxdim;++ind){
      O_11[OpKeyType<1>(ind,ind)]=Complex(1.0);
      if(ind<(maxdim-1))
	bd[OpKeyType<1>(ind+1,ind)] = Complex(sqrt(ind+1));
    }
    //you may wonder about this, but if the coefficients get too small (much smaller than machine precision), then, after application of this operator to
    //the GS at halffilling (which btw works fine), doing an orthonormalization gives pretty strange results. I THINK that this is an degeneracy issue.
    //this works now: if the non diagonal entries are very very small, the mpo matrix is 1x1. if the coeff is larger than a threshold, the usual
    //representation is used ... i think that before, i had multiple superpositions of the state with itself, which led to problems with the svd in 
    //orthonormalize
    bool first=true,last=true;
    for(uint s=0;s<(N);++s){
      bd.setSite(s);O_11.setSite(s);
      if(abs(coeff(s))>1e-15){
	if(first){
	  mat[ukey<2>(0,0)] = coeff(s)*bd;
	  mat[ukey<2>(0,1)] = O_11;
	  mat.Dl=1;mat.Dr=2;
	  (*this)[s]=mat;
	  mat.clear();
	  first=false;
	}else{
	  mat[ukey<2>(0,0)] = O_11;
	  mat[ukey<2>(1,0)] = coeff(s)*bd;
	  mat[ukey<2>(1,1)] = O_11;
	  mat.Dl=2;mat.Dr=2;
	  (*this)[s]=mat;
	  mat.clear();
	}
      }else{
	if(first){
	  mat[ukey<2>(0,0)] = O_11;
	  mat.Dl=1;mat.Dr=1;
	  (*this)[s]=mat;
	  mat.clear();
	}
	if(!first&&last){
	  mat[ukey<2>(0,0)] = O_11;
	  mat[ukey<2>(1,0)] = coeff(s-1)*bd;
	  mat.Dl=2;mat.Dr=1;
	  (*this)[s-1]=mat;
	  mat.clear();
	  last=false;
	}
	if(!first&&!last){
	  mat[ukey<2>(0,0)] = O_11;
	  mat.Dl=1;mat.Dr=1;
	  (*this)[s]=mat;
	  mat.clear();
	}
      }
    }
  }
  virtual ~gaussMPOBoseHub(){};
};




//ordering of local states at the measuring sites (particle states and spin states): |0,down>,|1,down>,|0,up>,|1,up>
class gaussianparticleDetectorMPO: public MPO<Complex> {
public:
  gaussianparticleDetectorMPO():MPO<Complex>(){};
  gaussianparticleDetectorMPO(uint N,Real x0,Real k0,Real sigma,VectorType positions,uint maxdim){
    CVectorType coeff(N);
    for(uint n=0;n<N;++n){
      Complex alpha=exp(Complex(-2*(pow(sigma,2)*pow(M_PI,2)*pow(n- x0,2)),M_PI*k0/N*(n - x0)));
      coeff(n)=abs(alpha)>1e-15?alpha:Complex(0.0,0.0);
      
    }

    MPOMat<1,Complex> mat;

    SparseLocalOperator<Complex> bd(maxdim,maxdim),id(maxdim,maxdim),idS(2,2);
    idS[OpKeyType<1>(0,0)]=1.0;idS[OpKeyType<1>(1,1)]=1.0;
    for(int ind = 0;ind<maxdim;++ind){
      id[OpKeyType<1>(ind,ind)]=Complex(1.0);
      if(ind<(maxdim-1))
	bd[OpKeyType<1>(ind+1,ind)] = Complex(sqrt(ind+1));
    }


    SparseLocalOperator<Complex> id_id=kronecker(id,idS);
    SparseLocalOperator<Complex> bd_id=kronecker(bd,idS);

    //you may wonder about this, but if the coefficients get too small (much smaller than machine precision), then, after application of this operator to
    //the GS at halffilling (which btw works fine), doing an orthonormalization gives pretty strange results. I THINK that this is an degeneracy issue.
    //this works now: if the non diagonal entries are very very small, the mpo matrix is 1x1. if the coeff is larger than a threshold, the usual
    //representation is used ... i think that before, i had multiple superpositions of the state with itself, which led to problems with the svd in 
    //orthonormalize
    bool first=true,last=true;
    uint pos=0;
    uint currentSite=positions(pos);
    for(uint s=0;s<(N);++s){
      if (s!=currentSite){
	bd.setSite(s);id.setSite(s);
	if(abs(coeff(s))>1e-15){
	  if(first){
	    mat[ukey<2>(0,0)] = coeff(s)*bd;
	    mat[ukey<2>(0,1)] = id;
	    mat.Dl=1;mat.Dr=2;
	    (*this)[s]=mat;
	    mat.clear();
	    first=false;
	  }else{
	    mat[ukey<2>(0,0)] = id;
	    mat[ukey<2>(1,0)] = coeff(s)*bd;
	    mat[ukey<2>(1,1)] = id;
	    mat.Dl=2;mat.Dr=2;
	    (*this)[s]=mat;
	    mat.clear();
	  }
	}else{
	  if(first){
	    mat[ukey<2>(0,0)] = id;
	    mat.Dl=1;mat.Dr=1;
	    (*this)[s]=mat;
	    mat.clear();
	  }
	  if(!first&&last){
	    mat[ukey<2>(0,0)] = id;
	    mat[ukey<2>(1,0)] = coeff(s-1)*bd;
	    mat.Dl=2;mat.Dr=1;
	    (*this)[s-1]=mat;
	    mat.clear();
	    last=false;
	  }
	  if(!first&&!last){
	    mat[ukey<2>(0,0)] = id;
	    mat.Dl=1;mat.Dr=1;
	    (*this)[s]=mat;
	    mat.clear();
	  }
	}
      }
      if (s==currentSite){
       	bd_id.setSite(s);id_id.setSite(s);
       	if(abs(coeff(s))>1e-15){
       	  if(first==true){
       	    mat[ukey<2>(0,0)] = coeff(s)*bd_id;
       	    mat[ukey<2>(0,1)] = id_id;
       	    mat.Dl=1;mat.Dr=2;
       	    (*this)[s]=mat;
       	    mat.clear();
       	    first=false;
       	  }else{
       	    mat[ukey<2>(0,0)] = id_id;
       	    mat[ukey<2>(1,0)] = coeff(s)*bd_id;
       	    mat[ukey<2>(1,1)] = id_id;
       	    mat.Dl=2;mat.Dr=2;
       	    (*this)[s]=mat;
       	    mat.clear();
       	  }
       	}else{
       	  if(first){
       	    mat[ukey<2>(0,0)] = id_id;
       	    mat.Dl=1;mat.Dr=1;
       	    (*this)[s]=mat;
       	    mat.clear();
       	  }
       	  if(!first&&last){
       	    mat[ukey<2>(0,0)] = id_id;
       	    mat[ukey<2>(1,0)] = coeff(s-1)*bd_id;
       	    mat.Dl=2;mat.Dr=1;
       	    (*this)[s-1]=mat;
       	    mat.clear();
       	    last=false;
       	  }
       	  if(!first&&!last){
       	    mat[ukey<2>(0,0)] = id_id;
       	    mat.Dl=1;mat.Dr=1;
       	    (*this)[s]=mat;
       	    mat.clear();
       	  }
       	}
       	if(pos<positions.size()-1){
       	  pos++;
       	  currentSite=positions(pos);
       	}
      }
    }
  }
  virtual ~gaussianparticleDetectorMPO(){};
};


//this is a fermi hubbard ladder system, string order operators could be wrong ...
//for completeness,here is the state numbering I use; first one is the lower chain, the secon one the upper
//0  0,0
//1  d,0
//2  u,0
//3  ud,0
//4  0,d
//5  d,d
//6  u,d
//7  ud,d
//8  0,u
//9  d,u
//10 u,u
//11 ud,u
//12 0,ud
//13 d,ud
//14 u,ud
//15 ud,ud
//sign convention: H=-t\sum_i (c_i^{\dagger}c_{i+1} +h.c.)+U\sum_i n_{i,up}n_{i,down})

template<typename T>
class BenzeneMPO:public MPO<T>
{
public:
  BenzeneMPO(VectorType hopbath,VectorType mubathUp,VectorType mubathDown,Real Ubenzene,Real hopbenzene,Real mubenzene);

  ~BenzeneMPO(){};
  virtual T measure(const MPS<AbKeyType,T>&mps)const {return MPO<T>::measure(mps);}
  virtual T measure(const CanonizedMPS<AbKeyType,T>&mps)const {return MPO<T>::measure(mps);}

};



template<typename T>
BenzeneMPO<T>::BenzeneMPO(VectorType hopbath,VectorType mubathUp,VectorType mubathDown,Real Ubenzene,Real hopbenzene,Real mubenzene){

 // cout<<"hopbath: "<<hopbath<<endl;
 // cout<<"mubathUp: "<<mubathUp<<endl;
 // cout<<"mubathDown: "<<mubathDown<<endl;
 // cout<<"Ubenzene: "<<Ubenzene<<endl;
 // cout<<"hopbenzene: "<<hopbenzene<<endl;
 // cout<<"mubenzene: "<<mubenzene<<endl;
  const uint N =hopbath.size()+1;//total length
  const uint Nl=(N-3)/2;//lenght of the bath chains
  // cout<<"Nl="<<Nl<<endl;
  
  if((N-3)%2!=0){
    cout<<N<<endl;
    cout<<"#N-bath="<<Nl<<": bath chains with different lenght found use N with (N-3) mod 2 = 0; aborting ... "<<endl;
    abort();
  }
  //START WITH THE BATH: LOCAL OPERATOR DEFINITION FOR THE BATH CHAINS
  SparseLocalOperator<T> nu(4,4),nd(4,4),cuD(4,4),cdD(4,4),cu(4,4),cd(4,4),idbath(4,4),pbath(4,4);
  nu[OpKeyType<1>(2,2)]=T(1.0);nu[OpKeyType<1>(3,3)]=T(1.0);
  nd[OpKeyType<1>(1,1)]=T(1.0);nd[OpKeyType<1>(3,3)]=T(1.0);

  cuD[OpKeyType<1>(2,0)] = T(1.0);cuD[OpKeyType<1>(3,1)] = T(1.0); 
  cu[OpKeyType<1>(0,2)] = T(1.0);cu[OpKeyType<1>(1,3)] = T(1.0);   
								   
  cdD[OpKeyType<1>(1,0)] = T(1.0);cdD[OpKeyType<1>(3,2)] = T(-1.0);
  cd[OpKeyType<1>(0,1)] = T(1.0);cd[OpKeyType<1>(2,3)] = T(-1.0);  

  idbath[OpKeyType<1>(0,0)]=T(1);idbath[OpKeyType<1>(1,1)]=T(1);idbath[OpKeyType<1>(2,2)]=T(1);idbath[OpKeyType<1>(3,3)]=T(1);
  pbath[OpKeyType<1>(0,0)]=T(1);pbath[OpKeyType<1>(1,1)]=T(-1);pbath[OpKeyType<1>(2,2)]=T(-1);pbath[OpKeyType<1>(3,3)]=T(1);

  nu.setSite(0);nd.setSite(0);cu.setSite(0);cd.setSite(0);cuD.setSite(0);cdD.setSite(0);idbath.setSite(0);pbath.setSite(0);


  MPOMat<1,T> M;
  //START FROM LEFT, MOVE TO THE RIGHT UNTIL YOU REACH THE BENZENE MOLEKULE
  M[ukey<2>(0,0)] = mubathUp(0)*nu+mubathDown(0)*nd;
  M[ukey<2>(0,1)] = hopbath(0)*(cu*pbath);
  M[ukey<2>(0,2)] = -hopbath(0)*(cuD*pbath);
  M[ukey<2>(0,3)] = hopbath(0)*(cd*pbath);
  M[ukey<2>(0,4)] = -hopbath(0)*(cdD*pbath);
  M[ukey<2>(0,5)] = idbath;

  M.Dl=1;M.Dr=6;
  M.squeeze(1e-15);
  (*this)[0]=M;
  M.clear();

  for(uint site = 1;site<Nl;++site){
    nu.setSite(site);nd.setSite(site);cu.setSite(site);cd.setSite(site);cuD.setSite(site);cdD.setSite(site);idbath.setSite(site);pbath.setSite(site);
    M[ukey<2>(0,0)] = idbath;
    M[ukey<2>(1,0)] = cuD;
    M[ukey<2>(2,0)] = cu;
    M[ukey<2>(3,0)] = cdD;
    M[ukey<2>(4,0)] = cd;
    M[ukey<2>(5,0)] = mubathUp(site)*nu+mubathDown(site)*nd;
     
    M[ukey<2>(5,1)] = hopbath(site)*(cu*pbath);
    M[ukey<2>(5,2)] = -hopbath(site)*(cuD*pbath);
    M[ukey<2>(5,3)] = hopbath(site)*(cd*pbath);
    M[ukey<2>(5,4)] = -hopbath(site)*(cdD*pbath);
    M[ukey<2>(5,5)] = idbath;
    M.Dl=6;M.Dr=6;
    M.squeeze(1e-15);
    (*this)[site]=M;
    M.clear();
  }

  //DEFINITION OF LOCAL OPERATORS ON THE DOUBLE-SITE SYSTEM BENZENE
  //LOWER SITES

  SparseLocalOperator<T> nuU(16,16),ndU(16,16),cuDU(16,16),cdDU(16,16),cuU(16,16),cdU(16,16),id(16,16),p(16,16),pU(16,16),pL(16,16);
  SparseLocalOperator<T> nuL(16,16),ndL(16,16),cuDL(16,16),cdDL(16,16),cuL(16,16),cdL(16,16);

  cuDL[OpKeyType<1>(2,0)] = T(1.0);cuDL[OpKeyType<1>(3,1)] = T(1.0);cuDL[OpKeyType<1>(6,4)] = T(1.0);cuDL[OpKeyType<1>(7,5)] = T(1.0);
  cuDL[OpKeyType<1>(10,8)] = T(1.0);cuDL[OpKeyType<1>(11,9)] = T(1.0);cuDL[OpKeyType<1>(14,12)] = T(1.0);cuDL[OpKeyType<1>(15,13)] = T(1.0);
 
  cuL[OpKeyType<1>(0,2)] = T(1.0);cuL[OpKeyType<1>(1,3)] = T(1.0);cuL[OpKeyType<1>(4,6)] = T(1.0);cuL[OpKeyType<1>(5,7)] = T(1.0);
  cuL[OpKeyType<1>(8,10)] = T(1.0);cuL[OpKeyType<1>(9,11)] = T(1.0);cuL[OpKeyType<1>(12,14)] = T(1.0);cuL[OpKeyType<1>(13,15)] = T(1.0);
								    
  cdL[OpKeyType<1>(0,1)] = T(1.0);cdL[OpKeyType<1>(2,3)] = T(-1.0);cdL[OpKeyType<1>(4,5)] = T(1.0);cdL[OpKeyType<1>(6,7)] = T(-1.0);
  cdL[OpKeyType<1>(8,9)] = T(1.0);cdL[OpKeyType<1>(10,11)] = T(-1.0);cdL[OpKeyType<1>(12,13)] = T(1.0);cdL[OpKeyType<1>(14,15)] = T(-1.0);

  cdDL[OpKeyType<1>(1,0)] = T(1.0);cdDL[OpKeyType<1>(3,2)] = T(-1.0);cdDL[OpKeyType<1>(5,4)] = T(1.0);cdDL[OpKeyType<1>(7,6)] = T(-1.0);
  cdDL[OpKeyType<1>(9,8)] = T(1.0);cdDL[OpKeyType<1>(11,10)] = T(-1.0);cdDL[OpKeyType<1>(13,12)] = T(1.0);cdDL[OpKeyType<1>(15,14)] = T(-1.0);

  //UPPER SITES
  cuDU[OpKeyType<1>(8,0)] = T(1.0);cuDU[OpKeyType<1>(9,1)] = T(1.0);cuDU[OpKeyType<1>(10,2)] = T(1.0);cuDU[OpKeyType<1>(11,3)] = T(1.0);
  cuDU[OpKeyType<1>(12,4)] = T(1.0);cuDU[OpKeyType<1>(13,5)] = T(1.0);cuDU[OpKeyType<1>(14,6)] = T(1.0);cuDU[OpKeyType<1>(15,7)] = T(1.0);
     
  cuU[OpKeyType<1>(0,8)] = T(1.0);cuU[OpKeyType<1>(1,9)] = T(1.0);cuU[OpKeyType<1>(2,10)] = T(1.0);cuU[OpKeyType<1>(3,11)] = T(1.0);
  cuU[OpKeyType<1>(4,12)] = T(1.0);cuU[OpKeyType<1>(5,13)] = T(1.0);cuU[OpKeyType<1>(6,14)] = T(1.0);cuU[OpKeyType<1>(7,15)] = T(1.0);
     
  cdDU[OpKeyType<1>(4,0)] = T(1.0);cdDU[OpKeyType<1>(5,1)] = T(1.0);cdDU[OpKeyType<1>(6,2)] = T(1.0);cdDU[OpKeyType<1>(7,3)] = T(1.0);
  cdDU[OpKeyType<1>(12,8)] = T(-1.0);cdDU[OpKeyType<1>(13,9)] = T(-1.0);cdDU[OpKeyType<1>(14,10)] = T(-1.0);cdDU[OpKeyType<1>(15,11)] = T(-1.0);

  cdU[OpKeyType<1>(0,4)] = T(1.0);cdU[OpKeyType<1>(1,5)] = T(1.0);cdU[OpKeyType<1>(2,6)] = T(1.0);cdU[OpKeyType<1>(3,7)] = T(1.0);
  cdU[OpKeyType<1>(8,12)] = T(-1.0);cdU[OpKeyType<1>(9,13)] = T(-1.0);cdU[OpKeyType<1>(10,14)] = T(-1.0);cdU[OpKeyType<1>(11,15)] = T(-1.0);

  SparseLocalOperator<T> bla(16,16);
  
  for(uint ind=0;ind<16;++ind)
    id[OpKeyType<1>(ind,ind)]=T(1);

  //FOLLOWING TAKES CARE OF FERMIONIC SIGNS
  p[OpKeyType<1>(0,0)]=T(1);p[OpKeyType<1>(1,1)]=T(-1);p[OpKeyType<1>(2,2)]=T(-1);p[OpKeyType<1>(3,3)]=T(1);
  p[OpKeyType<1>(4,4)]=T(-1);p[OpKeyType<1>(5,5)]=T(1);p[OpKeyType<1>(6,6)]=T(1);p[OpKeyType<1>(7,7)]=T(-1);
  p[OpKeyType<1>(8,8)]=T(-1);p[OpKeyType<1>(9,9)]=T(1);p[OpKeyType<1>(10,10)]=T(1);p[OpKeyType<1>(11,11)]=T(-1);
  p[OpKeyType<1>(12,12)]=T(1);p[OpKeyType<1>(13,13)]=T(-1);p[OpKeyType<1>(14,14)]=T(-1);p[OpKeyType<1>(15,15)]=T(1);

  pL[OpKeyType<1>(0,0)]=T(1);pL[OpKeyType<1>(1,1)]=T(-1);pL[OpKeyType<1>(2,2)]=T(-1);pL[OpKeyType<1>(3,3)]=T(1);
  pL[OpKeyType<1>(4,4)]=T(1);pL[OpKeyType<1>(5,5)]=T(-1);pL[OpKeyType<1>(6,6)]=T(-1);pL[OpKeyType<1>(7,7)]=T(1);
  pL[OpKeyType<1>(8,8)]=T(1);pL[OpKeyType<1>(9,9)]=T(-1);pL[OpKeyType<1>(10,10)]=T(-1);pL[OpKeyType<1>(11,11)]=T(1);
  pL[OpKeyType<1>(12,12)]=T(1);pL[OpKeyType<1>(13,13)]=T(-1);pL[OpKeyType<1>(14,14)]=T(-1);pL[OpKeyType<1>(15,15)]=T(1);
 
  pU[OpKeyType<1>(0,0)]=T(1);pU[OpKeyType<1>(1,1)]=T(1);pU[OpKeyType<1>(2,2)]=T(1);pU[OpKeyType<1>(3,3)]=T(1);
  pU[OpKeyType<1>(4,4)]=T(-1);pU[OpKeyType<1>(5,5)]=T(-1);pU[OpKeyType<1>(6,6)]=T(-1);pU[OpKeyType<1>(7,7)]=T(-1);
  pU[OpKeyType<1>(8,8)]=T(1);pU[OpKeyType<1>(9,9)]=T(1);pU[OpKeyType<1>(10,10)]=T(1);pU[OpKeyType<1>(11,11)]=T(1);
  pU[OpKeyType<1>(12,12)]=T(-1);pU[OpKeyType<1>(13,13)]=T(-1);pU[OpKeyType<1>(14,14)]=T(-1);pU[OpKeyType<1>(15,15)]=T(-1);
  
  //THAT'S CRUCIAL!!!
  bla=pL*cuU;cuU=bla;
  bla=pL*cdU;cdU=bla;
  bla=pL*cuDU;cuDU=bla;
  bla=pL*cdDU;cdDU=bla;

  nuL[OpKeyType<1>(2,2)]=T(1);nuL[OpKeyType<1>(3,3)]=T(1);nuL[OpKeyType<1>(6,6)]=T(1);nuL[OpKeyType<1>(7,7)]=T(1);
  nuL[OpKeyType<1>(10,10)]=T(1);nuL[OpKeyType<1>(11,11)]=T(1);nuL[OpKeyType<1>(14,14)]=T(1);nuL[OpKeyType<1>(15,15)]=T(1);

  ndL[OpKeyType<1>(1,1)]=T(1);ndL[OpKeyType<1>(3,3)]=T(1);ndL[OpKeyType<1>(5,5)]=T(1);ndL[OpKeyType<1>(7,7)]=T(1);
  ndL[OpKeyType<1>(9,9)]=T(1);ndL[OpKeyType<1>(11,11)]=T(1);ndL[OpKeyType<1>(13,13)]=T(1);ndL[OpKeyType<1>(15,15)]=T(1);
  
  for(uint ind=8;ind<16;++ind)
    nuU[OpKeyType<1>(ind,ind)]=T(1);
  for(uint ind=4;ind<8;++ind)
    ndU[OpKeyType<1>(ind,ind)]=T(1);
  for(uint ind=12;ind<16;++ind)
    ndU[OpKeyType<1>(ind,ind)]=T(1);

  nuU.setSite(Nl);ndU.setSite(Nl);cuU.setSite(Nl);cdU.setSite(Nl);cuDU.setSite(Nl);cdDU.setSite(Nl);id.setSite(Nl);p.setSite(Nl);
  nuL.setSite(Nl);ndL.setSite(Nl);cuL.setSite(Nl);cdL.setSite(Nl);cuDL.setSite(Nl);cdDL.setSite(Nl);pU.setSite(Nl);pL.setSite(Nl);
  
  M[ukey<2>(0,0)] = id;
  M[ukey<2>(1,0)] = cuDL;
  M[ukey<2>(2,0)] = cuL;
  M[ukey<2>(3,0)] = cdDL;
  M[ukey<2>(4,0)] = cdL;
  M[ukey<2>(5,0)] = Ubenzene*((nuL*ndL)+(nuU*ndU))+mubenzene*((nuL+ndL)+(nuU+ndU))-hopbenzene*((cuDL*cuU)+(cuDU*cuL)+(cdDL*cdU)+(cdDU*cdL));
  //M[ukey<2>(5,0)] = Ubenzene*(nuL*ndL+nuU*ndU)+mubenzene*(nuL+ndL+nuU+ndU);

  M[ukey<2>(5,1)] = -hopbenzene*(cuL*p);
  M[ukey<2>(5,2)] = hopbenzene*(cuDL*p);//
  M[ukey<2>(5,3)] = -hopbenzene*(cuU*p);
  M[ukey<2>(5,4)] = hopbenzene*(cuDU*p);//
  M[ukey<2>(5,5)] = -hopbenzene*(cdL*p);
  M[ukey<2>(5,6)] = hopbenzene*(cdDL*p);//
  M[ukey<2>(5,7)] = -hopbenzene*(cdU*p);
  M[ukey<2>(5,8)] = hopbenzene*(cdDU*p);//
  M[ukey<2>(5,9)] = id;

  M.Dl=6,M.Dr=10;
  M.squeeze(1e-15);
  (*this)[Nl]=M;
  M.clear();
  
  nuU.setSite(Nl+1);ndU.setSite(Nl+1);cuU.setSite(Nl+1);cdU.setSite(Nl+1);cuDU.setSite(Nl+1);cdDU.setSite(Nl+1);id.setSite(Nl+1);p.setSite(Nl+1);
  nuL.setSite(Nl+1);ndL.setSite(Nl+1);cuL.setSite(Nl+1);cdL.setSite(Nl+1);cuDL.setSite(Nl+1);cdDL.setSite(Nl+1);pU.setSite(Nl+1);pL.setSite(Nl+1);

  M[ukey<2>(0,0)] = id;
  M[ukey<2>(1,0)] = cuDL;
  M[ukey<2>(2,0)] = cuL;
  M[ukey<2>(3,0)] = cuDU;
  M[ukey<2>(4,0)] = cuU;
  M[ukey<2>(5,0)] = cdDL;
  M[ukey<2>(6,0)] = cdL;
  M[ukey<2>(7,0)] = cdDU;
  M[ukey<2>(8,0)] = cdU;
  M[ukey<2>(9,0)] = Ubenzene*((nuL*ndL)+(nuU*ndU))+mubenzene*((nuL+ndL)+(nuU+ndU));

  M[ukey<2>(9,1)] = -hopbenzene*(cuL*p);
  M[ukey<2>(9,2)] = hopbenzene*(cuDL*p);//
  M[ukey<2>(9,3)] = -hopbenzene*(cuU*p);
  M[ukey<2>(9,4)] = hopbenzene*(cuDU*p);//
  M[ukey<2>(9,5)] = -hopbenzene*(cdL*p);
  M[ukey<2>(9,6)] = hopbenzene*(cdDL*p);//
  M[ukey<2>(9,7)] = -hopbenzene*(cdU*p);
  M[ukey<2>(9,8)] = hopbenzene*(cdDU*p);//
  M[ukey<2>(9,9)] = id;

  M.Dl=10,M.Dr=10;
  M.squeeze(1e-15);
  (*this)[Nl+1]=M;
  M.clear();

  nuU.setSite(Nl+2);ndU.setSite(Nl+2);cuU.setSite(Nl+2);cdU.setSite(Nl+2);cuDU.setSite(Nl+2);cdDU.setSite(Nl+2);id.setSite(Nl+2);p.setSite(Nl+2);
  nuL.setSite(Nl+2);ndL.setSite(Nl+2);cuL.setSite(Nl+2);cdL.setSite(Nl+2);cuDL.setSite(Nl+2);cdDL.setSite(Nl+2);pU.setSite(Nl+2);pL.setSite(Nl+2);

  M[ukey<2>(0,0)] = id;
  M[ukey<2>(1,0)] = cuDL;
  M[ukey<2>(2,0)] = cuL;
  M[ukey<2>(3,0)] = cuDU;
  M[ukey<2>(4,0)] = cuU;
  M[ukey<2>(5,0)] = cdDL;
  M[ukey<2>(6,0)] = cdL;
  M[ukey<2>(7,0)] = cdDU;
  M[ukey<2>(8,0)] = cdU;
  M[ukey<2>(9,0)] = Ubenzene*(nuL*ndL+nuU*ndU)+mubenzene*(nuL+ndL+nuU+ndU)-hopbenzene*(cuDL*cuU+cuDU*cuL+cdDL*cdU+cdDU*cdL);
  //M[ukey<2>(9,0)] = Ubenzene*(nuL*ndL+nuU*ndU)+mubenzene*(nuL+ndL+nuU+ndU);

  M[ukey<2>(9,1)] = hopbath(Nl+2)*(cuL*p);
  M[ukey<2>(9,2)] = -hopbath(Nl+2)*(cuDL*p);
  M[ukey<2>(9,3)] = hopbath(Nl+2)*(cdL*p);
  M[ukey<2>(9,4)] = -hopbath(Nl+2)*(cdDL*p);
  M[ukey<2>(9,5)] = id;

  M.Dl=10,M.Dr=6;
  M.squeeze(1e-15);
  (*this)[Nl+2]=M;
  M.clear();
  for(uint site = Nl+3;site<N-1;++site){
    nu.setSite(site);nd.setSite(site);cu.setSite(site);cd.setSite(site);cuD.setSite(site);cdD.setSite(site);idbath.setSite(site);pbath.setSite(site);
    M[ukey<2>(0,0)] = idbath;
    M[ukey<2>(1,0)] = cuD;
    M[ukey<2>(2,0)] = cu;
    M[ukey<2>(3,0)] = cdD;
    M[ukey<2>(4,0)] = cd;
    M[ukey<2>(5,0)] = mubathUp(site)*nu+mubathDown(site)*nd;
     
    M[ukey<2>(5,1)] = hopbath(site)*(cu*pbath);
    M[ukey<2>(5,2)] = -hopbath(site)*(cuD*pbath);
    M[ukey<2>(5,3)] = hopbath(site)*(cd*pbath);
    M[ukey<2>(5,4)] = -hopbath(site)*(cdD*pbath);
    M[ukey<2>(5,5)] = idbath;

    M.Dl=6;M.Dr=6;
    M.squeeze(1e-15);
    (*this)[site]=M;
    M.clear();
  }
  nu.setSite(N-1);nd.setSite(N-1);cu.setSite(N-1);cd.setSite(N-1);cuD.setSite(N-1);cdD.setSite(N-1);id.setSite(N-1);p.setSite(N-1);
  M[ukey<2>(0,0)] = idbath;
  M[ukey<2>(1,0)] = cuD;
  M[ukey<2>(2,0)] = cu;
  M[ukey<2>(3,0)] = cdD;
  M[ukey<2>(4,0)] = cd;
  M[ukey<2>(5,0)] = mubathUp(N-1)*nu+mubathDown(N-1)*nd;
  M.Dl=6;M.Dr=1;
  M.squeeze(1e-15);
  (*this)[N-1]=M;
}

template<typename T>
class twoOrbitalSIAMMPO:public MPO<T>
{
public:
  twoOrbitalSIAMMPO(std::string bath0,std::string bath1,std::string H_local);
  ~twoOrbitalSIAMMPO(){};
  virtual T measure(const MPS<AbKeyType,T>&mps)const {return MPO<T>::measure(mps);}
  virtual T measure(const CanonizedMPS<AbKeyType,T>&mps)const {return MPO<T>::measure(mps);}
  uint N0() const{return _N0;}
  uint N1() const{return _N1;}
  void shiftandrescale(const Real E0,const Real a, const Real W){
    SparseLocalOperator<Real> id(4,4);
    const int N=this->size();
    const Real E=E0*1.0/N;
    id[OpKeyType<1>(0,0)]=E;id[OpKeyType<1>(1,1)]=E;id[OpKeyType<1>(2,2)]=E;id[OpKeyType<1>(3,3)]=E;
    id.setSite(0);
    (*this)[0][ukey<2>(0,0)]-=id;
    for(uint site=1;site<N;++site){
      const uint lend=(*this)[site].Dl;
      id.setSite(site);
      (*this)[site][ukey<2>(lend-1,0)]-=id;
    }
    if(abs(W)>1e-14){
      id[OpKeyType<1>(0,0)]=Real(1.);id[OpKeyType<1>(1,1)]=Real(1.);
      const Real _w=W*1.0/N;
      id*=_w;
      id.setSite(0);
      (*this)[0][ukey<2>(0,0)]-=id;
      for(uint site=1;site<N;++site){
	const uint lend=(*this)[site].Dl;
	id.setSite(site);
	(*this)[site][ukey<2>(lend-1,0)]-=id ;
      }
    }
  }
  void shiftandrescale(const Complex E0,const Complex a, const Complex W){
    SparseLocalOperator<Complex> id(4,4);
    const int N=this->size();
    const Complex E=Complex(E0.real()*1.0/N,E0.imag()*1.0/N);
    id[OpKeyType<1>(0,0)]=Complex(E);id[OpKeyType<1>(1,1)]=Complex(E);id[OpKeyType<1>(2,2)]=Complex(E);id[OpKeyType<1>(3,3)]=Complex(E);
    id.setSite(0);
    (*this)[0][ukey<2>(0,0)]-=id;
    for(uint site=1;site<N;++site){
      const uint lend=(*this)[site].Dl;
      id.setSite(site);
      (*this)[site][ukey<2>(lend-1,0)]-=id;
    }
    if(abs(W)>1e-14){
      id[OpKeyType<1>(0,0)]=Complex(1.);id[OpKeyType<1>(1,1)]=Complex(1.);
      const Complex _w=Complex(W.real()*1.0/N,W.imag()*1.0/N);
      id*=_w;
      id.setSite(0);
      (*this)[0][ukey<2>(0,0)]-=id;
      for(uint site=1;site<N;++site){
	const uint lend=(*this)[site].Dl;
	id.setSite(site);
	(*this)[site][ukey<2>(lend-1,0)]-=id ;
      }
    }
  }

  template<class Archive>
  void serialize(Archive & ar, const unsigned int version){
    ar & boost::serialization::base_object<MPO<T> >(*this);
    ar &_hop;
    ar &_N0;
    ar &_N1;
  }

  VectorType _hop;
private:
  uint _N0,_N1;


};



template<typename T>
twoOrbitalSIAMMPO<T>::twoOrbitalSIAMMPO(std::string bath_0,std::string bath_1,std::string H0){
  std::fstream infile;
  if(!fexist(bath_0)){
    cout<<"#error reading "<<bath_0<<": no such file or directory"<<endl;
    abort();
  }
  infile.open(bath_0.c_str(),ios_base::in);
  if(!infile.good()){
    cout<<"#error reading "<<bath_0<<": file is broken"<<endl;
    abort();
  }
  cout<<"#reading data from "<<bath_0<<endl;
  int lbath;
  infile>>lbath;
  _N0=lbath+1;
  VectorType hop0(_N0-1),mu0(_N0);
  boost_setTo(mu0,0.0);boost_setTo(hop0,0.0);
  //================read hoppings from file ==========================
  Real hyb0,D0;  
  infile>>D0;
  infile>>hyb0;
  hop0(_N0-2)= hyb0;
  for(uint site=3;site<=_N0;++site){
    Real tn,en;
    infile>>tn;
    infile>>en;
    hop0(_N0-site)=tn;
    mu0(_N0-site+1)=en;
  }
  Real en;
  infile>>en;
  mu0(0)=en;
  infile.close();
  cout<<"#done..."<<endl;

  //================read the hopping parameters from infile========================
  if(!fexist(bath_1)){
    cout<<"#error reading "<<bath_1<<": no such file or directory"<<endl;
    abort();
  }
  infile.open(bath_1.c_str(),ios_base::in);
  if(!infile.good()){
    cout<<"#error reading "<<bath_1<<": file is broken"<<endl;
    abort();
  }
  cout<<"#reading data from "<<bath_1<<endl;
  infile>>lbath;
  _N1=lbath+1;
  VectorType hop1(_N1-1),mu1(_N1);
  boost_setTo(mu1,0.0);boost_setTo(hop1,0.0);
  //================read hoppings from file ==========================
  Real hyb1,D1;  
  infile>>D1;
  infile>>hyb1;
  hop1(0)= hyb1;
  for(uint site=1;site<_N1-1;++site){
    Real tn,en;
    infile>>tn;
    infile>>en;
    hop1(site)=tn;
    mu1(site)=en;
  }
  infile>>en;
  mu1(_N1-1)=en;
  infile.close();
  cout<<"#done..."<<endl;
  //####################################that has still to be done!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if(!fexist(H0)){
    cout<<"#error reading "<<H0<<": no such file or directory"<<endl;
    abort();
  }
  infile.open(H0.c_str(),ios_base::in);
  if(!infile.good()){
    cout<<"#error reading "<<H0<<": file is broken"<<endl;
    abort();
  }
  cout<<"#reading data from "<<H0<<endl;
  Real e0d,e0u,e1d,e1u,U,J;
  infile>>e0d;
  infile>>e0u;
  infile>>e1d;
  infile>>e1u;
  infile>>U;
  infile>>J;
  infile.close();
  cout<<"#done..."<<endl;
  uint N=_N0+_N1;
  
  VectorType hop(N-1),muu(N),mud(N);;
  boost_setTo(hop,0.0);boost_setTo(muu,0.0);boost_setTo(mud,0.0);

  for(uint n=0;n<N;++n){
    if(n<_N0-1){
      muu(n)=mu0(n);
      mud(n)=mu0(n);
      hop(n)=hop0(n);
    }else if(n>=_N0){
      if(n<N-1)
	hop(n)=hop1(n-_N0);
      if(n>_N0){
	muu(n)=mu1(n-_N0-1);
	mud(n)=mu1(n-_N0-1);
      }
    }
  }
  //dont forget the local mu's!!!!!
  muu(_N0-1)=e0u;
  mud(_N0-1)=e0d;
  muu(_N0)=e1u;
  mud(_N0)=e1d;
  cout<<_N0<<endl;
  cout<<_N1<<endl;
  cout<<"#following parameters will be used in the solver:"<<endl;
  cout<<"#t: "<<hop<<endl;
  cout<<"#muu: "<<muu<<endl;
  cout<<"#mud: "<<mud<<endl;
  cout<<"#U: "<<U<<endl;
  cout<<"#J: "<<J<<endl;
  _hop=hop;
  MPOMat<1,T> M;
  SparseLocalOperator<T> nu(4,4),nd(4,4),cuD(4,4),cdD(4,4),cu(4,4),cd(4,4),id(4,4),p(4,4);
  nu[OpKeyType<1>(2,2)]=T(1.0);nu[OpKeyType<1>(3,3)]=T(1.0);
  nd[OpKeyType<1>(1,1)]=T(1.0);nd[OpKeyType<1>(3,3)]=T(1.0);

  cuD[OpKeyType<1>(2,0)] = T(1.0);cuD[OpKeyType<1>(3,1)] = T(1.0); 
  cu[OpKeyType<1>(0,2)] = T(1.0);cu[OpKeyType<1>(1,3)] = T(1.0);   
								   
  cdD[OpKeyType<1>(1,0)] = T(1.0);cdD[OpKeyType<1>(3,2)] = T(-1.0);
  cd[OpKeyType<1>(0,1)] = T(1.0);cd[OpKeyType<1>(2,3)] = T(-1.0);  

  id[OpKeyType<1>(0,0)]=T(1);id[OpKeyType<1>(1,1)]=T(1);id[OpKeyType<1>(2,2)]=T(1);id[OpKeyType<1>(3,3)]=T(1);
  p[OpKeyType<1>(0,0)]=T(1);p[OpKeyType<1>(1,1)]=T(-1);p[OpKeyType<1>(2,2)]=T(-1);p[OpKeyType<1>(3,3)]=T(1);

  nu.setSite(0);nd.setSite(0);cu.setSite(0);cd.setSite(0);cuD.setSite(0);cdD.setSite(0);id.setSite(0);p.setSite(0);

  M[ukey<2>(0,0)] = muu(0)*nu+mud(0)*nd;
  M[ukey<2>(0,1)] = -hop(0)*(cuD*p);
  M[ukey<2>(0,2)] = -hop(0)*(cdD*p);
  M[ukey<2>(0,3)] = hop(0)*(cu*p);
  M[ukey<2>(0,4)] = hop(0)*(cd*p);
  M[ukey<2>(0,5)] = id;

  M.Dl=1;M.Dr=6;
  M.squeeze(1e-15);
  (*this)[0]=M;
  M.clear();


  for(uint site = 1;site<(_N0-1);++site){
    MPOMat<1,T> Mb;
    nu.setSite(site);nd.setSite(site);cu.setSite(site);cd.setSite(site);cuD.setSite(site);cdD.setSite(site);id.setSite(site);p.setSite(site);
    Mb[ukey<2>(0,0)] = id;
    Mb[ukey<2>(1,0)] = cu;
    Mb[ukey<2>(2,0)] = cd;
    Mb[ukey<2>(3,0)] = cuD;
    Mb[ukey<2>(4,0)] = cdD;
    Mb[ukey<2>(5,0)] = muu(site)*nu+mud(site)*nd;
     
    Mb[ukey<2>(5,1)] = -hop(site)*cuD*p;
    Mb[ukey<2>(5,2)] = -hop(site)*cdD*p;
    Mb[ukey<2>(5,3)] = hop(site)*cu*p;
    Mb[ukey<2>(5,4)] = hop(site)*cd*p;
    Mb[ukey<2>(5,5)] = id;
    Mb.Dl=6;Mb.Dr=6;
    Mb.squeeze(1e-15);
    (*this)[site]=Mb;
  }
  uint site=_N0-1;
  MPOMat<1,T> Mb;
  nu.setSite(site);nd.setSite(site);cu.setSite(site);cd.setSite(site);cuD.setSite(site);cdD.setSite(site);id.setSite(site);p.setSite(site);
  Mb[ukey<2>(0,0)] = id;
  Mb[ukey<2>(1,0)] = cu;
  Mb[ukey<2>(2,0)] = cd;
  Mb[ukey<2>(3,0)] = cuD;
  Mb[ukey<2>(4,0)] = cdD;
  Mb[ukey<2>(5,0)] = muu(site)*nu+mud(site)*nd+U*(nu*nd);
  
  //fill in the appropriate local impurity hamiltonian!!
  Mb[ukey<2>(5,1)] = -J*(cuD*cd);
  Mb[ukey<2>(5,2)] = -J*(cdD*cu);
  Mb[ukey<2>(5,3)] = -J*(cdD*cuD);
  Mb[ukey<2>(5,4)] = -J*(cd*cu);
  Mb[ukey<2>(5,5)] = nd;
  Mb[ukey<2>(5,6)] = nu;
  Mb[ukey<2>(5,7)] = id;
  Mb.Dl=6;Mb.Dr=8;
  Mb.squeeze(1e-15);
  (*this)[site]=Mb;
  Mb.clear();

  site=_N0;
  Real Up=U-2.*J;
  Real Upp=U-3.*J;

  nu.setSite(site);nd.setSite(site);cu.setSite(site);cd.setSite(site);cuD.setSite(site);cdD.setSite(site);id.setSite(site);p.setSite(site);
  Mb[ukey<2>(0,0)] = id;
  Mb[ukey<2>(1,0)] = cdD*cu;
  Mb[ukey<2>(2,0)] = cuD*cd;
  Mb[ukey<2>(3,0)] = cu*cd;
  Mb[ukey<2>(4,0)] = cuD*cdD;
  Mb[ukey<2>(5,0)] = Up*nu+Upp*nd;
  Mb[ukey<2>(6,0)] = Up*nd+Upp*nu;
  Mb[ukey<2>(7,0)] = muu(site)*nu+mud(site)*nd+U*(nu*nd);
  
  Mb[ukey<2>(7,1)] = -hop(site)*cuD*p;
  Mb[ukey<2>(7,2)] = -hop(site)*cdD*p;
  Mb[ukey<2>(7,3)] = hop(site)*cu*p;
  Mb[ukey<2>(7,4)] = hop(site)*cd*p;
  Mb[ukey<2>(7,5)] = id;
  Mb.Dl=8;Mb.Dr=6;
  Mb.squeeze(1e-15);
  (*this)[site]=Mb;
  Mb.clear();

  for(uint site = _N0+1;site<(N-1);++site){
    MPOMat<1,T> Mb;
    nu.setSite(site);nd.setSite(site);cu.setSite(site);cd.setSite(site);cuD.setSite(site);cdD.setSite(site);id.setSite(site);p.setSite(site);
    Mb[ukey<2>(0,0)] = id;
    Mb[ukey<2>(1,0)] = cu;
    Mb[ukey<2>(2,0)] = cd;
    Mb[ukey<2>(3,0)] = cuD;
    Mb[ukey<2>(4,0)] = cdD;
    Mb[ukey<2>(5,0)] = muu(site)*nu+mud(site)*nd;
     
    Mb[ukey<2>(5,1)] = -hop(site)*cuD*p;
    Mb[ukey<2>(5,2)] = -hop(site)*cdD*p;
    Mb[ukey<2>(5,3)] = hop(site)*cu*p;
    Mb[ukey<2>(5,4)] = hop(site)*cd*p;
    Mb[ukey<2>(5,5)] = id;
    Mb.Dl=6;Mb.Dr=6;
    Mb.squeeze(1e-15);
    (*this)[site]=Mb;
  }

  nu.setSite(N-1);nd.setSite(N-1);cu.setSite(N-1);cd.setSite(N-1);cuD.setSite(N-1);cdD.setSite(N-1);id.setSite(N-1);p.setSite(N-1);
  M[ukey<2>(0,0)] = id;
  M[ukey<2>(1,0)] = cu;
  M[ukey<2>(2,0)] = cd;
  M[ukey<2>(3,0)] = cuD;
  M[ukey<2>(4,0)] = cdD;
  M[ukey<2>(5,0)] = muu(N-1)*nu+mud(N-1)*nd;
  M.Dl=6;M.Dr=1;
  M.squeeze(1e-15);
  (*this)[N-1]=M;
  M.clear();
}


template<typename T>
class twoOrbitalSIAMMPO2:public MPO<T>
{
public:
  twoOrbitalSIAMMPO2(std::string bath,std::string H_local);
  ~twoOrbitalSIAMMPO2(){};
  virtual T measure(const MPS<AbKeyType,T>&mps)const {return MPO<T>::measure(mps);}
  virtual T measure(const CanonizedMPS<AbKeyType,T>&mps)const {return MPO<T>::measure(mps);}
  uint N0() const{return _N0;}
  uint N1() const{return _N1;}
  void shiftandrescale(const Real E0,const Real a, const Real W){
    SparseLocalOperator<Real> id(4,4);
    const int N=this->size();
    const Real E=E0*1.0/N;
    id[OpKeyType<1>(0,0)]=E;id[OpKeyType<1>(1,1)]=E;id[OpKeyType<1>(2,2)]=E;id[OpKeyType<1>(3,3)]=E;
    id.setSite(0);
    (*this)[0][ukey<2>(0,0)]-=id;
    for(uint site=1;site<N;++site){
      const uint lend=(*this)[site].Dl;
      id.setSite(site);
      (*this)[site][ukey<2>(lend-1,0)]-=id;
    }
    if(abs(W)>1e-14){
      id[OpKeyType<1>(0,0)]=Real(1.);id[OpKeyType<1>(1,1)]=Real(1.);
      const Real _w=W*1.0/N;
      id*=_w;
      id.setSite(0);
      (*this)[0][ukey<2>(0,0)]-=id;
      for(uint site=1;site<N;++site){
	const uint lend=(*this)[site].Dl;
	id.setSite(site);
	(*this)[site][ukey<2>(lend-1,0)]-=id ;
      }
    }
  }
  void shiftandrescale(const Complex E0,const Complex a, const Complex W){
    SparseLocalOperator<Complex> id(4,4);
    const int N=this->size();
    const Complex E=Complex(E0.real()*1.0/N,E0.imag()*1.0/N);
    id[OpKeyType<1>(0,0)]=Complex(E);id[OpKeyType<1>(1,1)]=Complex(E);id[OpKeyType<1>(2,2)]=Complex(E);id[OpKeyType<1>(3,3)]=Complex(E);
    id.setSite(0);
    (*this)[0][ukey<2>(0,0)]-=id;
    for(uint site=1;site<N;++site){
      const uint lend=(*this)[site].Dl;
      id.setSite(site);
      (*this)[site][ukey<2>(lend-1,0)]-=id;
    }
    if(abs(W)>1e-14){
      id[OpKeyType<1>(0,0)]=Complex(1.);id[OpKeyType<1>(1,1)]=Complex(1.);
      const Complex _w=Complex(W.real()*1.0/N,W.imag()*1.0/N);
      id*=_w;
      id.setSite(0);
      (*this)[0][ukey<2>(0,0)]-=id;
      for(uint site=1;site<N;++site){
	const uint lend=(*this)[site].Dl;
	id.setSite(site);
	(*this)[site][ukey<2>(lend-1,0)]-=id ;
      }
    }
  }

  template<class Archive>
  void serialize(Archive & ar, const unsigned int version){
    ar & boost::serialization::base_object<MPO<T> >(*this);
    ar &_hopu;
    ar &_hopd;
    ar &_muu;
    ar &_mud;
    ar &_N0;
    ar &_N1;
  }



  uint _N0,_N1;
  VectorType _hopu,_hopd,_muu,_mud;
};



template<typename T>
twoOrbitalSIAMMPO2<T>::twoOrbitalSIAMMPO2(std::string bath,std::string H0){
  cout<<"#reading data from "<<bath<<endl;
  MatrixType data=readMatrix(bath);
  if (data.size2()!=8){
    cerr<<"wrong number of elements in row "<<data.size2()<<" of bath-parameter file; expecting 8 columns"<<endl;
    abort();
  }
  //get the chain lenght 
  uint N0_0, N0_1, N1_0, N1_1,x=0;

  while((fabs(data(x,0))>1E-10)&&(x<data.size1())){
    x++;
  }
  N0_0=x;
  cout<<N0_0<<endl;
  x=0;
  while((fabs(data(x,4))>1E-10)&&(x<data.size1())){
    x++;
  }
  N0_1=x;
  cout<<N0_1<<endl;
  this->_N0=std::max(N0_0,N0_1);

  x=0;
  while((fabs(data(x,2))>1E-10)&&(x<data.size1())){
    x++;
  }
  N1_0=x;

  x=0;
  while((fabs(data(x,6))>1E-10)&&(x<data.size1())){
    x++;
  }
  N1_1=x;
  
  this->_N1=std::max(N1_0,N1_1);
  

  VectorType hopU0(_N0),muU0(_N0),hopD0(_N0),muD0(_N0), hopU1(_N1),muU1(_N1),hopD1(_N1),muD1(_N1);
  boost_setTo(muU0,0.0);boost_setTo(hopU0,0.0);boost_setTo(muD0,0.0);boost_setTo(hopD0,0.0);
  boost_setTo(muU1,0.0);boost_setTo(hopU1,0.0);boost_setTo(muD1,0.0);boost_setTo(hopD1,0.0);
  //================read hoppings from file ==========================
  cout<<"at line 5105 in Hamiltonians: check if the hybs are associated with the correct orbitals"<<endl;
  for(uint site=0;site<_N0;++site){
    hopD0(_N0-1-site)=data(site,0);
    hopU0(_N0-1-site)=data(site,4);
    muD0(_N0-1-site)=data(site,1);
    muU0(_N0-1-site)=data(site,5);
  }
  for(uint site=0;site<_N1;++site){
    hopD1(site)=data(site,2);
    hopU1(site)=data(site,6);
    muD1(site)=data(site,3);
    muU1(site)=data(site,7);
  }


  if(!fexist(H0)){
    cout<<"#error reading "<<H0<<": no such file or directory"<<endl;
    abort();
  }

  std::ifstream infile;
  infile.open(H0.c_str(),ios_base::in);
  if(!infile.good()){
    cout<<"#error reading "<<H0<<": file is broken"<<endl;
    abort();
  }
  cout<<"#reading local Hamiltonian from "<<H0<<endl;
  Real e0d,e0u,e1d,e1u,U,J;
  infile>>e0d;
  infile>>e1d;
  infile>>e0u;
  infile>>e1u;
  infile>>U;
  infile>>J;
  infile.close();

  cout<<"#done..."<<endl;
  uint N=_N0+_N1+2;

  _hopu.resize(N-1),_hopd.resize(N-1),_muu.resize(N),_mud.resize(N);
  boost_setTo(_hopu,0.0);boost_setTo(_hopd,0.0);boost_setTo(_muu,0.0);boost_setTo(_mud,0.0);

  for(uint n=0;n<_N0;++n){
    _muu(n)=muU0(n);
    _mud(n)=muD0(n);
    _hopu(n)=hopU0(n);
    _hopd(n)=hopD0(n);
  }
  for(uint n=_N0+1;n<(N-1);++n){
    _muu(n+1)=muU1(n-(_N0+1));
    _mud(n+1)=muD1(n-(_N0+1));
    _hopu(n)=hopU1(n-(_N0+1));
    _hopd(n)=hopD1(n-(_N0+1));

  }
  //dont forget the local mu's
  _muu(_N0)=e0u;
  _mud(_N0)=e0d;
  _muu(_N0+1)=e1u;
  _mud(_N0+1)=e1d;
  //cout<<_N0<<endl;
  //cout<<_N1<<endl;
  //
  //cout<<"#following parameters will be used in the solver:"<<endl;
  //cout<<"hopping up-spins"<<endl;
  //cout<<_hopu<<endl;
  //cout<<"hopping down-spins"<<endl;
  //cout<<_hopd<<endl;
  //cout<<"mu up-spins"<<endl;
  //cout<<_muu<<endl;
  //cout<<"mu down-spins"<<endl;
  //cout<<_mud<<endl;
  //cout<<"#U: "<<U<<endl;
  //cout<<"#J: "<<J<<endl;

  MPOMat<1,T> M;
  SparseLocalOperator<T> nu(4,4),nd(4,4),cuD(4,4),cdD(4,4),cu(4,4),cd(4,4),id(4,4),p(4,4);
  nu[OpKeyType<1>(2,2)]=T(1.0);nu[OpKeyType<1>(3,3)]=T(1.0);
  nd[OpKeyType<1>(1,1)]=T(1.0);nd[OpKeyType<1>(3,3)]=T(1.0);
  
  cuD[OpKeyType<1>(2,0)] = T(1.0);cuD[OpKeyType<1>(3,1)] = T(1.0); 
  cu[OpKeyType<1>(0,2)] = T(1.0);cu[OpKeyType<1>(1,3)] = T(1.0);   
  								   
  cdD[OpKeyType<1>(1,0)] = T(1.0);cdD[OpKeyType<1>(3,2)] = T(-1.0);
  cd[OpKeyType<1>(0,1)] = T(1.0);cd[OpKeyType<1>(2,3)] = T(-1.0);  
  
  id[OpKeyType<1>(0,0)]=T(1);id[OpKeyType<1>(1,1)]=T(1);id[OpKeyType<1>(2,2)]=T(1);id[OpKeyType<1>(3,3)]=T(1);
  p[OpKeyType<1>(0,0)]=T(1);p[OpKeyType<1>(1,1)]=T(-1);p[OpKeyType<1>(2,2)]=T(-1);p[OpKeyType<1>(3,3)]=T(1);
  
  nu.setSite(0);nd.setSite(0);cu.setSite(0);cd.setSite(0);cuD.setSite(0);cdD.setSite(0);id.setSite(0);p.setSite(0);
  
  M[ukey<2>(0,0)] = _muu(0)*nu+_mud(0)*nd;
  M[ukey<2>(0,1)] = -_hopu(0)*(cuD*p);
  M[ukey<2>(0,2)] = -_hopd(0)*(cdD*p);
  M[ukey<2>(0,3)] = _hopu(0)*(cu*p);
  M[ukey<2>(0,4)] = _hopd(0)*(cd*p);
  M[ukey<2>(0,5)] = id;
  
  M.Dl=1;M.Dr=6;
  M.squeeze(1e-15);
  (*this)[0]=M;
  M.clear();
  
  
  for(uint site = 1;site<_N0;++site){
    MPOMat<1,T> Mb;
    nu.setSite(site);nd.setSite(site);cu.setSite(site);cd.setSite(site);cuD.setSite(site);cdD.setSite(site);id.setSite(site);p.setSite(site);
    Mb[ukey<2>(0,0)] = id;
    Mb[ukey<2>(1,0)] = cu;
    Mb[ukey<2>(2,0)] = cd;
    Mb[ukey<2>(3,0)] = cuD;
    Mb[ukey<2>(4,0)] = cdD;
    Mb[ukey<2>(5,0)] = _muu(site)*nu+_mud(site)*nd;
   
    Mb[ukey<2>(5,1)] = -_hopu(site)*cuD*p;
    Mb[ukey<2>(5,2)] = -_hopd(site)*cdD*p;
    Mb[ukey<2>(5,3)] = _hopu(site)*cu*p;
    Mb[ukey<2>(5,4)] = _hopd(site)*cd*p;
    Mb[ukey<2>(5,5)] = id;
    Mb.Dl=6;Mb.Dr=6;
    Mb.squeeze(1e-15);
    (*this)[site]=Mb;
  }
  uint site=_N0;
  MPOMat<1,T> Mb;
  nu.setSite(site);nd.setSite(site);cu.setSite(site);cd.setSite(site);cuD.setSite(site);cdD.setSite(site);id.setSite(site);p.setSite(site);
  Mb[ukey<2>(0,0)] = id;
  Mb[ukey<2>(1,0)] = cu;
  Mb[ukey<2>(2,0)] = cd;
  Mb[ukey<2>(3,0)] = cuD;
  Mb[ukey<2>(4,0)] = cdD;
  Mb[ukey<2>(5,0)] = _muu(site)*nu+_mud(site)*nd+U*(nu*nd);
 
  //fill in the appropriate local impurity hamiltonian!!
  Mb[ukey<2>(5,1)] = -J*(cuD*cd);
  Mb[ukey<2>(5,2)] = -J*(cdD*cu);
  Mb[ukey<2>(5,3)] = -J*(cdD*cuD);
  Mb[ukey<2>(5,4)] = -J*(cd*cu);
  Mb[ukey<2>(5,5)] = nd;
  Mb[ukey<2>(5,6)] = nu;
  Mb[ukey<2>(5,7)] = id;
  Mb.Dl=6;Mb.Dr=8;
  Mb.squeeze(1e-15);
  (*this)[site]=Mb;
  Mb.clear();

  site=_N0+1;
  Real Up=U-2.*J;
  Real Upp=U-3.*J;
  nu.setSite(site);nd.setSite(site);cu.setSite(site);cd.setSite(site);cuD.setSite(site);cdD.setSite(site);id.setSite(site);p.setSite(site);
  Mb[ukey<2>(0,0)] = id;
  Mb[ukey<2>(1,0)] = cdD*cu;
  Mb[ukey<2>(2,0)] = cuD*cd;
  Mb[ukey<2>(3,0)] = cu*cd;
  Mb[ukey<2>(4,0)] = cuD*cdD;
  Mb[ukey<2>(5,0)] = Up*nu+Upp*nd;
  Mb[ukey<2>(6,0)] = Up*nd+Upp*nu;
  Mb[ukey<2>(7,0)] = _muu(site)*nu+_mud(site)*nd+U*(nu*nd);
  
  Mb[ukey<2>(7,1)] = -_hopu(site)*cuD*p;
  Mb[ukey<2>(7,2)] = -_hopd(site)*cdD*p;
  Mb[ukey<2>(7,3)] = _hopu(site)*cu*p;
  Mb[ukey<2>(7,4)] = _hopd(site)*cd*p;
  Mb[ukey<2>(7,5)] = id;
  Mb.Dl=8;Mb.Dr=6;
  Mb.squeeze(1e-15);
  (*this)[site]=Mb;
  Mb.clear();
  
  for(uint site = _N0+2;site<(N-1);++site){
    MPOMat<1,T> Mb;
    nu.setSite(site);nd.setSite(site);cu.setSite(site);cd.setSite(site);cuD.setSite(site);cdD.setSite(site);id.setSite(site);p.setSite(site);
    Mb[ukey<2>(0,0)] = id;
    Mb[ukey<2>(1,0)] = cu;
    Mb[ukey<2>(2,0)] = cd;
    Mb[ukey<2>(3,0)] = cuD;
    Mb[ukey<2>(4,0)] = cdD;
    Mb[ukey<2>(5,0)] = _muu(site)*nu+_mud(site)*nd;
   
    Mb[ukey<2>(5,1)] = -_hopu(site)*cuD*p;
    Mb[ukey<2>(5,2)] = -_hopd(site)*cdD*p;
    Mb[ukey<2>(5,3)] = _hopu(site)*cu*p;
    Mb[ukey<2>(5,4)] = _hopd(site)*cd*p;
    Mb[ukey<2>(5,5)] = id;
    Mb.Dl=6;Mb.Dr=6;
    Mb.squeeze(1e-15);
    (*this)[site]=Mb;
  }
  
  nu.setSite(N-1);nd.setSite(N-1);cu.setSite(N-1);cd.setSite(N-1);cuD.setSite(N-1);cdD.setSite(N-1);id.setSite(N-1);p.setSite(N-1);
  M[ukey<2>(0,0)] = id;
  M[ukey<2>(1,0)] = cu;
  M[ukey<2>(2,0)] = cd;
  M[ukey<2>(3,0)] = cuD;
  M[ukey<2>(4,0)] = cdD;
  M[ukey<2>(5,0)] = _muu(N-1)*nu+_mud(N-1)*nd;
  M.Dl=6;M.Dr=1;
  M.squeeze(1e-15);
  (*this)[N-1]=M;
  M.clear();
}





//this is an mpo for decoupled spinless fermions, which is equivalent to two chains of spinful fermions which interact only with its own type (i.e. only up with up and down with down)
template<typename T>
class DecoupledFermionsMPO:public MPO<T>
{
public:
  //hop: hopping, U: local interaction; chUp(chDown: chemical potential for up/down-spin electrons;
  DecoupledFermionsMPO(VectorType hop,VectorType V,VectorType chempot);
  //V: next-nearest neighbor interaction
  ~DecoupledFermionsMPO(){};
  virtual T measure(const MPS<AbKeyType,T>&mps)const {return MPO<T>::measure(mps);}
  virtual T measure(const CanonizedMPS<AbKeyType,T>&mps)const {return MPO<T>::measure(mps);}
};

template<typename T>
DecoupledFermionsMPO<T>::DecoupledFermionsMPO(VectorType hop,VectorType V,VectorType chempot){
  uint N =chempot.size();
  MPOMat<1,T> M;
  //  this->resize(N);
  SparseLocalOperator<T> nu(4,4),nd(4,4),cuD(4,4),cdD(4,4),cu(4,4),cd(4,4),id(4,4),p(4,4);//all statistics are set to bosonic (_statistics=1)
  nu[OpKeyType<1>(2,2)]=T(1.0);nu[OpKeyType<1>(3,3)]=T(1.0);
  nd[OpKeyType<1>(1,1)]=T(1.0);nd[OpKeyType<1>(3,3)]=T(1.0);

  cuD[OpKeyType<1>(2,0)] = T(1.0);cuD[OpKeyType<1>(3,1)] = T(1.0); 
  cu[OpKeyType<1>(0,2)] = T(1.0);cu[OpKeyType<1>(1,3)] = T(1.0);   

  cdD[OpKeyType<1>(1,0)] = T(1.0);cdD[OpKeyType<1>(3,2)] = T(-1.0);
  cd[OpKeyType<1>(0,1)] = T(1.0);cd[OpKeyType<1>(2,3)] = T(-1.0);  

  id[OpKeyType<1>(0,0)]=T(1);id[OpKeyType<1>(1,1)]=T(1);id[OpKeyType<1>(2,2)]=T(1);id[OpKeyType<1>(3,3)]=T(1);
  p[OpKeyType<1>(0,0)]=T(1);p[OpKeyType<1>(1,1)]=T(-1);p[OpKeyType<1>(2,2)]=T(-1);p[OpKeyType<1>(3,3)]=T(1);
//p[OpKeyType<1>(0,0)]=T(1);p[OpKeyType<1>(1,1)]=T(1);p[OpKeyType<1>(2,2)]=T(-1);p[OpKeyType<1>(3,3)]=T(1);

  nu.setSite(0);nd.setSite(0);cu.setSite(0);cd.setSite(0);cuD.setSite(0);cdD.setSite(0);id.setSite(0);p.setSite(0);

  M[ukey<2>(0,0)] = chempot(0)*nu+chempot(0)*nd;
  M[ukey<2>(0,1)] = -hop(0)*(cuD*p);
  M[ukey<2>(0,2)] = hop(0)*(cdD*p);
  M[ukey<2>(0,3)] = hop(0)*(cu*p);//.
  M[ukey<2>(0,4)] = -hop(0)*(cd*p);
  M[ukey<2>(0,5)] = V(0)*nu;
  M[ukey<2>(0,6)] = -V(0)*nd;
  M[ukey<2>(0,7)] = id;

  M.Dl=1;M.Dr=8;
  M.squeeze(1e-15);
  (*this)[0]=M;
  M.clear();

  nu.setSite(N-1);nd.setSite(N-1);cu.setSite(N-1);cd.setSite(N-1);cuD.setSite(N-1);cdD.setSite(N-1);id.setSite(N-1);p.setSite(N-1);

  M[ukey<2>(0,0)] = id;
  M[ukey<2>(1,0)] = cu;
  M[ukey<2>(2,0)] = cd;
  M[ukey<2>(3,0)] = cuD;
  M[ukey<2>(4,0)] = cdD;
  M[ukey<2>(5,0)] = nu;
  M[ukey<2>(6,0)] = nd;
  M[ukey<2>(7,0)] = chempot(N-1)*nu+chempot(N-1)*nd;
  M.Dl=8;M.Dr=1;
  M.squeeze(1e-15);
  (*this)[N-1]=M;

  for(uint site = 1;site<(N-1);++site){
    MPOMat<1,T> Mb;
    nu.setSite(site);nd.setSite(site);cu.setSite(site);cd.setSite(site);cuD.setSite(site);cdD.setSite(site);id.setSite(site);p.setSite(site);
    Mb[ukey<2>(0,0)] = id;
    Mb[ukey<2>(1,0)] = cu;
    Mb[ukey<2>(2,0)] = cd;
    Mb[ukey<2>(3,0)] = cuD;
    Mb[ukey<2>(4,0)] = cdD;
    Mb[ukey<2>(5,0)] = nu;
    Mb[ukey<2>(6,0)] = nd;
    Mb[ukey<2>(7,0)] = chempot(site)*nu+chempot(site)*nd;
     
    Mb[ukey<2>(7,1)] = -hop(site)*cuD*p;
    Mb[ukey<2>(7,2)] = hop(site)*cdD*p;
    Mb[ukey<2>(7,3)] = hop(site)*cu*p;
    Mb[ukey<2>(7,4)] = -hop(site)*cd*p;
    Mb[ukey<2>(7,5)] = V(site)*nu;
    Mb[ukey<2>(7,6)] = -V(site)*nd;
    Mb[ukey<2>(7,7)] = id;
    Mb.Dl=8;Mb.Dr=8;
    Mb.squeeze(1e-15);
    (*this)[site]=Mb;
  }
}



template<typename T>
class TransverseIsingNNNMPO:public MPO<T>{
public:
  TransverseIsingNNNMPO(VectorType Jz,VectorType Jz2,VectorType B);
  ~TransverseIsingNNNMPO(){};
  void shiftandrescale(const T E0,const T a, const T W);
};


template<typename T>
void TransverseIsingNNNMPO<T>::shiftandrescale(const T E0,const T a, const T W){
  SparseLocalOperator<T> id(2,2);
  id[OpKeyType<1>(0,0)]=T(1.);id[OpKeyType<1>(1,1)]=T(1.);
  id*=E0;
  cout<<(*this)[0]<<endl;
  (*this)[0][ukey<2>(0,0)]-=id;
  typename MPOMat<1,T>::iterator o;
  for(o=(*this)[0].begin();o!=(*this)[0].end();++o)
    o->second/=a;

  id[OpKeyType<1>(0,0)]=T(1.);id[OpKeyType<1>(1,1)]=T(1.);
  id*=W;
  (*this)[0][ukey<2>(0,0)]-=id;
  cout<<(*this)[0]<<endl;
  //that's it ...
}

template<typename T>
TransverseIsingNNNMPO<T>::TransverseIsingNNNMPO(VectorType Jz1,VectorType Jz2,VectorType B) {
  uint N = B.size();
  MPOMat<1,T> M;
  SparseLocalOperator<T> Sz(2,2),Sx(2,2),id(2,2);
  Sz[OpKeyType<1>(0,0)]=T(-1.0);Sz[OpKeyType<1>(1,1)]=T(1.0);
  Sx[OpKeyType<1>(0,1)]=T(1.0);Sx[OpKeyType<1>(1,0)]=T(1.0);
  id[OpKeyType<1>(0,0)]=T(1);id[OpKeyType<1>(1,1)]=T(1);

  Sz.setSite(0);Sx.setSite(0);id.setSite(N-1);
  M[ukey<2>(0,0)] = B(0)*Sx;
  M[ukey<2>(0,1)] = Jz1(0)*Sz;
  M[ukey<2>(0,2)] = Jz2(0)*Sz;
  M[ukey<2>(0,3)] = id;
  M.Dl=1,M.Dr=4;
  M.squeeze(1e-16);
  (*this)[0]=M;

  M.clear();
  Sz.setSite(N-1);Sx.setSite(N-1);id.setSite(N-1);
  M[ukey<2>(0,0)] = id;
  M[ukey<2>(1,0)] = Sz;
  M[ukey<2>(3,0)] = B(N-1)*Sx;
  M.Dl=4,M.Dr=1;
  M.squeeze(1e-16);
  (*this)[N-1]=M;

  for(uint site = 1;site<(N-1);++site){
    MPOMat<1,T> Mb;
    Sz.setSite(site);Sx.setSite(site);id.setSite(site);
    Mb[ukey<2>(0,0)] = id;
    Mb[ukey<2>(1,0)] = Sz;
    Mb[ukey<2>(2,1)] = id;
    Mb[ukey<2>(3,0)] = B(site)*Sx;
    Mb[ukey<2>(3,1)] = Jz1(site)*Sz;
    Mb[ukey<2>(3,2)] = Jz2(site)*Sz;
    Mb[ukey<2>(3,3)] = id;

    Mb.Dl=4,Mb.Dr=4;
    Mb.squeeze(1e-16);
    (*this)[site]=Mb;
  }
}

//definition of bose hubbard model:
//-\sum_i hop(i)b_i^{\dagger}b_{i+1}+h.c. + \sum_i U(i)/2 n_i(n_i-1)
template<typename T>
class LiebLinigerMPO:public MPO<T>
{
public:
  LiebLinigerMPO(const Real mass, const Real g, const VectorType mu, const uint maxdim,const Real dx);
  ~LiebLinigerMPO(){};
  virtual T measure(const MPS<AbKeyType,T>&mps)const {return MPO<T>::measure(mps);}
  virtual T measure(const CanonizedMPS<AbKeyType,T>&mps)const {return MPO<T>::measure(mps);}

};

template<typename T>
LiebLinigerMPO<T>::LiebLinigerMPO(const Real mass, const Real g, const VectorType mu, const uint maxdim,const Real dx){
  uint N =mu.size();
  MPOMat<1,T> M;

  SparseLocalOperator<T> c(maxdim,maxdim),cd(maxdim,maxdim),id(maxdim,maxdim),n(maxdim,maxdim),cdcdcc(maxdim,maxdim);
  for(int ind = 0;ind<maxdim;++ind){
    n[OpKeyType<1>(ind,ind)]=T(ind*1.0);
    id[OpKeyType<1>(ind,ind)]=T(1.0);
    if(ind!=0)
      c[OpKeyType<1>(ind-1,ind)] = T(sqrt(ind));
    if(ind<(maxdim-1))
      cd[OpKeyType<1>(ind+1,ind)] = T(sqrt(ind+1));
  }
  cdcdcc=(((cd*cd)*c)*c);
  cout<<cdcdcc;
  cin.get();
  c.setSite(0);cd.setSite(0);id.setSite(0);n.setSite(0);cdcdcc.setSite(0);

  M[ukey<2>(0,0)] = ((mu(0)+1.0/(mass*(dx*dx)))*n)+(g/dx)*cdcdcc;
  M[ukey<2>(0,1)] = -1.0/(2.0*mass*(dx*dx))*c;
  M[ukey<2>(0,2)] = -1.0/(2.0*mass*(dx*dx))*cd;
  M[ukey<2>(0,3)] = id;
  M.Dl=1,M.Dr=4;
  M.squeeze(1e-15);
  (*this)[0]=M;
  M.clear();


  for(uint site = 1;site<(N-1);++site){
    MPOMat<1,T> M;
    c.setSite(site);cd.setSite(site);id.setSite(site);n.setSite(site),cdcdcc.setSite(site);

    M[ukey<2>(0,0)]= id;
    M[ukey<2>(1,0)]= cd;
    M[ukey<2>(2,0)] = c;
    M[ukey<2>(3,0)] = ((mu(site)+1.0/(mass*(dx*dx)))*n)+(g/dx)*cdcdcc;
     
    M[ukey<2>(3,1)] = -1.0/(2.0*mass*(dx*dx))*c;
    M[ukey<2>(3,2)] = -1.0/(2.0*mass*(dx*dx))*cd;
    M[ukey<2>(3,3)] = id;
    M.Dl=4,M.Dr=4;
    M.squeeze(1e-15);
    (*this)[site]=M;
    M.clear();
  }
  
  c.setSite(N-1);cd.setSite(N-1);id.setSite(N-1);n.setSite(N-1);cdcdcc.setSite(N-1);
  M[ukey<2>(0,0)]= id;
  M[ukey<2>(1,0)]= cd;
  M[ukey<2>(2,0)] = c;
  M[ukey<2>(3,0)] = ((mu(N-1)+1.0/(mass*(dx*dx)))*n)+(g/dx)*cdcdcc;

  M.Dl=4,M.Dr=1;
  M.squeeze(1e-15);
  (*this)[N-1]=M;

}





//ordering of local states at the measuring sites (particle states and spin states): |0,down>,|1,down>,|0,up>,|1,up>
template<typename T>
class ParticleDetectorMPO:public MPO<T>
{
public:
  ParticleDetectorMPO(const VectorType hop,const VectorType mu,const VectorType coupling,const UintVectorType couplingPositions,\
		      const uint maxdim);
  
  ~ParticleDetectorMPO(){};
  virtual T measure(const MPS<AbKeyType,T>&mps)const {return MPO<T>::measure(mps);}
  virtual T measure(const CanonizedMPS<AbKeyType,T>&mps)const {return MPO<T>::measure(mps);}

};

template<typename T>
ParticleDetectorMPO<T>::ParticleDetectorMPO(const VectorType hop,const VectorType mu,const VectorType coupling,const UintVectorType couplingPositions,const uint maxdim){
  uint N =mu.size();
  MPOMat<1,T> M;
  assert(maxdim>1);
  SparseLocalOperator<T> b(maxdim,maxdim),bd(maxdim,maxdim),id(maxdim,maxdim),n(maxdim,maxdim);
  SparseLocalOperator<T> Sp(2,2),Sm(2,2),idS(2,2);
  Sp[OpKeyType<1>(1,0)]=T(1.0);Sm[OpKeyType<1>(0,1)]=T(1.0);
  idS[OpKeyType<1>(0,0)]=T(1.0);idS[OpKeyType<1>(1,1)]=T(1.0);


  for(int ind = 0;ind<maxdim;++ind){
    n[OpKeyType<1>(ind,ind)]=T(ind*1.0);
    id[OpKeyType<1>(ind,ind)]=T(1.0);
    if(ind!=0)
      b[OpKeyType<1>(ind-1,ind)] = T(sqrt(ind));
    if(ind<(maxdim-1))
      bd[OpKeyType<1>(ind+1,ind)] = T(sqrt(ind+1));
  }


  SparseLocalOperator<T> n_sp=kronecker(n,Sp);
  SparseLocalOperator<T> n_sm=kronecker(n,Sm);
  SparseLocalOperator<T> id_id=kronecker(id,idS);
  SparseLocalOperator<T> n_id=kronecker(n,idS);
  SparseLocalOperator<T> b_id=kronecker(b,idS);
  SparseLocalOperator<T> bd_id=kronecker(bd,idS);

  b.setSite(0);bd.setSite(0);id.setSite(0);n.setSite(0);
  M[ukey<2>(0,0)] = -mu(0)*n;
  M[ukey<2>(0,1)] = -hop(0)*b;
  M[ukey<2>(0,2)] = -hop(0)*bd;
  M[ukey<2>(0,3)] = id;
  M.Dl=1,M.Dr=4;
  M.squeeze(1e-15);
  (*this)[0]=M;
  M.clear();

  b.setSite(N-1);bd.setSite(N-1);id.setSite(N-1);n.setSite(N-1);

  M[ukey<2>(0,0)] = id;
  M[ukey<2>(1,0)] = bd;
  M[ukey<2>(2,0)] = b;
  M[ukey<2>(3,0)] = -mu(N-1)*n;

  M.Dl=4,M.Dr=1;
  M.squeeze(1e-15);
  (*this)[N-1]=M;
  
  uint position=0;
  uint currentSite=couplingPositions(position);
  assert(currentSite>0);
  assert(couplingPositions(couplingPositions.size()-1)!=N);
  for(uint site = 1;site<(N-1);++site){
    MPOMat<1,T> Mb;
    if (site!=currentSite){
      b.setSite(site);bd.setSite(site);id.setSite(site);n.setSite(site);
      Mb[ukey<2>(0,0)]= id;
      Mb[ukey<2>(1,0)]= bd;
      Mb[ukey<2>(2,0)] = b;
      Mb[ukey<2>(3,0)] = mu(site)*n;
     
      Mb[ukey<2>(3,1)] = hop(site)*b;
      Mb[ukey<2>(3,2)] = hop(site)*bd;
      Mb[ukey<2>(3,3)] = id;
      Mb.Dl=4,Mb.Dr=4;
      Mb.squeeze(1e-15);
      (*this)[site]=Mb;
    }


    if (site==currentSite){
      n_sp.setSite(site),n_sm.setSite(site),id_id.setSite(site),b_id.setSite(site),bd_id.setSite(site);
      Mb[ukey<2>(0,0)]=id_id;
      Mb[ukey<2>(1,0)]=bd_id;
      Mb[ukey<2>(2,0)]=b_id;
      Mb[ukey<2>(3,0)]=(mu(site)*n_id)+(coupling(position)*n_sp)+(conj(coupling(position))*n_sm);
     
      Mb[ukey<2>(3,1)] = hop(site)*b_id;
      Mb[ukey<2>(3,2)] = hop(site)*bd_id;
      Mb[ukey<2>(3,3)] = id_id;
      Mb.Dl=4,Mb.Dr=4;
      Mb.squeeze(1e-15);
      (*this)[site]=Mb;
      if(position<couplingPositions.size()-1){
	position++;
	currentSite=couplingPositions(position);
      }
    }

  }
}

//ordering of local states at the measuring sites (particle states and spin states): |0,down>,|1,down>,|0,up>,|1,up>
template<typename T>
class NMRMPO:public MPO<T>
{
public:
  NMRMPO(const VectorType hop,const VectorType mu,const VectorType coupling,const VectorType Bz,const UintVectorType couplingPositions,const VectorType x,const Real omega,
		      const uint maxdim);
  
  ~NMRMPO(){};
  virtual T measure(const MPS<AbKeyType,T>&mps)const {return MPO<T>::measure(mps);}
  virtual T measure(const CanonizedMPS<AbKeyType,T>&mps)const {return MPO<T>::measure(mps);}

};

template<typename T>
NMRMPO<T>::NMRMPO(const VectorType hop,const VectorType mu,const VectorType coupling,const VectorType Bz,const UintVectorType couplingPositions,const VectorType x,const Real omega,const uint maxdim){
  uint N =mu.size();
  MPOMat<1,T> M;
  assert(maxdim>1);
  SparseLocalOperator<T> b(maxdim,maxdim),bd(maxdim,maxdim),id(maxdim,maxdim),n(maxdim,maxdim);
  SparseLocalOperator<T> Sp(2,2),Sm(2,2),idS(2,2),Sz(2,2),Sx(2,2);
  Sz[OpKeyType<1>(0,0)]=T(-0.5);Sz[OpKeyType<1>(1,1)]=T(0.5);
  Sx[OpKeyType<1>(0,1)]=T(0.5);Sx[OpKeyType<1>(1,0)]=T(0.5);
  Sp[OpKeyType<1>(1,0)]=T(1.0);Sm[OpKeyType<1>(0,1)]=T(1.0);
  idS[OpKeyType<1>(0,0)]=T(1.0);idS[OpKeyType<1>(1,1)]=T(1.0);


  for(int ind = 0;ind<maxdim;++ind){
    n[OpKeyType<1>(ind,ind)]=T(ind*1.0);
    id[OpKeyType<1>(ind,ind)]=T(1.0);
    if(ind!=0)
      b[OpKeyType<1>(ind-1,ind)] = T(sqrt(ind));
    if(ind<(maxdim-1))
      bd[OpKeyType<1>(ind+1,ind)] = T(sqrt(ind+1));
  }


  SparseLocalOperator<T> n_sp=kronecker(n,Sp);
  SparseLocalOperator<T> n_sm=kronecker(n,Sm);
  SparseLocalOperator<T> id_sz=kronecker(id,Sz);
  SparseLocalOperator<T> id_sp=kronecker(id,Sp);
  SparseLocalOperator<T> id_sm=kronecker(id,Sm);
  //SparseLocalOperator<T> id_sx=kronecker(id,Sx);
  SparseLocalOperator<T> id_id=kronecker(id,idS);
  SparseLocalOperator<T> n_id=kronecker(n,idS);
  SparseLocalOperator<T> b_id=kronecker(b,idS);
  SparseLocalOperator<T> bd_id=kronecker(bd,idS);

  b.setSite(0);bd.setSite(0);id.setSite(0);n.setSite(0);
  M[ukey<2>(0,0)] = -mu(0)*n;
  M[ukey<2>(0,1)] = -hop(0)*b;
  M[ukey<2>(0,2)] = -hop(0)*bd;
  M[ukey<2>(0,3)] = id;
  M.Dl=1,M.Dr=4;
  M.squeeze(1e-15);
  (*this)[0]=M;
  M.clear();

  b.setSite(N-1);bd.setSite(N-1);id.setSite(N-1);n.setSite(N-1);

  M[ukey<2>(0,0)] = id;
  M[ukey<2>(1,0)] = bd;
  M[ukey<2>(2,0)] = b;
  M[ukey<2>(3,0)] = -mu(N-1)*n;

  M.Dl=4,M.Dr=1;
  M.squeeze(1e-15);
  (*this)[N-1]=M;
  
  uint position=0;
  uint currentSite=couplingPositions(position);
  assert(currentSite>0);
  assert(couplingPositions(couplingPositions.size()-1)!=N);
  for(uint site = 1;site<(N-1);++site){

    MPOMat<1,T> Mb;
    if (site!=currentSite){
      b.setSite(site);bd.setSite(site);id.setSite(site);n.setSite(site);
      Mb[ukey<2>(0,0)]= id;
      Mb[ukey<2>(1,0)]= bd;
      Mb[ukey<2>(2,0)] = b;
      Mb[ukey<2>(3,0)] = mu(site)*n;
     
      Mb[ukey<2>(3,1)] = hop(site)*b;
      Mb[ukey<2>(3,2)] = hop(site)*bd;
      Mb[ukey<2>(3,3)] = id;
      Mb.Dl=4,Mb.Dr=4;
      Mb.squeeze(1e-15);
      (*this)[site]=Mb;
    }


    if (site==currentSite){
      n_sp.setSite(site),n_sm.setSite(site),id_id.setSite(site),b_id.setSite(site),bd_id.setSite(site),id_sz.setSite(site);
      Mb[ukey<2>(0,0)]=id_id;
      Mb[ukey<2>(1,0)]=bd_id;
      Mb[ukey<2>(2,0)]=b_id;

      Mb[ukey<2>(3,0)]=(mu(site)*n_id)+(coupling(position)*exp(Complex(0,-omega*x(site)))*n_sp)+(conj(coupling(position)*exp(Complex(0,-omega*x(site))))*n_sm)+Bz(position)*id_sz;
      //+Bz(position)*(id_sp+id_sm);
     
      Mb[ukey<2>(3,1)] = hop(site)*b_id;
      Mb[ukey<2>(3,2)] = hop(site)*bd_id;
      Mb[ukey<2>(3,3)] = id_id;
      Mb.Dl=4,Mb.Dr=4;
      Mb.squeeze(1e-15);
      (*this)[site]=Mb;
      if(position<couplingPositions.size()-1){
	position++;
	currentSite=couplingPositions(position);
      }
    }

  }
}



#endif
