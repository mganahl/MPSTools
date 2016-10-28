#ifndef LANZCOS_HPP__
#define LANZCOS_HPP__

#include "BlockTensor.hpp"
#include "MPSfunctions.hpp"
#include "MPSTypes.hpp"
#include "MPOTypes.hpp"
#include "blasroutines.hpp"
#include "lapack_routines.hpp"
#include "utilities.hpp"

using boost::numeric::ublas::matrix;
using boost::numeric::ublas::column_major;
using namespace BoostTypedefs;

template<int Ns,class KeyType,typename T>
class LanzcosEngine{
public:
  LanzcosEngine(){_verbose=0;}
  // LanzcosEngine(const BlockMatrix<1,KeyType,T>&L,const BlockMatrix<1,KeyType,T>&R,const BlockMatrix<1,KeyType,T>&mps);
  // LanzcosEngine(const BlockMatrix<1,KeyType,T>&L,const BlockMatrix<1,KeyType,T>&R,const BlockMatrix<1,KeyType,T>&mpsl,const BlockMatrix<1,KeyType,T>&mpsr);


  virtual int simulate(const int maxsteps,const Real delta,const bool testenergy){return 0;}
  BlockMatrix<Ns,KeyType,T> getGroundState();
  BlockMatrix<Ns,KeyType,T> getState(const int j);
  T E(const int i)const {assert(i<_es.size());return _es[i];}
  uint Elargerthan(const T E){
    uint s=0;
    while((s<_es.size())&&(abs(_es[s])<abs(E))){
      ++s;
    }
    // cout<<s<<endl;
    // cout<<_es<<endl;
    // cout<<endl;
    return s;
  }
  void setAinit(const BlockMatrix<Ns,KeyType,T>&A){_Ainit=A;}
protected:
  BlockMatrix<1,KeyType,T>_R,_L;
  BlockMatrix<Ns,KeyType,T>_Ainit;
  BlockMatrix<Ns,KeyType,T>_GS;
  //  MPOMat<1,T>_mpol,_mpor;
  //  std::vector<i_to_key<KeyType> >_itokey;
  int _Ndiag;
  std::vector<T>_es,_ks;
  std::vector<BlockMatrix<Ns,KeyType,T> >_As;
  matrix<Real,column_major>_eigs;
  uint _verbose;
private:

};

// template<int Ns,class KeyType,typename T>
// LanzcosEngine<Ns,KeyType,T>::LanzcosEngine(){}

// template<int Ns,class KeyType,typename T>
// LanzcosEngine<1,KeyType,T>::LanzcosEngine(const BlockMatrix<1,KeyType,T>&L,const BlockMatrix<1,KeyType,T>&R,const BlockMatrix<1,KeyType,T>&mps):
//   _L(L),_R(R){
//   _Ndiag=5;
//   _Ainit=mps;
// }

// template<int Ns,class KeyType,typename T>
// LanzcosEngine<2,KeyType,T>::LanzcosEngine(const BlockMatrix<1,KeyType,T>&L,const BlockMatrix<1,KeyType,T>&R,const BlockMatrix<1,KeyType,T>&mpsl,const BlockMatrix<1,KeyType,T>&mpsr):
//   _L(L),_R(R),_Ainit(mpsl,mpsr){
//   _Ndiag=5;
// }

template<int Ns,class KeyType,typename T>
BlockMatrix<Ns,KeyType,T>LanzcosEngine<Ns,KeyType,T>::getGroundState(){
  //  BlockMatrix<Ns,KeyType,T>A=_As[0]*T(_eigs.data()[0]);
  _GS.clear();
  _GS=_As[0]*T(_eigs.data()[0]);
  for(uint i = 1;i<_es.size();i++){
    _GS +=(_As[i]*=T(_eigs.data()[i]));
  }
  return _GS;
}

template<int Ns,class KeyType,typename T>
BlockMatrix<Ns,KeyType,T>LanzcosEngine<Ns,KeyType,T>::getState(const int j){
  if(j>=_eigs.size2()){
    cout<<"LanzcosEngine::getState(const int j): j<_eigs.size2() failed"<<endl;
    abort();
  }
  BlockMatrix<Ns,KeyType,T>A=_As[0]*T(_eigs(0,j));
  for(uint i = 1;i<_es.size();i++){
    A +=(_As[i]*T(_eigs(i,j)));
  }
  return A;
}



template<class KeyType,typename T>
class TSLanzcosEngine:public LanzcosEngine<2,KeyType,T>{
public:
  TSLanzcosEngine():LanzcosEngine<2,KeyType,T>(){};
  TSLanzcosEngine(const BlockMatrix<1,KeyType,T>&L,const BlockMatrix<1,KeyType,T>&R,const BlockMatrix<1,KeyType,T>&mpsl,
		  const BlockMatrix<1,KeyType,T>&mpsr,const MPOMat<1,T> &mpol,const MPOMat<1,T> &mpor,i_to_key<KeyType> itokeyl,
		  i_to_key<KeyType> itokeyr,const uint verbose){
    _mpol=mpol;
    _mpor=mpor;
    this->_L=L;
    this->_R=R;
    this->_Ndiag=5;
    this->_Ainit=BlockMatrix<2,KeyType,T>(mpsl,mpsr);
    _itokey.clear();
    _itokey.push_back(itokeyl);
    _itokey.push_back(itokeyr);
    this->_verbose=verbose;
  }
  virtual int simulate(const int maxsteps,const Real delta,const bool testenergy,const bool bla=false);

  //does a restricted lanzcos with Gram-Schmidt orthonormalization which projects out parts parallel to state
  virtual int simulate(const int maxsteps,const Real delta,const bool testenergy,BlockMatrix<2,KeyType,T>&state,const T delta_ortho);
  //for Operator, a correction vector is returned
  map_container<Real,BlockMatrix<2,KeyType,T> >
  getCorrectionVector(const BlockMatrix<2,KeyType,T>&GS,const BlockMatrix<1,KeyType,T> &Lo,const MPOMat<1,T> &opl,const MPOMat<1,T> &opr,
		      const BlockMatrix<1,KeyType,T>&Ro,const std::vector<Real> om,const std::vector<Real> eta,const int maxsteps,
		      const Real delta,const bool testenergy);

protected:

private:
  MPOMat<1,T>_mpol,_mpor;
  std::vector<i_to_key<KeyType> >_itokey;
};


template<class KeyType,typename T>
map_container<Real,BlockMatrix<2,KeyType,T> >TSLanzcosEngine<KeyType,T>::
getCorrectionVector(const BlockMatrix<2,KeyType,T>&GS,const BlockMatrix<1,KeyType,T> &Lo,const MPOMat<1,T> &opl,const MPOMat<1,T> &opr,
		    const BlockMatrix<1,KeyType,T>&Ro,const std::vector<Real> om,const std::vector<Real> eta,const int maxsteps,const Real delta,
		    const bool testenergy){
  //first get the initial state ...
  map_container<Real,BlockMatrix<2,KeyType,T> > out;
  //  if(this->_GS.size()==0)
  //  BlockMatrix<2,KeyType,T>bla=  this->getGroundState();//check if the Groundstate has really been computed ...
  this->_Ainit.clear();
  contractTwosite(this->_Ainit,Lo,Ro,GS,opl,opr,_itokey);
  int N=TSLanzcosEngine<KeyType,T>::simulate(maxsteps,delta,testenergy);
  //now expand the state in the krylov basis ...e
  if(N!=0){
    for(uint ind=0;ind<om.size();++ind){
      BlockMatrix<2,KeyType,T> A,state;
      state=this->getState(0);
      T c=state*this->_Ainit;
      c*=1.0/(om[ind]+Complex(0.0,eta[ind]));
      A=state*c;
      for(int i=1;i<N;++i){
	state=this->getState(i);
	T c=state*this->_Ainit;
	c*=1.0/(om[ind]+this->E(0)-this->E(i)+Complex(0.0,eta[ind]));
	A+=state*c;
      }
      out[om[ind]]=A;
    }
  }
  return out;
  //ole ole ...
}

template<class KeyType,typename T>
int TSLanzcosEngine<KeyType,T>::simulate(const int itmax,const Real delta,const bool testenergy,const bool debug){
  // if(debug)
  //   cout<<"in lanczos inside the class"<<endl;
  //bring everything into right form:
  int it = 0;
  this->_As.resize(itmax);
  T k = 0,e,eold = T(1e100);
  double *eig = NULL;
  bool converged = false, first = true;

  vector<T>estmp,kstmp;
  BlockMatrix<2,KeyType,T>A0,A1=this->_Ainit,A2,A3;
  this->_es.clear();
  estmp.clear();
  this->_ks.clear();
  kstmp.clear();
  Real resi=1e5;
  A0.clear();A2.clear();A3.clear();
  while(!converged){
    // if(debug)
    //cout<<"doing a lanczos step"<<endl;
    k=(sqrt(A1*A1));
    if(abs(k)<delta){converged = true;break;}
    if(!first){this->_ks.push_back(k);}
    A1/=k;
    // if(debug){
    //   cout<<"these are the matrices wich are contracted"<<endl;
    //   cout<<this->_L<<endl;
    //   cout<<A1<<endl;
    //   cout<<this->_R<<endl;
    // }
    contractTwosite(A2,this->_L,this->_R,A1,_mpol,_mpor,_itokey,8,debug);

    // if(debug){
    //   cout<<"that's the result of the contraction"<<endl;
    //   cout<<A2<<endl;
    //   cin.get();
    // }
    if(it<itmax)
      this->_As[it]=A1;
    else if(it>=itmax)
      this->_As.push_back(A1);
    e = A2*A1;
    this->_es.push_back(e);
    if(it%this->_Ndiag==0&&it!=0){
      if(testenergy){
	kstmp = this->_ks;
	estmp = this->_es;
	tridiag(&estmp[0],estmp.size(),&kstmp[0],eig,'N');
	cout.precision(18);
	if(abs(estmp[0]-eold)<delta){/*cout<<"E="<<estmp[0]<<" after "<<it<<" iterations; ";*/converged = true;break;}
	eold = estmp[0];
      }else{
	kstmp = this->_ks;
	estmp = this->_es;
	BlockMatrix<2,KeyType,T>A4,A5,A6;
	double*eig2 = new double[estmp.size()*estmp.size()];
	tridiag(&estmp[0],estmp.size(),&kstmp[0],eig2,'V');
	A4 = this->_As[0]*T(eig2[0]);
	for(uint i = 1;i<estmp.size();i++){
	  A5 = (A4+(this->_As[i]*T(eig2[i])));
	  A4 = A5;
	}
	delete[] eig2;
	//now test for convergence;
	contractTwosite(A5,this->_L,this->_R,A4,_mpol,_mpor,_itokey);
	A6=A5+A4*(-estmp[0]);
	resi = abs(sqrt(A6*A6));
	if(abs(resi)<delta){/*cout<<"state converged at "<<estmp[0]<<" after "<<it<<" iterations at a residuum of "<<resi<<"; ";*/converged = true;break;}
      }
    }
    if(!first){
      A3 = (((A1*(-e))+=A2)+=(A0*=(-k)));
      A0 = A1;
      A1=A3;
    }
    if(first){
      A0 = A1;
      A3= ((A1*(-e))+A2);
      A1=A3;
      first = false;
    }
    it++;
    if(it>=itmax){
      if(testenergy){
	cout<<"LanzcosEngine:: no convergence after "<<it<<" iteration steps; exit Lanzcos loop at energy "<<estmp[0]<<endl;
	converged = true;break;}
      else{cout<<"LanzcosEngine:: no convergence after "<<it<<" iteration steps; exit Lanzcos loop at energy "<<estmp[0]
	       <<" and a residuum of "<<resi<<endl;converged = true;break;}
    }

  }
  if(converged&&this->_es.size()!=0){
    const int size = this->_es.size();
    double*eig2 = new double[size*size];

    tridiag(&this->_es[0],this->_es.size(),&this->_ks[0],eig2,'V');
    this->_eigs.resize(size,size);
    memcpy(&this->_eigs.data()[0],eig2,size*size*sizeof(double));
    if((it<itmax)&&(this->_verbose==1))
      cout<<"E="<<this->_es[0]<<" after "<<it<<" iterations; ";
    delete[] eig2;
  }
  if(abs(k)<delta&&this->_es.size()!=0){if(this->_verbose==1)
      cout<<"k < "<<delta<<": E="<<this->_es[0]<<"; ";
  }
  if(abs(k)<delta&&this->_es.size()==0){
    if(this->_verbose==1)
      cout<<"k < "<<delta<<"; no energies found ";
  }
  
  return it;
}


template<class KeyType,typename T>
int TSLanzcosEngine<KeyType,T>::simulate(const int itmax,const Real delta,const bool testenergy,BlockMatrix<2,KeyType,T>&state,const T delta_ortho){
  //bring everything into right form:
  int it = 0;
  this->_As.resize(itmax);
  T k = 0,e,eold = T(1e100);
  double *eig = NULL;
  bool converged = false, first = true;

  vector<T>estmp,kstmp;
  BlockMatrix<2,KeyType,T>A0,A1=this->_Ainit,A2,A3,A3_,A1_;
  this->_es.clear();
  estmp.clear();
  this->_ks.clear();
  kstmp.clear();
  Real resi=1e5;
  state.squeeze(1e-15);
  T state_norm=state*state;
  A0.clear();A2.clear();A3.clear();A3_.clear();
  A1_=A1-state*((state*A1)/state_norm);
  A1=A1_;
  //  A1_.clear();
      
  while(!converged){
    k=(sqrt(A1*A1));
    if(abs(k)<delta){converged = true;break;}
    if(!first){this->_ks.push_back(k);}
    A1/=k;
    contractTwosite(A2,this->_L,this->_R,A1,_mpol,_mpor,_itokey);
    if(it<itmax)
      this->_As[it]=A1;
    else if(it>=itmax)
      this->_As.push_back(A1);
    e = A2*A1;
    this->_es.push_back(e);
    if(it%this->_Ndiag==0&&it!=0){
      if(testenergy){
	kstmp = this->_ks;
	estmp = this->_es;
	tridiag(&estmp[0],estmp.size(),&kstmp[0],eig,'N');
	cout.precision(18);
	if(abs(estmp[0]-eold)<delta){/*cout<<"E="<<estmp[0]<<" after "<<it<<" iterations; ";*/converged = true;break;}
	eold = estmp[0];
      }else{
	kstmp = this->_ks;
	estmp = this->_es;
	BlockMatrix<2,KeyType,T>A4,A5,A6;
	double*eig2 = new double[estmp.size()*estmp.size()];
	tridiag(&estmp[0],estmp.size(),&kstmp[0],eig2,'V');
	A4 = this->_As[0]*T(eig2[0]);
	for(uint i = 1;i<estmp.size();i++){
	  A5 = (A4+(this->_As[i]*T(eig2[i])));
	  A4 = A5;
	}
	delete[] eig2;
	//now test for convergence;
	contractTwosite(A5,this->_L,this->_R,A4,_mpol,_mpor,_itokey);
	A6=A5+A4*(-estmp[0]);
	resi = abs(sqrt(A6*A6));
	if(abs(resi)<delta){/*cout<<"state converged at "<<estmp[0]<<" after "<<it<<" iterations at a residuum of "<<resi<<"; ";*/converged = true;break;}
      }
    }

    if(!first){
      A3_ = (((A1*(-e))+=A2)+=(A0*=(-k)));
      A3=A3_-state*((state*A3_)/state_norm);
  //      //A3 = (((A1*(-e))+=A2));
//
//      if(A1_*A3>delta_ortho){
//	A3_=A3;
//	cout<<"warning in TSlanczos: loss of orthogonality detected; starting Gram-Schmidt orthonormalization: "<<A1_*A3<<endl;
//	cout<<"at step "<<it<<endl;
//	cin.get();
//	typename std::vector<BlockMatrix<2,KeyType,T> >::iterator iA;
//	uint i=0;
//	for(iA=this->_As.begin();iA!=this->_As.end();++iA){
//	  if(i>=it)
//	    break;
//	  ++i;
//	  //BlockMatrix<2,KeyType,T>temp=(*it);
//	  A3_-=((*iA)*(((*iA)*A3)/((*iA)*(*iA))));
//	  //A3-=(temp*(temp*A3)/(temp*temp));
//	}
//	A3=A3_;
//	cout<<A1_*A3<<endl;
//      }
//


      //      cout<<A1_*A3<<endl;
      //      cin.get();
      //}
      A0 = A1;
      A1=A3;
    }
    if(first){
      A0 = A1;
      A3_= ((A1*(-e))+A2);
      A3=A3_-state*((state*A3_)/state_norm);
      //      A3_=A3-state*((state*A3)/state_norm);
      //A3=A3_;
      //A3_.clear();
      A1=A3;
      first = false;
    }
    it++;
    if(it>=itmax){
      if(testenergy){
	cout<<"LanzcosEngine:: no convergence after "<<it<<" iteration steps; exit Lanzcos loop at energy "<<estmp[0]<<endl;
	converged = true;break;}
      else{cout<<"LanzcosEngine:: no convergence after "<<it<<" iteration steps; exit Lanzcos loop at energy "<<estmp[0]
	       <<" and a residuum of "<<resi<<endl;converged = true;break;}
    }
  }
  if(converged&&this->_es.size()!=0){
    const int size = this->_es.size();
    double*eig2 = new double[size*size];

    tridiag(&this->_es[0],this->_es.size(),&this->_ks[0],eig2,'V');
    this->_eigs.resize(size,size);
    memcpy(&this->_eigs.data()[0],eig2,size*size*sizeof(double));
    if((it<itmax)&&(this->_verbose==1))
      cout<<"E="<<this->_es[0]<<" after "<<it<<" iterations; ";
    delete[] eig2;
  }
  if(abs(k)<delta&&this->_es.size()!=0){if(this->_verbose==1)
      cout<<"k < "<<delta<<": E="<<this->_es[0]<<"; ";
  }
  if(abs(k)<delta&&this->_es.size()==0){
    if(this->_verbose==1)
      cout<<"k < "<<delta<<"; no energies found ";
  }
  return it;
}




template<class KeyType,typename T>
class SSLanzcosEngine:public LanzcosEngine<1,KeyType,T>{
public:
  SSLanzcosEngine():LanzcosEngine<1,KeyType,T>(){}
  SSLanzcosEngine(const BlockMatrix<1,KeyType,T>&L,const BlockMatrix<1,KeyType,T>&R,const BlockMatrix<1,KeyType,T>&mps,
		  const MPOMat<1,T> &mpo,i_to_key<KeyType> itokey,const uint verbose){
    this->_L=L;
    this->_R=R;
    this->_Ndiag=5;
    this->_Ainit=mps;
    _mpo=mpo;
    _itokey.clear();
    _itokey=itokey;
    this->_verbose=verbose;
  }

  virtual int simulate(const int maxsteps,const Real delta,const bool testenergy);
  virtual int simulate(const int itmax,const Real delta,const bool testenergy,BlockMatrix<1,KeyType,T>&state,const T delta_ortho);
protected:

private:
  MPOMat<1,T>_mpo;
  i_to_key<KeyType>_itokey;
};



// template<class KeyType,typename T>
// class SSLanzcosEngine{
// public:
//   SSLanzcosEngine(const BlockMatrix<1,KeyType,T>&L,const BlockMatrix<1,KeyType,T>&R,const BlockMatrix<1,KeyType,T>&mps,
// 		const MPOMat<1,T> &mpo,i_to_key<KeyType> itokey):_L(L),_R(R),_Ainit(mps){
//     _mpo=mpo;
//     _itokey.clear();
//     _itokey=itokey;
//     _Ndiag=5;
//   }
//   int simulate(const int maxsteps,const Real delta,const bool testenergy);
//   BlockMatrix<1,KeyType,T> getGroundState();
//   BlockMatrix<1,KeyType,T> getState(const int j);
//   T E(const int i)const {assert(i<_es.size());return _es[i];}
//   uint Elargerthan(const T E){
//     uint s=0;
//     while((s<_es.size())&&(_es[s]<E)){
//       ++s;
//     }
//     return s;
//   }
// protected:

// private:
//   BlockMatrix<1,KeyType,T>_R,_L;
//   BlockMatrix<1,KeyType,T>_Ainit;
//   MPOMat<1,T>_mpo;
//   i_to_key<KeyType>_itokey;
//   int _Ndiag;
//   std::vector<T>_es,_ks;
//   std::vector<BlockMatrix<1,KeyType,T> >_As;
//   matrix<Real,column_major>_eigs;
// };

// template<class KeyType,typename T>
// BlockMatrix<1,KeyType,T>SSLanzcosEngine<KeyType,T>::getGroundState(){
//   BlockMatrix<1,KeyType,T>A=_As[0]*T(_eigs.data()[0]);
//   for(uint i = 1;i<_es.size();i++){
//     A +=(_As[i]*T(_eigs.data()[i]));
//   }
//   return A;
// }

// template<class KeyType,typename T>
// BlockMatrix<1,KeyType,T>SSLanzcosEngine<KeyType,T>::getState(const int j){
//   if(j>=_eigs.size2()){
//     cout<<"SSLanzcosEngine::getState(const int j): j<_eigs.size2() failed"<<endl;
//     abort();
//   }
//   BlockMatrix<1,KeyType,T>A=_As[0]*T(_eigs(0,j));
//   for(uint i = 1;i<_es.size();i++){
//     A +=(_As[i]*T(_eigs(i,j)));
//   }
//   return A;
// }

template<class KeyType,typename T>
int SSLanzcosEngine<KeyType,T>::simulate(const int itmax,const Real delta,const bool testenergy){
  //bring everything into right form:
  int it = 0;
  this->_As.resize(itmax);
  T k = 0,e,eold = T(1e100);
  double *eig = NULL;
  bool converged = false, first = true;

  vector<T>estmp,kstmp;
  BlockMatrix<1,KeyType,T>A0,A1=this->_Ainit,A2,A3;
  this->_es.clear();
  estmp.clear();
  this->_ks.clear();
  kstmp.clear();
  Real resi=1e5;
  A0.clear();A2.clear();A3.clear();
  while(!converged){
    k=(sqrt(A1*A1));
    if(abs(k)<delta){converged = true;break;}
    if(!first){
      this->_ks.push_back(k);
    }
    A1/=k;
    contractSingle(A2,this->_L,A1,_mpo,this->_R,_itokey);
    if(it<itmax)
      this->_As[it]=A1;
    else if(it>=itmax)
      this->_As.push_back(A1);
    e = A2*A1;
    this->_es.push_back(e);
    if(it%this->_Ndiag==0&&it!=0){
      if(testenergy){
	kstmp = this->_ks;
	estmp = this->_es;
	tridiag(&estmp[0],estmp.size(),&kstmp[0],eig,'N');
	cout.precision(18);
	if(abs(estmp[0]-eold)<delta){/*cout<<"E="<<estmp[0]<<" after "<<it<<" iterations; ";*/converged = true;break;}
	eold = estmp[0];
      }else{
	kstmp = this->_ks;
	estmp = this->_es;
	BlockMatrix<1,KeyType,T>A4,A5,A6;
	double*eig2 = new double[estmp.size()*estmp.size()];
	tridiag(&estmp[0],estmp.size(),&kstmp[0],eig2,'V');
	A4 = this->_As[0]*T(eig2[0]);
	for(uint i = 1;i<estmp.size();i++){
	  A5 = (A4+(this->_As[i]*T(eig2[i])));
	  A4 = A5;
	}
	delete[] eig2;
	//now test for convergence;
	contractSingle(A5,this->_L,A4,_mpo,this->_R,_itokey);
	A6=A5+A4*(-estmp[0]);
	resi = abs(sqrt(A6*A6));
	if(abs(resi)<delta){
	  //cout<<"state converged at "<<estmp[0]<<" after "<<it<<" iterations at a residuum of "<<resi<<"; ";
	  converged = true;break;
	}
      }
    }
    if(!first){
      A3 = (((A1*(-e))+=A2)+=(A0*=(-k)));
      A0 = A1;
      A1=A3;
    }
    if(first){
      A0 = A1;
      A3= ((A1*(-e))+A2);
      A1=A3;
      first = false;
    }
    it++;
    if(it>=itmax){
      if(testenergy){
	//cout<<"SSLanzcosEngine:: no convergence after "<<it<<" iteration steps; exit Lanzcos loop at energy "<<estmp[0]<<"; ";
	converged = true;break;
      }
      else{
	//	cout<<"SSLanzcosEngine:: no convergence after "<<it<<" iteration steps; exit Lanzcos loop at energy "<<estmp[0]<<" and a residuum of "<<resi<<endl;
	converged = true;break;}
    }
  }
  if(converged&&this->_es.size()!=0){
    const int size = this->_es.size();
    double*eig2 = new double[size*size];
    
    tridiag(&this->_es[0],this->_es.size(),&this->_ks[0],eig2,'V');
    this->_eigs.resize(size,size);
    memcpy(&this->_eigs.data()[0],eig2,size*size*sizeof(double));
    // if(it<itmax)
    //   cout<<"E="<<this->_es[0]<<" after "<<it<<" iterations; "<<endl;
    delete[] eig2;
  }
  //  if(abs(k)<delta&&this->_es.size()!=0){cout<<"k < "<<delta<<": E="<<this->_es[0]<<"; "<<endl;}
  return it;
}




template<class KeyType,typename T>
int SSLanzcosEngine<KeyType,T>::simulate(const int itmax,const Real delta,const bool testenergy,BlockMatrix<1,KeyType,T>&state,const T delta_ortho){
  //bring everything into right form:
  int it = 0;
  this->_As.resize(itmax);
  T k = 0,e,eold = T(1e100);
  double *eig = NULL;
  bool converged = false, first = true;

  vector<T>estmp,kstmp;
  BlockMatrix<1,KeyType,T>A0,A1=this->_Ainit,A2,A3,A3_,A1_;
  this->_es.clear();
  estmp.clear();
  this->_ks.clear();
  kstmp.clear();
  Real resi=1e5;
  A0.clear();A2.clear();A3.clear();
  T state_norm=state*state;
  A1_=A1-state*((state*A1)/state_norm);
  A1=A1_;
  while(!converged){
    k=(sqrt(A1*A1));
    if(abs(k)<delta){converged = true;break;}
    if(!first){
      this->_ks.push_back(k);
    }
    A1/=k;
    contractSingle(A2,this->_L,A1,_mpo,this->_R,_itokey);
    if(it<itmax)
      this->_As[it]=A1;
    else if(it>=itmax)
      this->_As.push_back(A1);
    e = A2*A1;
    this->_es.push_back(e);
    if(it%this->_Ndiag==0&&it!=0){
      if(testenergy){
	kstmp = this->_ks;
	estmp = this->_es;
	tridiag(&estmp[0],estmp.size(),&kstmp[0],eig,'N');
	cout.precision(18);
	if(abs(estmp[0]-eold)<delta){/*cout<<"E="<<estmp[0]<<" after "<<it<<" iterations; ";*/converged = true;break;}
	eold = estmp[0];
      }else{
	kstmp = this->_ks;
	estmp = this->_es;
	BlockMatrix<1,KeyType,T>A4,A5,A6;
	double*eig2 = new double[estmp.size()*estmp.size()];
	tridiag(&estmp[0],estmp.size(),&kstmp[0],eig2,'V');
	A4 = this->_As[0]*T(eig2[0]);
	for(uint i = 1;i<estmp.size();i++){
	  A5 = (A4+(this->_As[i]*T(eig2[i])));
	  A4 = A5;
	}
	delete[] eig2;
	//now test for convergence;
	contractSingle(A5,this->_L,A4,_mpo,this->_R,_itokey);
	A6=A5+A4*(-estmp[0]);
	resi = abs(sqrt(A6*A6));
	if(abs(resi)<delta){
	  //cout<<"state converged at "<<estmp[0]<<" after "<<it<<" iterations at a residuum of "<<resi<<"; ";
	  converged = true;break;
	}
      }
    }
    if(!first){
      A3_ = (((A1*(-e))+=A2)+=(A0*=(-k)));
      A3=A3_-state*((state*A3_)/state_norm);
      A0 = A1;
      A1=A3;
    }
    if(first){
      A0 = A1;
      A3_= ((A1*(-e))+A2);
      A3=A3_-state*((state*A3_)/state_norm);
      A1=A3;
      first = false;
    }
    it++;
    if(it>=itmax){
      if(testenergy){
	//cout<<"SSLanzcosEngine:: no convergence after "<<it<<" iteration steps; exit Lanzcos loop at energy "<<estmp[0]<<"; ";
	converged = true;break;
      }
      else{
	//	cout<<"SSLanzcosEngine:: no convergence after "<<it<<" iteration steps; exit Lanzcos loop at energy "<<estmp[0]<<" and a residuum of "<<resi<<endl;
	converged = true;break;}
    }
  }
  if(converged&&this->_es.size()!=0){
    const int size = this->_es.size();
    double*eig2 = new double[size*size];
    
    tridiag(&this->_es[0],this->_es.size(),&this->_ks[0],eig2,'V');
    this->_eigs.resize(size,size);
    memcpy(&this->_eigs.data()[0],eig2,size*size*sizeof(double));
    // if(it<itmax)
    //   cout<<"E="<<this->_es[0]<<" after "<<it<<" iterations; "<<endl;
    delete[] eig2;
  }
  //  if(abs(k)<delta&&this->_es.size()!=0){cout<<"k < "<<delta<<": E="<<this->_es[0]<<"; "<<endl;}
  return it;



}


#endif
