#ifndef BICGSTAB_HPP__
#define BICGSTAB_HPP__

#include "BlockTensor.hpp"
#include "MPSfunctions.hpp"
#include "MPSTypes.hpp"
#include "MPOTypes.hpp"
#include "blasroutines.hpp"
#include "lapack_routines.hpp"
#include "utilities.hpp"



/*
Biconjugate gradient engine for MPS
author: martin ganahl (martin.ganahl@tugraz.at)


 */
using boost::numeric::ublas::matrix;
using boost::numeric::ublas::column_major;
using namespace BoostTypedefs;
template<class KeyType,typename T>
class BiCGStabEngine{
public:
  BiCGStabEngine(){};
  BiCGStabEngine(const BlockMatrix<1,KeyType,T>&L,const BlockMatrix<1,KeyType,T>&R,const BlockMatrix<1,KeyType,T>&mpsl,const BlockMatrix<1,KeyType,T>&mpsr,
		 const BlockMatrix<2,KeyType,T>&mps0,const MPOMat<1,T> &mpol,const MPOMat<1,T> &mpor,i_to_key<KeyType> itokeyl,i_to_key<KeyType> itokeyr){
    this->_x0=BlockMatrix<2,KeyType,T>(mpsl,mpsr);
    this->_b=mps0;
    _R=R;
    _L=L;
    _itokey.clear();
    _mpol=mpol;_mpor=mpor;
    _itokey.push_back(itokeyl);
    _itokey.push_back(itokeyr);

  };
  quadruple<BlockMatrix<2,KeyType,T>,uint,T,int> simulate(const int maxsteps,const Real tol,const Real TOL);
  void setX0(BlockMatrix<2,KeyType,T> const x0){_x0=x0;}
protected:
private:
  MPOMat<1,T>_mpol,_mpor;
  BlockMatrix<1,KeyType,T>_R,_L;
  BlockMatrix<2,KeyType,T>_x0,_b;
  std::vector<i_to_key<KeyType> >_itokey;


};

template<class KeyType,typename T>
quadruple<BlockMatrix<2,KeyType,T>,uint,T,int> BiCGStabEngine<KeyType,T>::simulate(const int maxsteps,const Real tol,const Real TOL) 
{
  int info=0;
  BlockMatrix<2,KeyType,T>r_i,r_im1,BM1,rh_0,v_i,v_im1,p_i,p_im1,s,t,x_i,x_im1=_x0,test;
  contractTwosite(BM1,this->_L,this->_R,_x0,_mpol,_mpor,_itokey);
  
  
  //im1=i-1=0
  r_im1=_b-BM1;
  rh_0=r_im1;
  Real norm_s_sq=abs(r_im1*r_im1),norm_b_sq=abs(_b*_b);
  if(norm_b_sq<1e-15)
    norm_b_sq=1.0;

  T rho_im1=T(1.0),alpha=T(1.0),omega_im1=T(1.0),rho_i=T(0.0),omega_i=T(0.0),beta_i=T(0.0);
  

  uint it=0;

  if((norm_s_sq/norm_b_sq)<=tol)
    return quadruple<BlockMatrix<2,KeyType,T>,uint,T,int>(x_im1,it,norm_s_sq,0);


  bool converged=false;
  while(!converged){
    rho_i=rh_0*r_im1;beta_i=(rho_i/rho_im1)*(alpha/omega_im1);

    if(abs(rho_i)<5e-14){
      if(it==0){
	return quadruple<BlockMatrix<2,KeyType,T>,uint,T,int>(x_im1,it,norm_s_sq,-2);
      }
      break;
    }

    if(it==0){
      p_i=r_im1;
    }
    else if(it!=0){
      p_i=r_im1+((p_im1-(v_im1*omega_im1))*beta_i);
    }
    contractTwosite(v_i,this->_L,this->_R,p_i,_mpol,_mpor,_itokey);
    T temp=rh_0*v_i;

    alpha=rho_i/temp;
    s=r_im1-(v_i*alpha);
    norm_s_sq=abs(s*s);

    if((norm_s_sq/norm_b_sq)<=tol){
      x_i=x_im1+((p_i*alpha));
      converged=true;
      info=0;
      break;
    }
    if(norm_s_sq>TOL){
      cout<<"bicgstab: residual > "<<TOL<<"; stopping bicgstab"<<endl;
      // T norm_x=sqrt(x_im1*x_im1);
      // x_im1*=(1.0/norm_x);
      return quadruple<BlockMatrix<2,KeyType,T>,uint,T,int>(x_im1,it,norm_s_sq,-1);
    }
    contractTwosite(t,this->_L,this->_R,s,_mpol,_mpor,_itokey);    
    temp=t*t;
    omega_i=(t*s)/temp;    
    x_i=x_im1+((p_i*alpha)+(s*omega_i));
    //check if x_i is accurate enough
    // contractTwosite(test,this->_L,this->_R,x_i,_mpol,_mpor,_itokey);
    // cout<<"true solution"<<endl;
    // cout<<_b<<endl;
    // cout<<"approximate solution"<<endl;
    // cout<<test<<endl;
    r_i=s-t*omega_i;
    T norm_r=r_i*r_i;
    //cout<<"residual norms at step "<<it<<": ("<<abs(norm_s)<<","<<abs(norm_r)<<")"<<endl;
    //now rename all variables for next iteration
    x_im1=x_i;
    r_im1=r_i;
    p_im1=p_i;
    v_im1=v_i;
    rho_im1=rho_i;
    omega_im1=omega_i;
    //    cout<<"at step "<<it<<" : rho0="<<rho_im1<<endl;
    ++it;
    if(it>maxsteps)
      converged=true;
  }
  return quadruple<BlockMatrix<2,KeyType,T>,uint,T,int>(x_im1,it,norm_s_sq,info);
}



#endif
