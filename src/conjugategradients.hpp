//work in progress; this file implements a conjugate gradient procedure for MPS (to be used similar to the lanzcos engine in dmrg) but
//it is currently not working and it is not needed by the standard container.hpp-engines.
#ifndef CGRADIENTS_HPP__
#define CGRADIENTS_HPP__

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

template<class KeyType,typename T>
class TSBiCGengine{
public:
  TSBiCGengine(const BlockMatrix<1,KeyType,T>&L,const BlockMatrix<1,KeyType,T>&R,const BlockMatrix<1,KeyType,T>&mpsl,const BlockMatrix<1,KeyType,T>&mpsr,
	       const BlockMatrix<1,KeyType,T>&mpsl_,const BlockMatrix<1,KeyType,T>&mpsr_,const MPOMat<1,T> &mpol,const MPOMat<1,T> &mpor,
	       i_to_key<KeyType> itokeyl,i_to_key<KeyType> itokeyr){
    _mpol=mpol;
    _mpor=mpor;
    this->_L=L;
    this->_R=R;
    this->_Ndiag=5;
    this->_r=BlockMatrix<2,KeyType,T>(mpsl,mpsr);
    _p=_r;
    this->_r_=BlockMatrix<2,KeyType,T>(mpsl_,mpsr_);
    _p_=_r_;
    _itokey.clear();
    _itokey.push_back(itokeyl);
    _itokey.push_back(itokeyr);
  }
  int simulate(const int maxsteps,const Real delta);
protected:

private:
  BlockMatrix<1,KeyType,T>_R,_L;
  BlockMatrix<Ns,KeyType,T>_r,_r_,_p,_p_;
  MPOMat<1,T>_mpol,_mpor;
  std::vector<i_to_key<KeyType> >_itokey;
};
template<class KeyType,typename T>
int TSBiCGengine<KeyType,T>::simulate(const int maxsteps,const Real delta){
  T alpha_k,number;

  BlockMatrix<Ns,KeyType,T>temp;
  while(!converged){
    contractTwosite(temp,_L,_R,_p,_mpol,_mpor,_itokey);
    alpha_k=_r_*_r/(_p_*temp);
    _r=_r-alpha_k*temp;
    
  }

  
}







#endif
