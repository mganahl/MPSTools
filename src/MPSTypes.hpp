#ifndef MPSTYPES_HPP__
#define MPSTYPES_HPP__
#include "../boost_typedefs.hpp"
#include "BlockTensor.hpp"
#include <iostream>
#include <fstream>
#include "typedefs.hpp"
#include "serialization.hpp"
#include <omp.h>
#include <fstream>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

using boost::numeric::ublas::column_major;
using boost::numeric::ublas::identity_matrix;
using boost::numeric::ublas::vector_range;

//this exception is not used;
class mps_exception: public std::exception
{
public:
  mps_exception(): std::exception() {}
  mps_exception(const char* message): std::exception() { std::cerr << message << std::endl; abort();}
};

//Base class of all MPS; it basically contains all feature of MPS
template<class KeyType,typename T>
class MPSBase : public std::vector<BlockMatrix<1,KeyType,T> >{
public:
  //constructors ...
  MPSBase():_lOrtho(false),_rOrtho(false){};
  MPSBase(const uint N):_lOrtho(false),_rOrtho(false){resize(N);}
  MPSBase(const UintVectorType dimP,std::vector<i_to_key<KeyType> > itokey){
    resize(dimP.size());_dimP=dimP;_itokey=itokey;_lOrtho=false;_rOrtho=false;
  }
  MPSBase(const UintVectorType dim_site_space,const UintVectorType local_state,std::vector<i_to_key<KeyType> > itokey);
  MPSBase(const MPSBase<KeyType,T>&mps):std::vector<BlockMatrix<1,KeyType,T> >(mps),_lambdas(mps.getLambdas()){
	  initialize(mps);
  };
  virtual ~MPSBase(){};

  //squeezes out all values in the MPS matrices which are smaller than delta
  uint squeeze(const Real delta=1e-15){
    uint e=0;
    for(uint s=0;s<this->size();++s)
      e+=(*this)[s].squeeze(delta);
    return e;
  }
  //saves the MPS in a file with the name filename.mps
  void save(string filename);
  
  //loads an MPS from a file with the name filename.mps. If cast=true, then the MPS is cast from real to complex; the file then has to be of typy Real, and the MPS of type Complex
  void load(string filename,bool cast = false);
  
  //changes the length of the MPS
  void resize(const uint N){
    std::vector<BlockMatrix<1,KeyType,T> >::resize(N);
    _itokey.resize(N);
    _dimP.resize(N,true);
    _dimAux.resize(N,true);
    _numSites=N;
  }

  //for measurements: calculates the entanglement entropy at site 
  Real getEntropy(int site) const;

  //clearing the lambdas of MPS:
  void clearLambdas(){
    //cout<<"clearing lambdas ..."<<endl;
    const uint s=_lambdas.size();
    for(uint i=0;i<s;++i){
      _lambdas[i].clear();
    }
    _lambdas.clear();
  }
  //clears the MPS matrices
  void clearMatrices(){
    //    cout<<"clearing matrices ..."<<endl;
    const int s = this->size();
    for(uint i=0;i<s;++i){
      (*this)[i].clear();
    }
    std::vector<BlockMatrix<1,KeyType,T> >::clear();
    this->resize(s);
    _numSites=s;
  }
  //clears both lambdas and the mps matrices
  void clear(){clearLambdas();clearMatrices();};

  //get information:
  map_container<int,DiagBlockMatrix<KeyType,Real > > getLambdas()const{return _lambdas;} //return the full set of lambdas for all sites
  void setLambdas(const map_container<int,DiagBlockMatrix<KeyType,Real > >&lam){_lambdas=lam;}//set lambdas of the mps to lam
  DiagBlockMatrix<KeyType,Real>&lambda(const int site){return _lambdas[site];}//return lambdas at a certain site; should really be int, should NOT be const, is used by prepareSite
  DiagBlockMatrix<KeyType,Real>lambda(const int site)const{return _lambdas.find(site)->second;}//return lambdas at a certain site, key hs to be of type int
  UintVectorType getPSpace()const {return _dimP;}//returns a boost vector containing the physical dimensions of the MPS
  uint getPSpace(const uint s)const {return _dimP(s);} //return physical dimension of MPS at site s;
  UintVectorType getAuxSpace();//returns a vector containing the auxilliary dimensions of the mps
  UintVectorType getAuxSpace()const;//returns a vector containing the auxilliary dimensions of the mps
  uint chi(const int site)const;//returns matrix-dimension at site
  void normalize(const bool store=false){orthonormalize(*this,"lr",store);}
  DiagBlockMatrix<KeyType,Real> norm()const{return _norm;};
  DiagBlockMatrix<KeyType,Real> &norm(){return _norm;};
  void setPSpace(const UintVectorType &dimP){_dimP=dimP;}
  void setAuxSpace(const UintVectorType &dimAux){ _dimAux=dimAux;}
  std::set<KeyType> getQN()const;
  void setLOrtho(const bool &o){_lOrtho=o;};
  void setROrtho(const bool &o){ _rOrtho=o;};
  bool getLOrtho()const {return _lOrtho;};
  bool getROrtho()const {return _rOrtho;};
  bool check()const;
  bool dropDisco(){
    bool stop=true;
    for(int site = 0;site<(this->size()-1);++site){
      bool dd=dropDangling((*this)[site],(*this)[site+1]);
      if(!dd)
	stop = false;
    }
    while(!stop)
      stop=this->dropDisco();
    return stop;
  }
  std::vector<i_to_key<KeyType> > getItoKey(){return _itokey;}
  i_to_key<KeyType> getItoKey(const int site)const {return _itokey[site];};
  void setItoKey(const std::vector<i_to_key<KeyType> > itokey){ _itokey=itokey;};
  void setItoKey(const int site,i_to_key<KeyType>itokey){ _itokey[site]=itokey;};
  void generateBoundary();
  BlockMatrix<1,KeyType,T>_leftDummy,_rightDummy;  //these are boundary BlockMatrices
  void pop_back(){
    _numSites=this->size()-1;
    std::vector<BlockMatrix<1,KeyType,T> >::pop_back();
    _dimAux.resize(_numSites,true);
    _dimP.resize(_numSites,true);
    _itokey.pop_back();
  }
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version){
    ar & boost::serialization::base_object<std::vector<BlockMatrix<1,KeyType,T> > >(*this);
    ar & _dimP;
    ar & _dimAux;
    ar & _lOrtho;
    ar & _rOrtho;
    ar & _conservedQN;
    ar & _itokey;
    ar & _lambdas;
    ar & _leftDummy;
    ar & _rightDummy;
  }

  //operations ...
  MPSBase<KeyType,T>& operator +=(const MPSBase<KeyType,T>&mps);
  MPSBase<KeyType,T>& operator -=(MPSBase<KeyType,T>mps); //do not pass ass reference, it is changed ...
  MPSBase<KeyType,T> operator -(const MPSBase<KeyType,T>&mps)const;
  MPSBase<KeyType,T>& operator *=(const T num);
  MPSBase<KeyType,T> operator *(const T num)const{MPSBase<KeyType,T>t=(*this);return t*=num;}
  MPSBase<KeyType,T> operator +(const MPSBase<KeyType,T>&mps)const;
  void conjugate();
protected:
  void initialize(const MPSBase<KeyType,T>&mps){
    _dimP = mps.getPSpace();
    _dimAux = mps.getAuxSpace();
    _lOrtho=mps.getLOrtho();
    _rOrtho=mps.getROrtho();
    _itokey.resize(mps.size());
    _numSites=mps.size();
    for(int site = 0;site<mps.size();++site)
      _itokey[site]=mps.getItoKey(site);
  }
  void fixSizes();
protected:
  //_lambdas contains Schmidt-eigenvalues; DiagBlockMatrix is a class derived from a boost vector (defined in BlockTensor.hpp) and stores the lambda of 
  //each particle number block
  map_container<int,DiagBlockMatrix<KeyType,Real > > _lambdas;//key should be int, cause there's a beyond-the-leftmost-site lambda
  DiagBlockMatrix<KeyType,Real> _norm;//should be int, cause there's a beyond-the-leftmost-site lambda
  UintVectorType _dimP;
  UintVectorType _dimAux;
  uint _numSites;
  std::vector<KeyType>  _conservedQN; //this could be removed but then older mps files can no longer be loaded
  bool _rOrtho;
  bool _lOrtho;
  std::vector<i_to_key<KeyType> > _itokey;
};

template<class KeyType,typename T>
void MPSBase<KeyType,T>::conjugate(){
  for(uint site=0;site<this->size();++site){
    (*this)[site].conjugate();
  }
}

template<class KeyType,typename T>
void MPSBase<KeyType,T>::fixSizes(){
  typename BlockMatrix<1,KeyType,T>::iterator b;
  typename std::map<KeyType,uint>::iterator si;
  for(uint s=0;s<(this->size()-1);++s){
    std::map<KeyType,uint> sizes=getlargestonbond((*this)[s],(*this)[s+1]);
    typename std::map<KeyType,uint>::iterator si;
    
    for(b=(*this)[s].begin();b!=(*this)[s].end();++b){
      si=sizes.find(b->first[1]);
      assert(si!=sizes.end());
      uint newsize=si->second;
      if(newsize==b->second.size2())
	continue;
      matrix<T,column_major> tmp(b->second.size1(),newsize);
      assert(newsize>=b->second.size2());
      boost_setTo(tmp,T(0.0));
      boost::numeric::ublas::matrix_range< matrix< T, column_major > > (tmp,Range(0,b->second.size1()),Range(0,b->second.size2()))=b->second;
      b->second=tmp;
    }
    for(b=(*this)[s+1].begin();b!=(*this)[s+1].end();++b){
      si=sizes.find(b->first[0]);
      assert(si!=sizes.end());
      uint newsize=si->second;
      if(newsize==b->second.size1())
	continue;
      matrix<T,column_major> tmp(newsize,b->second.size2());
      assert(newsize>=b->second.size1());
      boost_setTo(tmp,T(0.0));
      matrix_range< matrix < T, column_major > >(tmp,Range(0,b->second.size1()),Range(0,b->second.size2()))=b->second;
      b->second=tmp;
    }
  }
}

template<class KeyType,typename T>
MPSBase<KeyType,T>& MPSBase<KeyType,T>::operator *=(const T num){
  //OK, where do we put the number?
  if(this->size()==0){
    cout<<"in MPSBase<KeyType,T>& MPSBase<KeyType,T>::operator *=(const T num): trying to access (*this)[0], but this->size()=0;"<<endl;
    abort();
  }
  (*this)[0]*=num;
  return *this;
}

template<class KeyType,typename T>
MPSBase<KeyType,T>& MPSBase<KeyType,T>::operator-=(MPSBase<KeyType,T>mps){
  mps*=T(-1.0);
  (*this)+=mps;
  return (*this);
}

template<class KeyType,typename T>
MPSBase<KeyType,T>& MPSBase<KeyType,T>::operator +=(const MPSBase<KeyType,T>&mps){
  if(this==&mps){
    cout<<"in MPSBase<KeyType,T>::operator +=: self-assignment"<<endl;
    abort();
  }
  assert(mps.size()==this->size());
  //first and last sites are treated different from the bulk ...
  BlockMatrix<1,KeyType,T> A=(*this)[0];
  typename BlockMatrix<1,KeyType,T>::iterator a;
  typename BlockMatrix<1,KeyType,T>::const_iterator b;
  for(b=mps[0].begin();b!=mps[0].end();++b){
    a=(*this)[0].find(b->first);
    if(a==(*this)[0].end())
      A[b->first]=b->second;
    else{
      matrix<T,column_major>m,m1=a->second,m2=b->second;
      const int s11=m1.size1();
      const int s12=m1.size2();
      const int s21=m2.size1();
      const int s22=m2.size2();
      m.resize(1,s12+s22);
      boost_setTo(m,T(0.0));
      matrix_range<matrix<T, column_major> > (m,Range(0,1),Range(0,s12))=m1;
      matrix_range<matrix<T, column_major> > (m,Range(0,1),Range(s12,s12+s22))=m2;
      A[b->first]=m;
    }
  }
  (*this)[0]=A;
  A=(*this)[this->size()-1];
  for(b=mps[this->size()-1].begin();b!=mps[this->size()-1].end();++b){
    a=(*this)[this->size()-1].find(b->first);
    if(a==(*this)[this->size()-1].end())
      A[b->first]=b->second;
    else{
      matrix<T,column_major>m,m1=a->second,m2=b->second;
      const int s11=m1.size1();
      const int s12=m1.size2();
      const int s21=m2.size1();
      const int s22=m2.size2();
      m.resize(s11+s21,1);
      boost_setTo(m,T(0.0));
      matrix_range<matrix<T, column_major> > (m,Range(0,s11),Range(0,1))=m1;
      matrix_range<matrix<T, column_major> > (m,Range(s11,s11+s21),Range(0,1))=m2;
      A[b->first]=m;
    }
  }
  (*this)[this->size()-1]=A;
  for(int s=1;s<mps.size()-1;++s){
    (*this)[s]=outersum((*this)[s],mps[s]);
  }
  //OK; in principle we're done, but now there are two problems: a) we will have a size mismatch
  //of the matrices of BlockMatrices; and b) there may be dangling matrices. this has to be fixed ...
  //this->dropDisco();
  this->fixSizes();
  //now the mps surely is not orthonormal anymore ...
  this->_rOrtho=false;
  this->_lOrtho=false;
  boost_setTo(this->_dimAux,uint(0));
  this->clearLambdas();
  return *this;
}


template<class KeyType,typename T>
MPSBase<KeyType,T> MPSBase<KeyType,T>::operator +(const MPSBase<KeyType,T>&mps)const{
  MPSBase<KeyType,T> cpy=(*this);
  cpy+=mps;
  return cpy;
}

template<class KeyType,typename T>
MPSBase<KeyType,T> MPSBase<KeyType,T>::operator -(const MPSBase<KeyType,T>&mps)const{
  MPSBase<KeyType,T> cpy=(*this);
  return cpy-=mps;
}
//generates a boundary blockmatrix which absorbes the norm of the mps
template<class KeyType,typename T>
void MPSBase<KeyType,T>::generateBoundary(){
  if(!this->size()){
    cout<<"got an empty MPS ..." <<endl;
    return;
  }
  //first clear it ...
  _rightDummy.clear();
  typename BlockMatrix<1,KeyType,T>::iterator a;
  for(a=(*this)[this->size()-1].begin();a!=(*this)[this->size()-1].end();++a){
    TensorKeyType<2,1,KeyType>TKm(a->first[1],a->first[1],0);
    if(_rightDummy.find(TKm)==_rightDummy.end()){
      identity_matrix<T,column_major> mat(a->second.size2());
      //      mat.resize(a->second.size2(),a->second.size2());
      _rightDummy[TKm]=mat;
    }
  }
}

template<class KeyType,typename T>
UintVectorType MPSBase<KeyType,T>::getAuxSpace(){
  _dimAux.resize(this->size()-1);
  for(uint site=0;site<(this->size()-1);++site){
    _dimAux(site)=this->chi(site);
  }
  return _dimAux;
}

template<class KeyType,typename T>
UintVectorType MPSBase<KeyType,T>::getAuxSpace()const{
  UintVectorType dimAux(this->size()-1);
  for(int site=0;site<(this->size()-1);++site){
    dimAux(site)=this->chi(site);
  }
  return dimAux;
}


template<class KeyType,typename T>
uint MPSBase<KeyType,T>::chi(const int site)const{
  uint size=0;
  map_container<KeyType,uint> sizes;
  typename map_container<KeyType,uint>::iterator s;
  for(typename BlockMatrix<1,KeyType,T>::const_iterator a=(*this)[site].begin();a!=(*this)[site].end();++a){
    s=sizes.find(a->first[1]);
    if(s==sizes.end())
      sizes[a->first[1]]=a->second.size2();
  }
  uint chi=0;
  for(s=sizes.begin();s!=sizes.end();++s)
    chi+=s->second;
  return chi;
}


template<class KeyType,typename T>
void MPSBase<KeyType,T>::save(string filename){
  filename+=".mps";
  ofstream os(filename.c_str(), ios::binary);
  boost::archive::binary_oarchive oar(os);
  _numSites=this->size();
  oar << _numSites;
  for (int i=0; i<this->size(); i++){
    oar << (*this)[i];
  }
  oar << _dimP;
  oar << _dimAux;
  oar << _lOrtho;
  oar << _rOrtho;
  oar << _conservedQN;
  oar << _itokey;
  oar << _lambdas;
  oar << _leftDummy;
  oar << _rightDummy;
  os.close();
}


template<class KeyType,typename T>
void MPSBase<KeyType,T>::load(string filename,bool cast){
  // load mps from file
  this->clear();
  filename+=".mps";
  ifstream is(filename.c_str(), ios::binary);
  std::cout <<"#loading mps from file "<< filename<<": ";
  boost::archive::binary_iarchive iar(is);
  iar >> _numSites;
  this->resize(_numSites);
  //std::cout <<"#\t length of mps:"<< _numSites << std::endl;
  for (int i=0; i<_numSites; i++){
    iar >> (*this)[i];
  }
  iar >> _dimP;
  iar >> _dimAux;
  iar >> _lOrtho;
  iar >> _rOrtho;
  iar >> _conservedQN;
  iar >> _itokey;
  iar >> _lambdas;
  iar >> _leftDummy;
  iar >> _rightDummy;
  is.close();
  std::cout <<"#done..."<<std::endl;
  //initialize(*this);
}

template<class KeyType,typename T>
bool MPSBase<KeyType,T>::check()const{
  bool isOK = true;
  for(int site =0;site<this->size()-1;++site)
    if(!checkBMs((*this)[site],(*this)[site+1])) 
      isOK=false;
  return isOK;
}

template<class KeyType,typename T>
std::set<KeyType> MPSBase<KeyType,T>::getQN()const{
  std::set<KeyType> QN;
  typename BlockMatrix<1,KeyType,T>::const_iterator i;

  for(i=(*this)[this->size()-1].begin();i!=(*this)[this->size()-1].end();++i){
    QN.insert(i->first[1]);
  }
  return QN;
}

//constructor for product states
template<class KeyType,typename T>
MPSBase<KeyType,T>::MPSBase(const UintVectorType dimSiteSpace,const UintVectorType localState,std::vector<i_to_key<KeyType> > itokey):
  _dimP(dimSiteSpace),_lOrtho(false),_rOrtho(false)
{
  this->resize(localState.size());
  _numSites=localState.size();
  _itokey=itokey;
  matrix<T,column_major> mat(1,1);
  mat(0,0) = T(1);
  //  KeyType keyl(_itokey[0].begin()->second.size()),keyr;
  KeyType keyl=_itokey[0].begin()->second,keyr;
  reset(keyl);
  TensorKeyType<2,1,KeyType>TK(keyl,keyl,0);
  _leftDummy[TK] = mat;
  //cout<<keyl<<endl;
  for(int site=0;site<this->size();++site){
    assert(_itokey[site].size()==dimSiteSpace(site));
    assert(localState(site)<dimSiteSpace(site));
    (*this)[site].clear();
    if(site==0){reset(keyl);
    }else if(site<this->size()){keyl = keyr;}
    keyr = keyl+_itokey[site][localState[site]];
    TensorKeyType<2,1,KeyType>Tkey(keyl,keyr,localState(site));
    (*this)[site][Tkey]= mat;
  }
  TensorKeyType<2,1,KeyType>TKm(keyr,keyr,0);
  _rightDummy[TKm] = mat;
}


//not really necessary, but nicer ...
template<class KeyType,typename T>
class MPS : public MPSBase<KeyType,T>{
  template<class KT, typename VT>
  friend bool prepareSite(MPSBase<KT,VT> &mps,const int site,const string dir,const bool multiply,const bool storenorm,const bool truncate,const Real delta);
public:
  MPS():MPSBase<KeyType,T>(){};
  MPS(const uint size):MPSBase<KeyType,T>(size){};
  MPS(const UintVectorType dimP,std::vector<i_to_key<KeyType> > itokey):MPSBase<KeyType,T>(dimP,itokey){}
  MPS(const UintVectorType dim_site_space,const UintVectorType local_state,std::vector<i_to_key<KeyType> > itokey):
    MPSBase<KeyType,T>(dim_site_space,local_state,itokey){};
  MPS(const MPS<KeyType,T>&mps):MPSBase<KeyType,T>(mps){};
  MPS(const MPSBase<KeyType,T>&mps):MPSBase<KeyType,T>(mps){};
  virtual ~MPS(){};

  //operations ...
  MPS<KeyType,T>& operator +=(const MPS<KeyType,T>&mps){
    MPSBase<KeyType,T>::operator+=(mps);
    return *this;
  }
  MPS<KeyType,T> operator +(const MPS<KeyType,T>&mps)const{
    MPS<KeyType,T> cpy=*this;
    return cpy+=mps;
  }
};



template<class KeyType,typename T>
class CanonizedMPS : public MPSBase<KeyType,T>{
  //template<class KT, typename VT>
  //  friend bool prepareSite(CanonizedMPS<KT,VT> &mps,const int site,const string dir,const bool multiply,const bool truncate,const Real delta);
public:
  CanonizedMPS():MPSBase<KeyType,T>(){};
  CanonizedMPS(UintVectorType dim_site_space,UintVectorType local_state,std::vector<i_to_key<KeyType> > itokey):
    MPSBase<KeyType,T>(dim_site_space,local_state,itokey){};
  CanonizedMPS(const CanonizedMPS<KeyType,T>&mps):MPSBase<KeyType,T>(mps){this->_lOrtho=false;this->_rOrtho=false;};
  CanonizedMPS(const MPSBase<KeyType,T>&mps):MPSBase<KeyType,T>(mps){this->_lOrtho=false;this->_rOrtho=false;};
  CanonizedMPS(const MPS<KeyType,T>&mps):MPSBase<KeyType,T>(mps){this->_lOrtho=false;this->_rOrtho=false;};
  virtual ~CanonizedMPS(){};
  
  MPS<KeyType,T> toMPS()const 		//NUSS 09 05 2012 added this function
  {
    MPS<KeyType,T> mps(*this);
    mps[0] = (*this)[0];
    for(int siteInd = 1; siteInd < this->size(); siteInd++){
      if(this->_lambdas.find(siteInd-1)!=this->_lambdas.end()){
	BlockMatrix<1,KeyType,T>  bm = this->_lambdas.find(siteInd-1)->second * (*this)[siteInd];
	mps[siteInd] = bm;
      }else if(this->_lambdas.find(siteInd-1)==this->_lambdas.end()){
	cerr<<"in CanonizedMPS<KeyType,T>::toMPS(): no entry found in this->_lambdas at site "<<siteInd<<endl;
	abort();
      }
      //BlockMatrix<1,KeyType,T>  bm = this->_lambdas.at(siteInd-1) * (*this)[siteInd];
      //BlockMatrix<1,KeyType,T>  bm = this->_lambdas.find(siteInd-1) * this->find(siteInd);

    }
    mps.setLOrtho(true);
    mps.setROrtho(false);
    return mps;
  }
  
  //operations ...
  CanonizedMPS<KeyType,T>& operator +=(const CanonizedMPS<KeyType,T>&mps){
    MPSBase<KeyType,T>::operator+=(mps);
    for(uint site=0;site<this->size()-1;++site){
      typename DiagBlockMatrix<KeyType,Real>::iterator it;
      typename DiagBlockMatrix<KeyType,Real>::iterator itthis;
      for(it=mps.lambdas(site).begin();it!=mps.lambdas(site).end();++it){
	itthis=this->_lambas[site].find(it->first);
	if(itthis==this->_lambas[site].end())
	  this->_lambas[site][it->first]=it->second;
	else{
	  VectorType lam1=itthis->second,lam2=it->second;
	  lam1.resize(lam1.size()+lam2.size(),true);
	  boost::numeric::ublas::vector_range<boost::numeric::ublas::vector<double> > vr(lam1, boost::numeric::ublas::range(lam1.size(),lam1.size()+lam2.size()));
	  vr=lam2;
	}
      }
    }
    return *this;
  }
  CanonizedMPS<KeyType,T> operator +(const CanonizedMPS<KeyType,T>&mps)const{
    CanonizedMPS<KeyType,T> cpy=*this;
    return cpy+=mps;
  }

};
//=================== other operations ... ==========================

template<class KeyType,typename T>
MPSBase<KeyType,T> operator *(const MPSBase<KeyType,T>&mps,const T num){
  MPSBase<KeyType,T>cpy=mps;
  return cpy*=num;
}

template<class KeyType,typename T>
MPSBase<KeyType,T> operator *(const T num,const MPSBase<KeyType,T>&mps){
  MPSBase<KeyType,T>cpy=mps;
  return cpy*=num;
}

//this is a friend function of MPS
//always stores the lambas, so be sure that they originate from a valid preparation, where the mps is in an orthonormal state (either left or right,
//depending on dir)
//it also can store a lambda at the last+1 site, which usually doesn't belong to the system. if storenorm=true, then it stores the norm of the mps; 
template<class KeyType,typename T>
bool prepareSite(MPSBase<KeyType,T> &mps,const int site,const string dir, const bool multiply,const bool storenorm=false,
		 const bool truncate=false,const Real delta=0.0){
  bool erase=false;
  if(dir=="l"){
    assert(site<mps.size());
    assert(mps.size()>0);
    if(site<(mps.size()-1)){
      erase=prepareMatrix(mps[site], mps.lambda(site),mps[site+1],dir,multiply,truncate,delta);
    }
    else if(site==(mps.size()-1)){
      mps.generateBoundary();
      erase=prepareMatrix(mps[site], mps.lambda(site),mps._rightDummy,dir,multiply,truncate,delta);
      for(typename DiagBlockMatrix<KeyType,Real>::iterator l=mps.lambda(site).begin();l!=mps.lambda(site).end();++l){
	if(storenorm)
	  mps.norm()[l->first]=l->second;
      }
    }
  }
  if(dir=="r"){
    assert(site<mps.size());
    assert(mps.size()>0);
    if(site!=0){
      erase=prepareMatrix(mps[site-1], mps.lambda(site-1),mps[site],dir,multiply,truncate,delta);
    }
    else if(site==0){
      erase=prepareMatrix(mps._leftDummy, mps.lambda(site-1),mps[site],dir,multiply,truncate,delta);
      for(typename DiagBlockMatrix<KeyType,Real>::iterator l=mps.lambda(site-1).begin();l!=mps.lambda(site-1).end();++l){
	if(storenorm)
	  mps.norm()[l->first]=l->second;
      }
    }
  }
  if(dir!="l"&&dir!="r"){cout<<"prepareSite: wrong input string dir "<<dir<<"; use \"l\" or \"r\" instead"<<endl;}
  //TODO: normalize the lambdas; can be done inlyl after the FULL matrix is obtained!!!!
  return erase;
}


template<class KeyType,typename T1,typename T2>
MPS<KeyType,T1> cast(const MPS<KeyType,T2>&mps){
  MPS<KeyType,T1> mps_cast;
  mps_cast.resize(mps.size());
  for(uint site=0;site<mps.size();++site){
    mps_cast[site] = cast<1,KeyType,T1,T2>(mps[site]);
    mps_cast.setLambdas(mps.getLambdas());
    i_to_key<KeyType> itokey = mps.getItoKey(site);
    mps_cast.setItoKey(site,itokey);
  }
  mps_cast.setPSpace(mps.getPSpace());
  mps_cast.setAuxSpace(mps.getAuxSpace());
  mps_cast.setROrtho(mps.getROrtho());
  mps_cast.setLOrtho(mps.getLOrtho());
  return mps_cast;
}



#endif
