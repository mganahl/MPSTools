#ifndef MPOTYPES_HPP__
#define MPOTYPES_HPP__
#include "keyTypes.hpp"
#include "blasroutines.hpp"
#include "lapack_routines.hpp"
#include "MPSfunctions.hpp"
#include "serialization.hpp"
#include <map>
#include <omp.h>
#include <stdio.h>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>




template<typename T>
class Operator
{
public:
  Operator(){};
  //  Operator(const Operator &op):_dim1(op.size1()),_dim2(op.size2()){}
  virtual ~Operator(){_dim1=0;_dim2=0;};
  Operator(const Operator &o){
    _dim1 = o.size1();
    _dim2 = o.size2();
  }
  virtual T measure(const MPS<AbKeyType,T> &mps)const=0;
  virtual T measure(const CanonizedMPS<AbKeyType,T> &mps)const=0;

  // virtual int site()const{return -100000;};
  // virtual void setSite(const int site){};
  uint size1()const{return _dim1;}
  uint size2()const{return _dim2;}
  void resize(const uint s1,const uint s2){_dim1=s1;_dim2=s2;};

  template<class Archive>
  void serialize(Archive & ar, const unsigned int version){
    ar & _dim1;
    ar & _dim2;
  }
  
protected:
  uint _dim1;
  uint _dim2;
};
//====================================================SPARSE OPERATOR ==============================================
template<int Ns,typename T>
class SparseOperator:public map_container<OpKeyType<Ns>,T>,public Operator<T>
{
public:
  typedef typename map_container<OpKeyType<Ns>,T>::iterator iterator;
  typedef typename map_container<OpKeyType<Ns>,T>::const_iterator const_iterator;
  SparseOperator():Operator<T>(){_statistics=1;_Ns=Ns;_dims1.resize(Ns);_dims2.resize(Ns);_sites.resize(Ns);}
  SparseOperator(const int s1,const int s2){
    this->resize(s1,s2);   
    _statistics=1;_Ns=Ns;
    _sites.resize(_Ns);
  }
  SparseOperator(const SparseOperator<1,T>&o1,const SparseOperator<1,T>&o2){
    this->clear();
    assert(Ns==2);
    _Ns=Ns;
    this->resize(o1.size1()*o2.size1(),o1.size2()*o2.size2());
    (*this) = kron(o1,o2);
    _dims1.push_back(o1.size1());
    _dims1.push_back(o2.size1());
    _dims2.push_back(o1.size2());
    _dims2.push_back(o2.size2());
    _sites.resize(_Ns);
    _sites[0] = o1.sites()[0];
    _sites[1] = o2.sites()[0];
  }
  SparseOperator(const SparseOperator<1,T>&o1):map_container<OpKeyType<Ns>,T>(o1),Operator<T>(o1){
    _statistics = o1.getSign();
    _Ns = Ns;
    _dims1 = o1.dims1();
    _dims2 = o1.dims2();
    _sites=o1.sites();
  }
  int getSign()const{return _statistics;}
  void setSign(const int sign){assert(sign==-1||sign==1);_statistics=sign;}
  virtual T measure(const MPS<AbKeyType,T> &mps)const {return T(-1000);}
  virtual T measure(const CanonizedMPS<AbKeyType,T> &mps)const {return T(-1000);}
  std::vector<int> sites()const {return _sites;}
  void setSites(const std::vector<int> sites){_sites=sites;};
  uint length()const {return _Ns;}
  std::vector<uint>dims1()const{return _dims1;}
  uint &dims1(const int ind){return _dims1[ind];}
  void setDims1(const std::vector<uint>&d){_dims1 = d;}

  std::vector<uint>dims2()const{return _dims2;}
  uint &dims2(const int ind){return _dims2[ind];}
  void setDims2(const std::vector<uint>&d){_dims2 = d;}

  uint squeeze(const Real thresh);

  template<class Archive>
  void serialize(Archive & ar, const unsigned int version){
    ar & boost::serialization::base_object<Operator<T> >(*this);
    ar & boost::serialization::base_object<map_container<OpKeyType<Ns>,T > >(*this);

    ar &_statistics;
    ar &_Ns;
    ar &_dims1;
    ar &_dims2;
    ar &_sites;
  }

  
//operators;
  SparseOperator<Ns,T>&operator+=(const SparseOperator<Ns,T>&in);
  SparseOperator<Ns,T> operator+(const SparseOperator<Ns,T>&in)const;
  SparseOperator<Ns,T>&operator-=(const SparseOperator<Ns,T>&in);
  SparseOperator<Ns,T> operator-(const SparseOperator<Ns,T>&in)const;
  SparseOperator<Ns,T>&operator*=(const Real&in);
  SparseOperator<Ns,T>&operator*=(const Complex&in);
  SparseOperator<Ns,T>&operator/=(const Real&in){return (*this)*=(1./in);}
  SparseOperator<Ns,T>&operator/=(const Complex&in){return (*this)*=(1./in);}
  SparseOperator<Ns,T>&operator*=(const SparseOperator<Ns,T>&in);
  SparseOperator<Ns,T> operator*(const SparseOperator<Ns,T>&in)const;
protected:
  int _statistics;
  int _Ns;
  std::vector<uint>_dims1,_dims2;
  std::vector<int> _sites;
};

template<int Ns,typename T>
uint SparseOperator<Ns,T>::squeeze(const Real thresh){
  iterator i;
  uint out = 0;
  std::vector<iterator> erase;
  typename std::vector<iterator>::iterator e;
  for(i = this->begin();i!=this->end();++i){
    if(abs(i->second)<thresh)
      erase.push_back(i);
  }
  out = erase.size();
  for(e = erase.begin();e!=erase.end();++e){
    this->erase(*e);
  }
  return out;
}


template<int Ns,typename T>
SparseOperator<Ns,T>&SparseOperator<Ns,T>::operator+=(const SparseOperator<Ns,T>&in){
  if(this==&in){cout<<"in SparseOperator<Ns,T> operator +=: selfassignment"<<endl;abort();}
  assert(_Ns==in.length());
  assert(this->_sites==in.sites());
  assert(this->size1()==in.size1());
  assert(this->size2()==in.size2());
  for(const_iterator it=in.begin();it!=in.end();++it){
    //better save than sorry: don't know if T is initialized as T(0.0) on every platform
    if(this->find(it->first)==this->end()){
      (*this)[it->first]=it->second;
    }else{
      (*this)[it->first]+=it->second;
    }
  }
  _dims1=in.dims1();
  _dims2=in.dims2();
  return *this;
}

template<int Ns,typename T>
SparseOperator<Ns,T>&SparseOperator<Ns,T>::operator-=(const SparseOperator<Ns,T>&in){
  if(this==&in){cout<<"in SparseOperator<Ns,T> operator +=: selfassignment"<<endl;abort();}
  assert(_Ns==in.length());
  assert(this->_sites==in.sites());
  assert(this->size1()==in.size1());
  assert(this->size2()==in.size2());
  for(const_iterator it=in.begin();it!=in.end();++it){
    //better save than sorry: don't know if T is initialized as T(0.0) on every platform
    if(this->find(it->first)==this->end()){
      (*this)[it->first]=-it->second;
    }else{
      (*this)[it->first]-=it->second;
    }
  }
  _dims1=in.dims1();
  _dims2=in.dims2();
  return *this;
}

template<int Ns,typename T>
SparseOperator<Ns,T> SparseOperator<Ns,T>::operator+(const SparseOperator<Ns,T>&in)const{
  SparseOperator<Ns,T> cpy =(*this);
  return cpy+=in;
}

template<int Ns,typename T>
SparseOperator<Ns,T> SparseOperator<Ns,T>::operator-(const SparseOperator<Ns,T>&in)const{
  SparseOperator<Ns,T> cpy =(*this);
  return cpy-=in;
}

template<int Ns,typename T>
SparseOperator<Ns,T>&SparseOperator<Ns,T>::operator*=(const Real &in){
  for(iterator it=this->begin();it!=this->end();++it){
    (*this)[it->first]*=in;
  }
  return *this;
}

template<int Ns,typename T>
SparseOperator<Ns,T>&SparseOperator<Ns,T>::operator*=(const Complex &in){
  for(iterator it=this->begin();it!=this->end();++it){
    (*this)[it->first]*=in;
  }
  return *this;
}

//usual matrix product
template<int Ns,typename T>
SparseOperator<Ns,T>&SparseOperator<Ns,T>::operator*=(const SparseOperator<Ns,T>&in){
  assert(_Ns==in.length());
  assert(this->_sites==in.sites());
  assert(this->size1()==in.size1());
  assert(this->size2()==in.size2());
  SparseOperator<Ns,T> tmp(this->size1(),in.size2());
  iterator t,i;
  const_iterator j;
  for(i=this->begin();i!=this->end();++i){
    for(j=in.begin();j!=in.end();++j){
      if(i->first.getIn()==j->first.getOut()){
	OpKeyType<Ns>ky(i->first.getOut(),j->first.getIn());
	T val = i->second*j->second;
	t=tmp.find(ky);
	if(t!=tmp.end())
	  val+=t->second;
	tmp[ky] = val;
      }
    }
  }
  tmp.setDims1(this->dims1());
  tmp.setDims2(in.dims2());
  (*this) = tmp;
  return (*this);
}

template<int Ns,typename T>
SparseOperator<Ns,T>SparseOperator<Ns,T>:: operator*(const SparseOperator<Ns,T>&in)const {
  SparseOperator<Ns,T> cpy=(*this);
  return cpy*=in;
}

//====================================================END OF SPARSE OPERATOR ==============================================

template<int Ns,typename T1,typename T2>
SparseOperator<Ns,T1> operator*(const SparseOperator<Ns,T1>&op,const T2&num){
  SparseOperator<Ns,T1> cpy=op;
  return cpy*=num;
}

template<int Ns,typename T1,typename T2>
SparseOperator<Ns,T1> operator*(const T2&num,const SparseOperator<Ns,T1>&op){
  SparseOperator<Ns,T1> cpy=op;
  return cpy*=num;
}

template<int Ns,typename T1,typename T2>
SparseOperator<Ns,T1> operator/(const SparseOperator<Ns,T1>&op,const T2&num){
  SparseOperator<Ns,T1> cpy=op;
  return cpy/=num;
}

//I think we don't need this any more, at least not foer Fermi Hubbard; the signs of the gates are correct
//when we use the P-operators in the MPO
//template<typename T>
// SparseOperator<2,T> kron(const SparseOperator<1,T>&op1,const SparseOperator<1,T>&op2,i_to_key<AbKeyType> itokeyl,i_to_key<AbKeyType> itokeyr)
// {
//   SparseOperator<2,T>opout(op1.size1()*op2.size1(),op1.size2()*op2.size2());
//   typename SparseOperator<1,T>::const_iterator it1,it2;
//   for(it1 = op1.begin();it1!=op1.end();++it1){
//     for(it2 = op2.begin();it2!=op2.end();++it2){
//       std::vector<uint> in,out;
//       in.push_back(it1->first[0]);in.push_back(it2->first[0]);
//       out.push_back(it1->first(0));out.push_back(it2->first(0));
//       OpKeyType<2>k(out,in);
//       long unsigned int ql=0,deltaQ=0;
//       AbKeyType kyl=itokeyl[in[0]],kyr=itokeyr[in[1]],kyrp=itokeyr[out[1]];
//       assert(kyl.size());
//       assert(kyl.size()==kyr.size());
//       assert(kyl.size()==kyrp.size());
//       for(uint i=0;i<kyl.size();++i){
// 	ql+=kyl(i);
// 	deltaQ=abs(static_cast<int>(kyr(i))-static_cast<int>(kyrp(i)));
//       }
//       Real sig=pow(itokeyl.getSign(),ql*deltaQ);
//       //      opout[k]=it1->second*it2->second*pow(op2.getSign(),n1)*pow(op1.getSign(),n2p)*pow(nt,n1+n1p);
//       opout[k]=it1->second*it2->second*sig;
//     }
//   }
//   std::vector<uint>d1,d2;
//   d1.push_back(op1.size1());  d1.push_back(op2.size1());
//   d2.push_back(op1.size2());  d2.push_back(op2.size2());
//   opout.setDims1(d1);
//   opout.setDims2(d2);
//   assert(opout.getSign());
//   opout.setSign(1);
//   return opout;
// } 

//that's the usual kron operation for sparse operators living at different sites ...
//it's important that the _dims1 and _dims2 member is  set for both input operators!!!

template<typename T>
SparseOperator<2,T> kron(const SparseOperator<1,T>&op1,const SparseOperator<1,T>&op2)
{
  SparseOperator<2,T>opout(op1.size1()*op2.size1(),op1.size2()*op2.size2());
  typename SparseOperator<1,T>::const_iterator it1,it2;
  for(it1 = op1.begin();it1!=op1.end();++it1){
    for(it2 = op2.begin();it2!=op2.end();++it2){
      std::vector<uint> in,out;
      in.push_back(it1->first[0]);in.push_back(it2->first[0]);
      out.push_back(it1->first(0));out.push_back(it2->first(0));
      OpKeyType<2>k(out,in);
      opout[k]=it1->second*it2->second;
    }
  }
  std::vector<uint>d1,d2;
  d1.push_back(op1.size1());  d1.push_back(op2.size1());
  d2.push_back(op1.size2());  d2.push_back(op2.size2());
  opout.setDims1(d1);
  opout.setDims2(d2);
  assert(opout.getSign());
  opout.setSign(1);
  return opout;
} 


template<int N1,int N2,typename T>
SparseOperator<N1+N2,T> kron(const SparseOperator<N1,T>&op1,const SparseOperator<N2,T>&op2)
{
  SparseOperator<N1+N2,T>opout(op1.size1()*op2.size1(),op1.size2()*op2.size2());
  typename SparseOperator<N1,T>::const_iterator it1;
  typename SparseOperator<N2,T>::const_iterator it2;
  for(it1 = op1.begin();it1!=op1.end();++it1){
    for(it2 = op2.begin();it2!=op2.end();++it2){
      std::vector<uint> in,out;
      for(int l=0;l<N1;++l)
	in.push_back(it1->first[l]);
      for(int l=0;l<N2;++l)
	in.push_back(it2->first[l]);
      for(int l=0;l<N1;++l)
	out.push_back(it1->first(l));
      for(int l=0;l<N2;++l)
	out.push_back(it2->first(l));
      OpKeyType<N1+N2>k(out,in);
      opout[k]=it1->second*it2->second;
    }
  }
  std::vector<uint>d1,d2;
  for(int l=0;l<N1;++l)
    d1.push_back(op1.dims1()[l]);
  for(int l=0;l<N2;++l)
    d1.push_back(op2.dims1()[l]);

  for(int l=0;l<N1;++l)
    d2.push_back(op1.dims2()[l]);
  for(int l=0;l<N2;++l)
    d2.push_back(op2.dims2()[l]);

  opout.setDims1(d1);
  opout.setDims2(d2);
  opout.setSign(1);//that's a relict
  return opout;
} 


//that need's to be different name (kronecker instead of kron) because of overloading problems
//template<typename T>
//SparseOperator<1,T> kronecker(const SparseOperator<1,T>&op1,const SparseOperator<1,T>&op2){
//  SparseOperator<1,T>opout(op1.size1()*op2.size1(),op1.size2()*op2.size2());
//  typename SparseOperator<1,T>::const_iterator it1;
//  typename SparseOperator<1,T>::const_iterator it2;
//  uint i1,i2,j1,j2;
//  for(it1 = op1.begin();it1!=op1.end();++it1){
//    for(it2 = op2.begin();it2!=op2.end();++it2){
//      i1=it1->first(0);
//      i2=it2->first(0);
//
//      j1=it1->first[0];
//      j2=it2->first[0];
//
//      OpKeyType<1>k(i1+op1.size1()*i2,j1+op1.size2()*j2);
//      opout[k]=it1->second*it2->second;
//    }
//  }
//  std::vector<uint> d1,d2;
//  d1.push_back(op1.size1()*op2.size1());
//  d2.push_back(op1.size2()*op2.size2());
//  opout.setDims1(d1);
//  opout.setDims2(d2);
//  opout.setSign(1);//that's a relict
//  return opout;
//} 
//



template<int Ns,typename T>
ostream& operator<<(ostream& o,const SparseOperator<Ns,T>&in){
  typename SparseOperator<Ns,T>::const_iterator i;
  o<<"\t size of Operator: "<<in.size()<<endl;
  for(i=in.begin();i!=in.end();++i){
    cout<<"\t \t in"<<i->first.getK()<<" out"<<i->first.getP()<<"  "<<i->second<<endl;
  }
  return o;
}

//====================================================SPARSE LOCAL OPERATOR ==============================================

template<typename T>
class SparseLocalOperator:public SparseOperator<1,T>{
public:
  typedef typename SparseOperator<1,T>::const_iterator const_iterator;
  typedef typename SparseOperator<1,T>::iterator iterator;
  SparseLocalOperator():SparseOperator<1,T>(){this->_sites.resize(1);this->_sites[0] = -100000;}
  SparseLocalOperator(const SparseLocalOperator<T>&op):SparseOperator<1,T>(op){
    this->_sites.resize(1);
    this->_sites[0] = op.site();
  }
  SparseLocalOperator(const SparseOperator<1,T>&op,const int site):SparseOperator<1,T>(op){
    this->_sites.resize(1);
    this->_sites[0] = site;
  }

  SparseLocalOperator(const int s1,const int s2){
    this->resize(s1,s2);   
    this->_statistics=1;this->_Ns=1;
  }

  virtual T measure(const MPS<AbKeyType,T> &mps)const{return T(-1000);};
  virtual T measure(const CanonizedMPS<AbKeyType,T> &mps)const;
  //  virtual Complex measure(const MPS<AbKeyType,Complex> &mps)const;  
  //int site;
  int site()const {if(this->_sites[0]==-10000){cerr<<"warning in SparseLocalOperator::site(): requesting uninitialized value of _site"<<endl;} return this->_sites[0];}
  void setSite(const int site){this->_sites[0]=site;}

  SparseLocalOperator<T>&operator+=(const SparseLocalOperator<T>&in){
    assert(in.site()==this->_sites[0]);
    SparseOperator<1,T>::operator+=(in);
    return *this;
  }
  SparseLocalOperator<T> operator+(const SparseLocalOperator<T>&in)const{
    assert(in.site()==this->_sites[0]);
    SparseLocalOperator<T> cpy =(*this);
    return cpy+=in;
  }
  SparseLocalOperator<T>&operator-=(const SparseLocalOperator<T>&in){
    assert(in.site()==this->_sites[0]);
    SparseOperator<1,T>::operator-=(in);
    return *this;
  }
  SparseLocalOperator<T> operator-(const SparseLocalOperator<T>&in)const{
    assert(in.site()==this->_sites[0]);
    SparseLocalOperator<T> cpy =(*this);
    return cpy-=in;
  }
  SparseLocalOperator<T>&operator*=(const Real&in){
    SparseOperator<1,T>::operator*=(in);
    return *this;
  }
  SparseLocalOperator<T>&operator/=(const Real&in){
    SparseOperator<1,T>::operator/=(in);
    return *this;
  }
  SparseLocalOperator<T>&operator*(const Real&in)const{
    SparseLocalOperator<T> cpy=(*this);
    return cpy*=(in);
  }
  SparseLocalOperator<T>&operator*=(const Complex&in){
    SparseOperator<1,T>::operator*=(in);
    return *this;
  }
  SparseLocalOperator<T>&operator/=(const Complex&in){
    SparseOperator<1,T>::operator/=(in);
    return *this;
  }

  SparseLocalOperator<T>&operator*(const Complex&in)const{
    SparseLocalOperator<T> cpy=(*this);
    return cpy*=(in);
  }

  SparseLocalOperator<T>&operator*=(const SparseLocalOperator<T>&in){
    assert(in.site()==this->_sites[0]);
    SparseOperator<1,T>::operator*=(in);
    this->_sites[0]=in.site();
    return *this;
  }

  SparseLocalOperator<T>operator*(const SparseLocalOperator<T>&in)const{
    assert(in.site()==this->_sites[0]);
    SparseLocalOperator<T>cpy=(*this);
    return cpy*=in;
  }
};

template<typename T>
SparseLocalOperator<T> kronecker(const SparseLocalOperator<T>&op1,const SparseLocalOperator<T>&op2){
  SparseLocalOperator<T>opout(op1.size1()*op2.size1(),op1.size2()*op2.size2());
  typename SparseLocalOperator<T>::const_iterator it1;
  typename SparseLocalOperator<T>::const_iterator it2;
  uint i1,i2,j1,j2;
  for(it1 = op1.begin();it1!=op1.end();++it1){
    for(it2 = op2.begin();it2!=op2.end();++it2){
      i1=it1->first(0);
      i2=it2->first(0);

      j1=it1->first[0];
      j2=it2->first[0];

      OpKeyType<1>k(i1+op1.size1()*i2,j1+op1.size2()*j2);
      opout[k]=it1->second*it2->second;
    }
  }
  std::vector<uint> d1,d2;
  d1.push_back(op1.size1()*op2.size1());
  d2.push_back(op1.size2()*op2.size2());
  opout.setDims1(d1);
  opout.setDims2(d2);
  opout.setSign(1);//that's a relict
  return opout;
} 

template<typename T>
T SparseLocalOperator<T>:: measure(const CanonizedMPS<AbKeyType,T> &mps)const{

  if(mps.size()==0){cout<<"in compute Mean(SparseLocalOperator, MPS): size of CanonizedMPS<AbKeyType,T> = 0"<<endl;return-1000;}
  int site = this->site();
  assert(site>=-1);
  T mean = T(0.0);
  using boost::numeric::ublas::matrix;
  using boost::numeric::ublas::column_major;
  matrix<T,column_major>tmp;
  typename BlockMatrix<1,AbKeyType,T >::const_iterator it,it_p;
  SparseLocalOperator<T>::const_iterator o;
  BlockMatrix<1,AbKeyType,T> M;
  M = (mps.lambda(site-1)*mps[site])*mps.lambda(site);
  for(o=this->begin();o!=this->end();++o){
    if(abs(o->second)<1e-16)continue;
    for(it = M.begin();it!=M.end();++it){
      if(it->first(0)!=o->first[0])continue;
      AbKeyType ql = it->first[0];
      AbKeyType qr = it->first[1];
      AbKeyType qr_prime = ql+mps.getItoKey(site)[o->first(0)];
      if(qr==qr_prime){
	it_p = M.find(TensorKeyType<2,1,AbKeyType>(ql,qr_prime,o->first(0)));
	if(it_p!=M.end()){
	  MatMatProd(T(1.0),conj(it_p->second),it->second,tmp,'t','n');
	  mean += o->second*trace(tmp);
	}
      }
    }
  }
  return mean;
}

//====================================================END OF SPARSE LOCAL OPERATOR ==============================================

template<typename T1,typename T2>
SparseLocalOperator<T1> operator*(const SparseLocalOperator<T1>&op,const T2&num){
  SparseLocalOperator<T1> cpy=op;
  return cpy*=num;
}

template<typename T1,typename T2>
SparseLocalOperator<T1> operator*(const T2&num,const SparseLocalOperator<T1>&op){
  SparseLocalOperator<T1> cpy=op;
  return cpy*=num;
}
template<typename T1,typename T2>
SparseLocalOperator<T1> operator/(const SparseLocalOperator<T1>&op,const T2&num){
  SparseLocalOperator<T1> cpy=op;
  return cpy*=num;
}


template<typename T>
class Sparse2SiteOperator:public SparseOperator<2,T>{
public:
  typedef typename SparseOperator<2,T>::const_iterator const_iterator;
  typedef typename SparseOperator<2,T>::iterator iterator;
  Sparse2SiteOperator():SparseOperator<2,T>(){this->_sites.resize(2);this->_sites[0] = -100000;this->_sites[1] = -100000;}
  Sparse2SiteOperator(const Sparse2SiteOperator<T>&op):SparseOperator<2,T>(op){
    this->_sites.resize(2);
    this->_sites = op.sites();
  }
  Sparse2SiteOperator(const int s1,const int s2){
    this->resize(s1,s2);   
    this->_statistics=1;this->_Ns=1;
  }

  virtual T measure(const MPS<AbKeyType,T> &mps)const{return T(-1000);};
  virtual T measure(const CanonizedMPS<AbKeyType,T> &mps)const;
  //  virtual Complex measure(const MPS<AbKeyType,Complex> &mps)const;  
  //int site;
  std::vector<int> sites()const {if(this->_sites[0]==-10000){cerr<<"warning in Sparse2SiteOperator::site(): requesting uninitialized value of _site"<<endl;}
    return this->_sites;
  }
  void setSites(const int site1, const int site2){this->_sites.resize(2);this->_sites[0]=site1;this->_sites[1] = site2;}

  Sparse2SiteOperator<T>&operator+=(const Sparse2SiteOperator<T>&in){
    assert(in.sites()==this->_sites);
    SparseOperator<2,T>::operator+=(in);
    return *this;
  }
  Sparse2SiteOperator<T> operator+(const Sparse2SiteOperator<T>&in)const{
    assert(in.sites()==this->_sites);
    Sparse2SiteOperator<T> cpy =(*this);
    return cpy+=in;
  }

  Sparse2SiteOperator<T>&operator*=(const Real&in){
    SparseOperator<2,T>::operator*=(in);
    return *this;
  }
  Sparse2SiteOperator<T>&operator*=(const Complex&in){
    SparseOperator<2,T>::operator*=(in);
    return *this;
  }

  Sparse2SiteOperator<T>&operator*=(const Sparse2SiteOperator<T>&in){
    assert(in.sites()==this->_sites);
    SparseOperator<2,T>::operator*=(in);
    return *this;
  }

  Sparse2SiteOperator<T>operator*(const Sparse2SiteOperator<T>&in)const{
    assert(in.sites()==this->_sites);
    Sparse2SiteOperator<T>cpy=(*this);
    return cpy*=in;
  }
};


template<typename T>
T Sparse2SiteOperator<T>:: measure(const CanonizedMPS<AbKeyType,T> &mps)const{

  if(mps.size()==0){cout<<"in compute Mean(SparseLocalOperator, MPS): size of CanonizedMPS<AbKeyType,T> = 0:"<<endl;return T(-1000);}
  int site1=this->sites()[0];
  int site2=this->sites()[1];
  
  if(this->sites().size()!=2){cerr<<"in Sparse2SiteOperator<Ns,T>::measure(MPS<KeyType,T>): only two-site operators are implemented; leaving function ..."<<endl; 
    return T(-1000);
  }
  if(site2!=site1+1){
    cerr<<"in Sparse2SiteOperator<T>::measure(MPS<KeyType,T>): only nearest neighbor measurement implemented; leaving function ..."<<endl;return T(-1000);
  }
  T mean = T(0.0);
  using boost::numeric::ublas::matrix;
  using boost::numeric::ublas::column_major;
  matrix<T,column_major>tmp;
  typename BlockMatrix<2,AbKeyType,T >::const_iterator it,it_p;
  const_iterator o;
  BlockMatrix<1,AbKeyType,T> Ml,Mr;

  Ml = (mps.lambda(site1-1)*mps[site1])*mps.lambda(site1);
  Mr = mps[site2]*mps.lambda(site2);

  const BlockMatrix<2,AbKeyType,T> M(Ml,Mr);
  for(o=this->begin();o!=this->end();++o){
    if(abs(o->second)<1e-16)continue;
    for(it = M.begin();it!=M.end();++it){
      if(it->first(0)!=o->first[0])continue;
      if(it->first(1)!=o->first[1])continue;
      AbKeyType ql = it->first[0];
      AbKeyType qr = it->first[1];
      AbKeyType qr_prime = ql+mps.getItoKey(site1)[o->first(0)]+mps.getItoKey(site2)[o->first(1)];
      if(qr==qr_prime){
   	it_p = M.find(TensorKeyType<2,2,AbKeyType>(ql,qr_prime,o->first.getOut()));
   	if(it_p!=M.end()){
   	  MatMatProd(T(1.0),conj(it_p->second),it->second,tmp,'t','n');
   	  mean += o->second*trace(tmp);
   	}
      }
    }
  }
  return mean;
}

//==================================================== END OF SPARSE2SITEOPERATOR=============================================================0


template<typename T1,typename T2>
Sparse2SiteOperator<T1> operator*(const Sparse2SiteOperator<T1>&op,const T2&num){
  Sparse2SiteOperator<T1> cpy=op;
  return cpy*=num;
}

template<typename T1,typename T2>
Sparse2SiteOperator<T1> operator*(const T2&num,const Sparse2SiteOperator<T1>&op){
  Sparse2SiteOperator<T1> cpy=op;
  return cpy*=num;
}

template<typename T>
class Moment:public Operator<T>
{
public:
  Moment(uint site, int m){_site = site; this->m=m;}
  Moment(){};
  ~Moment(){};
  int site()const{return _site;}
  void setSite(const int site){_site=site;}
  virtual T measure(const MPS<AbKeyType,T> &mps)const{return compMom(mps);}
  virtual T measure(const CanonizedMPS<AbKeyType,T> &mps)const{return compMom(mps);}
private:
  int _site, m;
  T compMom(const MPSBase<AbKeyType,T>&mps)const;
};

template<typename T>
T Moment<T>::compMom(const MPSBase<AbKeyType,T>&mps)const{
  DiagBlockMatrix<AbKeyType,Real>::iterator it;
  DiagBlockMatrix<AbKeyType,Real> lam = mps.lambda(_site);
  T mom = T(0.0);
  for(it = lam.begin();it!=lam.end();++it){
      VectorType::iterator itv;
      T sl = T(0.0);
      for(itv = it->second.begin();itv!=it->second.end();++itv){
    	  sl+=(*itv)*(*itv);
      }
      int qn=it->first.getQ();
      mom+=pow(qn,m)*sl;
    }
  return mom;
}

template<typename T>
class Entropy:public Operator<T>
{
public:
  Entropy(uint site){_site = site;}
  Entropy(){};
  ~Entropy(){};
  int site()const{return _site;}
  void setSite(const int site){_site=site;}
  virtual T measure(const MPS<AbKeyType,T> &mps)const{return compEnt(mps);}
  virtual T measure(const CanonizedMPS<AbKeyType,T> &mps)const{return compEnt(mps);}
private:
  int _site;
  T compEnt(const MPSBase<AbKeyType,T>&mps)const;
};
//base-2 logarithm is used in the calculations
template<typename T>
T Entropy<T>::compEnt(const MPSBase<AbKeyType,T>&mps)const{
  DiagBlockMatrix<AbKeyType,Real>::iterator it;
  DiagBlockMatrix<AbKeyType,Real> lam = mps.lambda(_site);
  T ent = T(0.0);
  for(it = lam.begin();it!=lam.end();++it){
    VectorType::iterator itv;
    for(itv = it->second.begin();itv!=it->second.end();++itv){
      ent-=((*itv)*(*itv))*log((*itv)*(*itv))/log(2);
    }
  }
  return ent;
}

template<typename T>
class maxNorm:public Operator<T>{
public:
  maxNorm(const BlockMatrix<1,AbKeyType,T> &BM,const int site){_BM=BM;_site=site;};
  ~maxNorm(){};
  virtual T measure(const CanonizedMPS<AbKeyType,T> &mps)const;
  virtual T measure(const MPS<AbKeyType,T> &mps)const{return -1000;};
private:
  BlockMatrix<1,AbKeyType,T> _BM;
  int _site;
};
template<typename T>
T maxNorm<T>::measure(const CanonizedMPS<AbKeyType,T> &mps)const{
  typename BlockMatrix<1,AbKeyType,T>::iterator b1,it;
  typename BlockMatrix<1,AbKeyType,T>::const_iterator b2;
  Real val=0.0;
  if(_site>=mps.size()){
    cout<<"#in maxNorm: _site>=mps.size(); skipping measurement"<<endl;
    return Complex(1000.,1000.);
  }
  BlockMatrix<1,AbKeyType,T> t=mps[_site]*mps.lambda(_site);
  //cout<<_BM<<endl;
  it=t.begin();
  uint u=t.size();
  uint bla=0;
  while(bla<(u-u%2)/2){
    ++bla;
    ++it;
  }
  cout<<it->first<<endl;
  for(uint n=0;n<it->second.size1()*it->second.size2();++n){
    cout<<abs(it->second.data()[n])<<endl;
  }
    
  //  printboostmatrix(abs(it->second));

  for(b1=t.begin();b1!=t.end();++b1){
    b2=_BM.find(b1->first);
    if(b2!=_BM.end()){
      matrix<T,column_major> mat1=b1->second;
      matrix<T,column_major> mat2=b2->second;
      uint size1=mat1.size1()>mat2.size1()?mat1.size1():mat2.size1();
      uint size2=mat1.size2()>mat2.size2()?mat1.size2():mat2.size2();
      matrix<T,column_major> t1(size1,size2);
      matrix<T,column_major> t2(size1,size2);
      boost_setTo(t1,T(0.0));boost_setTo(t2,T(0.0));
      matrix_range<matrix<T, column_major> >(t1,Range(0,mat1.size1()),Range(0,mat1.size2()))= mat1;
      matrix_range<matrix<T, column_major> >(t2,Range(0,mat2.size1()),Range(0,mat2.size2()))= mat2;
      for(uint n=0;n<size1*size2;++n){
	val+=(abs(t1.data()[n])-abs(t2.data()[n]));
      }
    }
    if(b2==_BM.end()){
      for(uint n=0;n<b1->second.size1()*b1->second.size2();++n){
	val+=abs(b1->second.data()[n]);
      }
    }
      
  }
  for(b2=_BM.begin();b2!=_BM.end();++b2){
    if(t.find(b2->first)==t.end()){
      for(uint n=0;n<b2->second.size1()*b2->second.size2();++n){
	val+=abs(b2->second.data()[n]);
      }
    }
  }
  cout<<"the difference: "<<val<<endl;
  my_pause();

  return Complex(val,0.0);
}

//purely local terms should be at the left bottom of the MPO
template<int Ns,typename T>
class MPOMat:public map_container<ukey<2>,SparseOperator<Ns,T> >{
public:
  typedef typename  map_container<ukey<2>,SparseOperator<Ns,T> >::iterator iterator;
  typedef typename map_container<ukey<2>,SparseOperator<Ns,T> >::const_iterator const_iterator;

  template<class Archive>
  void serialize(Archive & ar, const unsigned int version){
    ar & boost::serialization::base_object<map_container<ukey<2>,SparseOperator<Ns,T> > >(*this);
    ar &Dl;
    ar &Dr;
    ar &_Ns;
  }

  MPOMat(){_Ns=Ns;Dl=0;Dr=0;};
  template<class KeyType>
  MPOMat(const MPOMat<1,T>&m1,const MPOMat<1,T>&m2,i_to_key<KeyType> itokeyl,i_to_key<KeyType> itokeyr){
    assert(Ns==2);
    _Ns=Ns;
    this->clear();
    typename MPOMat<1,T>::const_iterator it1,it2;
    for(it1 = m1.begin();it1!=m1.end();++it1){
      for(it2 = m2.begin();it2!=m2.end();++it2){
	if(it1->first[1]==it2->first[0]){
	  uint d1 =it1->second.size1()*it2->second.size1();
	  uint d2 =it1->second.size2()*it2->second.size2();
	  SparseOperator<2,T> lo(d1,d2);
	  ukey<2>k(it1->first[0],it2->first[1]);
	  lo = kron(it1->second,it2->second,itokeyl,itokeyr);
	  if(this->find(k)==this->end())
	    (*this)[k] = lo;
	  else
	    (*this)[k]+=lo;
	}
      }
    }
    assert(m1.Dl);
    assert(m2.Dr);
    Dl = m1.Dl;
    Dr = m2.Dr;
  }
  MPOMat(const MPOMat<1,T>&m1):map_container<ukey<2>,SparseOperator<Ns,T> >(m1){_Ns = Ns;Dl=m1.Dl;Dr = m1.Dr;}
  ~MPOMat(){};
  void clear(){map_container<ukey<2>,SparseOperator<Ns,T> >::clear();Dl=0;Dr=0;}
  uint squeeze(const Real thresh);
  uint Dl;
  uint Dr;

private:
  int _Ns;
};

template<int Ns,typename T>
uint MPOMat<Ns,T>::squeeze(const Real thresh){
  iterator i;
  uint out = 0;
  for(i = this->begin();i!=this->end();++i){
    out+=i->second.squeeze(thresh);
  }
  return out;
}

//========================================outstream stuff===============================
template<int Ns,typename T>
ostream& operator <<(ostream & o,const MPOMat<Ns,T>&M)
{
  typename MPOMat<Ns,T>::const_iterator m;
  o<<"size of MPOMat<"<<Ns<<",T>: "<<M.size()<<endl;
  for(m=M.begin();m!=M.end();++m){
    o<<"aux-key: ("<<m->first[0]<<","<<m->first[1]<<")"<<endl;
    cout<<m->second<<endl;
  }
  return o;
}

//watch out: Operator<T> and std::vector both have member resize overloaded
template<typename T>
class MPO:public map_container<int,MPOMat<1,T> >,public Operator<T>{
public:
  typedef typename map_container<int,MPOMat<1,T> >::iterator iterator;
  typedef typename map_container<int,MPOMat<1,T> >::const_iterator const_iterator;
  MPO(){};
  MPO(const MPO<T>&mpo):map_container<int,MPOMat<1,T> >(mpo){};
  MPO(const SparseLocalOperator<T> &o1,const SparseLocalOperator<T> &o2);
  ~MPO(){};
  SparseOperator<2,Complex> getNNGate(const MPS<AbKeyType,T>&mps,const int site1,const int site2,const Complex dt);		//NUSS 09 05 2012 dt to Complex
  SparseOperator<2,Real> getNNGateReal(const MPS<AbKeyType,T>&mps,const int site1,const int site2,const Real dt);		//NUSS 09 05 2012 dt to Complex
  SparseOperator<3,Complex> getNNNGate(const MPS<AbKeyType,T>&mps,const int site1,const int site2,const int site3,const Complex dt);//NUSS 09 05 2012 dt to Complex
  MPO<Real> computeNNGateMPOReal(const int even,const MPS<AbKeyType,Real>&mps,const Real dt);  
  void save(string filename);
  void load(string filename);
  virtual T measure(const MPS<AbKeyType,T>&mps)const;
  virtual T measure(const CanonizedMPS<AbKeyType,T>&mps)const;
  virtual void construct(const uint step){};//does nothing yet: on a higher level, construct is used in derived classes to construct time dependent hamiltonians
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version){
    ar & boost::serialization::base_object<Operator<T> >(*this);
    ar & boost::serialization::base_object<map_container<int,MPOMat<1,T> > >(*this);
  }
  //  void resize(const uint N){std::vector<MPOMat<1,T> >::resize(N);}
};


template<typename T>
MPO<T>::MPO(const SparseLocalOperator<T> &o1,const SparseLocalOperator<T> &o2){
  MPOMat<1,T> m1,m2;
  m1[ukey<2>(0,0)] = o1;  m2[ukey<2>(0,0)] = o2;
  int site1 = o1.site();
  int site2 = o2.site();
  assert(site1!=site2);
  this->insert(pair<int,MPOMat<1,T> >(site1,m1));
  (*this)[site1].Dl=1;
  (*this)[site1].Dr=1;
  this->insert(pair<int,MPOMat<1,T> >(site2,m2));
  (*this)[site2].Dl=1;
  (*this)[site2].Dr=1;
}


template<typename T>
void MPO<T>::save(string filename){
    ofstream os(filename.c_str(), ios::binary);
    boost::archive::binary_oarchive oar(os);
    int numSites=this->size();
    oar << numSites;
    for (int i=0; i<this->size(); i++){
    	oar << (*this)[i];
    }
    os.close();
}


template<typename T>
void MPO<T>::load(string filename){
  // load mps from file
  this->clear();
  ifstream is(filename.c_str(), ios::binary);
  std::cout << filename<< std::endl;
  boost::archive::binary_iarchive iar(is);
  int numSites;
  iar >> numSites;
  for (int i=0; i<numSites; i++){
    iar >>(*this)[i];
    //std::cout << (*this)[i].size() <<std::endl;
  }
  is.close();
  //initialize(*this);
}

//the full contraction is done... this is not optimal
template<typename T>
T MPO<T>::measure(const MPS<AbKeyType,T>&mps)const{
  BlockMatrix<1,AbKeyType,T> L;
  MPO<T>_mpo=(*this);
  //  computeL(mps,(*this),0, R);
  const int N=mps.size();
  assert(N==this->size());
  //if(!mps.getLOrtho()){cerr<<"computeL: mps is not l-orthonormal"<<endl;}
  L.clear();
  L= getLMin(mps);
  for(int s = 0;s<N;++s){
    L = addContractionLayer(L,mps[s],_mpo[s],mps.getItoKey(s),1);
  }

  typename BlockMatrix<1,AbKeyType,T>::iterator l;
  T out = T(0.0);
  for(l=L.begin();l!=L.end();++l){
    if(l->first[0]==l->first[1]){
      matrix<T,column_major>mat=l->second;
      int size=mat.size1()>mat.size2()?mat.size1():mat.size2();
      mat.resize(size,size);
      boost_setTo(mat,T(0.0));
      matrix_range<matrix<T,column_major> >(mat,Range(0,l->second.size1()),Range(0,l->second.size2()))=l->second;
      out+=trace(mat);
    }
  }
  // return (L[0].begin()->second(0,0));
  return out;
}



//does a PARTIAL contraction of the network; 
template<typename T>
T MPO<T>::measure(const CanonizedMPS<AbKeyType,T>&mps)const{
  if(!this->size()){
    cout<<"MPO<T>::measure(const CanonizedMPS<AbKeyType,T>): detected empty map"<<endl;
    return T(-1000);
  }
  if(!mps.size()){
    cout<<"MPO<T>::measure(const CanonizedMPS<AbKeyType,T>): detected size of CanonizedMPS<AbKeyType,T> = 0; skipping measurement"<<endl;
    return T(-1000);
  }
  MPO<T>cpy=(*this);
  map_container<int,BlockMatrix<1,AbKeyType,T> >L;
  typename BlockMatrix<1,AbKeyType,T>::iterator itL;
  iterator it;
  it= cpy.begin();
  T val = T(0.0);
  BlockMatrix<1,AbKeyType,T> A;
  int site1=it->first;
  int index = -1;
  A = mps.lambda(site1-1)*mps[site1]*mps.lambda(site1);
  L[index] = getBoundaryBL(mps[site1],0);
  ++index;
  L[index]=addContractionLayer(L[index-1],A,cpy.find(site1)->second,mps.getItoKey(site1),1);
  ++index;
  uint Dr=cpy.find(site1)->second.Dr;
  ++it;
  while(it!=cpy.end()){
    int site2=it->first;
    if(site1!=site2){
      ++site1;
      while(site1<(site2)){
	SparseLocalOperator<T> localId(mps.getPSpace()(site1),mps.getPSpace()(site1));
	for(uint ind = 0;ind<mps.getPSpace()(site1);++ind){
	  localId[OpKeyType<1>(ind,ind)] = T(1.0);      
	}
	MPOMat<1,T> id;
	for(uint auxi=0;auxi<Dr;++auxi)
	  id[ukey<2>(auxi,auxi)] = localId;
	A = mps[site1]*mps.lambda(site1);
	L[index]=addContractionLayer(L[index-1],A,id,mps.getItoKey(site1),1);
	++index;
	++site1;
      }
      A = mps[site2]*mps.lambda(site2);
      L[index]=addContractionLayer(L[index-1],A,cpy.find(site2)->second,mps.getItoKey(site2),1);
      index++;
      Dr=cpy.find(site2)->second.Dr;
      site1=site2;
      ++it;
    }else{
      cout<<"in MPO<T>::measure(CanonizedMPS<AbKeyType,T>): site1==site2 detected;"<<endl;
      return T(-1000);
    } 
  }
  for(itL=L[index-1].begin();itL!=L[index-1].end();itL++){
    if(itL->first[0]==itL->first[1])
      val+=trace(itL->second);
  }
  return val;
}


//IMPORTANT: assumes that the full two site interaction can be obtained from the last row of m1 * the first column of m2
//the B field term (and any other local terms) must be located at the left bottom of the Matrix
//all local terms involving something like chemical potentials or B-field really NEED to be put at the left-bottom corner of the MPO matrix;
//otherwise the gates at the boundaries are wrong;

template<typename T>
SparseOperator<2,Complex> MPO<T>::getNNGate(const MPS<AbKeyType,T>&mps,const int site1,const int site2,const Complex dt){//NUSS 09 05 2012 dt to Complex

  if(site2!=site1+1){cout<<"MPO<T>::getNNGate(): assertion site2=site1+1 failed"<<endl;abort();}
  MPOMat<1,T> m1=(*this)[site1],m2=(*this)[site2];
  i_to_key<AbKeyType>itokeyl = mps.getItoKey(site1),itokeyr = mps.getItoKey(site2);
  SparseOperator<2,T>gate(m1.begin()->second.size1()*m2.begin()->second.size1(),m1.begin()->second.size2()*m2.begin()->second.size2());
  assert(m1.Dr==m2.Dl);
  const uint lend=m1.Dl-1;
  for(uint i = 0;i<m1.Dr;++i){
    if(site1==0&&site2!=(this->size()-1)){
      if(i==0){
	if(m1.find(ukey<2>(0,i))==m1.end())
	  continue;
	if(m2.find(ukey<2>(i,0))==m2.end())
	  continue;
	gate+=kron(m1[ukey<2>(0,i)],m2[ukey<2>(i,0)]);
      }
      else if(i==(m1.Dr-1)){
	if(m1.find(ukey<2>(0,i))==m1.end())
	  continue;
	if(m2.find(ukey<2>(i,0))==m2.end())
	  continue;

	gate+=kron(m1[ukey<2>(0,i)],m2[ukey<2>(i,0)]*0.5);
      }
      else if(i!=0&&i!=(m1.Dr-1)){
	if(m1.find(ukey<2>(0,i))==m1.end())
	  continue;
	if(m2.find(ukey<2>(i,0))==m2.end())
	  continue;
      	gate+=kron(m1[ukey<2>(0,i)],m2[ukey<2>(i,0)]);
      }
    }
    if((site2==(this->size()-1))&&site1!=0){
      if(i==0){
	if(m1.find(ukey<2>(lend,i))==m1.end())
	  continue;
	if(m2.find(ukey<2>(i,0))==m2.end())
	  continue;
      	gate+=kron(m1[ukey<2>(lend,i)]*0.5,m2[ukey<2>(i,0)]);
      }
      else if(i==(m1.Dr-1)){
	if(m1.find(ukey<2>(lend,i))==m1.end())
	  continue;
	if(m2.find(ukey<2>(i,0))==m2.end())
	  continue;
	gate+=kron(m1[ukey<2>(lend,i)],m2[ukey<2>(i,0)]);
      }
      else if(i!=(m1.Dr-1)){
	if(m1.find(ukey<2>(lend,i))==m1.end())
	  continue;
	if(m2.find(ukey<2>(i,0))==m2.end())
	  continue;
	gate+=kron(m1[ukey<2>(lend,i)],m2[ukey<2>(i,0)]);
      }
    }
    if((site2!=(this->size()-1))&&site1!=0)
      if(i==0){
	if(m1.find(ukey<2>(lend,i))==m1.end())
	  continue;
	if(m2.find(ukey<2>(i,0))==m2.end())
	  continue;
	gate+=kron(m1[ukey<2>(lend,i)]*0.5,m2[ukey<2>(i,0)]);
      }
      else if(i==(m1.Dr-1)){
	if(m1.find(ukey<2>(lend,i))==m1.end())
	  continue;
	if(m2.find(ukey<2>(i,0))==m2.end())
	  continue;
	gate+=kron(m1[ukey<2>(lend,i)],m2[ukey<2>(i,0)]*0.5);
      }
      else if(i!=0&&i!=(m1.Dr-1)){
	if(m1.find(ukey<2>(lend,i))==m1.end())
	  continue;
	if(m2.find(ukey<2>(i,0))==m2.end())
	  continue;
	gate+=kron(m1[ukey<2>(lend,i)],m2[ukey<2>(i,0)]);
      }
  }
  // cout<<"$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"<<endl;
  // cout<<"that's the two site gate before taking the exponential"<<endl;
  // cout<<gate<<endl;
  // cout<<"$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"<<endl;
  gate.squeeze(1e-14);
  matrix<Complex,column_major> H(gate.size1(),gate.size2()),tmp;
  boost_setTo(H,Complex(0.0));
  typename SparseOperator<2,T>::iterator it;
  for(it=gate.begin();it!=gate.end();++it){
    std::vector<uint>out = it->first.getOut(),in=it->first.getIn();
    uint x = MuToSiIndex(out,gate.dims1());
    uint y = MuToSiIndex(in,gate.dims2());
    H(x,y)=Complex(it->second);
  }

  VectorType E;
  CMatrixTypeColMaj DiagE,tmp2;
  BoostDiag(H,E);
  DiagE.resize(E.size(),E.size());
  boost_setTo(DiagE,Complex(0.0));
  //for(uint ind = 0;ind<E.size();ind++)DiagE(ind,ind) = exp(-Complex(0.0,E(ind)*dt));
  for(uint ind = 0;ind<E.size();ind++)DiagE(ind,ind) = exp(-Complex(imag(dt),real(dt))*E(ind));		//NUSS 09 05 2012
  
  MatMatProd(T(1.0),H,DiagE,tmp,'n','n');
  MatMatProd(T(1.0),tmp,herm(H),tmp2,'n','n');

  SparseOperator<2,Complex>gateC(gate.size1(),gate.size2());
  gateC.setDims1(gate.dims1());
  gateC.setDims2(gate.dims2());
  for(uint x=0;x<gate.size1();++x){
    std::vector<uint>out = SiToMuIndex(x,gate.dims1());
    for(uint y=0;y<gate.size2();++y){
      std::vector<uint>in = SiToMuIndex(y,gate.dims2());
      gateC[OpKeyType<2>(out,in)]=tmp2(x,y);
    }
  }
  gateC.squeeze(1e-15);//remove everything below machine precision (1e-16 is too small, cause sometimes, the svd gives back entries> 1e-16 which actually 
                       //should be 0
  return gateC;
}

template<typename T>
MPO<Real> MPO<T>::computeNNGateMPOReal(const int even,const MPS<AbKeyType,Real>&mps,const Real dt){  
  MPO<Real> GateMPO;
  SparseOperator<2,Real> gate;
  SparseLocalOperator<Real> id1(mps.getItoKey(0).size(),mps.getItoKey(0).size());
  SparseLocalOperator<Real> idN(mps.getItoKey(mps.size()-1).size(),mps.getItoKey(mps.size()-1).size());
  id1.setSite(0);
  idN.setSite(mps.size()-1);
  for(uint index=0;index<mps.getItoKey(0).size();++index)
    id1[OpKeyType<1>(index,index)]=1.0;
  for(uint index=0;index<mps.getItoKey(mps.size()-1).size();++index)
    idN[OpKeyType<1>(index,index)]=1.0;
  for(uint site=0;site<(mps.size()-1);++site){

    if(!(site%2)&&even<0)
      continue;
    if(site%2&&even>=0)
      continue;
    gate=this->getNNGateReal(mps,site,site+1,dt);
    std::vector<uint>dims1(2),dims2(2);
    dims1[0]=gate.dims1()[0];dims1[1]=gate.dims2()[0];
    dims2[0]=gate.dims1()[1];dims2[1]=gate.dims2()[1];

    matrix<Real,column_major> H(dims1[0]*dims1[1],dims2[0]*dims2[1]),tmp,U,D,Vd;
    boost_setTo(H,0.0);
    typename SparseOperator<2,Real>::iterator it;
    for(it=gate.begin();it!=gate.end();++it){
      std::vector<uint>out(2),in(2);
      out[0]=it->first.getOut()[0];out[1]=it->first.getIn()[0];
      in[0]=it->first.getOut()[1];in[1]=it->first.getIn()[1];

      uint x = MuToSiIndex(out,dims1);
      uint y = MuToSiIndex(in,dims2);
      H(x,y)=it->second;
    }
    // cout<<"at site "<<site<<endl;
    // printboostmatrix(H);
    // cout<<endl;
    // cout<<gate<<endl;
    // my_pause();
    BoostSVD(H,U,D,Vd);
    tmp=prec_prod(D,Vd);
    // printboostmatrix(U);
    // cout<<endl;
    // printboostmatrix(D);
    // cout<<endl;
    // printboostmatrix(Vd);

    Vd=tmp;
    // cout<<endl;
    // printboostmatrix(Vd);
    // my_pause();
    // MatMatProd(T(1.0),H,DiagE,tmp,'n','n');
    // MatMatProd(T(1.0),tmp,herm(H),tmp2,'n','n');
    MPOMat<1,Real> matl,matr;
    for(uint ind=0;ind<U.size2();++ind){
      SparseLocalOperator<Real>Opl(dims1[0],dims2[0]);//,Opr(dims1[1],dims2[1]);
      Opl.setSite(site);
      for(uint x=0;x<U.size1();++x){
	std::vector<uint>left = SiToMuIndex(x,dims1);
	Opl[OpKeyType<1>(left[0],left[1])]=U(x,ind);
      }
      matl[ukey<2>(0,ind)]=Opl;
    }
    matl.Dl=1;matl.Dr=U.size2();
    matl.squeeze(1e-14);

    for(uint ind=0;ind<Vd.size1();++ind){
      SparseLocalOperator<Real>Opr(dims1[1],dims2[1]);
      Opr.setSite(site+1);
      for(uint x=0;x<Vd.size2();++x){
	std::vector<uint>right = SiToMuIndex(x,dims2);
	Opr[OpKeyType<1>(right[0],right[1])]=Vd(ind,x);
      }
      matr[ukey<2>(ind,0)]=Opr;
    }
    matr.Dl=Vd.size2();matr.Dr=1;
    matr.squeeze(1e-14);
    GateMPO[site]=matl;
    GateMPO[site+1]=matr;
    // SparseOperator<2,Real> gate1;
    // for(uint i=0;i<matl.Dr;++i)
    //   gate1+=kron(matl[ukey<2>(0,i)],matr[ukey<2>(i,0)]);
    // gate1.squeeze(1e-14);
    // cout<<"============================================================================"<<endl;
    // cout<<"a two site gate at sites "<<site<<" and "<<site+1<<" after SVD"<<endl;
    // cout<<gate1<<endl;
    // cout<<"a two site gate at sites "<<site<<" and "<<site+1<<endl;
    // cout<<gate<<endl;

    // cout<<"at site "<<site<<endl;
    // cout<<matl<<endl;
    // cout<<matr<<endl;
    // cout<<"============================================================================"<<endl;
    // my_pause();

  }
  MPOMat<1,Real> mat;
  mat.clear();
  if(even>=0&&(mps.size()%2)){
    mat[ukey<2>(0,0)]=idN;
    mat.Dr=1;mat.Dl=1;
    GateMPO[mps.size()-1]=mat;
  }
  if(even<0){
    mat[ukey<2>(0,0)]=id1;
    mat.Dr=1;mat.Dl=1;
    GateMPO[0]=mat;
    mat.clear();
    if(!(mps.size()%2)){
      mat[ukey<2>(0,0)]=idN;
      mat.Dr=1;mat.Dl=1;
      GateMPO[mps.size()-1]=mat;
    }
  }
  return GateMPO;
}


//computes exp(-tau H), with tau real.
template<typename T>
SparseOperator<2,Real> MPO<T>::getNNGateReal(const MPS<AbKeyType,T>&mps,const int site1,const int site2,const Real dt){

  if(site2!=site1+1){cout<<"MPO<T>::getNNGate(): assertion site2=site1+1 failed"<<endl;abort();}
  MPOMat<1,T> m1=(*this)[site1],m2=(*this)[site2];
  i_to_key<AbKeyType>itokeyl = mps.getItoKey(site1),itokeyr = mps.getItoKey(site2);
  SparseOperator<2,T>gate(m1.begin()->second.size1()*m2.begin()->second.size1(),m1.begin()->second.size2()*m2.begin()->second.size2());
  assert(m1.Dr==m2.Dl);
  const uint lend=m1.Dl-1;
  for(uint i = 0;i<m1.Dr;++i){
    if(site1==0&&site2!=(this->size()-1)){
      if(i==0){
	if(m1.find(ukey<2>(0,i))==m1.end())
	  continue;
	if(m2.find(ukey<2>(i,0))==m2.end())
	  continue;
	gate+=kron(m1[ukey<2>(0,i)],m2[ukey<2>(i,0)]);
      }
      else if(i==(m1.Dr-1)){
	if(m1.find(ukey<2>(0,i))==m1.end())
	  continue;
	if(m2.find(ukey<2>(i,0))==m2.end())
	  continue;
	gate+=kron(m1[ukey<2>(0,i)],m2[ukey<2>(i,0)]*0.5);
      }
      else if(i!=0&&i!=(m1.Dr-1)){
	if(m1.find(ukey<2>(0,i))==m1.end())
	  continue;
	if(m2.find(ukey<2>(i,0))==m2.end())
	  continue;
      	gate+=kron(m1[ukey<2>(0,i)],m2[ukey<2>(i,0)]);
      }
    }
    if((site2==(this->size()-1))&&site1!=0){
      if(i==0){
	if(m1.find(ukey<2>(lend,i))==m1.end())
	  continue;
	if(m2.find(ukey<2>(i,0))==m2.end())
	  continue;
      	gate+=kron(m1[ukey<2>(lend,i)]*0.5,m2[ukey<2>(i,0)]);
      }
      else if(i==(m1.Dr-1)){
	if(m1.find(ukey<2>(lend,i))==m1.end())
	  continue;
	if(m2.find(ukey<2>(i,0))==m2.end())
	  continue;
	gate+=kron(m1[ukey<2>(lend,i)],m2[ukey<2>(i,0)]);
      }
      else if(i!=(m1.Dr-1)){
	if(m1.find(ukey<2>(lend,i))==m1.end())
	  continue;
	if(m2.find(ukey<2>(i,0))==m2.end())
	  continue;
	gate+=kron(m1[ukey<2>(lend,i)],m2[ukey<2>(i,0)]);
      }
    }
    if((site2!=(this->size()-1))&&site1!=0)
      if(i==0){
	if(m1.find(ukey<2>(lend,i))==m1.end())
	  continue;
	if(m2.find(ukey<2>(i,0))==m2.end())
	  continue;
	gate+=kron(m1[ukey<2>(lend,i)]*0.5,m2[ukey<2>(i,0)]);
      }
      else if(i==(m1.Dr-1)){
	if(m1.find(ukey<2>(lend,i))==m1.end())
	  continue;
	if(m2.find(ukey<2>(i,0))==m2.end())
	  continue;
	gate+=kron(m1[ukey<2>(lend,i)],m2[ukey<2>(i,0)]*0.5);
      }
      else if(i!=0&&i!=(m1.Dr-1)){
	if(m1.find(ukey<2>(lend,i))==m1.end())
	  continue;
	if(m2.find(ukey<2>(i,0))==m2.end())
	  continue;
	gate+=kron(m1[ukey<2>(lend,i)],m2[ukey<2>(i,0)]);
      }
  }
  gate.squeeze(1e-14);
  // cout<<"$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"<<endl;
  // cout<<"that's the two site gate before taking the exponential"<<endl;
  // cout<<gate<<endl;
  // cout<<"$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"<<endl;
  matrix<Real,column_major> H(gate.size1(),gate.size2()),tmp;
  boost_setTo(H,0.0);
  typename SparseOperator<2,T>::iterator it;
  for(it=gate.begin();it!=gate.end();++it){
    std::vector<uint>out = it->first.getOut(),in=it->first.getIn();
    uint x = MuToSiIndex(out,gate.dims1());
    uint y = MuToSiIndex(in,gate.dims2());
    H(x,y)=it->second;
  }

  VectorType E;
  MatrixTypeColMaj DiagE,tmp2;
  BoostDiag(H,E);
  DiagE.resize(E.size(),E.size());
  boost_setTo(DiagE,0.0);
  //for(uint ind = 0;ind<E.size();ind++)DiagE(ind,ind) = exp(-Complex(0.0,E(ind)*dt));
  for(uint ind = 0;ind<E.size();++ind)DiagE(ind,ind) = exp(-dt*E(ind));		//NUSS 09 05 2012
  
  MatMatProd(T(1.0),H,DiagE,tmp,'n','n');
  MatMatProd(T(1.0),tmp,herm(H),tmp2,'n','n');

  gate.clear();
  for(uint x=0;x<gate.size1();++x){
    std::vector<uint>out = SiToMuIndex(x,gate.dims1());
    for(uint y=0;y<gate.size2();++y){
      std::vector<uint>in = SiToMuIndex(y,gate.dims2());
      gate[OpKeyType<2>(out,in)]=tmp2(x,y);
    }
  }
  gate.squeeze(1e-15);//remove everything below machine precision (1e-16 is too small, cause sometimes, the svd gives back entries> 1e-16 which actually 
                       //should be 0
  return gate;
}


//finds all next nearest neighbor gates which act on sites seperated by ONE!!! It assumes that the first entry of the MPO matrix, e.g. the element (0,0), is 
//the identity matrix!! If that's not the case, results are plain wrong!!!
//all local contributions, like chemical potentials, will be dismissed in the NNNgates
template<typename T>
SparseOperator<3,Complex> MPO<T>::getNNNGate(const MPS<AbKeyType,T>&mps,const int site1,const int site2,const int site3,const Complex dt){		//NUSS 09 05 2012 dt to Complex

  if(site2!=site1+1){cout<<"MPO<T>::getNNNGate(): assertion site2=site1+1 failed"<<endl;abort();}
  if(site3!=site2+1){cout<<"MPO<T>::getNNNGate(): assertion site3=site2+1 failed"<<endl;abort();}
  MPOMat<1,T> m1=(*this)[site1],m2=(*this)[site2],m3=(*this)[site3];
  i_to_key<AbKeyType>itokeyl = mps.getItoKey(site1),itokeyc = mps.getItoKey(site2),itokeyr = mps.getItoKey(site3);
  SparseOperator<3,T>gate(m1.begin()->second.size1()*m2.begin()->second.size1()*m3.begin()->second.size1(),
			  m1.begin()->second.size2()*m2.begin()->second.size2()*m3.begin()->second.size2());
  SparseOperator<2,T>tempGate(m3.begin()->second.size1()*m2.begin()->second.size1(),m3.begin()->second.size2()*m2.begin()->second.size2());
  assert(m1.Dr==m2.Dl);
  assert(m2.Dr==m3.Dl);
  typename MPOMat<1,T>::iterator mit;
  for(uint i=1;i<m2.Dl-1;++i){
    for(uint j=1;j<m2.Dr-1;++j){
      mit=m2.find(ukey<2>(i,j));
      if(mit==m2.end())
	continue;
      SparseOperator<1,T>opc=mit->second;
      mit=m1.find(ukey<2>(m1.Dl-1,i));
      if(mit==m1.end())
	continue;
      SparseOperator<1,T>opl=mit->second;
      mit=m3.find(ukey<2>(j,0));
      if(mit==m3.end())
	continue;
      SparseOperator<1,T>opr=mit->second;
      tempGate=kron(opc,opr);
      gate+=kron(opl,tempGate);
    }
  }
  if(gate.size()==0){
    cout<<"MPO contains no next-nearest-neighbor interactions; identity-gate will be returned for at sites "<<site1<<","<<site2<<","<<site3<<endl;
  }
  //the dims1 and dims2 of the gate have to be set because the second kron above assumes that opl has dims1 and dims2 both of length 0; 
  //that's because _dims1 and _dims1 members of SparseOperator are not initialized by the constructor (maybe that should be done ..)
  std::vector<uint> dims(3);
  dims[0]=m1.begin()->second.size1();
  dims[1]=m2.begin()->second.size1();
  dims[2]=m3.begin()->second.size1();
  gate.setDims1(dims);
  gate.setDims2(dims);
  matrix<Complex,column_major> H(gate.size1(),gate.size2()),tmp;
  boost_setTo(H,Complex(0.0));
  typename SparseOperator<3,T>::iterator it;
  for(it=gate.begin();it!=gate.end();++it){
    std::vector<uint>out = it->first.getOut(),in=it->first.getIn();
    uint x = MuToSiIndex(out,gate.dims1());
    uint y = MuToSiIndex(in,gate.dims2());
    H(x,y)=Complex(it->second);
  }
  VectorType E;
  CMatrixTypeColMaj DiagE,tmp2;
  BoostDiag(H,E);
  DiagE.resize(E.size(),E.size());
  boost_setTo(DiagE,Complex(0.0));
  //for(uint ind = 0;ind<E.size();ind++)DiagE(ind,ind) = exp(-Complex(0.0,E(ind)*dt));
  for(uint ind = 0;ind<E.size();ind++)DiagE(ind,ind) = exp(-Complex(imag(dt),real(dt))*E(ind));		//NUSS 09 05 2012
  MatMatProd(T(1.0),H,DiagE,tmp,'n','n');
  MatMatProd(T(1.0),tmp,herm(H),tmp2,'n','n');

  SparseOperator<3,Complex>gateC(gate.size1(),gate.size2());
  gateC.setDims1(gate.dims1());
  gateC.setDims2(gate.dims2());
  for(uint x=0;x<gate.size1();++x){
    std::vector<uint>out = SiToMuIndex(x,gate.dims1());
    for(uint y=0;y<gate.size2();++y){
      std::vector<uint>in = SiToMuIndex(y,gate.dims2());
      gateC[OpKeyType<3>(out,in)]=tmp2(x,y);
    }
  }
  gateC.squeeze(1e-15);//remove everything below machine precision (1e-16 is too small, cause sometimes, the diagonalization gives back entries > 1e-16 which actually 
                       //should be 0
  return gateC;
}


// template<class KeyType,typename T1,typename T2>
// MPO<KeyType,T1> cast(const MPO<KeyType,T2>&mpo){
//   MPO<KeyType,T1> mpo_cast;
//   mpo_cast.resize(mpo.size());
//   for(uint site=0;site<mpo.size();++site){
//     mpo_cast[site] = cast<1,KeyType,T1,T2>(mpo[site]);
//   }
//   return mpo_cast;
// }


template<typename T>
class MPO2Site:public std::multimap<int,MPOMat<1,T> >,public Operator<T>{
public:
  typedef typename std::multimap<int,MPOMat<1,T> >::iterator iterator;
  typedef typename std::multimap<int,MPOMat<1,T> >::const_iterator const_iterator;
  MPO2Site(){};
  MPO2Site(const SparseLocalOperator<T> &o1,const SparseLocalOperator<T> &o2);
  MPO2Site(const MPO2Site<T>&mpo):std::multimap<int,MPOMat<1,T> >(mpo){};
  ~MPO2Site(){};

  virtual T measure(const MPS<AbKeyType,T>&mps)const{return T(-1000);};
  virtual T measure(const CanonizedMPS<AbKeyType,T>&mps)const ;
};

template<typename T>
MPO2Site<T>::MPO2Site(const SparseLocalOperator<T> &o1,const SparseLocalOperator<T> &o2){
  MPOMat<1,T> m1,m2;
  m1[ukey<2>(0,0)] = o1;  m2[ukey<2>(0,0)] = o2;
  int site1 = o1.site();
  int site2 = o2.site();
  assert(site1!=site2);
  this->insert(pair<int,MPOMat<1,T> >(site1,m1));
  this->insert(pair<int,MPOMat<1,T> >(site2,m2));
}

template<typename T>
T MPO2Site<T>::measure(const CanonizedMPS<AbKeyType,T>&mps)const{

  if(!this->size()){
    cout<<"MPO2Site::measure(const MPS<AbKeyType>): detected empty map"<<endl;
    return T(-1000);
  }
  if(!mps.size()){
    cout<<"MPO2Site::measure(const MPS<AbKeyType>): detected size of CanonizedMPS<AbKeyType,T> = 0; skipping measurement"<<endl;
    return T(-1000);
  }
  MPO2Site<T>cpy=(*this);
  map_container<int,BlockMatrix<1,AbKeyType,T> >L;
  typename BlockMatrix<1,AbKeyType,T>::iterator itL;
  iterator it;
  it= cpy.begin();
  int site1 = it->first; 
  ++it;
  int site2 = it->first;
  T val = T(0.0);
  if(site1!=site2){
    L[-1] = getBoundaryBL(mps[site1],0);
    for(int site = 0;site<=(site2-site1);++site){
      SparseLocalOperator<T> localId(mps.getPSpace()(site1+site),mps.getPSpace()(site1+site));
      for(uint ind = 0;ind<mps.getPSpace()(site1+site);++ind){
	localId[OpKeyType<1>(ind,ind)] = T(1.0);      
      }
      MPOMat<1,T> id;
      id[ukey<2>(0,0)] = localId;
      BlockMatrix<1,AbKeyType,T> A;
      if(site==(site2-site1)||site==0)
	A = mps.lambda(site1+site-1)*mps[site1+site]*mps.lambda(site1+site);
      else if(site!=(site2-site2)&&site!=0)
	A = mps[site1+site]*mps.lambda(site1+site);
      if(site == 0)
	L[site]=addContractionLayer(L[site-1],A,cpy.find(site1+site)->second,mps.getItoKey(site1+site),1);
      else if((site!=0)&&site!=(site2-site1))
	L[site]=addContractionLayer(L[site-1],A,id,mps.getItoKey(site1+site),1);
      else if(site==(site2-site1))
	L[site]=addContractionLayer(L[site-1],A,cpy.find(site1+site)->second,mps.getItoKey(site1+site),1);
    }
    for(itL=L[site2-site1].begin();itL!=L[site2-site1].end();++itL){
      val+=trace(itL->second);
    }
  }
  else{
    it= cpy.begin();
    SparseLocalOperator<T> o1(it->second[ukey<2>(0,0)],site1);
    ++it;
    SparseLocalOperator<T> o2(it->second[ukey<2>(0,0)],site2);
    SparseLocalOperator<T>prod = o1*o2;
    prod.setSite(site1);
    val=prod.measure(mps);
  }
  return val;
 }

class CharacteristicFunction:public MPO<std::complex<double> >{
public:
	CharacteristicFunction(int dim, double lambda, int L, int LA){
		SparseLocalOperator<std::complex<double> > op(dim,dim);
		std::complex<double> il(0,lambda);
		for (int n=0; n<dim; ++n){
			op[OpKeyType<1>(n,n)]=exp(il*(double)n);
		}
		for (int i=0; i<LA; ++i){
			MPOMat<1,std::complex<double> > M;
			op.setSite(i);
			M[ukey<2>(0,0)] = op;
			M.Dl=1,M.Dr=1;
			M.squeeze(1e-16);
			(*this)[op.site()]=M;
		}
	}

	~CharacteristicFunction(){}
};

template<typename T>
class Correlator:public MPO<T>{
public:
	Correlator(SparseLocalOperator<T> &op1, SparseLocalOperator<T> &op2){
		MPOMat<1,T> M;
		M[ukey<2>(0,0)] = op1;
		M.Dl=1,M.Dr=1;
		M.squeeze(1e-16);
		(*this)[op1.site()]=M;
		M.clear();

		M[ukey<2>(0,0)] = op2;
		M.Dl=1,M.Dr=1;
		M.squeeze(1e-16);
		(*this)[op2.site()]=M;
	}

	~Correlator(){}
};


template<typename T>
class Overlap:public Operator<T>{
public:
  Overlap(){};
  Overlap(const MPS<AbKeyType,T>&mps);

  virtual T measure(const MPS<AbKeyType,T>&mps)const;
  virtual T measure(const CanonizedMPS<AbKeyType,T>&mps)const;
private:
  MPS<AbKeyType,T>_mps;
  MPO<T> _id;

};

template<typename T>
Overlap<T>::Overlap(const MPS<AbKeyType,T>&mps){
  _mps=mps;
  for(uint site = 0;site<_mps.size();++site){
    SparseOperator<1,T> id(_mps.getPSpace()(site),_mps.getPSpace()(site));
    for(uint ind=0;ind<_mps.getPSpace()(site);++ind)
      id[OpKeyType<1>(ind,ind)]= T(1.);
    _id[site][ukey<2>(0,0)]=id;
  }

}

template<typename T>
T Overlap<T>::measure(const MPS<AbKeyType,T>&mps)const{
  if(mps.size()!=_mps.size()){
    cerr<<"Overlap::measure: wrong size of input mps"<<endl;
    return -1000;
  }
  return computeOverlap(_mps,mps);
}


//computes Overlap of a CanonizedMPS with a usual MPS
template<typename T>
T Overlap<T>::measure(const CanonizedMPS<AbKeyType,T>&mps)const{
  if(mps.size()!=_mps.size()){
    cerr<<"Overlap::measure: wrong size of input mps"<<endl;
    return -1000;
  }
  BlockMatrix<1,AbKeyType,T>L;
  L= getLMin(_mps);
  for(uint site=0;site<mps.size();++site){
    BlockMatrix<1,AbKeyType,T> temp=mps.lambda(site-1)*mps[site];
    L=addContractionLayer(L,_mps[site],_id.find(site)->second,temp,_mps.getItoKey(site),1);
  }
  return L.size()>0?L.begin()->second(0,0):T(0.0);
}




//calculates the Overlap of input mps with its time reversed version
template<typename T>
class TRIOverlap:public Operator<T>{
public:
  TRIOverlap(){};
  virtual T measure(const CanonizedMPS<AbKeyType,T>&mps)const;
  virtual T measure(const MPS<AbKeyType,T>&mps)const;
private:

};



template<typename T>
T TRIOverlap<T>::measure(const MPS<AbKeyType,T>&mps)const {

  BlockMatrix<1,AbKeyType,T>L;
  L= getLMin(mps);
  for(uint site=0;site<mps.size();++site){
    MPOMat<1,T> idMPO;
    SparseOperator<1,T> id(mps.getPSpace()(site),mps.getPSpace()(site));
    for(uint ind=0;ind<mps.getPSpace()(site);++ind)
      id[OpKeyType<1>(ind,ind)]= T(1.);
    idMPO[ukey<2>(0,0)]=id;
    BlockMatrix<1,AbKeyType,T>mat=mps[site];
    L=addContractionLayer(L,mps[site],idMPO,mat.conjugate(),mps.getItoKey(site),1);
  }
  return L.size()>0?L.begin()->second(0,0):T(0.0);

}

template<typename T>
T TRIOverlap<T>::measure(const CanonizedMPS<AbKeyType,T>&cmps)const{
  MPS<AbKeyType,T> _mps=cmps.toMPS();
  _mps.conjugate();
  BlockMatrix<1,AbKeyType,T>L;
  L= getLMin(_mps);
  for(uint site=0;site<cmps.size();++site){
    MPOMat<1,T> idMPO;
    BlockMatrix<1,AbKeyType,T> temp=cmps.lambda(site-1)*cmps[site];
    SparseOperator<1,T> id(_mps.getPSpace()(site),_mps.getPSpace()(site));
    for(uint ind=0;ind<_mps.getPSpace()(site);++ind)
      id[OpKeyType<1>(ind,ind)]= T(1.);
    idMPO[ukey<2>(0,0)]=id;
    L=addContractionLayer(L,_mps[site],idMPO,temp,_mps.getItoKey(site),1);
  }

  return L.size()>0?L.begin()->second(0,0):T(0.0);
}


template<typename T>
class OverlapWithStored:public Operator<T>{
public:
  OverlapWithStored(std::string fname,const int maxstep){
    _fname=fname;
    _maxstep=maxstep;
    std::fstream outfile;
    outfile.open("lasterased.dat",ios_base::out);
    outfile<<-1;
    outfile.close();
  }

  virtual T measure(const MPS<AbKeyType,T>&mps)const{return T(-1000);}
  virtual T measure(const CanonizedMPS<AbKeyType,T>&mps)const;

private:
  int _maxstep;
  std::string _fname;
};

//computes Overlap of a CanonizedMPS with a usual MPS
//works only with abelian keys
template<typename T>
T OverlapWithStored<T>::measure(const CanonizedMPS<AbKeyType,T>&mps)const{
  //first read in the correct data file; the correct one is the one with the smallest index number; 
  //the index number is computed by incrementing it, starting from zero. The routine has to check for each possible
  //data filename wether it exists. after measuring, the data-file HAS TO BE DELETED!!!
  //the files produces an output file "lasterased.dat" which contains the index of the last erased .mps file; DO NOT DELETE IT!!!!
  //if the routine can't find an .mps file with a certain index and assumes that there will appear one in the future it will wait for this file;
  //it knows that there will appear one if the Numfiles argument is is larger than the index in lasterased but there is no .mps file to load 
  //anymore.

  std::fstream file;
  file.open("lasterased.dat",ios_base::in);
  int lasterased;
  file>>lasterased;
  file.close();
  CanonizedMPS<AbKeyType,T>cmps;
  MPS<AbKeyType,T>localmps;
  MPO<T> _id;
  std::ostringstream convert;
  int n=lasterased;
  convert<<n;
  std::string name=_fname+"_CP"+convert.str()+"c";
  while(!fexist(name+".mps")&&(n<_maxstep)){
    convert.str("");//clear the ostringstream object
    n++;
    convert<<n;
    name=_fname+"_CP"+convert.str()+"c";
  }
  if(!fexist(name+".mps")&&(n==_maxstep)){
    //there wasn't any file; either the routine itself has already erased it, or it hasn't been produced yet:
    if(lasterased==_maxstep){
      //in fact, the routine should never be here...
      cout<<"#routine has already reached the iteration "<<_maxstep<<" in a previous iteration and erased the .mps"<<endl;
      return T(-2000);
    }else if(lasterased<_maxstep){
      //now we have to wait for the other routine to produce the file;
      n=lasterased+1;
      convert.str("");//clear the ostringstream object
      convert<<n;
      name=_fname+"_CP"+convert.str()+"c";
      while(!fexist(name+".mps")){
	cout<<"#cannot find "<<name<<".mps. Will try again in 60 seconds) ..."<<endl;
	cout<<"#please check that Numfiles is not larger than the largest iteration index of the .mps-production run"<<endl;
	sleep(60);
      }
      cout<<"#found file "<<name<<".mps; Continuing simulation"<<endl;
    }
    
  }
  //now check: _maxstep is the highest possible index. if n=_maxstep and last erased!= 

  
  //now name contains the filename of an existing file; load it!!
  cmps.load(name);
  localmps=cmps.toMPS();
  cmps.clear();
  if(mps.size()!=localmps.size()){
    cerr<<"OverlapWithStored::measure: wrong size of input mps"<<endl;
    return -1000;
  }
  for(uint site = 0;site<localmps.size();++site){
    SparseOperator<1,T> id(localmps.getPSpace()(site),localmps.getPSpace()(site));
    for(uint ind=0;ind<localmps.getPSpace()(site);++ind)
      id[OpKeyType<1>(ind,ind)]= T(1.);
    _id[site][ukey<2>(0,0)]=id;
  }

  BlockMatrix<1,AbKeyType,T>L;
  L= getLMin(localmps);
  for(uint site=0;site<mps.size();++site){
    BlockMatrix<1,AbKeyType,T> temp=mps.lambda(site-1)*mps[site];
    L=addContractionLayer(L,localmps[site],_id.find(site)->second,temp,localmps.getItoKey(site),1);
  }

  //delete the file!!!!
  std::string tempstring=name+".mps";
  remove(tempstring.c_str());

  file.open("lasterased.dat",ios_base::out);
  file<<n;
  file.close();

  return L.size()>0?L.begin()->second(0,0):T(0.0);
}

template<typename T>
class Norm:public Operator<T>{
public:
  Norm(){};
  virtual T measure(const MPS<AbKeyType,T>&mps)const;
  virtual T measure(const CanonizedMPS<AbKeyType,T>&mps)const;
};


template<typename T>
T Norm<T>::measure(const MPS<AbKeyType,T>&mps)const{
  return sqrt(computeOverlap(mps,mps));
}


template<typename T>
//there are problems using const here, because operator[] is not const for maps, and toMPS() need that operator
T Norm<T>::measure(const CanonizedMPS<AbKeyType,T>&cmps)const{
  CanonizedMPS<AbKeyType,T> _cmps=cmps;
  MPS<AbKeyType,T> _mps=_cmps.toMPS();
  return sqrt(computeOverlap(_mps,_mps));
}


#endif
