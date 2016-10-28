/*
 * BlockMatrix.h
 *
 *  Created on: 09.03.2011
 *      Author: martinjg
 */

#ifndef BLOCKTENSOR_H_
#define BLOCKTENSOR_H_

#include "../boost_typedefs.hpp"
#include <map>
#include <iostream>
#include "keyTypes.hpp"
#include "tensormap.hpp"
#include "utilities.hpp"
#include "../blas/blasroutines.hpp"
#include <random>
#include <chrono>

using namespace BoostTypedefs;
using boost::numeric::ublas::matrix;
using boost::numeric::ublas::column_major;
using boost::numeric::ublas::matrix_range;

template<typename KeyType,typename T>
class DiagBlockMatrix:public map_container<KeyType,boost::numeric::ublas::vector<T> > {
public:
  typedef typename map_container<KeyType,boost::numeric::ublas::vector<T> >::iterator iterator;
  typedef typename map_container<KeyType,boost::numeric::ublas::vector<T> >::const_iterator const_iterator;
  DiagBlockMatrix(){};
  virtual ~DiagBlockMatrix(){};
  Real normalize();
  Real norm();
  DiagBlockMatrix<KeyType,T>&operator*=(const T num);
  DiagBlockMatrix<KeyType,T>&operator/=(const T num){return (*this)*=(1/num);}
  pair<T,uint> truncate(const Real tr_weight,const uint chimax,const bool normalize=true);
  pair<T,uint> truncate(const uint chimax,const bool normalize=true);

  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)    {
      ar & boost::serialization::base_object<map_container<KeyType,boost::numeric::ublas::vector<T> > >(*this);
  }
};

template<class KeyType,typename T>
Real DiagBlockMatrix<KeyType,T>::normalize(){
  iterator it;
  T norm=T(0.0);
  for(it=this->begin();it!=this->end();++it){
    boost::numeric::ublas::vector<T>d;
    d=it->second;
    norm+=dot_prod(d.size(),it->second.data().begin(),d.data().begin());
  }
  (*this)/=sqrt(norm);
  return sqrt(norm);
}


//returns the sum of the squares of the diagonal elements;
template<class KeyType,typename T>
Real DiagBlockMatrix<KeyType,T>::norm(){
  iterator it;
  T norm=T(0.0);
  for(it=this->begin();it!=this->end();++it){
    boost::numeric::ublas::vector<T>d;
    d=it->second;
    norm+=dot_prod(d.size(),it->second.data().begin(),d.data().begin());
  }
  return sqrt(norm);
}

template<class KeyType,typename T>
DiagBlockMatrix<KeyType,T>&DiagBlockMatrix<KeyType,T>::operator*=(const T num){
  iterator it;
  for(it=this->begin();it!=this->end();++it){
    it->second*=num;
  }
  return (*this);
}



//truncates the block-matrix; hard limit is chi; if tr_weight can be kept at smaller chi, then the smaller chi is used
template<class KeyType,typename T>
pair<T,uint> DiagBlockMatrix<KeyType,T>::truncate(const Real tr_weight,const uint chimax,const bool normalize){
  std::multimap<T,KeyType>m;
  std::multimap<KeyType,T>im;
  iterator it;
  T _tw=T(this->norm());
  T norm=T(this->norm());
  //T _tw=T(this->norm());//There is a problem when using imaginary time; the truncated weight gets negative
  //fill all entries into a map
  // if(abs(this->norm()-1)>1e-4){
  //   cout<<"got a lambda with norm "<<this->norm()<<endl;
  //   cin.get();
  // }
  typename boost::numeric::ublas::vector<T>::iterator i;
  for(it=this->begin();it!=this->end();++it){
    for(i=it->second.begin();i!=it->second.end();++i)
      m.insert(pair<T,KeyType>(*i,it->first));
  }
  
  uint chi=0;
  typename multimap<T,KeyType>::reverse_iterator r=m.rbegin();;
  while((r!=m.rend())&&(fabs(_tw)>norm*tr_weight)&&chi<chimax){
    _tw-=((r->first*r->first)/norm);
    im.insert(pair<KeyType,T>(r->second,r->first));
    ++r;
    ++chi;
  }
  if(_tw<0&&fabs(_tw)<1e-14)
    _tw=0.0;
  else if(_tw<0&&fabs(_tw)>1e-14){
    cout<<"warning: in DiagBlockMarix<KeyType,T>::truncate(const Real, const uint, const bool): negative _tw ("<<_tw<<") detected!"<<endl;
  }

  this->clear();
  typename std::multimap<KeyType,T>::iterator mit;
  while(im.size()!=0){
    boost::numeric::ublas::vector<T> v(im.count(im.begin()->first));
    pair<typename std::multimap<KeyType,T>::iterator,typename std::multimap<KeyType,T>::iterator>  p ;
    p  = im.equal_range(im.begin()->first);
    uint k=0;
    for(mit=p.first;mit!=p.second;++mit){
      v(k)=mit->second;
      ++k;
    }
    (*this)[im.begin()->first]=v;
    im.erase(p.first,p.second);
  }

  if(normalize)
    this->normalize();
  return pair<T,uint>(_tw,chi);
}

//this is a different truncation function, but it's not save (it crashes from time to time) ...
// template<class KeyType,typename T>
// uint DiagBlockMatrix<KeyType,T>::truncate(const uint chimax){
//   this->normalize();
//   std::vector<triple<typename map<KeyType,VectorType>::iterator,VectorType::iterator,uint> > position(this->size());
//   typename std::vector<triple<typename map<KeyType,VectorType>::iterator,VectorType::iterator,uint> >::iterator position_it = position.begin();
//   pair<typename std::vector<triple<typename map<KeyType,VectorType>::iterator,VectorType::iterator,uint> >::iterator,Real>
//     largest(position.begin(),-1.0),
//     secondlargest(position.begin(),-1.0);
//   uint cnt=0;
//   uint chi = 0;
//   Real nm = 0;
//   typename DiagBlockMatrix<KeyType,T>::iterator itlam;
//   Real truncweight;
//   VectorType D;
//   Real Norm=0;
//   for(itlam = this->begin();itlam!=this->end();itlam++){
//     D = itlam->second;
//     Norm += dot_prod(D.size(),D.data().begin(),itlam->second.data().begin());
//     //		itlam->second=cut(D,1e-16);
//     chi+=itlam->second.size();
//     position[cnt] = triple<typename map<KeyType,VectorType>::iterator,VectorType::iterator,uint> (itlam,itlam->second.begin(),uint (0));
//     if(itlam->second.size()!=0&&itlam->second(0)>=largest.second){
//       largest.first = position_it;
//       largest.second = itlam->second(0);
//     }
//     cnt++;
//     position_it++;
//   }
//   cnt = 0;
//   position_it = position.begin();
//   for(position_it = position.begin();position_it != position.end();position_it++){
//     if((position_it->second!=position_it->first->second.end())&&*position_it->second>secondlargest.second&&*position_it->second<=largest.second){
//       if(position_it!=largest.first){
// 	secondlargest.first = position_it;
// 	secondlargest.second = *position_it->second;
//       }
//     }
//   }
//   cnt = 0;
//   uint cnt2=0;
//   while(cnt<chimax&&cnt<chi){
//     if(cnt2>(chimax<chi?chimax:chi)){
//       cout<<"in normalize: stuck in while loop!"<<endl;
//     }
//     while((largest.first->second!=largest.first->first->second.end())&&(*largest.first->second>=secondlargest.second)){
//       largest.first->second++;
//       largest.first->third++;
//       cnt++;
//       if(cnt>=chimax)break;
//     }
//     largest=secondlargest;
//     secondlargest.second=-1.0;
//     secondlargest.first = position.begin();
//     for(position_it = position.begin();position_it != position.end();position_it++){
//       if((position_it->second!=position_it->first->second.end())&&*position_it->second>secondlargest.second&&*position_it->second<=largest.second&&
// 	 (position_it!=largest.first)){
// 	secondlargest.first = position_it;
// 	secondlargest.second = *position_it->second;
//       }
//     }
//     cnt2++;
//   }
    
//   for(position_it = position.begin();position_it != position.end();position_it++){
//     position_it->first->second.resize(position_it->third,true);
//   }

//   std::vector<typename map<KeyType,VectorType>::iterator>erase_lam;
//   typename std::vector<typename map<KeyType,VectorType>::iterator>::iterator erase_lam_it;

//   Real newNorm = 0.0;
//   D.resize(0);
//   chi=0;
//   uint start = 0;
//   for(itlam = this->begin();itlam!=this->end();itlam++){
//     chi=itlam->second.size();
//     start+=chi;
//     if(itlam->second.size()==0){
//       erase_lam.push_back(itlam);
//     }
//     else{
//       D= itlam->second;
//       nm = dot_prod(D.size(),D.data().begin(),itlam->second.data().begin());
//       newNorm +=nm;
//     }
//   }
//   for(erase_lam_it = erase_lam.begin();erase_lam_it != erase_lam.end();erase_lam_it++){
//     this->erase(*erase_lam_it);
//   }

//   return start;
// }


//this does only a truncation and no normalization
template<class KeyType,typename T>
pair<T,uint> DiagBlockMatrix<KeyType,T>::truncate(const uint chi,const bool normalize){
  std::multimap<T,KeyType>m;
  std::multimap<KeyType,T>im;
  iterator it;
  typename boost::numeric::ublas::vector<T>::iterator i;
  for(it=this->begin();it!=this->end();++it){
    for(i=it->second.begin();i!=it->second.end();++i)
      m.insert(pair<T,KeyType>(*i,it->first));
  }

  typename multimap<T,KeyType>::reverse_iterator r=m.rbegin();
  uint cnt=0;
  T _tw=T(1.0);//T(this->norm());There is a problem when using imaginary time; the truncated weight gets negative
  // if(abs(this->norm()-1)>1e-4){
  //   cout<<"got a lambda with norm "<<this->norm()<<endl;
  //   cin.get();
  // }

  while((r!=m.rend())&&(cnt!=chi)&&(abs(r->first)>1e-15)){
    _tw-=(r->first*r->first);
    im.insert(pair<KeyType,T>(r->second,r->first));
    ++r;
    ++cnt;
  }
  this->clear();
  typename std::multimap<KeyType,T>::iterator mit;
  while(im.size()!=0){
    boost::numeric::ublas::vector<T> v(im.count(im.begin()->first));
    pair<typename std::multimap<KeyType,T>::iterator,typename std::multimap<KeyType,T>::iterator>  p ;
    p  = im.equal_range(im.begin()->first);
    uint k=0;
    for(mit=p.first;mit!=p.second;++mit){
      v(k)=mit->second;
      ++k;
    }
    (*this)[im.begin()->first]=v;
    im.erase(p.first,p.second);
  }
  if(normalize)
    this->normalize();
  return pair<T,uint>(_tw,cnt);
}


template<int Ns,class KeyType,typename T>
class BlockMatrix :public map_container<TensorKeyType<2,Ns,KeyType>,matrix<T,column_major> > {
public:
  typedef typename map_container<TensorKeyType<2,Ns,KeyType>,matrix<T,column_major> >::iterator iterator;
  typedef typename map_container<TensorKeyType<2,Ns,KeyType>,matrix<T,column_major> >::const_iterator const_iterator;
  BlockMatrix(){};
  BlockMatrix(const BlockMatrix<1,KeyType,T>&B1,const BlockMatrix<1,KeyType,T>&B2){
    assert(Ns==2);
    typename BlockMatrix<1,KeyType,T>::const_iterator it1,it2;
    for(it1 = B1.begin();it1!=B1.end();++it1){
      for(it2 = B2.begin();it2!=B2.end();++it2){
	if(it1->first[1]==it2->first[0]){
	  std::vector<uint> p(2);
	  p[0]=it1->first(0);
	  p[1]=it2->first(0);
	  TensorKeyType<2,Ns,KeyType> k(it1->first[0],it2->first[1],p);
	  matrix<T,column_major> Tmp;
	  MatMatProd(T(1.0),it1->second,it2->second,Tmp,'n','n');
	  (*this)[k]=Tmp;
	}
      }
    }
  }
  BlockMatrix(const BlockMatrix<1,KeyType,T>&B1,const BlockMatrix<1,KeyType,T>&B2,const BlockMatrix<1,KeyType,T>&B3){
    assert(Ns==3);
    typename BlockMatrix<1,KeyType,T>::const_iterator it1,it2,it3;
    for(it1 = B1.begin();it1!=B1.end();++it1){
      for(it2 = B2.begin();it2!=B2.end();++it2){
	for(it3 = B3.begin();it3!=B3.end();++it3){
	  if((it1->first[1]==it2->first[0])&&(it2->first[1]==it3->first[0])){
	    std::vector<uint> p(3);
	    p[0]=it1->first(0);
	    p[1]=it2->first(0);
	    p[2]=it3->first(0);
	    TensorKeyType<2,Ns,KeyType> k(it1->first[0],it3->first[1],p);
	    matrix<T,column_major> Tmp,Tmp2;
	    MatMatProd(T(1.0),it1->second,it2->second,Tmp,'n','n');
	    MatMatProd(T(1.0),Tmp,it3->second,Tmp2,'n','n');
	    (*this)[k]=Tmp2;
	  }
	}
      }
    }
  }
  virtual ~BlockMatrix(){};

  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)    {
      ar & boost::serialization::base_object<map_container<TensorKeyType<2,Ns,KeyType>,matrix<T,column_major> > >(*this);
  }
  
  std::set<TensorKeyType<2,Ns,KeyType> >keys()const{
    std::set<TensorKeyType<2,Ns,KeyType> > m;
    for(typename BlockMatrix<Ns,KeyType,T>::const_iterator it=this->begin();it!=this->end();++it){
      m.insert(it->first);
    }
    return m;
  }
  uint size1()const{
    map_container<KeyType,uint> sizes;
    typename map_container<KeyType,uint>::iterator s;
    for(typename BlockMatrix<1,KeyType,T>::const_iterator a=this->begin();a!=this->end();++a){
      s=sizes.find(a->first[0]);
      if(s==sizes.end())
	sizes[a->first[0]]=a->second.size2();
    }
    uint chi=0;
    for(s=sizes.begin();s!=sizes.end();++s)
      chi+=s->second;
    return chi;
  }
  uint size2()const{
    map_container<KeyType,uint> sizes;
    typename map_container<KeyType,uint>::iterator s;
    for(typename BlockMatrix<1,KeyType,T>::const_iterator a=this->begin();a!=this->end();++a){
      s=sizes.find(a->first[1]);
      if(s==sizes.end())
	sizes[a->first[1]]=a->second.size2();
    }
    uint chi=0;
    for(s=sizes.begin();s!=sizes.end();++s)
      chi+=s->second;
    return chi;
  }
  void randomize(const T a=T(0.0), const T b=T(1.0));
  void nrandomize(const T mean=T(0.0), const T var=T(1.0));
  uint squeeze(const Real thresh);
  BlockMatrix<Ns,KeyType,T> conjugate();
  //operators
  T operator *(const BlockMatrix<Ns,KeyType,T>&In) const;
  BlockMatrix<Ns,KeyType,T>& operator+=(const BlockMatrix<Ns,KeyType,T>&In);
  BlockMatrix<Ns,KeyType,T> operator +(const BlockMatrix<Ns,KeyType,T>&In) const;
  BlockMatrix<Ns,KeyType,T>& operator-=(const BlockMatrix<Ns,KeyType,T>&In);
  BlockMatrix<Ns,KeyType,T> operator -(const BlockMatrix<Ns,KeyType,T>&In) const;
  BlockMatrix<Ns,KeyType,T>& operator*=(const T num);
  //multiplies DiagBlockMatrix from the left onto BlockMatrix:
  BlockMatrix<Ns,KeyType,T>& operator*=(const DiagBlockMatrix<KeyType,Real>&lam);
  BlockMatrix<Ns,KeyType,T> operator*(const T num) const{BlockMatrix<Ns,KeyType,T> cpy=(*this);return cpy*=num;};
  BlockMatrix<Ns,KeyType,T>& operator/=(const T num){return (*this)*=(T(1.0)/num);}

};




template<int Ns,class KeyType,typename T>
BlockMatrix<Ns,KeyType,T> BlockMatrix<Ns,KeyType,T>::conjugate(){
  typename BlockMatrix<Ns,KeyType,T>::iterator b;
  for(b=this->begin();b!=this->end();++b){
    b->second=conj(b->second);
  }
  return *this;
}

template<int Ns,class KeyType,typename T>
void BlockMatrix<Ns,KeyType,T>::nrandomize(const T mean, const T var){
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::normal_distribution<double> nd(mean,var);
  std::default_random_engine generator (seed);

   iterator i;
   for(i=this->begin();i!=this->end();++i){
     for(uint n=0;n<i->second.size1()*i->second.size2();++n)
       i->second.data()[n]+=nd(generator);
   }
}

template<int Ns,class KeyType,typename T>
void BlockMatrix<Ns,KeyType,T>::randomize(const T a, const T b){
   iterator i;
   for(i=this->begin();i!=this->end();++i){
     for(uint n=0;n<i->second.size1()*i->second.size2();++n)
       i->second.data()[n]+=c_rand(a,b);
   }
}


template<int Ns,class KeyType,typename T>
uint BlockMatrix<Ns,KeyType,T>::squeeze(const Real thresh){
  uint out=0;
  iterator i;
  std::vector<iterator> erase;
  typename std::vector<iterator>::iterator e;
  
  for(i=this->begin();i!=this->end();++i){
    uint N=i->second.size1()*i->second.size2();
    Real max=0.0;
    bool er=true;
    for(uint n=0;n<N;++n){
      if(abs(i->second.data()[n])>thresh){
	er=false;
	break;
      }
    }
    if(er)
      erase.push_back(i);
  }
  out=erase.size();
  for(e = erase.begin();e!=erase.end();++e){
     this->erase(*e);
  }
  return out;
}

template<int Ns,class KeyType,typename T>
BlockMatrix<Ns,KeyType,T>& BlockMatrix<Ns,KeyType,T>::operator+=(const BlockMatrix<Ns,KeyType,T>&In){
  if(this==&In){cout<<"int BlockMatrix::operator +=: selfassignment"<<endl;abort();}
  iterator it2;
  for(const_iterator it=In.begin();it!=In.end();++it){
    it2=this->find(it->first);
    if(it2==this->end()){
      (*this)[it->first]=it->second;
    }else{
      it2->second+=it->second;
    }
    // if(!(this->count(it->first))){
    //   (*this)[it->first]=it->second;
    // }else{
    //   (*this)[it->first]+=it->second;
    // }
  }
  return *this;
}

template<int Ns,class KeyType,typename T>
BlockMatrix<Ns,KeyType,T>& BlockMatrix<Ns,KeyType,T>::operator-=(const BlockMatrix<Ns,KeyType,T>&In){
  if(this==&In){cout<<"int BlockMatrix::operator -=: selfassignment"<<endl;abort();}
  for(const_iterator it=In.begin();it!=In.end();++it){
    if(!(this->count(it->first))){
      (*this)[it->first]=T(-1.0)*it->second;
    }else{
      (*this)[it->first]-=it->second;
    }
  }
  return *this;
}


template<int Ns,class KeyType,typename T>
BlockMatrix<Ns,KeyType,T>& BlockMatrix<Ns,KeyType,T>::operator*=(const T num){
  for(iterator it=this->begin();it!=this->end();++it){
    it->second*=num;
  }
  return (*this);
}



//multiplies DiagBlockMatrix from the left onto BlockMatrix:
template<int Ns,class KeyType,typename T>
BlockMatrix<Ns,KeyType,T>& BlockMatrix<Ns,KeyType,T>::operator*=(const DiagBlockMatrix<KeyType,Real>&lam){
  typename BlockMatrix<Ns,KeyType,T>::iterator b;
  typename DiagBlockMatrix<KeyType,Real>::const_iterator l;
  BlockMatrix<Ns,KeyType,T> out;
  matrix<T,column_major>tmp;
  for(b=this->begin();b!=this->end();++b){
    l=lam.find(b->first[0]);
    if(l!=lam.end()){
      MatMatProd(T(1.0),l->second,b->second,tmp,'n','n');
      (*this)[b->first] = tmp;
    }else{
      cout<<"in operator*(const DiagBlockMatrix<KeyType,Real>&lam,const BlockMatrix<Ns,KeyType,T>&B): no lambda found for key "<<b->first[0]<<
	" of BlockMatrix"<<endl;
    }
  }
  return (*this);
}


template<int Ns,class KeyType,typename T>
BlockMatrix<Ns,KeyType,T> BlockMatrix<Ns,KeyType,T>::operator+(const BlockMatrix<Ns,KeyType,T>&In)const{
  BlockMatrix<Ns,KeyType,T>cpy = (*this);
  return cpy+=In;
}

template<int Ns,class KeyType,typename T>
BlockMatrix<Ns,KeyType,T> BlockMatrix<Ns,KeyType,T>::operator-(const BlockMatrix<Ns,KeyType,T>&In)const{
  BlockMatrix<Ns,KeyType,T>cpy = (*this);
  return cpy-=In;
}

//scalar product for BlockMatrices
template<int Ns,typename KeyType,typename T>
T  BlockMatrix<Ns,KeyType,T>::operator *(const BlockMatrix<Ns,KeyType,T>&In) const
{
  //  assert(Ns==2);
  typename BlockMatrix<Ns,KeyType,T>::iterator itthis;
  typename BlockMatrix<Ns,KeyType,T>::iterator it;
  BlockMatrix<Ns,KeyType,T> copy = *this;
  BlockMatrix<Ns,KeyType,T> m = In;
  T val = T(0.0);
  for(itthis = copy.begin();itthis!=copy.end();++itthis){
    it=m.find(itthis->first);
    if(it!=m.end()){
      uint size1 = itthis->second.size1()*itthis->second.size2();
      uint size2 = it->second.size1()*it->second.size2();
      if(size1==size2){
	val=val+dot_prod(size1,it->second.data().begin(),itthis->second.data().begin());
      }else{
	cout<<"in BlockMatrix::operator *(BlockMatrix<Ns,KeyType,T>): size1==size2 failed"<<endl;
	abort();
      }
    }
  }
  return val;
}
//==================         other operations ...         =======================



template<int Ns,class KeyType,typename T1,typename T2>
BlockMatrix<Ns,KeyType,T1> cast(const BlockMatrix<Ns,KeyType,T2>&BM){
  BlockMatrix<Ns,KeyType,T1> BM_cast;
  typename BlockMatrix<Ns,KeyType,T2>::const_iterator bm;
  for(bm=BM.begin();bm!=BM.end();++bm){
    BM_cast[bm->first] = cast<T1,T2>(bm->second);
  }
  return BM_cast;
}



//for every key on the bond connecting the two B's, it finds the largest block size and stores them in a map; needed for mps addition
template<int Ns,class KeyType,typename T>
std::map<KeyType,uint> getlargestonbond(const BlockMatrix<Ns,KeyType,T>&B1,const BlockMatrix<Ns,KeyType,T>&B2){
  typename BlockMatrix<Ns,KeyType,T>::const_iterator b;
  std::map<KeyType,uint> largest;
  typename std::map<KeyType,uint>::iterator l;
  for(b=B1.begin();b!=B1.end();++b){
    l=largest.find(b->first[1]);
    if(l==largest.end())
      largest[b->first[1]]=b->second.size2();
    else if(l!=largest.end()){
      largest[b->first[1]]=largest[b->first[1]]>b->second.size2()?largest[b->first[1]]:b->second.size2();
    }
  }
  for(b=B2.begin();b!=B2.end();++b){
    l=largest.find(b->first[0]);
    if(l==largest.end())
      largest[b->first[0]]=b->second.size1();
    else if(l!=largest.end()){
      largest[b->first[0]]=largest[b->first[0]]>b->second.size1()?largest[b->first[0]]:b->second.size1();
    }
  }
  return largest;
}


//gets the blocksizes of the BlockMatrix B, on left (side<0) or the right (side>0) side.
template<int Ns,class KeyType,typename T>
std::map<KeyType,uint> getBlockSizes(const BlockMatrix<Ns,KeyType,T>&B,const int side){
  typename BlockMatrix<Ns,KeyType,T>::const_iterator b;
  std::map<KeyType,uint> sizes;
  typename std::map<KeyType,uint>::iterator l;
  if(side>0){
    for(b=B.begin();b!=B.end();++b){
      l=sizes.find(b->first[1]);
      if(l==sizes.end())
	sizes[b->first[1]]=b->second.size2();
    }
  }
  if(side<=0){
    for(b=B.begin();b!=B.end();++b){
      l=sizes.find(b->first[0]);
      if(l==sizes.end())
	sizes[b->first[0]]=b->second.size1();
    }
  }
  return sizes;
}

template<int Ns,class KeyType,typename T>
BlockMatrix<Ns,KeyType,T> outersum(const BlockMatrix<Ns,KeyType,T>&B1,const BlockMatrix<Ns,KeyType,T>&B2){
  BlockMatrix<Ns,KeyType,T> out=B1;//that's important ...
  
  typename BlockMatrix<Ns,KeyType,T>::const_iterator b1,b2;
  std::map<KeyType,uint> sizeright=getBlockSizes(B1,1);
  std::map<KeyType,uint> sizeleft=getBlockSizes(B1,-1);
  typename std::map<KeyType,uint> ::iterator il,ir;
  for(b2=B2.begin();b2!=B2.end();++b2){
    b1=B1.find(b2->first);
    //two cases:
    //a: the block b2 is present in B1; in this case, we can just insert it with a regular outer sum operation
    if(b1!=B1.end()){
      out[b2->first]=outersum(b1->second,b2->second);
    }else{
      //b: the block b2 is NOT present in B1; this is a more delicate problem. When copying it over into the other MPS, we
      //have to insert it in a consistent way: b2 in principle has to be put in as a second block on the lower left part of b1, if b1 is present. Now b2 is not yet present
      //but other blocks with EITHER the left OR the right quantum number might be here. If that's the case then b2 has to be inserted consistent with these blocks.
      uint size1=0,size2=0;
      ir=sizeright.find(b2->first[1]);
      il=sizeleft.find(b2->first[0]);
      if(ir!=sizeright.end()){
	size2=ir->second;
      }
      if(il!=sizeleft.end()){
	size1=il->second;
      }
      matrix<T,column_major>mat(size1,size2);
      boost_setTo(mat,T(0.0));
      out[b2->first]=outersum(mat,b2->second);
    }

  }
  return out;
}

//template<int Ns,class KeyType,typename T>
//BlockMatrix<Ns,KeyType,T> operator*(const T num,const BlockMatrix<Ns,KeyType,T>&B){
//  BlockMatrix<Ns,KeyType,T>cpy=B;
//  return B*=num;
//}


template<int Ns,class KeyType,typename T>
BlockMatrix<Ns,KeyType,T> operator*(const DiagBlockMatrix<KeyType,Real>&lam,const BlockMatrix<Ns,KeyType,T>&B){
  typename BlockMatrix<Ns,KeyType,T>::const_iterator b;
  typename DiagBlockMatrix<KeyType,Real>::const_iterator l;
  BlockMatrix<Ns,KeyType,T> out;
  matrix<T,column_major>tmp;
  for(b=B.begin();b!=B.end();++b){
    l=lam.find(b->first[0]);
    if(l!=lam.end()){
      MatMatProd(T(1.0),l->second,b->second,tmp,'n','n');
      out[b->first] = tmp;
    }else{
      cout<<"in operator*(const DiagBlockMatrix<KeyType,Real>&lam,const BlockMatrix<Ns,KeyType,T>&B): no lambda found for key "<<b->first[0]<<
	" of BlockMatrix"<<endl;
    }
  }
  return out;
}

template<int Ns,class KeyType,typename T>
BlockMatrix<Ns,KeyType,T> operator*(const BlockMatrix<Ns,KeyType,T>&B,const DiagBlockMatrix<KeyType,Real>&lam){
  typename BlockMatrix<Ns,KeyType,T>::const_iterator b;
  typename DiagBlockMatrix<KeyType,Real>::const_iterator l;
  BlockMatrix<Ns,KeyType,T> out;
  matrix<T,column_major>tmp;
  for(b=B.begin();b!=B.end();++b){
    l=lam.find(b->first[1]);
    if(l!=lam.end()){
      MatMatProd(T(1.0),b->second,l->second,tmp,'n','n');
      out[b->first] = tmp;
    }else{
      cout<<"in operator*(const BlockMatrix<Ns,KeyType,T>&B,const DiagBlockMatrix<KeyType,Real>&lam): no lambda found for key "<<b->first[0]<<
	" of BlockMatrix"<<endl;
    }
  }
  return out;
}

template<int Ns,class KeyType,typename T>
BlockMatrix<Ns,KeyType,T> operator/(const BlockMatrix<Ns,KeyType,T>&B,const DiagBlockMatrix<KeyType,Real>&lam){
  typename BlockMatrix<Ns,KeyType,T>::const_iterator b;
  typename DiagBlockMatrix<KeyType,Real>::const_iterator l;
  BlockMatrix<Ns,KeyType,T> out;
  matrix<T,column_major>tmp;
  for(b=B.begin();b!=B.end();++b){
    l=lam.find(b->first[1]);
    if(l!=lam.end()){
      VectorType t(l->second.size());
      for(uint i=0;i<l->second.size();++i)
	t(i) = 1.0/l->second(i);
      MatMatProd(T(1.0),b->second,t,tmp,'n','n');
      out[b->first] = tmp;
    }else{
      cout<<"in operator/(const BlockMatrix<Ns,KeyType,T>&B,const DiagBlockMatrix<KeyType,Real>&lam): no lambda found for key "<<b->first[1]<<
	" of BlockMatrix"<<endl;
    }
  }
  return out;
}

template<int Ns,class KeyType,typename T>
BlockMatrix<Ns,KeyType,T> operator/(const DiagBlockMatrix<KeyType,Real>&lam,const BlockMatrix<Ns,KeyType,T>&B){
  typename BlockMatrix<Ns,KeyType,T>::const_iterator b;
  typename DiagBlockMatrix<KeyType,Real>::const_iterator l;
  BlockMatrix<Ns,KeyType,T> out;
  matrix<T,column_major>tmp;
  for(b=B.begin();b!=B.end();++b){
    l=lam.find(b->first[0]);
    if(l!=lam.end()){
      VectorType t(l->second.size());
      for(uint i=0;i<l->second.size();++i)
	t(i) = 1.0/l->second(i);
      MatMatProd(T(1.0),t,b->second,tmp,'n','n');
      out[b->first] = tmp;
    }else{
      cout<<"in operator/(const DiagBlockMatrix<KeyType,Real>&lam,const BlockMatrix<Ns,KeyType,T>&B): no lambda found for key "<<b->first[0]<<
	" of BlockMatrix"<<endl;
    }
  }
  return out;
}




template<int Ns,class KeyType,typename T>
ostream& operator <<(std::ostream &o,const BlockMatrix<Ns,KeyType,T> &B)
{
  typename BlockMatrix<Ns,KeyType,T>::const_iterator b;
  o<<"BlockMatrix<"<<Ns<<",KeyType,T> of size: "<<B.size()<<"\n";
  for(b=B.begin();b!=B.end();++b){
    o<<b->first<<endl;
    printboostmatrix(b->second);
    cout<<endl;
  }
  return o;
} 

template<int Ns,class KeyType,typename T>
T checkConvergence(const BlockMatrix<Ns,KeyType,T>&B,const string dir){
  typename BlockMatrix<Ns,KeyType,T>::const_iterator b;
  boost::numeric::ublas::matrix<T,column_major> out,m;
  map_container<KeyType,boost::numeric::ublas::matrix<T,column_major> > OUTS;
  typename map_container<KeyType,boost::numeric::ublas::matrix<T,column_major> >::iterator o;
  for(b=B.begin();b!=B.end();++b){
    m = b->second;
    if(dir=="l"){
      MatMatProd(T(1.0),conj(b->second),m,out,'t','n');
      if(OUTS.find(b->first[1])!=OUTS.end())
	 out+=OUTS[b->first[1]];
      OUTS[b->first[1]]=out;
    }else if(dir=="r"){
      MatMatProd(T(1.0),b->second,conj(m),out,'n','t');
      if(OUTS.find(b->first[0])!=OUTS.end())
	 out+=OUTS[b->first[0]];
      OUTS[b->first[0]]=out;
    }
  }  
  T tr=T(0.0);
  for(o=OUTS.begin();o!=OUTS.end();++o){
    tr+=trace(o->second);
  }
  return tr;
}



template<int Ns,class KeyType,typename T>
Real checkO(const BlockMatrix<Ns,KeyType,T>&B,const string dir){
  typename BlockMatrix<Ns,KeyType,T>::const_iterator b;
  boost::numeric::ublas::matrix<T,column_major> out,m;
  map_container<KeyType,boost::numeric::ublas::matrix<T,column_major> > OUTS;
  typename map_container<KeyType,boost::numeric::ublas::matrix<T,column_major> >::iterator o;
  Real max = -1e100;
  for(b=B.begin();b!=B.end();++b){
    m = b->second;
    if(dir=="l"){
      MatMatProd(T(1.0),conj(b->second),m,out,'t','n');
      if(OUTS.find(b->first[1])!=OUTS.end())
	 out+=OUTS[b->first[1]];
      OUTS[b->first[1]]=out;
    }else if(dir=="r"){
      MatMatProd(T(1.0),b->second,conj(m),out,'n','t');
      if(OUTS.find(b->first[0])!=OUTS.end())
	 out+=OUTS[b->first[0]];
      OUTS[b->first[0]]=out;
    }
  }  
  for(o=OUTS.begin();o!=OUTS.end();++o){
    out=o->second;
    assert(out.size1()==out.size2());
    boost::numeric::ublas::identity_matrix<T> id(out.size1());
    out-=id;
    for(uint x=0;x<out.size1();++x){
      for(uint y=0;y<out.size2();++y){
	max = abs(out(x,y))>max?abs(out(x,y)):max;
      }
    }
  }
  return max;
}

template<int Ns,class KeyType,typename T>
bool checkBMs(BlockMatrix<Ns,KeyType,T>&BM1,BlockMatrix<Ns,KeyType,T>&BM2){
  typename BlockMatrix<Ns,KeyType,T>::iterator i1,i2;
  std::vector<typename BlockMatrix<Ns,KeyType,T>::iterator>toremove;
  bool unconected = true;
  bool isOK=true;
  for(i1=BM1.begin();i1!=BM1.end();++i1){
    unconected=true;
    for(i2=BM2.begin();i2!=BM2.end();++i2){
      if(i1->first[1]==i2->first[0]){
	if (i1->second.size2()!=i2->second.size1()){
	  cout<<"size mismatch for center key "<< i1->first[1]<<endl;
	  isOK=false;
	}
	unconected = false;
      }
    }
    if(unconected){
      isOK=false;
      //   cout<<"unconnected outgoing key of BlockMatrix<Ns,KeyType,T> B1 found:\n \t "<<i1->first[1]<<endl;
      toremove.push_back(i1);
    }
  }
  for(uint i=0;i<toremove.size();++i){
    BM1.erase(toremove[i]);
  }

  unconected = true;
  toremove.clear();
  for(i2=BM2.begin();i2!=BM2.end();++i2){
    unconected=true;
    for(i1=BM1.begin();i1!=BM1.end();++i1){
      if(i1->first[1]==i2->first[0]){
	unconected = false;
      }
    }
    if(unconected){
      isOK=false;
      //      cout<<"unconected ingoing key of BlockMatrix<Ns,KeyType,T> B2 found:\n \t "<<i2->first[0]<<endl;
      toremove.push_back(i2);
    }
  }
  for(uint i=0;i<toremove.size();++i){
    BM2.erase(toremove[i]);
  }
  return isOK;
}

template<int Ns,class KeyType,typename T>
bool dropDangling(BlockMatrix<Ns,KeyType,T>&BM1,BlockMatrix<Ns,KeyType,T>&BM2){
  typename BlockMatrix<Ns,KeyType,T>::iterator i1,i2;
  std::vector<typename BlockMatrix<Ns,KeyType,T>::iterator>toremove;
  bool unconected = true;
  bool isOK=true;
  for(i1=BM1.begin();i1!=BM1.end();++i1){
    unconected=true;
    for(i2=BM2.begin();i2!=BM2.end();++i2){
      if(i1->first[1]==i2->first[0]){
	unconected = false;
      }
    }
    if(unconected){
      isOK=false;
      toremove.push_back(i1);
    }
  }
  for(uint i=0;i<toremove.size();++i)
    BM1.erase(toremove[i]);

  unconected = true;
  toremove.clear();
  for(i2=BM2.begin();i2!=BM2.end();++i2){
    unconected=true;
    for(i1=BM1.begin();i1!=BM1.end();++i1){
      if(i1->first[1]==i2->first[0]){
	unconected = false;
      }
    }
    if(unconected){
      isOK=false;
      toremove.push_back(i2);
    }
  }
  for(uint i=0;i<toremove.size();++i)
    BM2.erase(toremove[i]);
  return isOK;
}


template<class KeyType,typename T>
ostream& operator <<(std::ostream &o,const DiagBlockMatrix<KeyType,T> &B){
  typename DiagBlockMatrix<KeyType,T>::const_iterator b;
  o<<"size of DiagBlockMatrix: "<<B.size()<<"\n";
  for(b=B.begin();b!=B.end();++b){
    o<<b->first<<"\t lam: "<<b->second<<endl;;
  }
  return o;
} 

#endif 
