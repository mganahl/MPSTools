#ifndef KEYTYPES_HPP__
#define KEYTYPES_HPP__

#include <vector>
#include <utility>
#include "../boost_typedefs.hpp"
#include "helperfunctions.hpp"
#include "utilities.hpp"
#include <istream>
using namespace BoostTypedefs;
//Every new Key needs to have operator< be defined here:

//this class needs all routines like <,> == .... as virtual functions;

void reset(long unsigned int &k){
  k=0;
}


class AbKeyType
{
public:
  AbKeyType(){_qn=0;_N=0;}
  AbKeyType(const long unsigned int &in,const int &N){assert(N==1);_qn=in;_N=N;}
  AbKeyType(const long unsigned int &q1,const long unsigned int &q2,const int &N){
    assert(N==2);
    _qn = (q2<<16);
    _qn+=q1;
    _N=N;
  };
  long unsigned int getQ()const {return _qn;}
  uint size()const{return _N;}
  void resize(const uint N){_N=N;}
  void reset(){_qn=0;}
  long int operator()(const uint &ind)const{
    long unsigned int pat=0;
    uint n=0;
    while(n<16){
      pat+=(1<<n);
      ++n;
    }
    assert(ind<_N);
    long unsigned int o=_qn;
    return (o>>(ind*16))&pat;
  }
  AbKeyType &operator+=(const AbKeyType &in){
    if(this==&in){cout<<"in AbKeyType::operator += (): selfassignemt"<<endl;abort();}
    _qn += in.getQ();
    _N = in.size();
    return (*this);
  }
  AbKeyType operator+(const AbKeyType &in)const{
    AbKeyType cpy = (*this);
    return cpy+=in;
  }
  AbKeyType &operator-=(const AbKeyType &in){
    if(this==&in){cout<<"in AbKeyType::operator += (): selfassignemt"<<endl;abort();}
    _qn -= in.getQ();
    _N = in.size();
    return (*this);
  }
  AbKeyType operator-(const AbKeyType &in)const{
    AbKeyType cpy = (*this);
    return cpy-=in;
  }

  template<class Archive>
  void serialize(Archive & ar, const int version){
	  ar & _qn;
	  ar & _N;
  }
  AbKeyType operator ++(){
    ++_qn;
    return (*this);
  }

private:
  long unsigned int _qn;
  uint _N;

};



void reset(AbKeyType &k)
{k.reset();}


bool operator<(const AbKeyType &k1,const AbKeyType &k2){
  return k1.getQ()<k2.getQ();
}


bool operator==(const AbKeyType &k1,const AbKeyType &k2){
  return k1.getQ()==k2.getQ();
}



bool operator>(const AbKeyType &k1,const AbKeyType &k2){
  return k1.getQ()>k2.getQ();
}


bool operator>=(const AbKeyType &k1,const AbKeyType &k2){
  return k1.getQ()>=k2.getQ();
}


bool operator<=(const AbKeyType &k1,const AbKeyType &k2){
  return k1.getQ()<=k2.getQ();
}

bool operator!=(const AbKeyType &k1,const AbKeyType &k2){
  return k1.getQ()!=k2.getQ();
}

std::ostream&operator<<(std::ostream &o,const AbKeyType &ky){
  o<<"(";
  for(uint n=0;n<ky.size();++n)
    o<<ky(n)<<",";
  o<<")";
  return o;
}

template<typename K1, typename K2>
class doubleKey{
public:
  K1 first;
  K2 second;
  doubleKey():first(K1()),second(K2()) {}
  doubleKey( const K1& y,const K2& x) : first(x), second(y){}
  //  template <typename K1,typename K2>
  doubleKey(const doubleKey<K1,K2> &p) : first(p.first), second(p.second){}
};


template<typename K1,typename K2>
bool operator<(const doubleKey<K1,K2> &k1,const doubleKey<K1,K2> &k2)
{
  if(k1.first==k2.first)
    return k1.second<k2.second;
  else if(k1.first!=k2.first)
    return k1.first<k2.first;
  return 0;
}

template<typename K1,typename K2>
bool operator==(const doubleKey<K1,K2> &k1,const doubleKey<K1,K2> &k2)
{
  if(k1.first==k2.first)
    return k1.second==k2.second;
  else if(k1.first!=k2.first)
    return false;
  return false;
}

template<typename K1,typename K2>
ostream& operator <<(std::ostream &o,const doubleKey<K1,K2> &k)
{
  o<<"first: "<<k.first<<"; second: "<<k.second;
  return o;
  
} 

//==============================================================TENSORKEYTYPE
template<size_t Nb,size_t Ns,typename KeyType>
class TensorKeyType
{
public:
//  __attribute__((noinline))  TensorKeyType();
//  __attribute__((noinline))  TensorKeyType(const KeyType &key0);
//  __attribute__((noinline))  TensorKeyType(const KeyType &key0,const KeyType &key1);
//  __attribute__((noinline))  TensorKeyType(const KeyType &key0,const KeyType &key1,const uint &p);
  //  __attribute__((noinline))  TensorKeyType(const KeyType &key0,const KeyType &key1,const std::vector<uint> &p);
  //  __attribute__((noinline))  TensorKeyType(const std::vector<KeyType> &k,const std::vector<uint> &p);

  TensorKeyType();
  TensorKeyType(const KeyType &key0);
  TensorKeyType(const KeyType &key0,const KeyType &key1);
  TensorKeyType(const KeyType &key0,const KeyType &key1,const uint &p);
  TensorKeyType(const KeyType &key0,const KeyType &key1,const std::vector<uint> &p);
  TensorKeyType(const std::vector<KeyType> &k,const std::vector<uint> &p);

  //~TensorKeyType(){delete[] _key;delete[] _p;}

  //  void resize(const uint size1,const uint size2){_key.resize(size1);_p.resize(size2);_Nb = size1;_Ns=size2;}
  //  void resize(const uint size1,const uint size2){_key=new KeyType[size1];_p=new uint[size2];}
  
  //[] acccesses the aux index  of type KeyType
  KeyType operator[](const int &i) const {assert(i<_Nb);return _key[i];}
  KeyType &operator[](const int &i) {assert(i<_Nb);return _key[i];}
  KeyType operator[](const uint &i) const {assert(i<_Nb);return _key[i];}
  KeyType &operator[](const uint &i) {assert(i<_Nb);return _key[i];}

  //() accesses the physical index of type uint
  uint operator()(const int &i) const {assert(i<_Ns);return _p[i];}
  uint &operator()(const int &i) {assert(i<_Ns);return _p[i];}
  uint operator()(const uint &i) const {assert(i<_Ns);return _p[i];}
  uint &operator()(const uint &i) {assert(i<_Ns);return _p[i];}
  


  TensorKeyType<Nb,Ns,KeyType> & operator=(const TensorKeyType<Nb,Ns,KeyType>&ksrc);
  const std::vector<uint> getP()const {std::vector<uint>p(_Ns);for(uint n=0;n<_Ns;++n)p[n]=_p[n];return p;}
  const std::vector<KeyType> getK()const {std::vector<KeyType>k(_Nb);for(uint n=0;n<_Nb;++n)k[n]=_key[n];return k;}
  uint size() const{return _Nb;}
  uint length() const{return _Ns;}//_p.size();}


  template<class Archive>
  void serialize(Archive & ar, const unsigned int version){
	  ar & _p;
	  ar & _key;
	  ar & _Nb;
	  ar & _Ns;
  }
  

protected:
  // std::vector<uint> _p;
  // std::vector<KeyType> _key;
  KeyType _key[Nb];
  uint _p[Ns];
  uint _Nb,_Ns;
};

template<size_t Nb, size_t Ns, class KeyType>
//__attribute__((noinline))TensorKeyType<Nb,Ns,KeyType>:: TensorKeyType(){_Nb = Nb;_Ns=Ns;/*_key.resize(_Nb);_Ns = Ns;_p.resize(_Ns);*/}
TensorKeyType<Nb,Ns,KeyType>:: TensorKeyType(){_Nb = Nb;_Ns=Ns;/*_key.resize(_Nb);_Ns = Ns;_p.resize(_Ns);*/}

template<size_t Nb, size_t Ns, class KeyType>
TensorKeyType<Nb,Ns,KeyType>:: TensorKeyType(const KeyType &key0){
  //__attribute__((noinline))TensorKeyType<Nb,Ns,KeyType>:: TensorKeyType(const KeyType &key0){
  assert(Ns==1);
  assert(Nb==1);
  _Ns = Ns;
  _Nb = Nb;
  //  _key.resize(_Nb);
  _key[0] = key0;
  //_p.resize(_Ns);
  _p[0] = 1000;
}

template<size_t Nb, size_t Ns, class KeyType>
TensorKeyType<Nb,Ns,KeyType>:: TensorKeyType(const KeyType &key0,const KeyType &key1){
  //__attribute__((noinline))TensorKeyType<Nb,Ns,KeyType>:: TensorKeyType(const KeyType &key0,const KeyType &key1){
  assert(Nb==2);
  _Nb = Nb;
  //_key.resize(_Nb);
  _key[0] = key0;_key[1] = key1;
  _Ns=Ns;
  //_p.resize(_Ns);
  for(uint i=0;i<_Ns;++i)_p[i]=1000;
}

template<size_t Nb, size_t Ns, class KeyType>
TensorKeyType<Nb,Ns,KeyType>:: TensorKeyType(const KeyType &key0,const KeyType &key1,const uint &p){
  //__attribute__((noinline))TensorKeyType<Nb,Ns,KeyType>:: TensorKeyType(const KeyType &key0,const KeyType &key1,const uint &p){
  assert(Ns==1);
  assert(Nb==2);
  _Ns=Ns;
  //_p.resize(Ns);
  _p[0]=p;
  _Nb = Nb;
  //_key.resize(_Nb);
  _key[0] = key0;_key[1] = key1;
}

template<size_t Nb, size_t Ns, class KeyType>
TensorKeyType<Nb,Ns,KeyType>:: TensorKeyType(const KeyType &key0,const KeyType &key1,const std::vector<uint> &p)/*:_p(p)*/{
  //__attribute__((noinline))TensorKeyType<Nb,Ns,KeyType>:: TensorKeyType(const KeyType &key0,const KeyType &key1,const std::vector<uint> &p)/*:_p(p)*/{
  assert(p.size()==Ns);
  assert(Nb==2);
  _Ns=Ns;
  _Nb=Nb;
  //_key.resize(_Nb);
  _key[0] = key0;_key[1] = key1;
  for(uint n=0;n<p.size();++n)
    _p[n]=p[n];
}

template<size_t Nb, size_t Ns, class KeyType>
TensorKeyType<Nb,Ns,KeyType>:: TensorKeyType(const std::vector<KeyType> &k,const std::vector<uint> &p)/*:_key(k),_p(p)*/{
  //__attribute__((noinline))TensorKeyType<Nb,Ns,KeyType>:: TensorKeyType(const std::vector<KeyType> &k,const std::vector<uint> &p)/*:_key(k),_p(p)*/{
  assert(p.size()==Ns);
  assert(k.size()==Nb);
  _Ns = Ns;//_p.size();
  _Nb = Nb;// _key.size();
  for(uint n=0;n<p.size();++n)
    _p[n]=p[n];
  for(uint n=0;n<k.size();++n)
    _key[n]=k[n];
}


 template<size_t Nb, size_t Ns,typename KeyType>
 TensorKeyType<Nb,Ns,KeyType>&TensorKeyType<Nb,Ns,KeyType>:: operator=(const TensorKeyType<Nb,Ns,KeyType>&ksrc){
   if (this == &ksrc){std::cout<<"in TensorKeyType::operator=: selfassignment"<<std::endl;abort();}
   _Ns = ksrc.length();
   assert(_Ns==Ns);
   _Nb=ksrc.size();
   assert(_Nb==Nb);
   std::vector<uint> p=ksrc.getP();
   std::vector<KeyType> k=ksrc.getK();
   for(uint n=0;n<_Ns;++n)
     _p[n]=p[n];
   for(uint n=0;n<_Nb;++n)
     _key[n]=k[n];
   return *this;
 }

template<size_t Nb, size_t Ns,typename KeyType>
bool operator<(const TensorKeyType<Nb,Ns,KeyType> &key1,const TensorKeyType<Nb,Ns,KeyType> &key2){
  assert(key1.size()==key2.size());
  assert(key1.length()==key2.length());
  uint i = 0,j=0;
  while(i<key1.size()){
    if(key1[i]==key2[i]){
      ++i;
    }else{
      break;
    }
  }
  if(i==key1.size()){
    while(j<key1.length()){
      if(key1(j)==key2(j)){
	++j;
      }else{
	break;
      }
    }
  }
  bool r=false;
  if((j==key1.length())&&(i==key1.size()))
    r=false;
  if((i==key1.size())&&(j<key1.length()))
    r=key1(j)<key2(j);
  if(i<key1.size())
    r=key1[i]<key2[i];;

  return r;//(i==key1.size())?(key1[j]<key2[j]):(key1[i]<key2[i]);
}

template<size_t Nb, size_t Ns,typename KeyType>
bool operator<=(const TensorKeyType<Nb,Ns,KeyType> &key1,const TensorKeyType<Nb,Ns,KeyType> &key2){
  return (key1<key2)?true:((key2<key1)?false:true);
}

template<size_t Nb, size_t Ns,typename KeyType>
bool operator>(const TensorKeyType<Nb,Ns,KeyType> &key1,const TensorKeyType<Nb,Ns,KeyType> &key2){
  return (key2<key1);
}

template<size_t Nb, size_t Ns,typename KeyType>
bool operator>=(const TensorKeyType<Nb,Ns,KeyType> &key1,const TensorKeyType<Nb,Ns,KeyType> &key2){
  return (key2<key1)?true:((key1<key2)?false:true);
}

template<size_t Nb, size_t Ns,typename KeyType>
bool operator==(const TensorKeyType<Nb,Ns,KeyType> &key1,const TensorKeyType<Nb,Ns,KeyType> &key2){
  return (!(key1<key2))&&(!(key2<key1));
}

template<size_t Nb, size_t Ns,typename KeyType>
bool operator!=(const TensorKeyType<Nb,Ns,KeyType> &key1,const TensorKeyType<Nb,Ns,KeyType> &key2){
  return !(key1==key2);
}


template<size_t Nb,size_t Ns,class KeyType>
ostream& operator <<(std::ostream &o,const TensorKeyType<Nb,Ns,KeyType> &key)
{
  o<<"a:"<<key.getK()<<" p:"<<key.getP();  
  return o;
  
} 





//Ns is the number of sites om which it acts
template<size_t Ns>
class OpKeyType:public TensorKeyType<Ns,Ns,uint>
{
public:
  OpKeyType (){};
//  __attribute__((noinline)) OpKeyType(const std::vector<uint> &out,const std::vector<uint> &in);
//  __attribute__((noinline)) OpKeyType(const uint &out,const uint &in);
//  __attribute__((noinline)) OpKeyType(const int &out,const int &in);

  OpKeyType(const std::vector<uint> &out,const std::vector<uint> &in);
  OpKeyType(const uint &out,const uint &in);
  OpKeyType(const int &out,const int &in);

  std::vector<uint> getIn()const {return this->getK();}
  std::vector<uint> getOut()const {return this->getP();}
  //operator[uint] gives back incoming index at site uint
  //operator(uint) gives back outgoing index at site uint
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version){
    ar & boost::serialization::base_object<TensorKeyType<Ns,Ns,uint> >(*this);
  }

};

template<size_t Ns>
OpKeyType<Ns>:: OpKeyType(const std::vector<uint> &out,const std::vector<uint> &in):TensorKeyType<Ns,Ns,uint>(in,out){};
//__attribute__((noinline))OpKeyType<Ns>:: OpKeyType(const std::vector<uint> &out,const std::vector<uint> &in):TensorKeyType<Ns,Ns,uint>(in,out){};

template<size_t Ns>
OpKeyType<Ns>:: OpKeyType(const uint &out,const uint &in):TensorKeyType<Ns,Ns,uint>(){
  //__attribute__((noinline))OpKeyType<Ns>:: OpKeyType(const uint &out,const uint &in):TensorKeyType<Ns,Ns,uint>(){
  /*this->resize(1,1);*/(*this)[0]=in;(*this)(0)=out;}

template<size_t Ns>
OpKeyType<Ns>:: OpKeyType(const int &out,const int &in):TensorKeyType<Ns,Ns,uint>(){
  //__attribute__((noinline))OpKeyType<Ns>:: OpKeyType(const int &out,const int &in):TensorKeyType<Ns,Ns,uint>(){
  /*this->resize(1,1);*/(*this)[0]=uint(in);(*this)(0)=uint(out);
}


template<size_t Nb>
class ukey:public TensorKeyType<Nb,1,uint>{
public:
  ukey():TensorKeyType<Nb,1,uint>(){};
  ukey(uint i1,uint i2);
  //  ukey(int i1,int i2);
  ~ukey(){};
  
  uint operator[](const int &i) const {assert(i<Nb);assert(i>=0);return TensorKeyType<Nb,1,uint>::operator[](i);}
  uint &operator[](const int &i) {assert(i<Nb);assert(i>=0);return TensorKeyType<Nb,1,uint>::operator[](i);}
  uint operator[](const uint &i) const {assert(i<Nb);return TensorKeyType<Nb,1,uint>::operator[](i);}
  uint &operator[](const uint &i) {assert(i<Nb);return TensorKeyType<Nb,1,uint>::operator[](i);}
  uint operator()(const int &i) const {assert(i<Nb);assert(i>=0);return TensorKeyType<Nb,1,uint>::operator[](i);}
  uint &operator()(const int &i) {assert(i<Nb);assert(i>=0);return TensorKeyType<Nb,1,uint>::operator[](i);}
  uint operator()(const uint &i) const {assert(i<Nb);return TensorKeyType<Nb,1,uint>::operator[](i);}
  uint &operator()(const uint &i) {assert(i<Nb);return TensorKeyType<Nb,1,uint>::operator[](i);}

  void frombuffer(std::istream &buf){
    buf.read((char*)&(this->_key[0]), Nb*sizeof(unsigned int));
  }
  void tobuffer(std::ostream &buf) const{
    buf.write((char*)&(this->_key[0]), Nb*sizeof(unsigned int));
  }
  const uint* raw() const{return &(this->_key[0]);}

  template<class Archive>
  void serialize(Archive & ar, const unsigned int version){
    ar & boost::serialization::base_object<TensorKeyType<Nb,1,uint> >(*this);
  }

};

template<size_t Nb>
ukey<Nb>:: ukey(uint i1,uint i2):TensorKeyType<Nb,1,uint>(){
  this->_Nb=Nb;this->_Ns=1;
  TensorKeyType<Nb,1,uint>::operator[](0)=i1;
  TensorKeyType<Nb,1,uint>::operator[](1)=i2;
};

// template<size_t Nb>
// ukey<Nb>:: ukey(int i1,int i2):TensorKeyType<Nb,1,uint>(){this->_Nb=Nb;this->_Ns=1;TensorKeyType<Nb,1,uint>::operator[](0)=i1;TensorKeyType<Nb,1,uint>::operator[](1) = i2;};

// template<size_t Nb>
// ukey<Nb>:: ukey(int i1,int i2):TensorKeyType<Nb,1,uint>(){this->_Nb=Nb;this->_Ns=1;TensorKeyType<Nb,1,uint>::operator[](0)=uint(i1);TensorKeyType<Nb,1,uint>::operator[](1) = uint(i2);};



template<size_t Nb>
bool operator<(const ukey<Nb> &key1,const ukey<Nb> &key2){
  assert(key1.size()==key2.size());
  uint i = 0;
  while(i<key1.size()){
    if(key1[i]==key2[i]){
      ++i;
    }else{
      break;
    }
  }
  return (i==key1.size())?false:(key1[i]<key2[i]);
}

#endif









