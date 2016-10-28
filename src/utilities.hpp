#ifndef UTILITIES_HPP__ 
#define UTILITIES_HPP__

#include <vector>
#include <set>
#include <sys/stat.h>
#include <sstream>
#include <iostream>
#include <fstream>
//#include "tensormap.hpp"

using namespace BoostTypedefs;
using boost::numeric::ublas::matrix;
using boost::numeric::ublas::matrix_range;
using boost::numeric::ublas::column_major;


UintVectorType createTwoorbitalfilling(uint N0down,uint N0up,uint N0,uint N1down,uint N1up,uint N1){
  UintVectorType localState(N0+N1+2);
  boost_setTo(localState,uint(0));
    uint n0d=0;
    if(N0down>(N0+1)){cout<<"too many down particles in orbital 0!!! choose N0d<=N0"<<endl;abort();}
    while (n0d<N0down){
      int site =(int)std::floor(std::rand()%(N0+1));

      if (localState(site)==0){localState(site)=1;n0d++;}
      if (localState(site)==2){localState(site)=3;n0d++;}
    }

    uint n0u=0;
    if(N0up>(N0+1)){cout<<"too many up particles in orbital 0!!! choose N0u<=N0"<<endl;abort();}
    while (n0u<N0up){
      int site =(int)std::floor(std::rand()%(N0+1));

      if (localState(site)==0){localState(site)=2;n0u++;}
      if (localState(site)==1){localState(site)=3;n0u++;}
    }

    uint n1d=0;
    if(N1down>(N1+1)){cout<<"too many down particles in orbital 1!!! choose N1d<=N1"<<endl;abort();}
    while (n1d<N1down){
      int site =N0+1+(int)std::floor(std::rand()%(N1+1));

      if (localState(site)==0){localState(site)=1;n1d++;}
      if (localState(site)==2){localState(site)=3;n1d++;}
    }

    uint n1u=0;
    if(N1up>(N1+1)){cout<<"too many up particles in orbital 1!!! choose N1u<=N1"<<endl;abort();}
    while (n1u<N1up){
      int site =N0+1+(int)std::floor(std::rand()%(N1+1));

      if (localState(site)==0){localState(site)=2;n1u++;}
      if (localState(site)==1){localState(site)=3;n1u++;}
    }
    return localState;

}


bool fexist(std::string filename){
  struct stat buf;
  if(stat(filename.c_str(),&buf)==0)
    return true;
  return false;
}


MatrixType readMatrix(std::string filename){

  std::ifstream infile;
  if(!fexist(filename)){
    cout<<"in utilities.hpp: readMatrix::error reading "<<filename<<": no such file or directory"<<endl;
    abort();
  }
  infile.open(filename.c_str(),ios_base::in);
  if(!infile.good()){
    cout<<"in utilities.hpp: readMatrix::error reading "<<filename<<": file is broken"<<endl;
    abort();
  }

  std::string line;
  uint N=0,M1,M2;
  std::vector<double>vec;
  while(std::getline(infile,line)){
    std::stringstream iss(line);
    M1=0;
    double value;
    while(iss>>value){
      vec.push_back(value);
      M1++;
    }
    if(N>0){
      if(M1!=M2){
	cerr<<"in utilities.hpp: readMatrix:: data in file "<<filename<<" is not matrix-shaped (see element ("<<N<<","<<M1<<"))"<<endl;
	cerr<<"aborting"<<endl;
	abort();
      }
    }
    M2=M1;
    N++;
  }
  uint M=M1;
  infile.close();

  MatrixType data(N,M);
  for(uint x=0;x<N;x++){
    for(uint y=0;y<M;y++){
      data(x,y)=vec[M*x+y];
    }
  }
  //infile.open(filename.c_str(),ios_base::in);
  //uint x=0,y;
  //while(std::getline(infile,line)){
  //  std::stringstream iss(line);
  //  double value;
  //  y=0;
  //  while(iss>>value){
  //    data(x,y)=value;
  //    y++;
  //  }
  //  x++;
  //}
  return data;
}

template<typename T,typename StorageOrder>
T trace(boost::numeric::ublas::matrix<T,StorageOrder> &M)
{
  if(M.size1()!=M.size2()){cout<<"in trace: matrix must be square"<<endl;abort();}
  uint size = M.size1();
  T val = T(0.0);
  for(uint i = 0;i<size;i++)val+=M(i,i);
  return val;
}


double c_rand(double min,double max){
  assert(min<=max);
  return min +fabs(max-min)*rand()/RAND_MAX;
}


Complex c_rand(Complex min,Complex max){
  return Complex(c_rand(min.real(),max.real()),c_rand(min.imag(),max.imag()));
}


template<class storage_type>
void boost_random_mat(boost::numeric::ublas::matrix<Complex,storage_type> &mat,Real min=0,Real max = 1)
{
	for (uint index = 0; index < mat.size1()*mat.size2(); index++) {
			mat.data()[index]= Complex(c_rand(min, max), c_rand(min, max));
	}
}
template<class storage_type>
void boost_random_mat(boost::numeric::ublas::matrix<Real,storage_type> &mat,Real min=0,Real max = 1)
{
	for (uint index = 0; index < mat.size1()*mat.size2(); index++) {
			mat.data()[index]= c_rand(min, max);
	}
}





UintVectorType generateNextState(const UintVectorType old, const UintVectorType maximum){
  UintVectorType out=old;
  for(uint s=0;s<out.size();++s){
    if(out(s)<(maximum(s)-1)){
      ++out(s);
      for(uint s2=0;s2<s;++s2)
	out(s2)=0;
      break;
    }
  }
  return out;
}

template<typename T>
matrix<T,column_major> outersum(const matrix<T,column_major>&m1,const matrix<T,column_major>&m2){
  matrix<T,column_major>m;
  const int s11=m1.size1();
  const int s12=m1.size2();
  const int s21=m2.size1();
  const int s22=m2.size2();
  m.resize(s11+s21,s12+s22);
  boost_setTo(m,T(0.0));
  matrix_range<matrix<T, column_major> > (m,Range(0,s11),Range(0,s12))=m1;
  matrix_range<matrix<T, column_major> > (m,Range(s11,s11+s21),Range(s12,s12+s22))=m2;
  return m;
}

template<typename T1,typename T2>
matrix<T1,column_major> cast(const matrix<T2,column_major>&mat){
  matrix<T1,column_major>mat_cast(mat.size1(),mat.size2());
  for(uint x=0;x<mat.size1();++x){
    for(uint y=0;y<mat.size2();++y){
      mat_cast(x,y) = static_cast<T1>(mat(x,y));
    }
  }
  return mat_cast;
}



template<class T1,class T2,class T3>
struct triple
{
  T1 first;
  T2 second;
  T3 third;
  triple() : first(T1()), second(T2()),third(T3()) {}
  triple(const T1& x, const T2& y,const T3& z) : first(x), second(y),third(z){}
  template <class U, class V,class W>
  triple (const triple<U,V,W> &p) : first(p.first), second(p.second),third(p.third){}
};



template<class T1,class T2,class T3,class T4>
struct quadruple
{
  T1 first;
  T2 second;
  T3 third;
  T4 fourth;
  quadruple() : first(T1()), second(T2()),third(T3()),fourth(T3()){}
  quadruple(const T1& x, const T2& y,const T3& z,const T4& a) : first(x), second(y),third(z),fourth(a){}
  template <class U, class V,class W,class X>
  quadruple(const quadruple<U,V,W,X> &p) : first(p.first), second(p.second),third(p.third),fourth(p.fourth){}
};

uint MuToSiIndex(const std::vector<uint> &ind,const std::vector<uint> &dims){
  assert(ind.size()==dims.size());
  uint out=0;
  uint d=1;
  for(uint i=0;i<ind.size();++i){
    assert(ind[i]<dims[i]);
    out+=ind[i]*d;
    d*=dims[i];
  }
  return out;
}

std::vector<uint> SiToMuIndex(const uint index,const std::vector<uint> &dims){
  std::vector<uint> out(dims.size());
  uint ind=index;
  for(uint i=0;i<dims.size();++i){
    out[i]=ind%dims[i];
    ind-=out[i];
    ind/=dims[i];
  }
  return out;
}

template<typename T>
ostream& operator <<(std::ostream &o,const std::vector<T>&v){
  o<<"["<<v.size()<<"]"<<"(";
  uint i=0;
  if(v.size()!=0){
    for(i = 0;i<(v.size()-1);i++)
    o<<v[i]<<",";
  }
  o<<v[i]<<")";
  return o;
  
} 

template<typename T>
ostream& operator <<(std::ostream &o,const std::set<T>&v){
  typename std::set<T>::const_iterator it;
  o<<"["<<v.size()<<"]"<<"(";
  uint i=0;
  if(v.size()!=0){
    for(it = v.begin();it!=v.end();++it)
      o<<(*it)<<",";
  }
  return o;
} 




template<typename T>
bool operator<(boost::numeric::ublas::vector<T>& a,boost::numeric::ublas::vector<T>&b){
  typename boost::numeric::ublas::vector<T>::iterator ita = a.begin();
  typename boost::numeric::ublas::vector<T>::iterator itb = b.begin();
  assert(a.size()==b.size());
  while(ita!=a.end())
    {
      if(*ita==*itb){ita++;itb++;}
      else{
	break;
      }
    }
  return  ita==a.end()?false:(*ita<*itb);
}


//====================================================

// bool operator<=(UintVectorType &key1,UintVectorType &key2){
//   return (key1<key2)?true:((key2<key1)?false:true);
// }

// // bool operator<(UintVectorType &a,UintVectorType &b){
// //   typename boost::numeric::ublas::vector<uint>::iterator ita = a.begin();
// //   typename boost::numeric::ublas::vector<uint>::iterator itb = b.begin();

// //   assert(a.size()==b.size());
// //   while(ita!=a.end())
// //     {
// //       if(*ita==*itb){ita++;itb++;}
// //       else{
// //  	break;
// //       }
// //     }
// //   return  ita==a.end()?false:(*ita<*itb);
// // }

// bool operator<(UintVectorType &a,UintVectorType &b){
//   UintVectorType::iterator ita = a.begin();
//   UintVectorType::iterator itb = b.begin();

//   assert(a.size()==b.size());
//   while(ita!=a.end())
//     {
//       if(*ita==*itb){ita++;itb++;}
//       else{
//  	break;
//       }
//     }
//   return  ita==a.end()?false:(*ita<*itb);
// }

// bool operator>(UintVectorType &key1,UintVectorType &key2){
//   return (key2<key1);
// }

// bool operator>=(UintVectorType &key1,UintVectorType &key2){
//   return (key2<key1)?true:((key1<key2)?false:true);
// }

// bool operator==(UintVectorType &key1,UintVectorType &key2){
//   return (!(key1<key2))&&(!(key2<key1));
// }

// bool operator!=(UintVectorType &key1,UintVectorType &key2){
//   return !(key1==key2);
// }


template<typename T>
bool  operator <(const std::vector<T>& a,const std::vector<T>&b){
  assert(a.size()==b.size());
  uint i=0;
  while(i<a.size()){
    if(a[i]==b[i])++i;
    else break;
  }
  return i==a.size()?false:(a[i]<b[i]);
}

template<typename T>
std::vector<T> set2array(const std::set<T>& s)
{
    typedef typename std::set<T>::const_iterator iter;
    std::vector<T> arr(s.size());
    int j=0;
    for(iter i=s.begin();i!=s.end();++i)
        arr[j++] = *i;
    return arr;
}







#endif
