#ifndef TENSOR_H__
#define TENSOR_H__
#include "keyTypes.hpp"
#include <iostream>
using boost::numeric::ublas::column_major;
using boost::numeric::ublas::matrix;
using boost::numeric::ublas::matrix_range;

//this is just to hide matrix class; 



template<uint Nb,typename T>
class tensor:public matrix<T,column_major>{
public:
  tensor():matrix<T,column_major>(){};
  tensor(const matrix<T,column_major>&mat):matrix<T,column_major>(mat){};
  tensor(const uint d1):matrix<T,column_major>(d1,1){};
  tensor(const uint d1, const uint d2):matrix<T,column_major>(d1,d2){};
  tensor<Nb,T> & operator =(const tensor <Nb,T> &src){
    assert(Nb==2);
    if(this==&src){std::cout<<"in tensor<Nb,T> =: selfassignment"<<std::endl;abort();}
    this->resize(src.size1(),src.size2());
    for(uint i = 0;i<src.size1()*src.size2();++i){
      this->data()[i]=src.data()[i];
    }
    return *this;
  };
  //tensor(const tensor<Nb,T,StorageOrder> t):matrix<T,StorageOrder>(t){};
};


#endif
