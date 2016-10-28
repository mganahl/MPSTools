/*
 * measurementRoutines.hpp
 *
 *  
 *      Author: martinjg
 */

#ifndef MEASUREMENTROUTINES_HPP_
#define MEASUREMENTROUTINES_HPP_
#include "MPSTypes.hpp"
#include "MPOTypes.hpp"


//this defines a single measurement; the operators in the list Ops are measured and the value is stored in a file filename
template<typename T>
class measurementObject
{
public:
  //TODO: less function for operator !!!!!
  measurementObject(std::vector<Operator<T>*> Ops_,string filename_,uint measureStep_)
    :Ops(Ops_),filename(filename_),measureStep(measureStep_),measureRealPart(true),measureImagPart(false){step=1;};
  measurementObject(std::vector<Operator<T>*> Ops_,string filename_,uint measureStep_,bool measurereal,bool measureimag)
    :Ops(Ops_),filename(filename_),measureStep(measureStep_),measureRealPart(true){step=1;measureRealPart=measurereal;measureImagPart=measureimag;};
  std::vector<Operator<T>*> Ops;
  std::string filename;

  uint measureStep;
  uint step;
  template<class KeyType>
  void measure(const MPS<KeyType,Real> &MPS);
  template<class KeyType>
  void measure(const CanonizedMPS<KeyType,Real> &MPS);
  template<class KeyType>
  void measure(const MPS<KeyType,Complex> &MPS);
  template<class KeyType>
  void measure(const CanonizedMPS<KeyType,Complex> &MPS);
  void writeMeasurements();
  void setMeasureRealPart(bool measureRealPart){ this->measureRealPart=measureRealPart; }
  void setMeasureImagPart(bool measureImagPart){ this->measureImagPart=measureImagPart; }
  void setMeasurementStep(uint measureStep){ this->measureStep=measureStep; }
  void resetStep(){this->step=1;}
protected:

  std::vector<std::vector<Real> > data;
  bool measureRealPart;
  bool measureImagPart;

};

template<typename T>
template<class KeyType>
void measurementObject<T>::measure(const CanonizedMPS<KeyType,Real> &MPS)
{
  typename std::vector<Operator<T>*>::iterator it;
  std::vector<Real> d;
  for(it=Ops.begin();it!=Ops.end();it++)
    {
      d.push_back((*it)->measure(MPS));
    }
  data.push_back(d);

}


template<typename T>
template<class KeyType>
void measurementObject<T>::measure(const CanonizedMPS<KeyType,Complex> &MPS)
{
  typename std::vector<Operator<T>*>::iterator it;
  std::vector<Real> d;
  for(it=Ops.begin();it!=Ops.end();it++)
    {
      Complex out =(*it)->measure(MPS);
      // if(abs(out.real())>1e-14){
      // 	 if(out.imag()>1e-14){
      // 	   cout<<"would return real of order "<<out.real()<<", but imag is of order "<<out.imag()<<endl;
      // 	   abort();
      // 	 }
      // 	 d.push_back(out.real());
      // }
      // else if(abs(out.imag())>1e-14){
      // 	 // if(out.imag()>1e-13){
      // 	 //   cout<<"would return imag of order "<<out.imag()<<", but real is of order "<<out.real()<<endl;
      // 	 //   abort();
      // 	 // }
      // 	 d.push_back(out.imag());
      // }
      if (measureRealPart&&(!measureImagPart)){
	d.push_back(out.real());
      }else if(measureImagPart&&(!measureRealPart)){
	d.push_back(out.imag());
      }else if(measureImagPart&&measureRealPart){
	d.push_back(out.real());
	d.push_back(out.imag());
      }
    }
  data.push_back(d);
}

template<typename T>
template<class KeyType>
void measurementObject<T>::measure(const MPS<KeyType,Real> &MPS)
{
  typename std::vector<Operator<T>*>::iterator it;
  std::vector<Real> d;
  for(it=Ops.begin();it!=Ops.end();it++)
    {
      d.push_back((*it)->measure(MPS));
    }
  data.push_back(d);

}


template<typename T>
template<class KeyType>
void measurementObject<T>::measure(const MPS<KeyType,Complex> &MPS)
{
  typename std::vector<Operator<T>*>::iterator it;
  std::vector<Real> d;
  for(it=Ops.begin();it!=Ops.end();it++)
    {
      Complex out =(*it)->measure(MPS);
      if (measureRealPart){
    	  d.push_back(out.real());
      }else{
    	  d.push_back(out.imag());
      }
    }
  data.push_back(d);
}


template<typename T>
void measurementObject<T>::writeMeasurements(){
  ofstream file; file.precision(12);
  file.open(filename.c_str(),ios_base::app);
  for (uint i=0; i<data.size(); i++){
    for (uint j=0; j<data[i].size(); j++){
      file<<data[i][j]<<" ";
    }
    file<<endl;
  }
  file.close();
  data.clear();
}

#endif /* MEASUREMENTROUTINES_HPP_ */

