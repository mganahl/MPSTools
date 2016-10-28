#ifndef INTERFACE_HPP__
#define INTERFACE_HPP__
#include "tensormap.hpp"
#include "array.hpp"
#include "boost_typedefs.hpp"
#include "keyTypes.hpp"
#include <iostream>
using namespace BoostTypedefs;



template<size_t rnk> unsigned int conv(const ukey<rnk>& a) { return (unsigned int)(-1); }
template<> unsigned int conv(const ukey<2>& a) { return (a[0] <<16) | a[1]; }


template<size_t rnk> ukey<rnk> iconv(const unsigned int i) { return ukey<rnk>(); }
template<> ukey<2> iconv<2>(const unsigned int i) { return ukey<2>(i>>16, i&( (0x1<<16)-1) ); }

typedef tensor_map<ukey<2>,2,Real> atype;
typedef arrayclass::array<atype> mpstype;

//interface between iztok nrg and my dmrg

class interface{
public:
  interface(const i_to_key<AbKeyType> itokey):_itokey(itokey){};
  ~interface(){};
  void readnrg(const string filename);
  void writenrg(const string filename);
  void readdmrg(const string filename);
  void writedmrg(const string filename);
  MPS<AbKeyType,Real> &mps_dmrg(){return _mps_dmrg;}
  mpstype mps_nrg()const {return _mps_nrg;}
  void transfernrgtodmrg();
  void transferdmrgtonrg();
  void setmpsDMRG(const MPS<AbKeyType,Real> &mps){
    _mps_dmrg=mps;
  }
private:
  MPS<AbKeyType,Real> _mps_dmrg;
  mpstype _mps_nrg;
  i_to_key<AbKeyType>  _itokey;
};

void interface::readnrg(const string filename){
  cout<<"#Reading nrg-mps ..." <<endl;
  if (filename.substr(filename.find_last_of(".") + 1) != "mps"){
    cout<<"readnrg: unkown filetype"<<endl;
    return ;
  }

  std::ifstream bfile(filename.c_str(), std::ios::in | std::ios::binary);
  _mps_nrg.read(bfile);
  cout<<"#Read successfull:"<<endl;
  cout<<"#\t Got an nrg-mps of n = "<<_mps_nrg.size()<<endl;

}

void interface::writenrg(const string filename){
  string fname=filename;
  if (filename.substr(filename.find_last_of(".") + 1) != "mps")
    fname=filename+".mps";
  std::ofstream bfile(fname.c_str(), std::ios::out | std::ios::binary);
  _mps_nrg.write(bfile);
  bfile.close();
}

void interface::readdmrg(const string filename){
  _mps_dmrg.load(filename,false);
}

void interface::writedmrg(const string filename){
  _mps_dmrg.save(filename);
}

void interface::transfernrgtodmrg(){
  cout<<"#starting transcription from nrg tp dmrg mps..."<<endl;
  if(!_mps_nrg.size()){
    cout<<"in interface: got an empty _mps_nrg"<<endl;
    return;
  }
  const int N=_mps_nrg.size();
  _mps_dmrg.resize(N);
  for(uint site=0;site<N;++site){
    for(atype::iterator a=_mps_nrg[site].begin();a!=_mps_nrg[site].end();++a){
      matrix<Real,column_major> mat;
      mat.resize(a->second.M(0),a->second.M(1));
      memcpy(&mat.data()[0],a->second.raw(),sizeof(Real)*a->second.effsize());

      ukey<2> qin = iconv<2>(a->first[0]);
      ukey<2> qout = iconv<2>(a->first[1]);
      ukey<2> diff=iconv<2>(a->first[1]-a->first[0]);
      AbKeyType loc(diff[0],diff[1],2);
      i_to_key<AbKeyType>::iterator it;
      uint locind;
      for(it = _itokey.begin();it!=_itokey.end();++it){
	if(it->second==loc){
	  locind=it->first;
	  break;
	}
      }
      if(it==_itokey.end())
	cout<<"difference of outgoing and ingoing qn not found in itokey"<<endl;
      AbKeyType kyin(qin[0],qin[1],2);
      AbKeyType kyout(qout[0],qout[1],2);
      TensorKeyType<2,1,AbKeyType> ky(kyin,kyout,it->first);
      
      if(_mps_dmrg[site].find(ky)==_mps_dmrg[site].end())
	_mps_dmrg[site][ky]=mat;
      else{
	cout<<"the key "<<ky<<" has already been inserted ..."<<endl;
	abort();
      }
    }
    _mps_dmrg.setItoKey(site,_itokey);
  }
  _mps_dmrg.generateBoundary();
  cout<<"#transcription successful:"<<endl;
  cout<<"#\t Got a dmrg-mps of n = "<<_mps_dmrg.size()<<endl;
  cout<<"#clearing nrg-mps ..."<<endl;
  for(uint site=0;site<N;++site){
    _mps_nrg[site].clear();
  }
}

void interface::transferdmrgtonrg(){
  cout<<"#starting transcription from dmrg tp nrg mps..."<<endl;
  if(!_mps_dmrg.size()){
    cout<<"in interface: got an empty _mps_dmrg"<<endl;
    return;
  }
  const int N=_mps_dmrg.size();
  _mps_nrg.resize(N);
  
  for(uint site=0;site<N;++site){
    for(BlockMatrix<1,AbKeyType,Real>::iterator a=_mps_dmrg[site].begin();a!=_mps_dmrg[site].end();++a){
      tensor<2,Real> mat;
      mat.resize(ukey<2>(a->second.size1(),a->second.size2()));
      memcpy(mat.raw(),&a->second.data()[0],sizeof(Real)*mat.effsize());
      ukey<2>kin(a->first[0](0),a->first[0](1));
      ukey<2>kout(a->first[1](0),a->first[1](1));
      ukey<2> ky(conv<2>(kin),conv<2>(kout));
      if(_mps_nrg[site].find(ky)==_mps_nrg[site].end()){
	_mps_nrg[site][ky]=mat;
      }else{
	cout<<"the key "<<ky<<" has already been inserted ..."<<endl;
	abort();
      }
    }
  }
  cout<<"#transcription successful:"<<endl;
  cout<<"#\t Got a nrg-mps of n = "<<_mps_dmrg.size()<<endl;
  cout<<"#clearing dmrg-mps ..."<<endl;
  _mps_dmrg.clear();
}


#endif
