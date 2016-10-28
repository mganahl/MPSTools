/*
This file contains all simulation container classes (dmrg, tebd, ...).
 */
#ifndef CONTAINERS_HPP__
#define CONTAINERS_HPP__
#include "MPSTypes.hpp"
#include "MPOTypes.hpp"
#include "MPSfunctions.hpp"
#include "utilities.hpp"
#include <list>
#include "measurementClasses.hpp"
#include "lanzcos.hpp"
#include "bicgstab.hpp"

template<class KeyType,typename T>
class Container{
public:
  Container(){_checkpointStep=50;_simName="Container_default_name";_verbose=1;};
  Container(std::string simname){_checkpointStep=50;_verbose=1;_simName=simname;};
  Container(const MPS<KeyType,T> &mps ,MPO<T> & mpo,const std::string simName):_mps(mps),_mpo(mpo),_simName(simName){
    _mpo.construct(0);//the mpo is constructed by default
    this->_verbose=1;
    this->_checkpointStep=50;
  }
  Container(const MPS<KeyType,T> &mps ,MPO<T> & mpo,const std::string simName,const uint verbose):_mps(mps),_mpo(mpo),_simName(simName){
    _mpo.construct(0);//the mpo is constructed by default
    this->_verbose=verbose;
    this->_checkpointStep=50;
  }
  Container(const Container<KeyType,T>&c){
    _mps=c.getMPS();
    _cmps=c.getCMPS();
    _mpo=c.getMPO();
    _simName=c.getName();
    _checkpointStep=c.getCheckpointStep();
    _verbose=c.verbose();
  };
  ~Container(){};
  
  //get stuff
  uint verbose()const{return _verbose;}
  void verbose(const uint v){_verbose=v;}
  void setMPS(const MPS<KeyType,T> &mps){_mps=mps;};
  MPS<KeyType,T> getMPS()const{return _mps;};
  MPS<KeyType,T>&getMPS(){return _mps;}
  CanonizedMPS<KeyType,T> getCMPS()const{return _cmps;};
  CanonizedMPS<KeyType,T>&getCMPS(){return _cmps;}

  void setMPO(const MPO<T> &mpo){_mpo=mpo;};
  MPO<T> getMPO()const{return _mpo;};
  MPO<T> &getMPO(){return _mpo;};
  void setName(const string name){_simName=name;}
  string getName()const {return _simName;}
  void setCheckpointStep(int checkpointStep){this->_checkpointStep=checkpointStep;}
  int getCheckpointStep()const{return this->_checkpointStep;}
  //simulation routines
  T simulate();
  void measureCMPS(list<measurementObject<T> > &m)const;
  void measureMPS(list<measurementObject<T> > &m)const;
  void writeMeasurements(list<measurementObject<T> > &m);
  void orthonormalizeMPS(const string dir){orthonormalize(_mps,dir);}
  void canonizeMPS(){_cmps = canonize(_mps);}
  virtual void clear(){this->_mps.clear();_cmps.clear();_mpo.clear();}
  void save();//save the _mps and _cmps
  void save(std::string name);//save the _mps and _cmps
  void load(bool cast = false);
  void read(std::string);
  void write(std::string);//save to a .container file
  void read();
  void write();

protected:
  virtual void toarchive(boost::archive::binary_oarchive &oar); //in derived classes, these two have to be modified
  virtual void fromarchive(boost::archive::binary_iarchive &iar);//

  void doCMPSmeasurement(list<measurementObject<T> > &m)const;
  void doMPSmeasurement(list<measurementObject<T> > &m)const;
  MPS<KeyType,T>_mps;
  CanonizedMPS<KeyType,T>_cmps;
  MPO<T> _mpo;
  int _checkpointStep;
  std::string _simName;
  uint _verbose;

};

template<class KeyType,typename T>
void Container<KeyType,T>::toarchive(boost::archive::binary_oarchive &oar){
  while(!checkMPS(this->_mps));
  oar<<_mps;
  oar<<_cmps;
  oar<<_mpo;
  oar<<_checkpointStep;
  oar<<_simName;
}

template<class KeyType,typename T>
void Container<KeyType,T>::fromarchive(boost::archive::binary_iarchive &iar){
  _mps.clear();_cmps.clear();_mpo.clear();
  iar>>_mps;
  iar>>_cmps;
  iar>>_mpo;
  iar>>_checkpointStep;
  iar>>_simName;
}

template<class KeyType,typename T>
void Container<KeyType,T>::write(std::string fname){
  std::string name=fname+".container";
  cout<<endl;
  cout<<"#writing data to "<<name<<endl;
  ofstream os(name.c_str(), ios::binary);
  boost::archive::binary_oarchive oar(os);
  toarchive(oar);
  os.close();
}

template<class KeyType,typename T>
void Container<KeyType,T>::read(std::string fname){
  std::string name=fname+".container";
  cout<<"#reading data from "<<name<<endl;
  ifstream is(name.c_str(), ios::binary);
  boost::archive::binary_iarchive iar(is);
  fromarchive(iar);
  is.close();
}

template<class KeyType,typename T>
void Container<KeyType,T>::write(){
  std::string name=this->_simName+".container";
  cout<<"#writing data to "<<name<<endl;
  ofstream os(name.c_str(), ios::binary);
  boost::archive::binary_oarchive oar(os);
  toarchive(oar);
  os.close();
}

template<class KeyType,typename T>
void Container<KeyType,T>::read(){
  std::string name=this->_simName+".container";
  cout<<"#reading data from "<<name<<endl;
  ifstream is(name.c_str(), ios::binary);
  boost::archive::binary_iarchive iar(is);
  fromarchive(iar);
  is.close();
}


template<class KeyType,typename T>
void Container<KeyType,T>::save(){
  _mps.save(_simName);
  _cmps.save(_simName+"c");
}

template<class KeyType,typename T>
void Container<KeyType,T>::save(std::string name){
  _cmps.save(name+"c");
}

template<class KeyType,typename T>
void Container<KeyType,T>::load(bool cast){
  _mps.load(_simName,cast);
  _cmps.load(_simName+"c",cast);
}

//this function does all measurements
template<class KeyType,typename T>
void Container<KeyType,T>::doCMPSmeasurement(list<measurementObject<T> > &m)const{
  if(!_cmps.size()){
    cout<<"Container<KeyType,T>::measure(list<measurementObject<T>): empty CanonizedMPS<KeyType,T> found"<<endl;
    return ;
    
  }
  typename list<measurementObject<T> >::iterator it;
  for(it = m.begin();it!=m.end();it++){
    if(it->step%it->measureStep==0){
      it->measure(this->_cmps);
    }
    it->step++;
  }

}


//this function does all measurements
template<class KeyType,typename T>
void Container<KeyType,T>::measureCMPS(list<measurementObject<T> > &m)const{
  if(!_cmps.size()){
    cout<<"Container<KeyType,T>::measure(list<measurementObject<T>): empty CanonizedMPS<KeyType,T> found"<<endl;
    return ;
  }
  typename list<measurementObject<T> >::iterator it;
  for(it = m.begin();it!=m.end();++it)
    it->measure(this->_cmps);
}


template<class KeyType,typename T>
void Container<KeyType,T>::measureMPS(list<measurementObject<T> > &m)const{
  if(!_mps.size()){
    cout<<"Container<KeyType,T>::measure(list<measurementObject<T>): empty MPS<KeyType,T> found"<<endl;
    return ;
  }
  typename list<measurementObject<T> >::iterator it;
  for(it = m.begin();it!=m.end();++it)
    it->measure(_mps);
}



template<class KeyType,typename T>
void Container<KeyType,T>::doMPSmeasurement(list<measurementObject<T> > &m)const{
  if(!_mps.size()){
    cout<<"Container<KeyType,T>::measure(list<measurementObject<T>): empty MPS<KeyType,T> found"<<endl;
    return ;
  }
  typename list<measurementObject<T> >::iterator it;
  for(it = m.begin();it!=m.end();it++){
    if(it->step%it->measureStep==0){
      it->measure(_mps);
    }
    it->step++;
  }

}

template<class KeyType,typename T>
void Container<KeyType,T>::writeMeasurements(list<measurementObject<T> > &m){
  typename list<measurementObject<T> >::iterator it;
  for(it=m.begin(); it!=m.end(); it++){
    it->writeMeasurements();
  }
}


template<class KeyType,typename T>
class DMRGengine:public Container<KeyType,T>{
public:
  DMRGengine():Container<KeyType,T>(){_lanmax=500;}
  DMRGengine(const Container<KeyType,T>&c):Container<KeyType,T>(c){_lanmax=500;};
  DMRGengine(std::string simName):Container<KeyType,T>(simName){_lanmax=500;};
  DMRGengine(const MPS<KeyType,T> &mps ,MPO<T> & mpo,const std::string simName):Container<KeyType,T>(mps,mpo,simName){
    _lanmax=500;}
  DMRGengine(const MPS<KeyType,T> &mps ,MPO<T> & mpo,const std::string simName,const uint verbose):Container<KeyType,T>(mps,mpo,simName,verbose){
    _lanmax=500;}
  //delta is the Lanzcos delta ...
  //targetting is NOT working, use only targetstate=0!!!!!
  /*parameters: 
  maxsteps= steps of DMRG
  EnConv: if energy aftert two successice dmrg (lr+rl) sweeps does not differ by more then EnConv, then DMRG is stopped
  chimax: clear
  doTrunc: if true, then only eigenvalues of the density-matrix which are larger than tr_weight are kept, if this is possible within the specified chimax; if false, then chimax eigenvalues are kept
  tr_weight: truncated weight, see doTrunc desciption
  delta: lanzcos delta; if the norm of the new lanzcos vector falls below delta, then lanzcos is stopped
  testenergy: if true, then lanzcos takes 


  */
  T simulate(const uint maxsteps=uint(5),const Real EnConv=1e-10,const uint chimax=100,const bool doTrunc=false, const Real tr_weight=1e-12,
	     const Real delta=1e-13,const bool testenergy=true,const bool leaveoutlast=false,const uint targetstate=0);

  T simulate(const MPS<KeyType,T> &mps,const T delta_ortho,const uint maxsteps=uint(5),const Real EnConv=1e-10,const uint chimax=100,const bool doTrunc=false,
	     const Real tr_weight=1e-12,const Real delta=1e-13,const bool testenergy=true,const bool leaveoutlast=false);

  virtual void clear(){Container<KeyType,T>::clear();this->_R.clear();this->_L.clear();}
protected:
  virtual void toarchive(boost::archive::binary_oarchive &oar){
    Container<KeyType,T>::toarchive(oar);
    oar<<_R;
    oar<<_L;
    oar<<_lanmax;
  }
  virtual void fromarchive(boost::archive::binary_iarchive &iar){
    Container<KeyType,T>::fromarchive(iar);
    iar>>_R;
    iar>>_L;
    iar>>_lanmax;
  }

  map_container<int,BlockMatrix<1,KeyType,T> > _R;
  map_container<int,BlockMatrix<1,KeyType,T> > _L;
  int _lanmax;
private:
  virtual triple<T,Real,uint> 
  doStep(const int site,const int dir,const uint chimax,const bool doTrunc, const Real tr_weight,const Real delta, const bool testenergy,const uint targetstate);

  virtual triple<T,Real,uint> 
  doStepExcited(const int site,const int dir,const uint chimax,const bool doTrunc, const Real tr_weight,const Real delta, const bool testenergy,
		BlockMatrix<2,KeyType,T> &state,const T delta_ortho);
};

template<class KeyType,typename T>
triple<T,Real,uint> DMRGengine<KeyType,T>:: doStep(const int site,const int dir,const uint chimax,const bool doTrunc, const Real tr_weight,
						   const Real delta, const bool testenergy,const uint targetstate){
  //initialize the engine per constructor ...
  TSLanzcosEngine<KeyType,T> 
    lan(_L[site-1],_R[site+2],this->_mps[site],this->_mps[site+1],this->_mpo[site],this->_mpo[site+1],this->_mps.getItoKey(site),
	this->_mps.getItoKey(site+1),this->_verbose);
  //and run it ...
//  cout<<"site: "<<site<<endl;
//  cout<<_L[site-1]<<endl;
//  cout<<_R[site+2]<<endl;
//  cout<<this->_mps.getItoKey(site+1).size()<<endl;
//  cout<<endl;  cout<<endl;
//  
  bool debug=false;
  lan.simulate(_lanmax,delta,testenergy,debug);
  //now get the groundstate (other states may be targeted as well) ...
  
  //BlockMatrix<2,KeyType,T>A=lan.getGroundState();
  BlockMatrix<2,KeyType,T>A=lan.getState(targetstate);
  std::vector<BlockMatrix<1,KeyType,T> >Aout;
  std::vector<DiagBlockMatrix<KeyType,Real> >lambda;
  std::vector<i_to_key<KeyType> > itokey;
  itokey.push_back(this->_mps.getItoKey(site));
  itokey.push_back(this->_mps.getItoKey(site+1));
  //  cout<<"the size of the 2-site-MPS matrix at site "<<site<<" is "<<A.size()<<endl;
  Aout.clear();
  lambda.clear();
  pair<std::vector<Real>,std::vector<uint> >sre=split(A,Aout,lambda, itokey,0,chimax,doTrunc,tr_weight);
  if(dir>0){
    this->_mps[site]=Aout[0];
    this->_mps[site+1]=lambda[0]*Aout[1];
  }
  if(dir<=0){
    this->_mps[site]=Aout[0]*lambda[0];
    this->_mps[site+1]=Aout[1];
  }

  if(dir>0)
    _L[site] = addContractionLayer(_L[site-1],this->_mps[site],this->_mpo[site],this->_mps.getItoKey(site),1);
  if(dir<=0)
    _R[site+1] = addContractionLayer(_R[site+2],this->_mps[site+1],this->_mpo[site+1],this->_mps.getItoKey(site+1),-1);

  triple<T,Real,uint> out(lan.E(0),sre.first[0],sre.second[0]);
  return out;
}

template<class KeyType,typename T>
triple<T,Real,uint> DMRGengine<KeyType,T>:: doStepExcited(const int site,const int dir,const uint chimax,const bool doTrunc, const Real tr_weight,
							  const Real delta, const bool testenergy,BlockMatrix<2,KeyType,T> &state,const T delta_ortho){
  //initialize the engine per constructor ...
  TSLanzcosEngine<KeyType,T> 
    lan(_L[site-1],_R[site+2],this->_mps[site],this->_mps[site+1],this->_mpo[site],this->_mpo[site+1],this->_mps.getItoKey(site),this->_mps.getItoKey(site+1),
	this->_verbose);
  //and run it ...
  lan.simulate(_lanmax,delta,testenergy,state,delta_ortho);
  
  //now get the groundstate (other states may be targeted as well) ...
  //BlockMatrix<2,KeyType,T>A=lan.getGroundState();
  BlockMatrix<2,KeyType,T>A=lan.getState(0);

  std::vector<BlockMatrix<1,KeyType,T> >Aout;
  std::vector<DiagBlockMatrix<KeyType,Real> >lambda;
  std::vector<i_to_key<KeyType> > itokey;
  itokey.push_back(this->_mps.getItoKey(site));
  itokey.push_back(this->_mps.getItoKey(site+1));
  //  cout<<"the size of the 2-site-MPS matrix at site "<<site<<" is "<<A.size()<<endl;
  Aout.clear();
  lambda.clear();
  pair<std::vector<Real>,std::vector<uint> >sre=split(A,Aout,lambda, itokey,0,chimax,doTrunc,tr_weight);
  if(dir>0){
    this->_mps[site]=Aout[0];
    this->_mps[site+1]=lambda[0]*Aout[1];
  }
  if(dir<=0){
    this->_mps[site]=Aout[0]*lambda[0];
    this->_mps[site+1]=Aout[1];
  }

  if(dir>0)
    _L[site] = addContractionLayer(_L[site-1],this->_mps[site],this->_mpo[site],this->_mps.getItoKey(site),1);
  if(dir<=0)
    _R[site+1] = addContractionLayer(_R[site+2],this->_mps[site+1],this->_mpo[site+1],this->_mps.getItoKey(site+1),-1);

  triple<T,Real,uint> out(lan.E(0),sre.first[0],sre.second[0]);
  return out;
}


//returns the (converged) energy
template<class KeyType,typename T>
T DMRGengine<KeyType,T>::simulate(const uint maxsteps,const Real EnConv,const uint chimax,const bool doTrunc, const Real tr_weight,
				  const Real delta, const bool testenergy,const bool leaveoutlast,const uint targetstate){
  assert(EnConv>0);
  orthonormalize(this->_mps,"rl");
  assert(this->_mps.getROrtho());
  computeR(this->_mps,this->_mpo,2,_R);

  _L.clear();
  _L[-1] = getLMin(this->_mps);
  triple<T,Real,uint> out;
  T globalE=T(1e100);
  uint step=0;
  bool converged = false;
  const uint last = this->_mps.size()-2;
  while(!converged){
    for(int site=0;site<this->_mps.size()-2;++site){
      if(this->_verbose==1)
	cout<<"lanczos at site "<<site<<" at step "<<step<<"; ";
      assert(site<this->_mps.size()-1);
      out=doStep(site,1,chimax,doTrunc,tr_weight,delta,testenergy,targetstate);
      if(this->_verbose==1)
	cout<<" tw = "<<out.second<<" at chi = "<<out.third<<endl;
    }
    for(int site=this->_mps.size()-2;site>0;--site){
      if(leaveoutlast)
	if(site==last)
	  continue;
      if(this->_verbose==1)
	cout<<"lanczos at site "<<site<<" at step "<<step<<"; ";
      assert(site>=0);
      out=doStep(site,-1,chimax,doTrunc,tr_weight,delta,testenergy,targetstate);
      if(this->_verbose==1)
	cout<<" tw = "<<out.second<<" at chi = "<<out.third<<endl;
    }
    ++step; 
    if(abs(globalE-out.first)<EnConv){
      cout<<"DMRG converged to "<<EnConv<<". leaving DMRG loop"<<endl;//that should always be printed ...
      converged=true;
    }
    if(step>=maxsteps){
      cout<<"in DMRG: reached maximum iteration number "<<maxsteps<<". Energy converged to "<<abs(globalE-out.first)<<". leaving DMRG loop"<<endl;
      converged=true;
    }
    globalE=out.first;
  }

  out=doStep(0,-1,chimax,doTrunc,tr_weight,delta,testenergy,targetstate);
  if(this->_verbose==1)
    cout<<" tw = "<<out.second<<" at chi = "<<out.third<<endl;

  this->_mps.setLOrtho(false);
  this->_mps.setROrtho(true);
  return out.first;
  //returns the GS-energy

}


//returns the (converged) energy
template<class KeyType,typename T>
T DMRGengine<KeyType,T>::simulate(const MPS<KeyType,T> &mps,const T delta_ortho,const uint maxsteps,const Real EnConv,const uint chimax,const bool doTrunc, 
				  const Real tr_weight,const Real delta, const bool testenergy,const bool leaveoutlast){
  assert(EnConv>0);

  const int N=mps.size();
  if(mps.size()!=this->_mps.size()){
    cout<<"in T DMRGengine<KeyType,T>::simulate(const MPS<KeyType,T> &mps,const uint maxsteps,const Real EnConv,const uint chimax,const bool doTrunc,const Real tr_weight,const Real delta, const bool testenergy,const bool leaveoutlast): mps.size()!=this->_mps.size()"<<endl;
    abort();
  }
  orthonormalize(this->_mps,"rl");
  assert(this->_mps.getROrtho());
  computeR(this->_mps,this->_mpo,2,_R);  
  _L.clear();
  _L[-1] = getLMin(this->_mps);
  triple<T,Real,uint> out;
  T globalE=T(1e100);
  uint step=0;
  bool converged = false;

  //$$$$$$$$$$$$$$$$$$$$$$$$$$first, transform mps into basis used in current DMRG simulation:$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  //collect all things needed to do basis transformation:

  IdentityMPO<T> idMPO(mps);
  std::vector<i_to_key<KeyType> >itokey;
  map_container<int,BlockMatrix<1,KeyType,T> > overlapR,overlapL;

  //check for consistency ...
  if(mps.getQN()!=this->_mps.getQN()){
    cout<<"in T DMRGengine<KeyType,T>::simulate(const MPS<KeyType,T> &mps,const uint maxsteps,const Real EnConv,const uint chimax,const bool doTrunc, const Real tr_weight,const Real delta, const bool testenergy,const bool leaveoutlast): mps.getQN()!=this->_mps.getQN()"<<endl;
    abort();
  }
  //compute the E-matrices ...
  overlapR.clear();
  overlapR[N]= getRMax(mps);
  for(int s = N-1;s>=0;--s){
    assert(s>=0);
    //overlapR[s] = addContractionLayer(overlapR[s+1],this->_mps[s],idMPO[s],mps[s],mps.getItoKey(s),0);
    overlapR[s] = addContractionLayer(overlapR[s+1],mps[s],idMPO[s],this->_mps[s],mps.getItoKey(s),0);
  }
  overlapL[-1] = getLMin(mps);

  const uint last = this->_mps.size()-2;
  while(!converged){
    for(int site=0;site<this->_mps.size()-2;++site){
      if(this->_verbose==1)
	cout<<"lanczos at site "<<site<<" at step "<<step<<"; ";
      assert(site<this->_mps.size()-1);

      //transform mps into current local basis
      itokey.clear();
      itokey.push_back(mps.getItoKey(site));
      itokey.push_back(mps.getItoKey(site+1));
      BlockMatrix<2,KeyType,T>state_oldbasis(mps[site],mps[site+1]),state_newbasis;
      BlockMatrix<2,KeyType,T>tempstate(this->_mps[site],this->_mps[site+1]);
      contractTwosite(state_newbasis,overlapL[site-1],overlapR[site+2],state_oldbasis,idMPO[site],idMPO[site+1],itokey);
//      cout<<"from contraction: "<<state_newbasis*tempstate<<endl;
//      cout<<"from overlap: "<<computeOverlap(mps,this->_mps)<<endl;
//      cin.get();
//      

      //basis-transformation is done 

      out=doStepExcited(site,1,chimax,doTrunc,tr_weight,delta,testenergy,state_newbasis,delta_ortho);
      if(this->_verbose==1)
	cout<<" tw = "<<out.second<<" at chi = "<<out.third<<endl;

      //add a layer to the L-blocks of the overlap-E-matrices.
      //overlapL[site] = addContractionLayer(overlapL[site-1],this->_mps[site],idMPO[site],mps[site],mps.getItoKey(site),1);
      overlapL[site] = addContractionLayer(overlapL[site-1],mps[site],idMPO[site],this->_mps[site],mps.getItoKey(site),1);
      
    }
    for(int site=this->_mps.size()-2;site>0;--site){
      if(leaveoutlast)
	if(site==last)
	  continue;
      if(this->_verbose==1)
	cout<<"lanczos at site "<<site<<" at step "<<step<<"; ";
      assert(site>=0);

      //transform mps into current local basis
      itokey.clear();
      itokey.push_back(mps.getItoKey(site));
      itokey.push_back(mps.getItoKey(site+1));
      BlockMatrix<2,KeyType,T>state_oldbasis(mps[site],mps[site+1]),state_newbasis;
      contractTwosite(state_newbasis,overlapL[site-1],overlapR[site+2],state_oldbasis,idMPO[site],idMPO[site+1],itokey);
      //basis-transformation is done 

      out=doStepExcited(site,-1,chimax,doTrunc,tr_weight,delta,testenergy,state_newbasis,delta_ortho);
      if(this->_verbose==1)
	cout<<" tw = "<<out.second<<" at chi = "<<out.third<<endl;

      //add a layer to the R-blocks of the overlap-E-matrices.

      //overlapR[site+1] = addContractionLayer(overlapR[site+2],this->_mps[site+1],idMPO[site+1],mps[site+1],mps.getItoKey(site+1),-1);
      overlapR[site+1] = addContractionLayer(overlapR[site+2],mps[site+1],idMPO[site+1],this->_mps[site+1],mps.getItoKey(site+1),-1);
    }
    
    ++step; 
    if(abs(globalE-out.first)<EnConv){
      cout<<"DMRG converged to "<<EnConv<<". leaving DMRG loop"<<endl;//that should always be printed ...
      converged=true;
    }
    if(step>=maxsteps){
      cout<<"in DMRG: reached maximum iteration number "<<maxsteps<<". Energy converged to "<<abs(globalE-out.first)<<". leaving DMRG loop"<<endl;
      converged=true;
    }
    globalE=out.first;
  }

  //now do the two leftmost sites:
  uint site=0;
  itokey.clear();
  itokey.push_back(mps.getItoKey(site));
  itokey.push_back(mps.getItoKey(site+1));
  BlockMatrix<2,KeyType,T>state_oldbasis(mps[site],mps[site]),state_newbasis;
  contractTwosite(state_newbasis,overlapL[site-1],overlapR[site+2],state_oldbasis,idMPO[site],idMPO[site+1],itokey);

  out=doStepExcited(site,-1,chimax,doTrunc,tr_weight,delta,testenergy,state_newbasis,delta_ortho);
  if(this->_verbose==1)
    cout<<" tw = "<<out.second<<" at chi = "<<out.third<<endl;

  this->_mps.setLOrtho(false);
  this->_mps.setROrtho(true);
  return out.first;
  
  //done: returns the GS-energy

}


//for now use only T=Complex
//it's not working ....
template<class KeyType,typename T>
class DMPSengine:public Container <KeyType,T>{
public:
  DMPSengine():Container<KeyType,T>(){_Nmaxbicg=1000;}
  DMPSengine(const string simname){_Nmaxbicg=1000;this->_simName=simname;this->_verbose=1;this->_checkpointStep=10;}
  DMPSengine(const string simname,const uint verbose){_Nmaxbicg=1000;this->_simName=simname;this->_verbose=verbose;this->_checkpointStep=10;}
  DMPSengine(const MPS<KeyType,T> &mps,const MPS<KeyType,T> &mps_ihg,MPO<T> & mpo,const std::string simName,const uint verbose)
    :Container<KeyType,T>(mps,mpo,simName,verbose){
    assert(mps.size()==mps_ihg.size());
    _mps_ihg=mps_ihg;
    orthonormalize(this->_mps,"rl");
    computeR(this->_mps,this->_mpo,2,_R);
    _L[-1] = getLMin(this->_mps);
    for(uint site=0;site<this->_mps.size();++site){
      MPOMat<1,T> ID;
      SparseOperator<1,T> id(_mps_ihg.getPSpace()(site),this->_mps.getPSpace()(site));
      for(uint ind=0;ind<_mps_ihg.getPSpace()(site);++ind)
	id[OpKeyType<1>(ind,ind)]= T(1.);
      ID[ukey<2>(0,0)]=id;
      ID.Dl=1;ID.Dr=1;
      _mpo_id[site]=ID;
    }
  }

  T simulate(const uint maxsweeps,const uint Nmaxbicg,const uint chimax,const bool doTrunc, const Real tr_weight,const Real bicgtol,const Real fidconv);
  MPS<KeyType,T> getMPSinhom()const {return _mps_ihg;}
protected:
  map_container<int,BlockMatrix<1,KeyType,T> > _R,_L,_R_ihg,_L_ihg;
  virtual void toarchive(boost::archive::binary_oarchive &oar); //in derived classes, these two have to be modified
  virtual void fromarchive(boost::archive::binary_iarchive &iar);//

private:
  MPS<KeyType,T> _mps_ihg;
  MPO<T>_mpo_id;
  virtual quadruple<Real,uint,uint,T> doStep(const BlockMatrix<2,KeyType,T> &mps0,const int site,const int dir,const uint chimax,const bool doTrunc, 
					     const Real tr_weight,const Real tol,const bool normalize);
  int _Nmaxbicg;
};


template<class KeyType,typename T>
void DMPSengine<KeyType,T>::toarchive(boost::archive::binary_oarchive &oar){
  Container<KeyType,T>::toarchive(oar);
  oar<<_R;
  oar<<_L;
  oar<<_R_ihg;
  oar<<_L_ihg;
  oar<<_mps_ihg;
  oar<<_mpo_id;
  oar<<_Nmaxbicg;
}


template<class KeyType,typename T>
void DMPSengine<KeyType,T>::fromarchive(boost::archive::binary_iarchive &iar){
  _R.clear();_L.clear();_R_ihg.clear();_L_ihg.clear();_mps_ihg.clear();_mpo_id.clear();
  Container<KeyType,T>::fromarchive(iar);
  iar>>_R;
  iar>>_L;
  iar>>_R_ihg;
  iar>>_L_ihg;
  iar>>_mps_ihg;
  iar>>_mpo_id;
  iar>>_Nmaxbicg;
}


template<class KeyType,typename T>
quadruple<Real,uint,uint,T> DMPSengine<KeyType,T>:: doStep(const BlockMatrix<2,KeyType,T> &mps0,const int site,const int dir,const uint chimax,const bool doTrunc, 
							   const Real tr_weight,const Real tol,const bool normalize){



  Real TOL=1e10;
  //initialize the engine per constructor ...
  cout<<"starting bicgstab"<<endl;
  BiCGStabEngine<KeyType,T> 
    bicgstab(this->_L[site-1],this->_R[site+2],this->_mps[site],this->_mps[site+1],mps0,this->_mpo[site],this->_mpo[site+1],this->_mps.getItoKey(site),
	     this->_mps.getItoKey(site+1));

  //and run it ...
  quadruple<BlockMatrix<2,KeyType,T>,uint,T,int>bicgstabout=bicgstab.simulate(_Nmaxbicg,tol,TOL);
  int info=bicgstabout.fourth;
  uint cnter=0;  
  cout<<"done "<<endl;

  while(info==-1){
    cout<<"######################## at site "<<site<<": bicgstab divergence detected; restarting with random vector #######################"<<endl;
    BlockMatrix<2,KeyType,T>out=bicgstabout.first;
    out.randomize(T(-1.0),T(1.0));
    T norm=sqrt(out*out);
    out*=(1./norm);
    bicgstab.setX0(out);
    bicgstabout=bicgstab.simulate(_Nmaxbicg,tol,TOL);
    info=bicgstabout.fourth;
    cnter++;
    if(cnter>10){
      cout<<"######################################## bicgstab not converging: moving to next site ########################################"<<endl;
      out=bicgstabout.first;
      norm=sqrt(out*out);
      out*=(1./norm);
      bicgstabout.first=out;
      break;
    }
  }
  BlockMatrix<2,KeyType,T> A=bicgstabout.first;

  std::vector<BlockMatrix<1,KeyType,T> >Aout;
  std::vector<DiagBlockMatrix<KeyType,Real> > lambda;
  std::vector<i_to_key<KeyType> > itokey;
  itokey.push_back(this->_mps.getItoKey(site));
  itokey.push_back(this->_mps.getItoKey(site+1));
  cout<<"before split"<<endl;

  pair<std::vector<Real>,std::vector<uint> >sre=split(A,Aout,lambda,itokey,0,chimax,doTrunc,tr_weight,normalize);
  if(dir>0){
    this->_mps[site]=Aout[0];
    this->_mps[site+1]=lambda[0]*Aout[1];
    this->_mps[site].squeeze(1e-14);this->_mps[site+1].squeeze(1e-14);
    //    cout<<"#norm of lambda after bicgstab at site "<<site<<": "<<lambda[0].norm()<<"; "<<endl;
    _L[site] = addContractionLayer(this->_L[site-1],this->_mps[site],this->_mpo[site],this->_mps.getItoKey(site),1);
    _L_ihg[site]=addContractionLayer(_L_ihg[site-1],this->_mps[site],_mpo_id[site],_mps_ihg[site],_mps_ihg.getItoKey(site),1);
    //_L_ihg[site]=addContractionLayer(_L_ihg[site-1],_mps_ihg[site],_mpo_id[site],this->_mps_ihg[site],_mps_ihg.getItoKey(site),1);
  }
  if(dir<=0){
    this->_mps[site]=Aout[0]*lambda[0];
    this->_mps[site+1]=Aout[1];
    this->_mps[site].squeeze(1e-14);this->_mps[site+1].squeeze(1e-14);
    //    cout<<"#norm of lambda after bicgstab at site "<<site<<": "<<lambda[0].norm()<<"; "<<endl;
    this->_R[site+1] = addContractionLayer(this->_R[site+2],this->_mps[site+1],this->_mpo[site+1],this->_mps.getItoKey(site+1),-1);
    _R_ihg[site+1]=addContractionLayer(_R_ihg[site+2],this->_mps[site+1],_mpo_id[site+1],_mps_ihg[site+1],_mps_ihg.getItoKey(site+1),-1);
    //_R_ihg[site+1]=addContractionLayer(_R_ihg[site+2],_mps_ihg[site+1],_mpo_id[site+1],this->_mps[site+1],_mps_ihg.getItoKey(site+1),-1);
  }

  quadruple<Real,uint,uint,T> out(sre.first[0],sre.second[0],bicgstabout.second,bicgstabout.third);
  return out;
}

//returns the (converged) energy

template<class KeyType,typename T>
T DMPSengine<KeyType,T>::simulate(const uint maxsweeps,const uint Nmaxbicg,const uint chimax,const bool doTrunc, const Real tr_weight,
				  const Real tol,const Real fidconv){
  // for(int site = this->_mps.size()-1;site>1;--site){
  //   prepareSite(this->_mps,site,"r",true);
  // }
  //the initial state can be of any norm
  orthonormalize(this->_mps,"lr");
  orthonormalize(this->_mps,"rl");

  assert(this->_mps.getROrtho());

  _Nmaxbicg=Nmaxbicg;
  //get all R-matrices 
  computeR(this->_mps,this->_mpo,2,_R);  
  _R_ihg[this->_mps.size()]=getRMax(this->_mps);
  //get the identity mpo
  for(uint site=this->_mps.size()-1;site>1;--site){
    _R_ihg[site]=addContractionLayer(_R_ihg[site+1],this->_mps[site],_mpo_id[site],_mps_ihg[site],_mps_ihg.getItoKey(site),-1);
    //_R_ihg[site]=addContractionLayer(_R_ihg[site+1],_mps_ihg[site],_mpo_id[site],_mps_ihg[site],this->_mps.getItoKey(site),-1);
  }

  //get all L-matrices 
  _L.clear();
  _L_ihg.clear();
  _L[-1] = getLMin(this->_mps);
  _L_ihg[-1]=getLMin(this->_mps);


  //the EFFECTIVE inhomogeneity vector:
  BlockMatrix<2,KeyType,T>ihg_eff;
  uint step=0;
  bool converged = false;
  quadruple<Real,uint,uint,T> out;
  bool normalize;

  while(!converged){
    MPS<KeyType,T>mps_old=this->_mps;
    for(int site=0;site<this->_mps.size()-2;++site){
      cout<<"site "<<site<<endl;
      BlockMatrix<2,KeyType,T>ihg(_mps_ihg[site],_mps_ihg[site+1]);
      std::vector<i_to_key<KeyType> >itokey;
      itokey.push_back(_mps_ihg.getItoKey(site));
      itokey.push_back(_mps_ihg.getItoKey(site+1));
      //get the effective inhomogeneity vector
      contractTwosite(ihg_eff,_L_ihg[site-1],_R_ihg[site+2],ihg,_mpo_id[site],_mpo_id[site+1],itokey);
      if(this->_verbose==1)
	cout<<"bicgstab at site "<<site<<" at step "<<step<<"; ";
      assert(site<this->_mps.size()-1);

      //       normalize=(step==0);
      //       normalize=false;
      normalize=true;
      out=doStep(ihg_eff,site,1,chimax,doTrunc,tr_weight,tol,normalize);
      if(this->_verbose==1)
	cout<<"tw = "<<out.first<<" at chi = "<<out.second<<"; bicgstab iterations "<<out.third<<" at residual "<<out.fourth<<endl;
    }
    for(int site=this->_mps.size()-2;site>0;--site){
      BlockMatrix<2,KeyType,T>ihg(_mps_ihg[site],_mps_ihg[site+1]);
      std::vector<i_to_key<KeyType> >itokey;
      itokey.push_back(_mps_ihg.getItoKey(site));
      itokey.push_back(_mps_ihg.getItoKey(site+1));
      //get the effective inhomogeneity vector
      contractTwosite(ihg_eff,_L_ihg[site-1],_R_ihg[site+2],ihg,_mpo_id[site],_mpo_id[site+1],itokey);
      if(this->_verbose==1)
	cout<<"bicgstab at site "<<site<<" at step "<<step<<"; ";
      assert(site>=0);
      //       normalize=false;
      normalize=true;
      //normalize=((step==0)&&(site==this->_mps.size()-2));

      out=doStep(ihg_eff,site,-1,chimax,doTrunc,tr_weight,tol,normalize);
      if(this->_verbose==1)
	cout<<"tw = "<<out.first<<" at chi = "<<out.second<<"; bicgstab iterations "<<out.third<<" at residual "<<out.fourth<<endl;
    }
    T fidelity=T(1.0)-T(computeOverlap(mps_old,this->_mps))/(T(computeOverlap(mps_old,mps_old)));
    if(this->_verbose==1)
      cout<<"#############################################fidelity at step "<<step<<": "<<fidelity<<"#############################################"<<endl;

    if (step%this->_checkpointStep==0)
      this->write();

    if(abs(fidelity)<fidconv){
      cout<<"###################################in DMPSengine: fidelity converged to "<<fidelity<<". exiting DMPSengine########################"<<endl;
      converged=true;
    }
    if(step>maxsweeps){
      cout<<"in DMPSengine: reached maximum sweep number "<<maxsweeps<<". exiting DMPSengine"<<endl;
      break;
    }
    ++step; 
  }

  BlockMatrix<2,KeyType,T>ihg(_mps_ihg[0],_mps_ihg[1]);
  std::vector<i_to_key<KeyType> >itokey;
  itokey.push_back(_mps_ihg.getItoKey(0));
  itokey.push_back(_mps_ihg.getItoKey(1));
  //get the effective inhomogeneity vector
  contractTwosite(ihg_eff,_L_ihg[-1],_R_ihg[2],ihg,_mpo_id[0],_mpo_id[1],itokey);
  out=doStep(ihg_eff,0,-1,chimax,doTrunc,tr_weight,tol,false);
  if(this->_verbose==1)
    cout<<"tw = "<<out.first<<" at chi = "<<out.second<<"; bicgstab iterations "<<out.third<<" at residual "<<out.fourth<<endl;

  this->_mps.setLOrtho(false);
  this->_mps.setROrtho(true);
  //returns the GS-energy
  return out.first;
}



template<class KeyType,typename T>
class SSDMRGengine:public Container<KeyType,T>{
public:
  SSDMRGengine():Container<KeyType,T>(){_lanmax=500;this->_simName="SingleSiteDMRG_default_name";}
  SSDMRGengine(std::string simname):Container<KeyType,T>(simname){_lanmax=500;}
  SSDMRGengine(const Container<KeyType,T>&c):Container<KeyType,T>(c){_lanmax=500;this->_simName="SingleSiteDMRG_default_name";};
  SSDMRGengine(const MPS<KeyType,T> &mps ,MPO<T> & mpo,const std::string simName):Container<KeyType,T>(mps,mpo,simName){
    this->_checkpointStep=50;_lanmax=500;}
  SSDMRGengine(const MPS<KeyType,T> &mps ,MPO<T> & mpo,const std::string simName,const uint verbose):Container<KeyType,T>(mps,mpo,simName,verbose){
    this->_checkpointStep=50;_lanmax=500;}

  //delta is the Lanzcos delta ...
  T simulate(const uint maxsteps=uint(5),const Real EnConv=1e-10,const Real delta=1e-12,const bool testenergy=true);
  virtual void clear(){Container<KeyType,T>::clear();this->_R.clear();this->_L.clear();}
  // virtual void read();
  // virtual void write();
  // virtual void read(std::string name);
  // virtual void write(std::string name);

protected:
  virtual void toarchive(boost::archive::binary_oarchive &oar){
    Container<KeyType,T>::toarchive(oar);
    oar<<_R;
    oar<<_L;
    oar<<_lanmax;
  }
  virtual void fromarchive(boost::archive::binary_iarchive &iar){
    Container<KeyType,T>::fromarchive(iar);
    iar>>_R;
    iar>>_L;
    iar>>_lanmax;
  }
  map_container<int,BlockMatrix<1,KeyType,T> > _R;
  map_container<int,BlockMatrix<1,KeyType,T> > _L;
private:
  T doStep(const int site,const int dir,const Real delta, const bool testenergy);
  int _lanmax;
};

template<class KeyType,typename T>
T SSDMRGengine<KeyType,T>:: doStep(const int site,const int dir,const Real delta, const bool testenergy){
  //initialize the engine per constructor ...
  SSLanzcosEngine<KeyType,T> 
    lan(_L[site-1],_R[site+1],this->_mps[site],this->_mpo[site],this->_mps.getItoKey(site),this->_verbose);
  //and run it ...
  lan.simulate(_lanmax,delta,testenergy);
  //now get the groundstate (other states may be targeted as well) ...
  this->_mps[site]=lan.getGroundState();
  //don't forget do normalize it ...
  string direction = dir>0?"l":"r";
  prepareSite(this->_mps,site,direction,true);//multiplies the lambda to the left/right for efficient initial state
  if(dir>0)
    _L[site] = addContractionLayer(_L[site-1],this->_mps[site],this->_mpo[site],this->_mps.getItoKey(site),1);
  if(dir<=0)
    _R[site] = addContractionLayer(_R[site+1],this->_mps[site],this->_mpo[site],this->_mps.getItoKey(site),-1);
  return lan.E(0);
}

//returns the (converged) energy
template<class KeyType,typename T>
T SSDMRGengine<KeyType,T>::simulate(const uint maxsteps,const Real EnConv, const Real delta, const bool testenergy){
  assert(EnConv>0);
  orthonormalize(this->_mps,"rl");
  assert(this->_mps.getROrtho());
  computeR(this->_mps,this->_mpo,1,_R);  
  _L.clear();
  _L[-1] = getLMin(this->_mps);
  T out;
  T globalE=T(1e100);
  uint step=0;
  bool converged = false;
  while(!converged){
    for(int site=0;site<this->_mps.size()-1;++site){
      //      cout<<"lanczos at site "<<site<<" at step "<<step<<"; ";
      assert(site<this->_mps.size()-1);
      out=doStep(site,1,delta,testenergy);
      if(this->_verbose==1)
	cout<<"lanczos at site "<<site<<" at step "<<step<<"; E0= "<<out<<endl;
    }
    
    for(int site=this->_mps.size()-1;site>0;--site){
      //cout<<"lanczos at site "<<site<<" at step "<<step<<"; ";
      assert(site>=0);
      out=doStep(site,-1,delta,testenergy);
      if(this->_verbose==1)
	cout<<"lanczos at site "<<site<<" at step "<<step<<"; E0= "<<out<<endl;
    }
    ++step; 
    if(abs(globalE-out)<EnConv){
      cout<<"DMRG converged to "<<EnConv<<". leaving DMRG loop"<<endl;
      converged=true;
    }
    if(step>maxsteps){
      cout<<"in DMRG: reached maximum iteration number "<<maxsteps<<". Energy converged to "<<abs(globalE-out)<<". leaving DMRG loop"<<endl;
      converged=true;
    }
    globalE=out;
  }
  out=doStep(0,-1,delta,testenergy);
  this->_mps.setLOrtho(false);
  this->_mps.setROrtho(true);
  //returns the GS-energy
  return out;
}


template<class KeyType,typename T>
class TEBDengine:public Container<KeyType,T>{
public:
  TEBDengine():Container<KeyType,T>(){_blockthreads=1;_tebdthreads=1;_twMeasurementStep=1;_nameit=1;_storeallcp=false;_TEBDstep=1;}
  TEBDengine(std::string simname):Container<KeyType,T>(simname)
  {_blockthreads=1;_tebdthreads=1;_twMeasurementStep=1;_nameit=1;_storeallcp=false;_TEBDstep=1;}
  TEBDengine(std::string simname,const bool storeallcp):Container<KeyType,T>(simname)
  {_blockthreads=1;_tebdthreads=1;_twMeasurementStep=1;_nameit=1;_storeallcp=storeallcp;_TEBDstep=1;}

  TEBDengine(const Container<KeyType,T>&c):Container<KeyType,T>(c){_blockthreads=1;_tebdthreads=1;
    _twMeasurementStep=1;_nameit=1;_storeallcp=false;_TEBDstep=1;}
  TEBDengine(const MPS<KeyType,T> &mps ,MPO<T> & mpo,const std::string simName,const bool storeallcp=false):Container<KeyType,T>(mps,mpo,simName){
    this->_checkpointStep=50;
    _blockthreads=1;
    _tebdthreads=1;
    _twMeasurementStep=1;
    _nameit=1;
    _storeallcp=storeallcp;
    _TEBDstep=1;
  }
  
  void simulate(list<measurementObject<T> > &mList,const std::string order,const uint MAXSTEPS=100,const T dt=T(0.05),
		const uint chimax=500, const bool doTrunc=false,const Real tr_weight=1e-12,const bool normalize=true,const bool recanonize=false,
		const int blockthreads=1,const int tebdthreads=1);

  void simulateNNN(list<measurementObject<T> > &mList,const std::string order,const uint MAXSTEPS=100,const T dt=T(0.05),
		   const uint chimax=500, const bool doTrunc=false,const Real tr_weight=1e-12,const bool normalize=true,const bool recanonize=false,
		   const int blockthreads=1,const int tebdthreads=1);

  void measureCMPS(list<measurementObject<T> > &m){Container<KeyType,T>::measureCMPS(m);}
  virtual void clear(){Container<KeyType,T>::clear();_NNGates.clear();_NNGatesHalf.clear();}
  void setTwMeasurementStep(const uint gTwMeasurementStep){this->_twMeasurementStep=gTwMeasurementStep;}
  uint getTwMeasurementStep(){return this->_twMeasurementStep;}

  void storecheckpoints(const bool store){_storeallcp=store;}
  // virtual void read();
  // virtual void write();
  // virtual void read(std::string name);
  // virtual void write(std::string name);
  void resetTEBDstep(uint step){_TEBDstep=step;}
protected:
  virtual void toarchive(boost::archive::binary_oarchive &oar){
    Container<KeyType,T>::toarchive(oar);
    oar<<_blockthreads;
    oar<<_tebdthreads;
    oar<<_twMeasurementStep;
    oar<<_nameit;
    oar<<_storeallcp;
    oar<<_TEBDstep;
    
  }
  virtual void fromarchive(boost::archive::binary_iarchive &iar){
    Container<KeyType,T>::fromarchive(iar);
    iar>>_blockthreads;
    iar>>_tebdthreads;
    iar>>_twMeasurementStep;
    iar>>_nameit;
    iar>>_storeallcp;
    iar>>_TEBDstep;
  }

  //the single gate application
  //two-site gate
  pair<std::vector<Real>,std::vector<uint> > doStep(const SparseOperator<2,T> &gate,const int site,const uint chimax,const bool doTrunc,const Real tr_weight,
  						    const bool normalize=true);
  //three-site gate
  pair<std::vector<Real>,std::vector<uint> > doStep(const SparseOperator<3,T> &gate,const int site,const uint chimax,const bool doTrunc,const Real tr_weight,
						    const bool normalize=true);
 
  //does a full trotter sweep of even or odd bonds (depending on start and inc (=increment) )
  VectorType doTimeSlice(const std::vector<SparseOperator<2,T> >&gates,const int start,const int inc,const uint chimax,const bool doTrunc,const Real tr_weight,
			 const bool normalize=true,const bool recanonize=false);
  //does a full trotter sweep of first, second or third bonds (depending on start and inc (=increment) ) with three-site gates
  VectorType doTimeSlice(const std::vector<SparseOperator<3,T> >&gates,const int start,const int inc,const uint chimax,const bool doTrunc,const Real tr_weight,
			 const bool normalize=true,const bool recanonize=false);
  //for convenience, the gates are stored in vectors ...
  //all two-site gates
  std::vector<SparseOperator<2,T> >_NNGates,_NNGatesHalf,_NNGatesQuart;
  //all three-site gates
  std::vector<SparseOperator<3,T> >_NNNGates,_NNNGatesHalf;

  //these functions implement the somewhat cumbersome way measurements have to be done if you use a second order trotter ...
  //for NN interactions
  void doCMPSmeasurement(list<measurementObject<T> > &m)const{Container<KeyType,T>::doCMPSmeasurement(m);}
  bool doCMPSmeasurement(list<measurementObject<T> > &m,const uint chimax,const bool doTrunc,const Real tr_weight, bool saveIteration,VectorType &tw,
			 const bool normalize=true,const bool recanonize=false,const bool skipsecondhalfstep=false);
  //for Next NN interactions
  bool doCMPSmeasurementNNN(list<measurementObject<T> > &m,const uint chimax,const bool doTrunc,const Real tr_weight, bool saveIteration,VectorType &twe,
			    VectorType &two,const bool normalize=true,const bool recanonize=false,const bool skipsecondhalfstep=false);


private:
  int _blockthreads;
  int _tebdthreads;
  uint _twMeasurementStep;
  uint _nameit;
  bool _storeallcp;
  uint _TEBDstep;//stores the current TEBD step of the simulate routine; when later restarted from a read() command, the routine starts from 
                 //the corresponding position
};



// template<class KeyType,typename T>
// void TEBDengine<KeyType,T>::write(std::string fname){
//   std::string name=fname+".container";
//   ofstream os(name.c_str(), ios::binary);
//   boost::archive::binary_oarchive oar(os);
//   toarchive(oar);
//   os.close();
// }

// template<class KeyType,typename T>
// void TEBDengine<KeyType,T>::read(std::string fname){
//   std::string name=fname+".container";
//   ifstream is(name.c_str(), ios::binary);
//   boost::archive::binary_iarchive iar(is);
//   fromarchive(iar);
//   is.close();
// }

// template<class KeyType,typename T>
// void TEBDengine<KeyType,T>::write(){
//   std::string name=this->_simName+".container";
//   ofstream os(name.c_str(), ios::binary);
//   boost::archive::binary_oarchive oar(os);
//   toarchive(oar);
//   os.close();
// }

// template<class KeyType,typename T>
// void TEBDengine<KeyType,T>::read(){
//   std::string name=this->_simName+".container";
//   ifstream is(name.c_str(), ios::binary);
//   boost::archive::binary_iarchive iar(is);
//   fromarchive(iar);
//   is.close();
// }

//measurement during time-evolution

template<class KeyType,typename T>
bool TEBDengine<KeyType,T>::doCMPSmeasurement(list<measurementObject<T> > &m,const uint chimax,const bool doTrunc,const Real tr_weight, 
					      bool saveIteration, VectorType &tw, const bool normalize, const bool recanonize,
					      const bool skipsecondhalfstep){
  bool doHStep = true;
  if(!this->_cmps.size()){
    cout<<"TEBDengine<KeyType,T>::measure(list<measurementObject<T>, uint ,bool, Real): empty CanonizedMPS<KeyType,T> found"<<endl;
    return false;
  }
  VectorType tw1(tw.size()),tw2(tw.size());
  boost_setTo(tw1,0.0);  boost_setTo(tw2,0.0);
  typename list<measurementObject<T> >::iterator it;
  for(it = m.begin();it!=m.end();it++){
    if(it->step%it->measureStep==0){//if any measurement is to be made, then do a trotter half step and then measure
      if(doHStep){
	doHStep=false;
	tw1=doTimeSlice(_NNGatesHalf,1,2,chimax,doTrunc, tr_weight,normalize,recanonize);
      }
      it->measure(this->_cmps);
    }
    it->step++;
  }
  if(!doHStep){
    if(saveIteration){//checkpointing if it's a checkpointstep
      std::string name=this->_simName+"_CP";
      this->write(name);
      if(_storeallcp){
	std::ostringstream convert;
	convert<<_nameit;
	name+=convert.str();
	this->_nameit++;
	cout<<name<<endl;
      }
      this->save(name);
    }
    //if a measurement has been done, do the second half step
    if(!skipsecondhalfstep)
      tw2=doTimeSlice(_NNGatesHalf,1,2,chimax,doTrunc, tr_weight,normalize,recanonize);
  }
  tw=tw1+tw2;
  return !doHStep;
}


template<class KeyType,typename T>
bool TEBDengine<KeyType,T>::doCMPSmeasurementNNN(list<measurementObject<T> > &m,const uint chimax,const bool doTrunc,const Real tr_weight, bool saveIteration,
						 VectorType &twe,VectorType &two, const bool normalize,const bool recanonize,const bool skipsecondhalfstep){
  assert(twe.size()==two.size());
  bool doHStep = true;
  if(!this->_cmps.size()){
    cout<<"TEBDengine<KeyType,T>::measure(list<measurementObject<T>, uint ,bool, Real): empty CanonizedMPS<KeyType,T> found"<<endl;
    return false;
  }
  VectorType tw1(twe.size()),tw2(twe.size()),tw3(twe.size()),tw1_(twe.size()),tw2_(twe.size()),tw3_(twe.size());
  
  boost_setTo(tw1,0.0);  boost_setTo(tw2,0.0);  boost_setTo(tw3,0.0);  
  boost_setTo(tw1_,0.0);  boost_setTo(tw2_,0.0);  boost_setTo(tw3_,0.0);  
  typename list<measurementObject<T> >::iterator it;
  for(it = m.begin();it!=m.end();it++){
    if(it->step%it->measureStep==0){
      if(doHStep){
	doHStep=false;
	tw1=doTimeSlice(_NNGatesQuart,1,2,chimax,doTrunc, tr_weight,normalize,recanonize);//odd
	tw2=doTimeSlice(_NNGatesHalf,0,2,chimax,doTrunc, tr_weight,normalize,recanonize);//even
	tw3=doTimeSlice(_NNGatesQuart,1,2,chimax,doTrunc, tr_weight,normalize,recanonize);//odd
      }
      it->measure(this->_cmps);
    }
    it->step++;
  }
  if(!doHStep){
    if(saveIteration){
      std::string name=this->_simName+"_CP";
      this->write(name);
      this->save(name);
    }
    if(!skipsecondhalfstep){
      tw1_=doTimeSlice(_NNGatesQuart,1,2,chimax,doTrunc, tr_weight,normalize,recanonize);//odd
      tw2_=doTimeSlice(_NNGatesHalf,0,2,chimax,doTrunc, tr_weight,normalize,recanonize);//even
      tw3_=doTimeSlice(_NNGatesQuart,1,2,chimax,doTrunc, tr_weight,normalize,recanonize);//odd
    }
    twe=tw2+tw2_;
    two=tw1+tw3+tw1_+tw3_;
  }
  return !doHStep;
}




//does the TEBD step
//one should template the length of the gate, and also applyGate, then simulation could be done for all interaction lengths
template<class KeyType,typename T>
pair<std::vector<Real>,std::vector<uint> > TEBDengine<KeyType,T>::doStep(const SparseOperator<2,T> &gate,const int site,const uint chimax,
									 const bool doTrunc,const Real tr_weight,const bool normalize){
  assert(site<(this->_cmps.size()-1));
  assert(this->_cmps.size());
  assert(this->_cmps.getPSpace()(site)*this->_cmps.getPSpace()(site+1)==gate.size1());
  assert(this->_cmps.getPSpace()(site)*this->_cmps.getPSpace()(site+1)==gate.size2());
  //cout<<"at site "<<site<<": ";
  pair<std::vector<Real>,std::vector<uint> >out=applyGate(this->_cmps.lambda(site-1),this->_cmps[site],this->_cmps.lambda(site),
							  this->_cmps[site+1],this->_cmps.lambda(site+1),gate,this->_cmps.getItoKey(site),
							  this->_cmps.getItoKey(site+1),chimax,doTrunc,tr_weight,normalize,_blockthreads);
  return out;  
};

template<class KeyType,typename T>
pair<std::vector<Real>,std::vector<uint> > TEBDengine<KeyType,T>::doStep(const SparseOperator<3,T> &gate,const int site,const uint chimax,
									 const bool doTrunc, const Real tr_weight,const bool normalize){
  assert(site<(this->_cmps.size()-1));
  assert(this->_cmps.size());
  assert(this->_cmps.getPSpace()(site)*this->_cmps.getPSpace()(site+1)*this->_cmps.getPSpace()(site+2)==gate.size1());
  assert(this->_cmps.getPSpace()(site)*this->_cmps.getPSpace()(site+1)*this->_cmps.getPSpace()(site+1)==gate.size2());
  pair<std::vector<Real>,std::vector<uint> > out=
    applyGate(this->_cmps.lambda(site-1),this->_cmps[site],this->_cmps.lambda(site),this->_cmps[site+1],this->_cmps.lambda(site+1),
	      this->_cmps[site+2],this->_cmps.lambda(site+2),gate,this->_cmps.getItoKey(site),this->_cmps.getItoKey(site+1),
	      this->_cmps.getItoKey(site+2),chimax,doTrunc,tr_weight,normalize);
  return out;  
};



//this could be overloaded if you template the SparseOperator length;
template<class KeyType,typename T>
VectorType TEBDengine<KeyType,T>::doTimeSlice(const std::vector<SparseOperator<2,T> >&gates,const int start,const int inc,const uint chimax,const bool doTrunc,
					      const Real tr_weight,const bool normalize,const bool recanonize)
{
  cout << "."; cout.flush();
  const int N=this->_cmps.size();
  VectorType t_weights(N-1);
  UintVectorType chi(N-1);
  boost_setTo(t_weights,0.0);
  boost_setTo(chi,uint(0));
  const int oldmaxthreads=1;
  //omp_set_num_threads(_tebdthreads);
#pragma omp parallel for
  for(int site = start;site<(N-inc+1);site+=inc){
    pair<std::vector<Real>,std::vector<uint> > out=doStep(gates[site],site,chimax,doTrunc,tr_weight,normalize);
    t_weights(site)= out.first[0];
    chi(site) = out.second[0];
  }
  //omp_set_num_threads(oldmaxthreads);
  //while(!checkMPS(this->_cmps));

  //NUSS 09 05 2012 - added renormalization of state
  if(recanonize){
    MPS<KeyType,T> mps= this->_cmps.toMPS();
    //    orthonormalize(mps,"lr");
    orthonormalize(mps,"rl");
    orthonormalize(mps,"lr");
    this->_cmps=canonize(mps);
    for(uint s=0;s<mps.size();++s)
      if(abs(this->_cmps.lambda(s).norm()-1)>1e-8){
	cout<<"non-normalized lambda at site "<<s<<": "<<this->_cmps.lambda(s).norm()<<endl;
      }
  }
  
  return t_weights;
}


template<class KeyType,typename T>
VectorType TEBDengine<KeyType,T>::doTimeSlice(const std::vector<SparseOperator<3,T> >&gates,const int start,const int inc,const uint chimax,
					      const bool doTrunc,const Real tr_weight,const bool normalize,const bool recanonize){
  cout << "."; cout.flush();
  const int N=this->_cmps.size();
  VectorType t_weights(N-1);
  UintVectorType chi(N-1);
  boost_setTo(t_weights,0.0);
  boost_setTo(chi,uint(0));
  const int oldmaxthreads=1;
  //omp_set_num_threads(_tebdthreads);
#pragma omp parallel for
  for(int site = start;site<(N-inc+1);site+=inc){
    pair<std::vector<Real>,std::vector<uint> >out=doStep(gates[site],site,chimax,doTrunc,tr_weight,normalize);
    t_weights(site)= out.first[0];
    t_weights(site+1)= out.first[1];
    chi(site) = out.second[0];
    chi(site+1) = out.second[1];
  }
  //omp_set_num_threads(oldmaxthreads);
  //while(!checkMPS(this->_cmps));
  
  //NUSS 09 05 2012 - added recanonization of state
  if(recanonize){
    MPS<KeyType,T> mps= this->_cmps.toMPS();
    //    orthonormalize(mps,"lr");
    orthonormalize(mps,"rl");
    orthonormalize(mps,"lr");
    this->_cmps=canonize(mps);
    for(uint s=0;s<mps.size();++s)
      if(abs(this->_cmps.lambda(s).norm()-1)>1e-8){
	cout<<"non-normalized lambda at site "<<s<<": "<<this->_cmps.lambda(s).norm()<<endl;
      }
  }
  
  return t_weights;
}

template<class KeyType,typename T>
void TEBDengine<KeyType,T>::simulate(list<measurementObject<T> > &mList,const std::string order,const uint MAXSTEPS,const T dt,
				     const uint chimax, const bool doTrunc,const Real tr_weight,const bool normalize,const bool recanonize,
				     const int blockthreads,const int tebdthreads){
  _tebdthreads=tebdthreads;
  _blockthreads=blockthreads;
  uint step=_TEBDstep;

  if(this->_cmps.size()==0){
    cout<<"_mps is not canonized; canonizing ..."<<endl;
    this->_cmps=canonize(this->_mps);
  }
  VectorType twe(this->_cmps.size()-1),two(this->_cmps.size()-1);
  boost_setTo(two,0.0);
  boost_setTo(twe,0.0);

  ofstream tw_outfile;
  std::string fname = "truncatedWeight" + this->_simName;
  list<std::vector<Real> > tr_weightMeasurementList;
  std::vector<Real> tr_weightMeasurementVec;

  if(order=="first"){
    cout<<"starting first order TEBD"<<endl;
    while(step<=MAXSTEPS){
      this->_mpo.construct(step-1);
      for(uint i=0;i<(this->_mpo.size()-1);++i){
	_NNGates.push_back(this->_mpo.getNNGate(this->_mps,i,i+1,dt));
	_NNGates[i].squeeze(1e-14);
      }
 
      if(_NNGates.size()==0){cout<<"no second order time evolution operators defined; aborting "<<endl;abort();}
      if(_NNGates.size()==this->_mps.size()){cout<<"in TEBD: length of Time evolution operator differs from that of MPS"<<endl;abort();}
      twe=doTimeSlice(_NNGates,0,2,chimax,doTrunc, tr_weight,normalize,recanonize);
      two=doTimeSlice(_NNGates,1,2,chimax,doTrunc, tr_weight,normalize,recanonize);
      // new version of saving trunc weight 09.02.2012 - mechanism is made more similar to actual measurement
      if(step%_twMeasurementStep==0){	//append tw -lists
	for(uint site = 0;site<(this->_cmps.size()-1);++site){
	    tr_weightMeasurementVec.push_back( (site%2?two(site):twe(site)) );
	}
	tr_weightMeasurementList.push_back(tr_weightMeasurementVec);
	tr_weightMeasurementVec.clear();
      }
      if(step%this->_checkpointStep==0){	//access hdd only at those steps where it is anyway done by measurement
	tw_outfile.open(fname.c_str(),ios_base::app);
	for (list<std::vector<Real> >::const_iterator ci = tr_weightMeasurementList.begin(); ci != tr_weightMeasurementList.end(); ++ci) {
             tr_weightMeasurementVec = *ci;
	     for(uint site = 0;site<(this->_cmps.size()-1);++site){
	      tw_outfile<<tr_weightMeasurementVec[site]<<"\t";
	    }
	    tw_outfile<<endl;
	}
	tw_outfile.close();
	tr_weightMeasurementList.clear();
	tr_weightMeasurementVec.clear();
      }
      doCMPSmeasurement(mList);
      //save the MPS for checkpointing
      if (step%this->_checkpointStep==0){
	std::string name=this->_simName+"_CP";
	this->write(name);
	if(_storeallcp){
	  std::ostringstream convert;
	  convert<<_nameit;
	  name+=convert.str();
	  this->_nameit++;
	}
	this->save(name);
	this->writeMeasurements(mList);
	cout << endl;
      }
      step++;
      _TEBDstep=step;
    }
  }else if(order=="second"){
    cout<<"starting second order TEBD"<<endl;
    cout<<step<<endl;
    this->_mpo.construct(step-1);
    for(uint i=0;i<(this->_mpo.size()-1);++i){
      _NNGates.push_back(this->_mpo.getNNGate(this->_mps,i,i+1,dt));
      _NNGates[i].squeeze(1e-14);
      _NNGatesHalf.push_back(this->_mpo.getNNGate(this->_mps,i,i+1,dt/2.0));
      _NNGatesHalf[i].squeeze(1e-14);
    }

    two=doTimeSlice(_NNGatesHalf,1,2,chimax,doTrunc, tr_weight,normalize,recanonize);
    while(step<=MAXSTEPS){
      twe=doTimeSlice(_NNGates,0,2,chimax,doTrunc, tr_weight,normalize,recanonize);
      bool skipOddStep=doCMPSmeasurement(mList,chimax,doTrunc,tr_weight,!(bool)(step%this->_checkpointStep),two,normalize,recanonize,step==MAXSTEPS);
      if((!skipOddStep)&&(step<MAXSTEPS)){
	two=doTimeSlice(_NNGates,1,2,chimax,doTrunc, tr_weight,normalize,recanonize);
      }
      if((!skipOddStep)&&(step==MAXSTEPS)){
	two=doTimeSlice(_NNGatesHalf,1,2,chimax,doTrunc, tr_weight,normalize,recanonize);
      }

      if(step%this->_twMeasurementStep==0){	//append tw -lists
	for(uint site = 0;site<(this->_cmps.size()-1);++site){
	    tr_weightMeasurementVec.push_back( (site%2?two(site):twe(site)) );
	}
	tr_weightMeasurementList.push_back(tr_weightMeasurementVec);
	tr_weightMeasurementVec.clear();
      }
      if(step%this->_checkpointStep==0){	//access hdd only at those steps where it is anyway done by measurement
	tw_outfile.open(fname.c_str(),ios_base::app);
	for (list<std::vector<Real> >::const_iterator ci = tr_weightMeasurementList.begin(); ci != tr_weightMeasurementList.end(); ++ci) {
             tr_weightMeasurementVec = *ci;
	     for(uint site = 0;site<(this->_cmps.size()-1);++site){
	      tw_outfile<<tr_weightMeasurementVec[site]<<"\t";
	    }
	    tw_outfile<<endl;
	}
	tw_outfile.close();
	tr_weightMeasurementList.clear();
	tr_weightMeasurementVec.clear();
      }
      if (step%this->_checkpointStep==0){
	//this->save();
	this->writeMeasurements(mList);
	cout << endl;
      }
      step++;	//step is incremented at last so that the first step will be written to disk before every checkpointStep data is written
      _TEBDstep=step;
    }
  }
  else{
    cout<<"#in TEBDengine::doTEBD(): wrong input string order; use \"first\" or \"second\" "<<endl;
  } 
  this->save();
  std::string name=this->_simName+"_CP";
  this->write(name);
  this->writeMeasurements(mList);
}


//note that the second order expansion is only a real second order expansion of the NNN-Gates commute with each other
//if that's not the case, it's just a somewhat better first order expansion. You can however try to implement a real second order ...
//new version
template<class KeyType,typename T>
void TEBDengine<KeyType,T>::simulateNNN(list<measurementObject<T> > &mList,const std::string order,const uint MAXSTEPS,const T dt,const uint chimax, 
					const bool doTrunc,const Real tr_weight,const bool normalize,const bool recanonize,const int blockthreads,
					const int tebdthreads){
  _tebdthreads=tebdthreads;
  _blockthreads=blockthreads;
  //Complex damped_t = Complex(dt,damping);		//NUSS 09 05 2012
  for(uint i=0;i<(this->_mpo.size()-1);++i){
    _NNGates.push_back(this->_mpo.getNNGate(this->_mps,i,i+1,dt));			//NUSS 09 05 2012 	dt->damped_t
    _NNGates[i].squeeze(1e-14);
    _NNGatesHalf.push_back(this->_mpo.getNNGate(this->_mps,i,i+1,dt/2.0));		//NUSS 09 05 2012 	dt->damped_t
    _NNGatesHalf[i].squeeze(1e-14);
    _NNGatesQuart.push_back(this->_mpo.getNNGate(this->_mps,i,i+1,dt/4.0));		//NUSS 09 05 2012 	dt->damped_t
    _NNGatesQuart[i].squeeze(1e-14);
  }
  for(uint i=0;i<(this->_mpo.size()-2);++i){
    _NNNGates.push_back(this->_mpo.getNNNGate(this->_mps,i,i+1,i+2,dt));		//NUSS 09 05 2012 	dt->damped_t
    _NNNGates[i].squeeze(1e-14);
    _NNNGatesHalf.push_back(this->_mpo.getNNNGate(this->_mps,i,i+1,i+2,dt/2.0));		//NUSS 09 05 2012 	dt->damped_t
    _NNNGatesHalf[i].squeeze(1e-14);
  }

  if(_NNGates.size()==0){cout<<"#no second order time evolution operators defined; aborting "<<endl;abort();}
  if(_NNNGates.size()==0){cout<<"#no second order time evolution operators defined; aborting "<<endl;abort();}
  if(_NNGates.size()==this->_mps.size()){cout<<"#in TEBD: length of time evolution operator differs from that of MPS"<<endl;abort();}

  if(this->_cmps.size()==0){
    cout<<"#in TEBDengine: _mps is not canonized; canonizing ..."<<endl;
    this->_cmps=canonize(this->_mps);
  }

  VectorType twe(this->_cmps.size()-1),two(this->_cmps.size()-1),two1(this->_cmps.size()-1),two2(this->_cmps.size()-1);
  VectorType twf(this->_cmps.size()-1),tws(this->_cmps.size()-1),twt(this->_cmps.size()-1);
  boost_setTo(two,0.0);
  boost_setTo(twe,0.0);
  boost_setTo(two1,0.0);
  boost_setTo(two2,0.0);
  boost_setTo(twf,0.0);
  boost_setTo(tws,0.0);
  boost_setTo(twt,0.0);

  ofstream tw_outfile;
  uint step=_TEBDstep;
  std::string fname = "truncatedWeight_"+this->_simName ;
  std::string fname1 = "truncatedWeightNNN1_"+this->_simName ;
  std::string fname2 = "truncatedWeightNNN2_"+this->_simName ;
  std::string fname3 = "truncatedWeightNNN3_"+this->_simName ;
  list<std::vector<Real> > tr_weightMeasurementList;
  list<std::vector<Real> > tr_weightMeasurementListNNN1;
  list<std::vector<Real> > tr_weightMeasurementListNNN2;
  list<std::vector<Real> > tr_weightMeasurementListNNN3;
  std::vector<Real> tr_weightMeasurementVec;
  std::vector<Real> tr_weightMeasurementVecNNN1;
  std::vector<Real> tr_weightMeasurementVecNNN2;
  std::vector<Real> tr_weightMeasurementVecNNN3;
  if(order=="first"){
    cout<<"starting first order TEBD for NNN interactions"<<endl;
    while(step<MAXSTEPS){
      twe=doTimeSlice(_NNGates,0,2,chimax,doTrunc, tr_weight,normalize,recanonize);
      two=doTimeSlice(_NNGates,1,2,chimax,doTrunc, tr_weight,normalize,recanonize);
      twf=doTimeSlice(_NNNGates,0,3,chimax,doTrunc, tr_weight,normalize,recanonize);
      tws=doTimeSlice(_NNNGates,1,3,chimax,doTrunc, tr_weight,normalize,recanonize);
      twt=doTimeSlice(_NNNGates,2,3,chimax,doTrunc, tr_weight,normalize,recanonize);
      //IF not working, comment the following and uncomment the lower part
      //=================================================================================
      // new version of saving trunc weight 09.02.2012 - mechanism is made more similar to actual measurement
      if(step%_twMeasurementStep==0){	//append tw -lists
	for(uint site = 0;site<(this->_cmps.size()-1);++site){
	  tr_weightMeasurementVec.push_back((site%2?two(site):twe(site)));
	  tr_weightMeasurementVecNNN1.push_back(twf(site));
	  tr_weightMeasurementVecNNN2.push_back(tws(site));
	  tr_weightMeasurementVecNNN3.push_back(twt(site));
	}
	tr_weightMeasurementList.push_back(tr_weightMeasurementVec);
	tr_weightMeasurementVec.clear();
	tr_weightMeasurementListNNN1.push_back(tr_weightMeasurementVecNNN1);
	tr_weightMeasurementVecNNN1.clear();
	tr_weightMeasurementListNNN2.push_back(tr_weightMeasurementVecNNN2);
	tr_weightMeasurementVecNNN2.clear();
	tr_weightMeasurementListNNN3.push_back(tr_weightMeasurementVecNNN3);
	tr_weightMeasurementVecNNN3.clear();

      }
      if(step%this->_checkpointStep==0){	//access hdd only at those steps where it is anyway done by measurement
	tw_outfile.open(fname.c_str(),ios_base::app);
	for (list<std::vector<Real> >::const_iterator ci = tr_weightMeasurementList.begin(); ci != tr_weightMeasurementList.end(); ++ci) {
             tr_weightMeasurementVec = *ci;
	     for(uint site = 0;site<(this->_cmps.size()-1);++site){
	      tw_outfile<<tr_weightMeasurementVec[site]<<"\t";
	    }
	    tw_outfile<<endl;
	}
	tw_outfile.close();
	tr_weightMeasurementList.clear();
	tr_weightMeasurementVec.clear();


	tw_outfile.open(fname1.c_str(),ios_base::app);
	for (list<std::vector<Real> >::const_iterator ci = tr_weightMeasurementListNNN1.begin(); ci != tr_weightMeasurementListNNN1.end(); ++ci) {
             tr_weightMeasurementVecNNN1 = *ci;
	     for(uint site = 0;site<(this->_cmps.size()-1);++site){
	      tw_outfile<<tr_weightMeasurementVecNNN1[site]<<"\t";
	    }
	    tw_outfile<<endl;
	}
	tw_outfile.close();
	tr_weightMeasurementListNNN1.clear();
	tr_weightMeasurementVecNNN1.clear();


	tw_outfile.open(fname2.c_str(),ios_base::app);
	for (list<std::vector<Real> >::const_iterator ci = tr_weightMeasurementListNNN2.begin(); ci != tr_weightMeasurementListNNN2.end(); ++ci) {
             tr_weightMeasurementVecNNN2 = *ci;
	     for(uint site = 0;site<(this->_cmps.size()-1);++site){
	      tw_outfile<<tr_weightMeasurementVecNNN2[site]<<"\t";
	    }
	    tw_outfile<<endl;
	}
	tw_outfile.close();
	tr_weightMeasurementListNNN2.clear();
	tr_weightMeasurementVecNNN2.clear();


	tw_outfile.open(fname3.c_str(),ios_base::app);
	for (list<std::vector<Real> >::const_iterator ci = tr_weightMeasurementListNNN3.begin(); ci != tr_weightMeasurementListNNN3.end(); ++ci) {
             tr_weightMeasurementVecNNN3 = *ci;
	     for(uint site = 0;site<(this->_cmps.size()-1);++site){
	      tw_outfile<<tr_weightMeasurementVecNNN3[site]<<"\t";
	    }
	    tw_outfile<<endl;
	}
	tw_outfile.close();
	tr_weightMeasurementListNNN3.clear();
	tr_weightMeasurementVecNNN3.clear();

      }
      //=============================================
      //      if(step%10*this->_checkpointStep==0){
      //	outfile.open(fname.c_str(),ios_base::app);
      //	for(uint site = 0;site<(this->_cmps.size()-1);++site){
      //	  outfile<<(site%2?two(site):twe(site))<<"\t";
      //	}
      //	outfile<<endl;
      //	outfile.close();
      //      }
      //=============================================
      doCMPSmeasurement(mList);
      //save the MPS for checkpointing
      if (step%this->_checkpointStep==0){
	std::string name=this->_simName+"_CP";
	this->save(name);
	this->write(name);
	writeMeasurements(mList);
	cout << endl;
      }
      step++;
      _TEBDstep=step;
    }
  }
  else if(order=="second"){
    cout<<"starting second order TEBD for NNN"<<endl;
    two1=doTimeSlice(_NNGatesQuart,1,2,chimax,doTrunc, tr_weight,normalize,recanonize);//odd
    twe=doTimeSlice(_NNGatesHalf,0,2,chimax,doTrunc, tr_weight,normalize,recanonize);//even
    two2=doTimeSlice(_NNGatesQuart,1,2,chimax,doTrunc, tr_weight,normalize,recanonize);//odd
    //truncated weight of first step not stored
    while(step<MAXSTEPS){
      twf=doTimeSlice(_NNNGates,0,3,chimax,doTrunc, tr_weight,normalize,recanonize);//first
      tws=doTimeSlice(_NNNGates,1,3,chimax,doTrunc, tr_weight,normalize,recanonize);//second
      twt=doTimeSlice(_NNNGates,2,3,chimax,doTrunc, tr_weight,normalize,recanonize);//third
      bool skipOddStep=doCMPSmeasurementNNN(mList,chimax,doTrunc,tr_weight,!(bool)(step%this->_checkpointStep),twe,two,normalize,recanonize,step==MAXSTEPS);
      if((!skipOddStep)&&(step<MAXSTEPS)){
	two1=doTimeSlice(_NNGatesHalf,1,2,chimax,doTrunc, tr_weight,normalize,recanonize);//odd
	twe=doTimeSlice(_NNGates,0,2,chimax,doTrunc, tr_weight,normalize,recanonize);//even
	two2=doTimeSlice(_NNGatesHalf,1,2,chimax,doTrunc, tr_weight,normalize,recanonize);//odd
	two=two1+two2;
      }
      if((!skipOddStep)&&(step==MAXSTEPS)){
	two1=doTimeSlice(_NNGatesQuart,1,2,chimax,doTrunc, tr_weight,normalize,recanonize);//odd
	twe=doTimeSlice(_NNGatesHalf,0,2,chimax,doTrunc, tr_weight,normalize,recanonize);//even
	two2=doTimeSlice(_NNGatesQuart,1,2,chimax,doTrunc, tr_weight,normalize,recanonize);//odd
	two=two1+two2;
      }

      //=====================================this is new; if not working, comment it out and uncomment lower part ===========================================
      if(step%_twMeasurementStep==0){	//append tw -lists
	for(uint site = 0;site<(this->_cmps.size()-1);++site){
	    tr_weightMeasurementVec.push_back( (site%2?two(site):twe(site)) );
	    tr_weightMeasurementVecNNN1.push_back(twf(site));
	    tr_weightMeasurementVecNNN2.push_back(tws(site));
	    tr_weightMeasurementVecNNN3.push_back(twt(site));
	}
	tr_weightMeasurementList.push_back(tr_weightMeasurementVec);
	tr_weightMeasurementVec.clear();
	tr_weightMeasurementListNNN1.push_back(tr_weightMeasurementVecNNN1);
	tr_weightMeasurementVecNNN1.clear();
	tr_weightMeasurementListNNN2.push_back(tr_weightMeasurementVecNNN2);
	tr_weightMeasurementVecNNN2.clear();
	tr_weightMeasurementListNNN3.push_back(tr_weightMeasurementVecNNN3);
	tr_weightMeasurementVecNNN3.clear();

      }

      if(step%this->_checkpointStep==0){	//access hdd only at those steps where it is anyway done by measurement
	tw_outfile.open(fname.c_str(),ios_base::app);
	for (list<std::vector<Real> >::const_iterator ci = tr_weightMeasurementList.begin(); ci != tr_weightMeasurementList.end(); ++ci) {
             tr_weightMeasurementVec = *ci;
	     for(uint site = 0;site<(this->_cmps.size()-1);++site){
	      tw_outfile<<tr_weightMeasurementVec[site]<<"\t";
	    }
	    tw_outfile<<endl;
	}
	tw_outfile.close();
	tr_weightMeasurementList.clear();
	tr_weightMeasurementVec.clear();


	tw_outfile.open(fname1.c_str(),ios_base::app);
	for (list<std::vector<Real> >::const_iterator ci = tr_weightMeasurementListNNN1.begin(); ci != tr_weightMeasurementListNNN1.end(); ++ci) {
             tr_weightMeasurementVecNNN1 = *ci;
	     for(uint site = 0;site<(this->_cmps.size()-1);++site){
	      tw_outfile<<tr_weightMeasurementVecNNN1[site]<<"\t";
	    }
	    tw_outfile<<endl;
	}
	tw_outfile.close();
	tr_weightMeasurementListNNN1.clear();
	tr_weightMeasurementVecNNN1.clear();


	tw_outfile.open(fname2.c_str(),ios_base::app);
	for (list<std::vector<Real> >::const_iterator ci = tr_weightMeasurementListNNN2.begin(); ci != tr_weightMeasurementListNNN2.end(); ++ci) {
             tr_weightMeasurementVecNNN2 = *ci;
	     for(uint site = 0;site<(this->_cmps.size()-1);++site){
	      tw_outfile<<tr_weightMeasurementVecNNN2[site]<<"\t";
	    }
	    tw_outfile<<endl;
	}
	tw_outfile.close();
	tr_weightMeasurementListNNN2.clear();
	tr_weightMeasurementVecNNN2.clear();


	tw_outfile.open(fname3.c_str(),ios_base::app);
	for (list<std::vector<Real> >::const_iterator ci = tr_weightMeasurementListNNN3.begin(); ci != tr_weightMeasurementListNNN3.end(); ++ci) {
             tr_weightMeasurementVecNNN3 = *ci;
	     for(uint site = 0;site<(this->_cmps.size()-1);++site){
	      tw_outfile<<tr_weightMeasurementVecNNN3[site]<<"\t";
	    }
	    tw_outfile<<endl;
	}
	tw_outfile.close();
	tr_weightMeasurementListNNN3.clear();
	tr_weightMeasurementVecNNN3.clear();

      }

      if (step%this->_checkpointStep==0){
	//this->save();
	writeMeasurements(mList);
	cout << endl;
      }
      step++;	//step is incremented at last so that the first step will be written to disk before every checkpointStep data is written
      _TEBDstep=step;
      //===================================================================================
      //    if(step%this->_checkpointStep*10==0){
      //	outfile.open(fname.c_str(),ios_base::app);
      //	for(uint site = 0;site<(this->_cmps.size()-1);++site){
      //	  outfile<<(site%2?two(site):twe(site))<<"\t";
      //	}
      //	outfile<<endl;
      //	outfile.close();
      //    }
      //    step++;
      //    if (step%this->_checkpointStep==0){
      //	//this->save();
      //	writeMeasurements(mList);
      //	cout << endl;
      //    }
      //===================================================================================
    }
  }
  else{
    cout<<"#in TEBDengine::doTEBD(): wrong input string order; use \"first\" or \"second\" "<<endl;
  } 
  this->save();
  std::string name=this->_simName+"_CP";
  this->write(name);
  writeMeasurements(mList);

}




template<class KeyType,typename T>
class ChebychevEngine:public Container<KeyType,T>{
public:
  ChebychevEngine():Container<KeyType,T>(){this->_checkpointStep=5;};
  ChebychevEngine(string simName){this->_checkpointStep=5;this->_simName=simName;this->_verbose=1;}
  ChebychevEngine(string simName,const uint verbose){this->_checkpointStep=5;this->_simName=simName;this->_verbose=verbose;}
  ChebychevEngine(const MPS<KeyType,T>&mps,MPO<T>&ham,const MPO<T>mpoi,const MPO<T>mpof,const std::string simname):
    Container<KeyType,T>(mps,ham,simname),_opi(mpoi),_opf(mpof){this->_checkpointStep=5;}
  ChebychevEngine(const MPS<KeyType,T>&mps,MPO<T>&ham,const MPO<T>mpoi,const MPO<T>mpof,const std::string simname,const uint verbose):
    Container<KeyType,T>(mps,ham,simname,verbose),_opi(mpoi),_opf(mpof){this->_checkpointStep=5;};

  ~ChebychevEngine(){};
  void simulate(const int Nmax,const int chimax,const int Esweeps,const int Dmax,list<measurementObject<T> > &mList, const int ncmax=10, 
		const Real delta=1e-12, const bool reconstruct=false,const bool resume=false,const uint truncstep=1,const Real landelta=1e-12,
		const bool cleanup=true,const bool PSS=false,const bool debug=false);
  void simulateSSCompress(const int Nmax,const int chimax,const int Esweeps,const int Dmax,list<measurementObject<T> > &mList, const int ncmax=10, 
		const Real delta=1e-12, const bool reconstruct=false,const bool resume=false,const uint truncstep=1,const Real landelta=1e-12,
		const bool cleanup=true);

  void simulateExponential(const int Nmax,const int chimax,const Real dt, const int ncmax=10,const Real delta=1e-12,const int trotterorder=1,
			   const bool reconstruct=false,const bool resume=false,const bool flip=false);
  std::vector<T> moments()const{return _moments;}
  std::vector<T> reconstructedMoments()const{return _recons_moments;}
  virtual void clear(){Container<KeyType,T>::clear();_opi.clear();_opf.clear();_t0.clear();_t1.clear();_t2.clear();_moments.clear();}
  // void read();//this overrides the definition in Container
  // void write();
  // void read(const string);
  // void write(const string);
protected:
  virtual void fromarchive(boost::archive::binary_iarchive &iar);
  virtual void toarchive(boost::archive::binary_oarchive &oar);
private:
  void measureMPS(const MPS<KeyType,T>&mps,list<measurementObject<T> > &m){
    if(!mps.size()){
      cout<<"Container<KeyType,T>::measure(list<measurementObject<T>): empty MPS<KeyType,T> found"<<endl;
      return ;
    }
    typename list<measurementObject<T> >::iterator it;
    for(it = m.begin();it!=m.end();++it){
      if(it->step%it->measureStep==0){
	it->measure(mps);
      }
      it->step++;
    }
  }

  void doEnergyTruncation(MPS<KeyType,T>&mps,MPO<T>&mpo,const int itmax,const int Dmax,const T Ethresh);
  void doSSEnergyTruncation(MPS<KeyType,T>&mps,MPO<T>&mpo,const int itmax,const int Dmax,const T Ethresh,const Real landelta);
  VectorType doTimeSlice(MPS<KeyType,T> &cmps,std::vector<SparseOperator<2,T> >&gates,const int start,const int inc,
			 const int blockthreads=1);

  pair<std::vector<Real>,std::vector<uint> > doStep(MPS<KeyType,T> &mps,const SparseOperator<2,T> &gate,const int site,
						    const int blockthreads=1);

  MPS<KeyType,T> applyOperator( MPS<KeyType,T>&mps,const uint chimax,const uint ncmax,const Real delta,const int trotterorder,const bool flip_w=false,
			       const int blockthreads=1);


  //use the following two with caution. Boundary-lambdas may be wrong ...
  CanonizedMPS<KeyType,T> rawCanonize(MPS<KeyType,T> &mps,const int dir);
  void rawOrtho(MPS<KeyType,T> &mps,const int dir);


  std::vector<SparseOperator<2,T> >_NNGates,_NNGatesHalf;
  MPO<T>_opi;
  MPO<T>_opf;
  std::vector<T>_moments;

  std::vector<T>_recons_moments;
  MPS<KeyType,T>_t0,_t1,_t2,_mpsf;
  CanonizedMPS<KeyType,T>_ct0,_ct1,_ct2;
  std::vector<T>_conv1_container,_conv2_container,_etrunc_container,_etrunc_step,_norm;
};


// template<class KeyType,typename T>
// void ChebychevEngine<KeyType,T>::write(){
//   std::string name=this->_simName+".container";
//   ofstream os(name.c_str(), ios::binary);
//   boost::archive::binary_oarchive oar(os);
//   Container<KeyType,T>::toarchive(oar);
//   toarchive(oar);
//   os.close();
// }


// template<class KeyType,typename T>
// void ChebychevEngine<KeyType,T>::read(){
//   std::string name=this->_simName+".container";
//   cout<<"#reading data from "<<name<<endl;
//   ifstream is(name.c_str(), ios::binary);
//   boost::archive::binary_iarchive iar(is);
//   Container<KeyType,T>::fromarchive(iar);
//   fromarchive(iar);
//   is.close();
// }

// template<class KeyType,typename T>
// void ChebychevEngine<KeyType,T>::write(const string fname){
//   std::string name=fname+".container";
//   ofstream os(name.c_str(), ios::binary);
//   boost::archive::binary_oarchive oar(os);
//   Container<KeyType,T>::toarchive(oar);
//   toarchive(oar);
//   os.close();
// }


// template<class KeyType,typename T>
// void ChebychevEngine<KeyType,T>::read(const string fname){
//   std::string name=fname+".container";
//   cout<<"#reading data from "<<name<<endl;
//   ifstream is(name.c_str(), ios::binary);
//   boost::archive::binary_iarchive iar(is);
//   Container<KeyType,T>::fromarchive(iar);
//   fromarchive(iar);
//   is.close();
// }


template<class KeyType,typename T>
void ChebychevEngine<KeyType,T>::toarchive(boost::archive::binary_oarchive &oar){
  Container<KeyType,T>::toarchive(oar);
  oar<<_t0;
  oar<<_t1;
  oar<<_t2;
  oar<<_mpsf;
  oar<<_opi;
  oar<<_opf;
  oar<<_moments;
  oar<<_recons_moments;
  oar<<_conv1_container;
  oar<<_conv2_container;
  oar<<_etrunc_container;
  oar<<_etrunc_step;
  oar<<_norm;
}


template<class KeyType,typename T>
void ChebychevEngine<KeyType,T>::fromarchive(boost::archive::binary_iarchive &iar){
  _t0.clear();_t1.clear();_t2.clear();_mpsf.clear();_opi.clear();_opf.clear();_moments.clear();_recons_moments.clear();
  Container<KeyType,T>::fromarchive(iar);
  iar>>_t0;
  iar>>_t1;
  iar>>_t2;
  iar>>_mpsf;
  iar>>_opi;
  iar>>_opf;
  iar>>_moments;
  iar>>_recons_moments;
  iar>>_conv1_container;
  iar>>_conv2_container;
  iar>>_etrunc_container;
  iar>>_etrunc_step;
  iar>>_norm;


}


template<class KeyType,typename T>
void ChebychevEngine<KeyType,T>::doEnergyTruncation(MPS<KeyType,T>&mps,MPO<T>&mpo,const int itmax,const int Dmax,const T Ethresh){
  for(int site = 0;site<mps.size()-1;++site){
    prepareSite(mps,site,"l",true);
  }
  for(int site = mps.size()-1;site>0;--site){
    prepareSite(mps,site,"r",true);
  }
  map_container<int,BlockMatrix<1,KeyType,T> >R,L;
  computeR(mps,mpo,2,R);
  L.clear();
  L[-1] = getLMin(mps);
  if(this->_verbose==1)
    cout<<"#norm of the state before energy truncation: "<<sqrt(computeOverlap(mps,mps))<<endl;
  bool normalize=true;
  bool converged = false;
  uint it = 0;
  while(!converged){
    if(it>=itmax){
      converged=true;
      break;
    }
    for(uint site=0;site<mps.size()-2;++site){
      TSLanzcosEngine<KeyType,T> 
	lan(L[site-1],R[site+2],mps[site],mps[site+1],mpo[site],mpo[site+1],mps.getItoKey(site),mps.getItoKey(site+1),this->_verbose);
      int N=lan.simulate(Dmax,1e-12,true);
      cout<<" at site "<<site<<" at step "<<it<<endl;
      //now get the projector ...
      BlockMatrix<2,KeyType,T> A1(mps[site],mps[site+1]),psi(mps[site],mps[site+1]),A2;
      int i0 = lan.Elargerthan(T(Ethresh));
      for(int i=i0;i<N;++i){
	A2=lan.getState(i);
	T c=A2*psi;
	A1-=(A2*c);
      }
      //now we should have the stripped state ...
      std::vector<BlockMatrix<1,KeyType,T> >Aout;
      std::vector<DiagBlockMatrix<KeyType,Real> > lambda;
      std::vector<i_to_key<KeyType> > itokey;
      itokey.push_back(mps.getItoKey(site));
      itokey.push_back(mps.getItoKey(site+1));
      uint chimax=mps.chi(site);

      pair<std::vector<Real>,std::vector<uint> >sre=split(A1,Aout,lambda,itokey,0,chimax,true,1e-16,normalize);
      mps[site]=Aout[0];
      mps[site+1]=lambda[0]*Aout[1];
      L[site] = addContractionLayer(L[site-1],mps[site],mpo[site],mps.getItoKey(site),1);
    }
    for(uint site=mps.size()-2;site>0;--site){
      TSLanzcosEngine<KeyType,T> 
	lan(L[site-1],R[site+2],mps[site],mps[site+1],mpo[site],mpo[site+1],mps.getItoKey(site),mps.getItoKey(site+1),this->_verbose);
      int N=lan.simulate(Dmax,1e-12,true);
      cout<<" at site "<<site<<" at step "<<it<<endl;
      //now get the projector ...
      BlockMatrix<2,KeyType,T> A1(mps[site],mps[site+1]),psi(mps[site],mps[site+1]),A2;
      int i0 = lan.Elargerthan(T(Ethresh));
      for(int i=i0;i<N;++i){
	A2=lan.getState(i);
	T c=A2*psi;
	A1-=(A2*c);
      }
      //now we should have the stripped state ...
      std::vector<BlockMatrix<1,KeyType,T> >Aout;
      std::vector<DiagBlockMatrix<KeyType,Real> > lambda;
      std::vector<i_to_key<KeyType> > itokey;
      itokey.push_back(mps.getItoKey(site));
      itokey.push_back(mps.getItoKey(site+1));
      uint chimax=mps.chi(site);
      //no truncation is made here, but the state is normalized ...
      pair<std::vector<Real>,std::vector<uint> >sre=split(A1,Aout,lambda,itokey,0,chimax,true,1e-14,normalize);
      
      mps[site]=Aout[0]*lambda[0];
      mps[site+1]=Aout[1];
      R[site+1] = addContractionLayer(R[site+2],mps[site+1],mpo[site+1],mps.getItoKey(site+1),-1);  
    }
    ++it;
  }
  // prepareSite(mps,0,"r",true);
  // mps.setLOrtho(false);
  // mps.setROrtho(true);
  cout<<"#norm of the state after energy truncation: "<<sqrt(computeOverlap(mps,mps))<<endl;
}


template<class KeyType,typename T>
void ChebychevEngine<KeyType,T>::doSSEnergyTruncation(MPS<KeyType,T>&mps,MPO<T>&mpo,const int itmax,const int Dmax,const T Ethresh,const Real landelta){
  // if(this->_verbose)
  //   cout<<"#truncating energy: ";
  //cout<<"entering doSSEnergyTruncation : "<<endl;
  if(itmax==0){
    return;
  }
  //  cout<<"                in doSSEnergyTruncation: ||mps|| before prepareSite: "<<sqrt(computeOverlap(mps,mps))<<endl;
  for(int site = 0;site<mps.size()-1;++site){
    prepareSite(mps,site,"l",true);
  }
  for(int site = mps.size()-1;site>0;--site){
    prepareSite(mps,site,"r",true);
  }
  //  cout<<"                in doSSEnergyTruncation: ||mps|| after prepareSite: "<<sqrt(computeOverlap(mps,mps))<<endl;
  map_container<int,BlockMatrix<1,KeyType,T> >R,L;
  computeR(mps,mpo,1,R);
  L.clear();
  L[-1] = getLMin(mps);
  bool converged = false;
  uint it = 0;
  while(!converged){
    if(it>=itmax){
      converged=true;
      break;
    }
    // if(this->_verbose==1){
    //   cout<<".";
    //   cout.flush();
    // }
    for(uint site=0;site<mps.size()-1;++site){
      SSLanzcosEngine<KeyType,T> 
	lan(L[site-1],R[site+1],mps[site],mpo[site],mps.getItoKey(site),this->_verbose);
      //      cout<<"SSlanzcos  at site "<<site<<" at step "<<it<<endl;
      int N=lan.simulate(Dmax,landelta,true);
      //now get the projector ...
      BlockMatrix<1,KeyType,T> A1=mps[site],A2;
      int i0 = lan.Elargerthan(T(Ethresh));
      for(int i=i0;i<N;++i){
	A2=lan.getState(i);
	T c=A2*mps[site];
	A1-=(A2*=c);
      }
      //now we should have the stripped state ...
      mps[site]=A1;
      prepareSite(mps,site,"l",true);
      L[site] = addContractionLayer(L[site-1],mps[site],mpo[site],mps.getItoKey(site),1);
    }
    for(uint site=mps.size()-1;site>0;--site){
      SSLanzcosEngine<KeyType,T> 
	lan(L[site-1],R[site+1],mps[site],mpo[site],mps.getItoKey(site),this->_verbose);
      //      cout<<"SSlanzcos  at site "<<site<<" at step "<<it<<endl;
      int N=lan.simulate(Dmax,landelta,true);
      //now get the projector ...
      BlockMatrix<1,KeyType,T> A1=mps[site],A2;
      int i0 = lan.Elargerthan(T(Ethresh));
      for(int i=i0;i<N;++i){
	A2=lan.getState(i);
	T c=A2*mps[site];
	A1-=(A2*c);
      }
      //now we should have the stripped state ...
      mps[site]=A1;
      prepareSite(mps,site,"r",true);
      R[site] = addContractionLayer(R[site+1],mps[site],mpo[site],mps.getItoKey(site),-1);  
    }
    ++it;
    if(it>itmax)
      converged=true;
  }
  if(itmax!=0){
    uint site =0;
    SSLanzcosEngine<KeyType,T> 
      lan(L[site-1],R[site+1],mps[site],mpo[site],mps.getItoKey(site),this->_verbose);
    //      cout<<"SSlanzcos  at site "<<site<<" at step "<<it<<endl;
    int N=lan.simulate(Dmax,landelta,true);
    //now get the projector ...
    BlockMatrix<1,KeyType,T> A1=mps[site],A2;
    int i0 = lan.Elargerthan(T(Ethresh));
    for(int i=i0;i<N;++i){
      A2=lan.getState(i);
      T c=A2*mps[site];
      A1-=(A2*c);
    }
    //now we should have the stripped state ...
    mps[site]=A1;
    prepareSite(mps,site,"l",true);
    R[site] = addContractionLayer(R[site+1],mps[site],mpo[site],mps.getItoKey(site),-1);  
  }
  // if(this->_verbose==1){
  //   cout<<endl;
  //   cout<<"#done ..."<<endl;
  // }
  //we don't do this anymore because i think that we loose norm ... 
  // prepareSite(mps,0,"r",true);
  // mps.setLOrtho(false);
  // mps.setROrtho(true);
}



template<class KeyType,typename T>
void ChebychevEngine<KeyType,T>::simulate(const int Nmax,const int chimax,const int Esweeps,const int Dmax,list<measurementObject<T> > &mList, 
					  const int ncmax,const Real delta, const bool reconstruct,const bool resume,const uint truncstep,
					  const Real landelta,const bool cleanup,const bool PSS,const bool debug){
  if(debug)
    cout<<"debug flag set: two different compress functions are used and results are compared to each other"<<endl;

  MPS<KeyType,T>tmps,t_test;
  uint etstep=1;
  ofstream file;
  std::string fname;
  Real normsq=0.0,conv1=0.0;
  if(!resume){
    _moments.resize(0);
    _recons_moments.resize(0);
    _conv1_container.resize(0);
    _conv2_container.resize(0);
    _etrunc_container.resize(0);
    _etrunc_step.resize(0);

    cout<<"#starting chebychev-engine ... "<<endl;
    _t0=applyMPO(_opi,this->_mps);  
    //    cout<<"#erased "<<_t0.squeeze(1e-15)<<" entries in _t0"<<endl;
    measureMPS(_t0,mList);
    etstep++;
    while(!checkMPS(_t0));
    _mpsf=applyMPO(_opf,this->_mps);
    while(!checkMPS(_mpsf));
    normsq=computeOverlap(_t0,_t0);

    _norm.push_back(sqrt(normsq));
    _moments.push_back(computeOverlap(_mpsf,_t0));
    _conv1_container.push_back(0.0);
    _conv2_container.push_back(0.0);
    if(this->_verbose==1)
      cout<<"#norm of _t0 = "<<sqrt(normsq)<<"; mu[0] = "<<_moments[0]<<"; rel.trunc.weights: "<<0.0<<", "<<0.0<<endl;
    _t1=applyMPO(this->_mpo,_t0);
    etstep++;
    while(!checkMPS(_t1));

    normsq=computeOverlap(_t1,_t1);
    _moments.push_back(computeOverlap(_mpsf,_t1));
    if(reconstruct){
      _recons_moments.push_back(2*computeOverlap(_t0,_t0)-_moments[0]);
      _recons_moments.push_back(2*computeOverlap(_t1,_t0)-_moments[1]);
      _recons_moments.push_back(2*computeOverlap(_t1,_t1)-_moments[0]);
    }

    if(debug){
      t_test=_t1;
      cout<<"||t_test|| BEFORE compressTest: "<<sqrt(computeOverlap(t_test,t_test))<<endl;
      compressTest(t_test,chimax,ncmax,delta,false);
      cout<<"||t_test|| after compressTest: "<<sqrt(computeOverlap(t_test,t_test))<<endl;
    }
    if(this->_verbose>1)
      cout<<"||t1|| BEFORE compression: "<<sqrt(computeOverlap(_t1,_t1))<<endl;

    if(!PSS){
      conv1=compress(_t1,chimax,ncmax,delta,false);
    }else{
      conv1=compressTest(_t1,chimax,ncmax,delta,false);
    }

    if(this->_verbose>1)
      cout<<"||t1|| after compression: "<<sqrt(computeOverlap(_t1,_t1))<<endl;

    if (debug){
      cout<<"Main test: sqrt(<t1|t_test>) ="<<sqrt(computeOverlap(_t1,t_test))<<endl;
    }
    _norm.push_back(sqrt(normsq));
    _conv1_container.push_back(conv1);
    _conv2_container.push_back(0.0);
    if(this->_verbose>=1)
      cout<<"#norm of _t1 = "<<sqrt(normsq)<<"; mu[1] = "<<_moments[1]<<"; rel.trunc.weights: "<<conv1<<", "<<0.0<<endl;
    measureMPS(_t1,mList);

    if(etstep%truncstep==0){
      tmps=_t1;
      if (debug){
	t_test=_t1;
      }
      if(this->_verbose>1)
	cout<<"||t1|| BEFORE energy truncation:     "<<sqrt(computeOverlap(_t1,_t1))<<endl;
      doSSEnergyTruncation(_t1,this->_mpo,Esweeps,Dmax,1.0,landelta);
      if(this->_verbose>1)
	cout<<"||t1|| after energy truncation:      "<<sqrt(computeOverlap(_t1,_t1))<<endl;
      _etrunc_container.push_back(computeOverlap(_t1-tmps,_t1-tmps));
      _etrunc_step.push_back(1);
      tmps.clear();
      if (debug){
	cout<<"||t_test|| BEFORE energy truncation: "<<sqrt(computeOverlap(t_test,t_test))<<endl;
	doSSEnergyTruncation(t_test,this->_mpo,Esweeps,Dmax,1.0,landelta);
	cout<<"||t_test|| after energy truncation:  "<<sqrt(computeOverlap(t_test,t_test))<<endl;
      }
    }
  }
  if(resume){
    this->read();
    cout<<"#resuming simulation ..."<<endl;
  }
  for(uint n=_moments.size();n<Nmax;++n){
    tmps=applyMPO(this->_mpo,_t1);
    if(debug)
      t_test=tmps;
    etstep++;
    while(!checkMPS(tmps));
    if(this->_verbose>1)
      cout<<"||H t"<<(n-1)<<"|| BEFORE compression: "<<sqrt(computeOverlap(tmps,tmps))<<endl;
    if(!PSS){
      conv1=compress(tmps,chimax,ncmax,delta,false);
    }else{
      conv1=compressTest(tmps,chimax,ncmax,delta,false);
    }
    if(this->_verbose>1)
      cout<<"||H t"<<(n-1)<<"|| after compression:  "<<sqrt(computeOverlap(tmps,tmps))<<endl;

    if (debug){
      cout<<"||(H t)_test|| BEFORE compressTest:    "<<sqrt(computeOverlap(t_test,t_test))<<endl;
      compressTest(t_test,chimax,ncmax,delta,false);
      cout<<"||(H t)_test|| after compressTest:     "<<sqrt(computeOverlap(t_test,t_test))<<endl;
      cout<<"Main test: sqrt(<H t|(H t)_test>)  =   "<<sqrt(computeOverlap(tmps,t_test))<<endl;
    }

    tmps*=T(2.0);
    _t2=tmps-_t0;
    tmps.clear();
    _t0.clear();
     normsq=computeOverlap(_t2,_t2);
    _moments.push_back(computeOverlap(_mpsf,_t2));
    if(reconstruct){
      _recons_moments.push_back(2*computeOverlap(_t2,_t1)-_moments[1]);
      _recons_moments.push_back(2*normsq-_moments[0]);
    }

    if (debug)
      t_test=_t2;
    if(this->_verbose>1)
      cout<<"||t"<<n<<"|| BEFORE compression:       "<<sqrt(computeOverlap(_t2,_t2))<<endl;
    Real conv2=0;
    if(!PSS){
      conv2=compress(_t2,chimax,ncmax,delta,false);
    }else{
      conv2=compressTest(_t2,chimax,ncmax,delta,false);
    }
    if(this->_verbose>1)
      cout<<"||t"<<n<<"|| after compression:        "<<sqrt(computeOverlap(_t2,_t2))<<endl;

    if (debug){
      cout<<"||t"<<n<<"_test|| BEFORE compressTest: "<<sqrt(computeOverlap(t_test,t_test))<<endl;
      compressTest(t_test,chimax,ncmax,delta,false);
      cout<<"||t"<<n<<"_test|| after compressTest:  "<<sqrt(computeOverlap(t_test,t_test))<<endl;

      cout<<"Main test: <t"<<n<<"|t_"<<n<<"_test>            =  "<<sqrt(computeOverlap(_t2,t_test))<<endl;
    }

    _norm.push_back(sqrt(normsq));
    _conv1_container.push_back(conv1);
    _conv2_container.push_back(conv2);

    if(this->_verbose>=1)
      cout<<"#norm of _t"<<n<<" = "<<sqrt(normsq)<<"; mu["<<n<<"] = "<<_moments[n]<<"; rel.trunc.weights: "<<conv1<<", "<<conv2<<endl;

    measureMPS(_t2,mList);
    if(etstep%truncstep==0){
      tmps=_t2;
      if(debug)
	t_test=_t2;
      if(this->_verbose>1)
	cout<<"||t"<<n<<"|| BEFORE energy truncation:     "<<sqrt(computeOverlap(_t2,_t2))<<endl;
      doSSEnergyTruncation(_t2,this->_mpo,Esweeps,Dmax,1.0,landelta);
      if(this->_verbose>1)
	cout<<"||t"<<n<<"|| after  energy truncation:     "<<sqrt(computeOverlap(_t2,_t2))<<endl;
      _etrunc_container.push_back(computeOverlap(_t2-tmps,_t2-tmps));
      _etrunc_step.push_back(n);
      tmps.clear();
      if (debug){
	cout<<"||t"<<n<<"_test|| BEFORE energy truncation:"<<sqrt(computeOverlap(t_test,t_test))<<endl;
	doSSEnergyTruncation(t_test,this->_mpo,Esweeps,Dmax,1.0,landelta);
	cout<<"||t"<<n<<"_test|| after energy truncation: "<<sqrt(computeOverlap(t_test,t_test))<<endl;
      }

    }
    _t0=_t1;
    _t1=_t2;
    _t2.clearMatrices();
    _t2.clearLambdas();
    if((n+1)%this->_checkpointStep==0){
      //write moments to file
      fname="moments_"+this->_simName+".dat";
      std::ifstream infile;
      Real a,W;
      infile.open(fname.c_str(),std::ios_base::in);
      infile>>a;
      infile>>W;
      infile.close();

      file.precision(16);
      file.open(fname.c_str(),std::ios_base::out);
      file<<a<<" "<<W<<endl;
      for(int s=0;s<_moments.size();++s)
	file<<_moments[s]<<" ";
      file<<endl;
      file.close();

      //write truncated weights to file
      fname="truncated_weights_"+this->_simName+".dat";
      file.precision(16);
      file.open(fname.c_str(),std::ios_base::out);
      for(int s=0;s<_conv1_container.size();++s)
	file<<_conv1_container[s]<<" ";
      file<<endl;
      for(int s=0;s<_conv2_container.size();++s)
	file<<_conv2_container[s]<<" ";
      file<<endl;
      file.close();
      
      //write etruncation to file
      fname="etruncation_overlap_"+this->_simName+".dat";
      file.precision(16);
      file.open(fname.c_str(),std::ios_base::out);
      for(int s=0;s<_etrunc_step.size();++s)
	file<<_etrunc_step[s]<<" "<<_etrunc_container[s]<<endl;
      file<<endl;
      file.close();

      //write norm to file
      fname="norm_"+this->_simName+".dat";
      file.precision(16);
      file.open(fname.c_str(),std::ios_base::out);
      for(int s=0;s<_norm.size();++s)
	file<<_norm[s]<<" ";
      file<<endl;
      file.close();

      
      if(reconstruct){
	fname="rmoments_"+this->_simName+".dat";
	infile.open(fname.c_str(),std::ios_base::in);
	infile>>a;
	infile>>W;
	infile.close();
      
	file.precision(16);
	file.open(fname.c_str(),std::ios_base::out);
	file<<a<<" "<<W<<endl;
	for(uint s=0;s<_recons_moments.size();++s){
	  file<<_recons_moments[s]<<" ";
	}
	file<<endl;
      	file.close();
      }
      this->write();
      this->writeMeasurements(mList);
    }
  }

  fname="moments_"+this->_simName+".dat";
  std::ifstream infile;
  Real a,W;
  infile.open(fname.c_str(),std::ios_base::in);
  infile>>a;
  infile>>W;
  infile.close();

  file.precision(16);
  file.open(fname.c_str(),std::ios_base::out);
  file<<a<<" "<<W<<endl;
  for(int s=0;s<_moments.size();++s)
    file<<_moments[s]<<" ";
  file<<endl;
  file.close();
  if(reconstruct){
    fname="rmoments_"+this->_simName+".dat";
    file.precision(16);
    file.open(fname.c_str(),std::ios_base::out);
    file<<a<<" "<<W<<endl;
    for(uint s=0;s<_recons_moments.size();++s){
      file<<_recons_moments[s]<<" ";
    }
    file<<endl;
    file.close();
  }
  if(cleanup){
    string name=this->_simName+".container";
    if(remove(name.c_str())!=0){
      cerr<<"#Error deleting "<<name<<endl;
    }else{
      cout<<"#deleting "<<name<<endl;
    }
  }
}

template<class KeyType,typename T>
void ChebychevEngine<KeyType,T>::simulateSSCompress(const int Nmax,const int chimax,const int Esweeps,const int Dmax,list<measurementObject<T> > &mList, 
						    const int ncmax,const Real delta, const bool reconstruct,const bool resume,const uint truncstep,
						    const Real landelta,const bool cleanup){

  MPS<KeyType,T>tmps;
  uint etstep=1;
  ofstream file;
  std::string fname;
  Real normsq=0.0,conv1=0.0;
  if(!resume){
    _moments.resize(0);
    _recons_moments.resize(0);
    _conv1_container.resize(0);
    _conv2_container.resize(0);
    _etrunc_container.resize(0);
    _etrunc_step.resize(0);

    cout<<"#starting chebychev-engine ... "<<endl;
    _t0=applyMPO(_opi,this->_mps);  
    //    cout<<"#erased "<<_t0.squeeze(1e-15)<<" entries in _t0"<<endl;
    measureMPS(_t0,mList);
    etstep++;
    while(!checkMPS(_t0));
    _mpsf=applyMPO(_opf,this->_mps);
    while(!checkMPS(_mpsf));
    normsq=computeOverlap(_t0,_t0);

    _norm.push_back(sqrt(normsq));
    _moments.push_back(computeOverlap(_mpsf,_t0));
    _conv1_container.push_back(0.0);
    _conv2_container.push_back(0.0);
    if(this->_verbose>=1)
      cout<<"#norm of _t0 = "<<sqrt(normsq)<<"; mu[0] = "<<_moments[0]<<"; rel.trunc.weights: "<<0.0<<", "<<0.0<<endl;
    tmps=_t0;
    conv1=compressSingleSite(tmps,this->_mpo,chimax,ncmax,delta,false);
    //conv1=compressTwoSite(tmps,this->_mpo,chimax,ncmax,delta,true,1e-12,false);
    _t1=tmps;
    tmps.clear();

    etstep++;
    while(!checkMPS(_t1));

    normsq=computeOverlap(_t1,_t1);
    _moments.push_back(computeOverlap(_mpsf,_t1));
    if(reconstruct){
      _recons_moments.push_back(2*computeOverlap(_t0,_t0)-_moments[0]);
      _recons_moments.push_back(2*computeOverlap(_t1,_t0)-_moments[1]);
      _recons_moments.push_back(2*computeOverlap(_t1,_t1)-_moments[0]);
    }
    _norm.push_back(sqrt(normsq));
    _conv1_container.push_back(conv1);
    _conv2_container.push_back(0.0);
    if(this->_verbose>=1)
      cout<<"#norm of _t1 = "<<sqrt(normsq)<<"; mu[1] = "<<_moments[1]<<"; rel.trunc.weights: "<<conv1<<", "<<0.0<<endl;
    measureMPS(_t1,mList);

    if(etstep%truncstep==0){
      tmps=_t1;
      doSSEnergyTruncation(_t1,this->_mpo,Esweeps,Dmax,1.0,landelta);
      _etrunc_container.push_back(computeOverlap(_t1-tmps,_t1-tmps));
      _etrunc_step.push_back(1);
      tmps.clear();
    }
  }
  if(resume){
    this->read();
    cout<<"#resuming simulation ..."<<endl;
  }
  for(uint n=_moments.size();n<Nmax;++n){
    tmps=_t1;
    conv1=compressSingleSite(tmps,this->_mpo,chimax,ncmax,delta,false);
    //Real conv1=compressTwoSite(tmps,this->_mpo,chimax,ncmax,delta,true,1e-12,false);
    etstep++;
    while(!checkMPS(tmps));
    tmps*=T(2.0);
    _t2=tmps-_t0;
    tmps.clear();
    _t0.clear();
     normsq=computeOverlap(_t2,_t2);
    _moments.push_back(computeOverlap(_mpsf,_t2));
    if(reconstruct){
      _recons_moments.push_back(2*computeOverlap(_t2,_t1)-_moments[1]);
      _recons_moments.push_back(2*normsq-_moments[0]);
    }
    Real conv2=compress(_t2,chimax,ncmax,delta,false);
    _norm.push_back(sqrt(normsq));
    _conv1_container.push_back(conv1);
    _conv2_container.push_back(conv2);

    if(this->_verbose>=1)
      cout<<"#norm of _t"<<n<<" = "<<sqrt(normsq)<<"; mu["<<n<<"] = "<<_moments[n]<<"; rel.trunc.weights: "<<conv1<<", "<<conv2<<endl;

    measureMPS(_t2,mList);

    if(etstep%truncstep==0){
      tmps=_t2;
      doSSEnergyTruncation(_t2,this->_mpo,Esweeps,Dmax,1.0,landelta);
      _etrunc_container.push_back(computeOverlap(_t2-tmps,_t2-tmps));
      _etrunc_step.push_back(n);
      tmps.clear();
    }
    _t0=_t1;
    _t1=_t2;
    _t2.clearMatrices();
    _t2.clearLambdas();
    if((n+1)%this->_checkpointStep==0){
      //write moments to file
      fname="moments_"+this->_simName+".dat";
      std::ifstream infile;
      Real a,W;
      infile.open(fname.c_str(),std::ios_base::in);
      infile>>a;
      infile>>W;
      infile.close();

      file.precision(16);
      file.open(fname.c_str(),std::ios_base::out);
      file<<a<<" "<<W<<endl;
      for(int s=0;s<_moments.size();++s)
	file<<_moments[s]<<" ";
      file<<endl;
      file.close();

      //write truncated weights to file
      fname="truncated_weights_"+this->_simName+".dat";
      file.precision(16);
      file.open(fname.c_str(),std::ios_base::out);
      for(int s=0;s<_conv1_container.size();++s)
	file<<_conv1_container[s]<<" ";
      file<<endl;
      for(int s=0;s<_conv2_container.size();++s)
	file<<_conv2_container[s]<<" ";
      file<<endl;
      file.close();
      
      //write etruncation to file
      fname="etruncation_overlap_"+this->_simName+".dat";
      file.precision(16);
      file.open(fname.c_str(),std::ios_base::out);
      for(int s=0;s<_etrunc_step.size();++s)
	file<<_etrunc_step[s]<<" "<<_etrunc_container[s]<<endl;
      file<<endl;
      file.close();

      //write norm to file
      fname="norm_"+this->_simName+".dat";
      file.precision(16);
      file.open(fname.c_str(),std::ios_base::out);
      for(int s=0;s<_norm.size();++s)
	file<<_norm[s]<<" ";
      file<<endl;
      file.close();

      
      if(reconstruct){
	fname="rmoments_"+this->_simName+".dat";
	infile.open(fname.c_str(),std::ios_base::in);
	infile>>a;
	infile>>W;
	infile.close();
      
	file.precision(16);
	file.open(fname.c_str(),std::ios_base::out);
	file<<a<<" "<<W<<endl;
	for(uint s=0;s<_recons_moments.size();++s){
	  file<<_recons_moments[s]<<" ";
	}
	file<<endl;
      	file.close();
      }
      this->write();
      this->writeMeasurements(mList);
    }
  }

  fname="moments_"+this->_simName+".dat";
  std::ifstream infile;
  Real a,W;
  infile.open(fname.c_str(),std::ios_base::in);
  infile>>a;
  infile>>W;
  infile.close();

  file.precision(16);
  file.open(fname.c_str(),std::ios_base::out);
  file<<a<<" "<<W<<endl;
  for(int s=0;s<_moments.size();++s)
    file<<_moments[s]<<" ";
  file<<endl;
  file.close();
  if(reconstruct){
    fname="rmoments_"+this->_simName+".dat";
    file.precision(16);
    file.open(fname.c_str(),std::ios_base::out);
    file<<a<<" "<<W<<endl;
    for(uint s=0;s<_recons_moments.size();++s){
      file<<_recons_moments[s]<<" ";
    }
    file<<endl;
    file.close();
  }
  if(cleanup){
    string name=this->_simName+".container";
    if(remove(name.c_str())!=0){
      cerr<<"#Error deleting "<<name<<endl;
    }else{
      cout<<"#deleting "<<name<<endl;
    }
  }
}




//T has to be Real!!!
template<class KeyType,typename T>
void ChebychevEngine<KeyType,T>::simulateExponential(const int Nmax,const int chimax,const Real dt, const int ncmax,const Real delta,
						     const int trotterorder,const bool reconstruct,const bool resume,const bool flip){

  if(trotterorder!=1&&trotterorder!=2){
    cout<<"#input parameter trotteroder has wrong value "<<trotterorder<<". use 1 or 2"<<endl;
    abort();
  }

  Real conv=0.0;
  MPS<KeyType,T>tmps,tmps1;
  // for(uint n=0;n<this->_mpo.size();++n)
  //   for(typename MPOMat<1,T>::iterator it=this->_mpo[n].begin();it!=this->_mpo[n].end();++it)
  //     it->second.resize(2,2);
  MPO<T>mpoeven=this->_mpo.computeNNGateMPOReal(1,this->_mps,dt),mpoodd=this->_mpo.computeNNGateMPOReal(-1,this->_mps,dt);
  MPO<T>mpoevenhalf=this->_mpo.computeNNGateMPOReal(1,this->_mps,dt/2),mpooddhalf=this->_mpo.computeNNGateMPOReal(-1,this->_mps,dt/2);

  for(uint i=0;i<(this->_mpo.size()-1);++i){
    _NNGates.push_back(this->_mpo.getNNGateReal(this->_mps,i,i+1,dt));
    if(trotterorder==2)
      _NNGatesHalf.push_back(this->_mpo.getNNGateReal(this->_mps,i,i+1,dt/2.0));
  }

  ofstream file;
  std::string fname;
  Real normsq=0.0;
  if(!resume){
    _moments.resize(0);
    _recons_moments.resize(0);
    cout<<"#starting chebychev-engine ... "<<endl;
    _t0=applyMPO(_opi,this->_mps);  
    _mpsf=applyMPO(_opf,this->_mps);
    this->_mps.clear();
    while(!checkMPS(_t0));
    _moments.push_back(computeOverlap(_mpsf,_t0));
    if(this->_verbose>=1)
      cout<<"#norm of _t0 = "<<sqrt(computeOverlap(_t0,_t0))<<"; mu[0] = "<<_moments[0]<<endl;
    conv=compress(_t0,chimax,ncmax,delta,false);
    while(!checkMPS(_t0));

    _t1=applyOperator(_t0,chimax,ncmax,delta,trotterorder,flip);  
    while(!checkMPS(_t1));
    _moments.push_back(computeOverlap(_mpsf,_t1));
    if(this->_verbose>=1)
      cout<<"#norm of _t1 = "<<sqrt(computeOverlap(_t1,_t1))<<"; mu[1] = "<<_moments[1]<<endl;
    if(reconstruct){
      cout<<"recons"<<endl;
      _recons_moments.push_back(2*computeOverlap(_t0,_t0)-_moments[0]);
      _recons_moments.push_back(2*computeOverlap(_t1,_t0)-_moments[1]);
      _recons_moments.push_back(2*computeOverlap(_t1,_t1)-_moments[0]);
    }
    compress(_t1,chimax,ncmax,delta,false);
  }
  for(uint n=_moments.size();n<Nmax;++n){
    // cout<<"# ======== doing step "<<n<<" in Chebychev engine ... ================================================"<<endl;
    _t2=applyOperator(_t1,chimax,ncmax,delta,trotterorder,flip);  
    while(!checkMPS(_t2));
    _t2*=T(2.0);
    _t2-=_t0;
    tmps.clear();
    _t0.clear();
    normsq=computeOverlap(_t2,_t2);
    _moments.push_back(computeOverlap(_mpsf,_t2));

    if(reconstruct){
      _recons_moments.push_back(2*computeOverlap(_t2,_t1)-_moments[1]);
      _recons_moments.push_back(2*normsq-_moments[0]);

    }
    conv=compress(_t2,chimax,ncmax,delta,false);
    if(this->_verbose>=1)
      cout<<"#norm of _t"<<n<<" = "<<sqrt(normsq)<<"; mu["<<n<<"] = "<<_moments[n]<<"; rel.trunc.weight: "<<conv<<endl;
    _t0=_t1;
    _t1=_t2;
    _t2.clearMatrices();
    _t2.clearLambdas();

    if(n%this->_checkpointStep==0){
      fname="ChebyshevMoments_"+this->_simName;
      file.precision(10);
      file.open(fname.c_str(),std::ios_base::out);
      for(int s=0;s<_moments.size();++s)
	file<<_moments[s]<<" ";
      file<<endl;
      file.close();
      if(reconstruct){
	fname="ReconstructedChebyshevmoments_"+this->_simName;
	file.precision(10);
	file.open(fname.c_str(),std::ios_base::out);
	for(uint s=0;s<_recons_moments.size();++s){
	  file<<_recons_moments[s]<<" ";
	}
	file<<endl;
      	file.close();
      }
      this->write();
    }
  }
  
  fname="ChebyshevMoments_"+this->_simName;
  file.open(fname.c_str(),ios_base::app);
  file.precision(10);
  file<<endl;
  file.close();
  if(reconstruct){
    fname="ReconstructedChebyshevmoments_"+this->_simName;
    file.open(fname.c_str(),ios_base::app);
    file.precision(10);
    file<<endl;
    file.close();
  }
}



template<class KeyType,typename T>
void ChebychevEngine<KeyType,T>::rawOrtho(MPS<KeyType,T> &mps,const int dir){
  if(dir>=0){
    for(uint site = 0;site<mps.size()-1;++site){
      prepareSite(mps,site,"l",true);
    }
    for(uint site = mps.size()-1;site>0;--site){
      prepareSite(mps,site,"r",true);
    }
  }
  if(dir<0){
    for(uint site = mps.size()-1;site>0;--site){
      prepareSite(mps,site,"r",true);
    }
    for(uint site = 0;site<mps.size()-1;++site){
      prepareSite(mps,site,"l",true);
    }
  }
}

template<class KeyType,typename T>
CanonizedMPS<KeyType,T> ChebychevEngine<KeyType,T>::rawCanonize(MPS<KeyType,T> &mps,const int dir){

  CanonizedMPS<KeyType,T> cmps;
  cmps.clearMatrices();
  if(dir<0){
    cmps[mps.size()-1] = mps[mps.size()-1];
    for(int site = (mps.size()-2);site>=0;--site){
      if(mps.lambda(site).size()==0){
	cout<<"#in ChebychevEngine::rawCanonize: input mps contains no lambdas at site "<<site<<endl;
	abort();
      }
      cmps[site]=mps[site]/mps.lambda(site);
    }
  }  
  if(dir>=0){
    cmps[0] = mps[0];
    for(int site = 1;site<mps.size();++site){
      if(mps.lambda(site-1).size()==0){
	cout<<"#in ChebychevEngine::rawCanonize: input mps contains no lambdas at site "<<site<<endl;
	abort();
      }
      cmps[site]=mps.lambda(site-1)/mps[site];
    }
  }

  return cmps;
}


template<class KeyType,typename T>
pair<std::vector<Real>,std::vector<uint> > ChebychevEngine<KeyType,T>::doStep(MPS<KeyType,T> &mps,const SparseOperator<2,T> &gate,const int site,
									      const int blockthreads){
  assert(site<(mps.size()-1));
  assert(mps.size());
  assert(mps.getPSpace()(site)*mps.getPSpace()(site+1)==gate.size1());
  assert(mps.getPSpace()(site)*mps.getPSpace()(site+1)==gate.size2());
  //there is no normalization done 
  pair<std::vector<Real>,std::vector<uint> >out=applyGate(mps[site],mps[site+1],gate,mps.getItoKey(site),mps.getItoKey(site+1),blockthreads);
  return out;  
};


template<class KeyType,typename T>
VectorType ChebychevEngine<KeyType,T>::doTimeSlice(MPS<KeyType,T> &mps,std::vector<SparseOperator<2,T> >&gates,const int start,const int inc,const int blockthreads){
  //  cout << "."; cout.flush();
  const int N=mps.size();
  VectorType t_weights(N-1);
  UintVectorType chi(N-1);
  boost_setTo(t_weights,0.0);
  boost_setTo(chi,uint(0));
#pragma omp parallel for
  for(int site = start;site<(N-inc+1);site+=inc){
    pair<std::vector<Real>,std::vector<uint> > out=doStep(mps,gates[site],site,blockthreads);
    t_weights(site)= out.first[0];
    chi(site) = out.second[0];
  }
  return chi;
}



//this routine is PRIVATE and can ONLY DEAL WITH HAMILTONIANS WITH NEAREST NEIGHBOR INTERACTION AT MOST. 

template<class KeyType,typename T>
MPS<KeyType,T> ChebychevEngine<KeyType,T>::applyOperator(MPS<KeyType,T>&mps,const uint chimax,const uint ncmax,const Real delta,
							 const int trotterorder,const bool flip_w,const int blockthreads){
  MPS<KeyType,T>tmps=mps,tmps2=mps;//the original state must not be changed, it will be used again later on
  Real conv1,conv2,conv3;
  if(trotterorder==1){
    //    cout<<"order 1"<<endl;
    doTimeSlice(tmps,_NNGates,1,2);//odd
    conv1=compress(tmps,chimax,ncmax,delta,false);
    doTimeSlice(tmps,_NNGates,0,2);//even
    conv2=compress(tmps,chimax,ncmax,delta,false);
    if(flip_w){
      tmps2-=tmps;
      //tmps2=mps-tmps;
      while(!checkMPS(tmps2));
      tmps=tmps2;
      tmps2.clear();
      compress(tmps,chimax,ncmax,delta,false);
    }

    cout<<"trunc. weights: "<<conv1<<", "<<conv2<<endl;;
  }else if(trotterorder==2){
    //    cout<<"order 2"<<endl;
    doTimeSlice(tmps,_NNGatesHalf,1,2);//odd
    conv1=compress(tmps,chimax,ncmax,delta,false);
    doTimeSlice(tmps,_NNGates,0,2);//even
    conv2=compress(tmps,chimax,ncmax,delta,false);
    doTimeSlice(tmps,_NNGatesHalf,1,2);//odd
    conv3=compress(tmps,chimax,ncmax,delta,false);
    if(flip_w){
      tmps2-=tmps;
      //tmps2=mps-tmps;
      while(!checkMPS(tmps2));
      tmps=tmps2;
      tmps2.clear();
      compress(tmps,chimax,ncmax,delta,false);
    }
    cout<<"trunc. weights: "<<conv1<<", "<<conv2<<", "<<conv3<<endl;

  }
  return tmps;
}



//NOT WORKING
template<class KeyType,typename T>
class KrylovEvolutionEngine:public Container<KeyType,T>{
public:
  KrylovEvolutionEngine():Container<KeyType,T>(){}
  KrylovEvolutionEngine(const Container<KeyType,T>&c):Container<KeyType,T>(c){};
  KrylovEvolutionEngine(const MPS<KeyType,T> &mps ,MPO<T> & mpo,const std::string simName):Container<KeyType,T>(mps,mpo,simName){this->_checkpointStep=50;};
  void simulate(list<measurementObject<T> >&mList,const int Nmax,const Real compdelta, const int chimax,const int ncmax);
protected:

private:
};


template<class KeyType,typename T>
void KrylovEvolutionEngine<KeyType,T>::simulate(list<measurementObject<T> >&mList,const int Nmax,const Real compdelta, const int chimax,const int ncmax){
  MPS<KeyType,T> phi_1=this->_mps,phi_2;
  std::vector<MPS<KeyType,T> >krylVecs;
  matrix<T,column_major>H(Nmax,Nmax);
  boost_setTo(H,T(0.0));
  orthonormalize(phi_1,"lr");
  krylVecs.push_back(phi_1);

  int n=1;
  while(n<Nmax){
    cout<<"at krylov step "<<n<<endl;
    phi_1=applyMPO(this->_mpo,krylVecs[n-1]);
    while(!checkMPS(phi_1));
    compress(phi_1,chimax,ncmax,compdelta);
    for(uint j=0;j<n;++j){
      H(j,n-1)=computeOverlap(krylVecs[j],phi_1);
      phi_1-=(krylVecs[j]*(H(j,n-1)));
      compress(phi_1,chimax,ncmax,compdelta);
    }
    H(n,n-1)=sqrt(computeOverlap(phi_1,phi_1));
    cout<<"that's the norm of phi_1: "<<H(n,n-1);
    compress(phi_1,chimax,ncmax,compdelta);
    orthonormalize(phi_1,"rl");
    orthonormalize(phi_1,"lr");
    cout<<"the norm of phi_1 "<<computeOverlap(phi_1,phi_1)<<endl;
    
    krylVecs.push_back(phi_1);
    ++n;
  }
  for(uint n=0;n<Nmax;++n)
    for(uint n_=0;n_<Nmax;++n_)
      cout<<"overlap of <"<<n<<"|"<<n_<<">: "<<computeOverlap(krylVecs[n],krylVecs[n_])<<endl;


}

// //this still is NOT working, ....maybe one really has to do arnoldi here...
// template<class KeyType,typename T>
// void KrylovEvolutionEngine<KeyType,T>::simulate(list<measurementObject<T> >&mList,const int tmax,const int Nmax,const Real delta, const int chimax,const Real dcomp,
// 						const int ncmax,const Real dt){
//   MPS<KeyType,T> a1=this->_mps;
//   uint it=0;
//   //  while(it<tmax){
//   MPS<KeyType,T> a0,a2,a3,t;
//   std::vector<MPS<KeyType,T> >kryl;
//   std::vector<T>ks,es;
//   orthonormalize(a1,"lr");
//   int n=0;
//   T k=0.0,e=0.0;
//   bool first=true;
//   cout<<"============"<<chimax<<endl;
//   cout<<"============"<<dcomp<<endl;
//   while(n<Nmax){
//     cout<<"at step "<<n<<endl;
//     k=sqrt(computeOverlap(a1,a1));
//     cout<<"                               that's k"<<k<<endl;
//     if(abs(k.real())<delta){break;}
//     if(!first){ks.push_back(k);}

//     orthonormalize(a1,"lr");
//     orthonormalize(a1,"rl");
//     a2=applyMPO(this->_mpo,a1);
//     e = computeOverlap(a2,a1);
//     cout<<"                               that's e before compress "<<e<<endl;
//     cout<<a2.getAuxSpace()<<endl;
//     compress(a2,chimax,ncmax,dcomp);
//     e = computeOverlap(a2,a1);
//     cout<<"                               that's e after compress "<<e<<endl;
//     cout<<a2.getAuxSpace()<<endl;
//     kryl.push_back(a1);//is this really clever? ....
//     es.push_back(e);
//     if(!first){
//       t=a2-(a1*e);
//       //	compress(t,chimax,ncmax,dcomp);
//       a3=t-(a0*k);
//       //	compress(a3,chimax,ncmax,dcomp);
//       a0=a1;
//       a1=a3;
//     }
//     if(first){
//       a0=a1;
//       a3=((a1*(-e))+a2);
//       //compress(a3,chimax,ncmax,dcomp);
//       a1=a3;
//       first = false;
//     }
//     ++n;
//   }
//   a1.clear();a2.clear();a3.clear();t.clear();
//   const int size = es.size();
//   double*eig2 = new double[size*size];//this really has to be a double, not a T or something like that ...
//   tridiag(&es[0],es.size(),&ks[0],eig2,'V');
//   matrix<Real,column_major> eigs(size,size);//this also has to be double ...
//   memcpy(&eigs.data()[0],eig2,size*size*sizeof(double));
//   delete[] eig2;
//   //now get the coefficients for the expansion
//   for(uint u=0;u!=size;++u){
//     for(uint v=0;v!=size;++v)
//       cout<<"overlap of state "<<u<<"with state "<<v<<" = "<<computeOverlap(kryl[u],kryl[v])<<endl;
//   }
//   abort();
//   for(uint u=0;u!=size;++u){
//     //run over all eigenvectors ...
//     Complex c(1.0,0.0);
//     for(uint v=0;v!=size;++v){
//       c+=(exp(Complex(0,-(es[v].real()-es[0].real())*dt))*eigs(u,v)*eigs(v,0));
//     }
//     a1+=(kryl[u]*=c);
//     compress(a1,chimax,ncmax,dcomp);
//   }
//   compress(a1,chimax,ncmax,dcomp);
//   cout<<a1.getAuxSpace()<<endl;
//   orthonormalize(a1,"lr");
//   orthonormalize(a1,"rl");
//   cout<<a1.getAuxSpace()<<endl;
//   this->_cmps=canonize(a1);
//   doCMPSmeasurement(mList);
//   writeMeasurements(mList);    
// }
//that's it ..
//}


//only for fermi-hubbard like systems
//template<class KeyType,typename T>
//class TEVEngine:public Container<KeyType,T>{
//public:
//  TEVEngine():Container<KeyType,T>(){}
//  TEVEngine(const Container<KeyType,T>&c):Container<KeyType,T>(c){};
//  //although the system has a star geometry, the gates can be computed using
//  //the mpo-routines for chain-like geometry; this facilitates coding; 
//  TEVEngine(const MPS<KeyType,T> &mps ,const std::string filename,const std::string simname){
//    this->_simName=simname;
//    this->_mps=mps;
//    std::fstream infile;
//    infile.open(filename.c_str(),ios_base::in);
//    //TODO: read in the parameters;
//    _swap.resize(4,4);
//    _swap[OpKeyType<2>(0,0)]=1.0;
//    _swap[OpKeyType<2>(1,2)]=1.0;
//    _swap[OpKeyType<2>(2,1)]=1.0;
//    _swap[OpKeyType<2>(3,3)]=-1.0;//this is the fermionic minus sign;
//
//  }
//
//  void simulate(list<measurementObject<T> >&mList,const int start,const int stop,const T dt, const int chimax,const Real tr_weight);
//  virtual void toarchive(boost::archive::binary_oarchive &oar); //in derived classes, these two have to be modified
//  virtual void fromarchive(boost::archive::binary_iarchive &iar);//
//
//protected:
//  //add stuff
//private:
//  MatrixType _bathparams;
//  std::vector<SparseOperator<2,T> >_Gates;
//  SparseOperator<2,T> _swap;  
//  CVectorType _hop,_U,_mu;
//  uint _N,_T;
//};
//


//only for fermi-hubbard like systems
template<class KeyType,typename T>
class TEVEngine:public Container<KeyType,T>{
public:
  TEVEngine():Container<KeyType,T>(){}
  TEVEngine(const Container<KeyType,T>&c):Container<KeyType,T>(c){};
  //although the system has a star geometry, the gates can be computed using
  //the mpo-routines for chain-like geometry; this facilitates coding; 
  TEVEngine(const MPS<KeyType,T> &mps,const std::string filename,const std::string simname,const Real U){
    this->_simName=simname;
    this->_mps=mps;
    _localu=U;
    _twMeasurementStep=10;
    //read in parameters
    std::ifstream infile;
    infile.open(filename.c_str());
    infile>>_T;infile>>_N;
    cout<<"#reading parameters from file "<<filename<<": number of columns (T)="<<_T<<", number of rows (N)="<<_N<<endl;
    if(this->_mps.size()!=_N+1){
      cout<<"#in TEVEngine<KeyType,T>::TEVEngine(const MPS<KeyType,T> &mps,const std::string filename,const std::string simname): this->_mps.size()!=_N+1; aborting ..."<<endl;
      abort();
    }
    _bathparams.resize(_T,_N);
    boost_setTo(_bathparams,0.0);
    for(uint t=0;t<_T;++t){
      for(uint n=0;n<_N;++n){
	infile>>_bathparams(t,n);
      }
    }
    cout<<"#done"<<endl;
  }
  
  void simulate(list<measurementObject<T> >&mList,const int start,const int stop,const Real dt, const int chimax,const Real tr_weight,
		const int ncmax,const Real thresh);
  virtual void toarchive(boost::archive::binary_oarchive &oar){
    Container<KeyType,T>::toarchive(oar);
    oar<<_twMeasurementStep;
    oar<<_localu;
    oar<<_bathparams;

  }
  virtual void fromarchive(boost::archive::binary_iarchive &iar){
    Container<KeyType,T>::fromarchive(iar);
    iar>>_twMeasurementStep;
    iar>>_localu;
    iar>>_bathparams;
  }

  void setTwMeasurementStep(const uint gTwMeasurementStep){this->_twMeasurementStep=gTwMeasurementStep;}
  uint getTwMeasurementStep(){return this->_twMeasurementStep;}
  void sort(std::pair<int,int> &p,const uint t,const Real thresh);

protected:
  void computeGates(const uint t,const Real dt);
  void doTimeSlice(std::pair<int,int>p,const uint t,const Real dt);
  //add stuff
private:
  MatrixType _bathparams;
  std::vector<SparseOperator<2,T> >_Gates;
  SparseOperator<2,T> _swap;  
  VectorType _hop,_U,_mu;
  uint _N,_T;
  Real _localu;
  uint _twMeasurementStep;
};

/*
  1) there is a fine grid 
  2) all parameters of the hamiltonian are stored in a file 
  3) the site ordering is fixed in the code
  -parameters are stored in a boost matrix
  -routine for computing the gates doesn't start from mpo-description
  -define the swap gates
  -fix the way they are applied
       
  4) first only test the usual time evolution
*/

//this function will be called at every time step;
//ordering; first gate is impurity-first site interaction, second gate is impurity-second site interaction, ...
//t is a time index
template<class KeyType,typename T>
void TEVEngine<KeyType,T>::computeGates(const uint t,const Real dt){
  //clear the old gates!
  _Gates.clear();
  //take parameters at time-step index "t" and create an mpo from it to get the gates
  if(t>_bathparams.size1()){
    cout<<"#in TEVengine::computeGates(uint t, Real dt): t>_bathparams.size1(); aborting"<<endl;
    abort();
  }
  _hop.resize(3);
  _mu.resize(4);
  _U.resize(3);
  boost_setTo(_mu,0.0);
  boost_setTo(_U,0.0);
  for(uint site=0;site<_N;++site){
    if(site==_N/2){
      _U(1)=_localu;
    }else if (site>0){
      _U(1)=0;
    }
    _hop(1)=_bathparams(t,site);
  
    SpinlessFermionsMPO<T> mpo(_U,_hop,_mu);
    SparseOperator<2,T>g;
    if(site<_N-1)
      g=mpo.getNNGate(this->_mps,1,2,dt);
    else if(site==_N-1)
      g=mpo.getNNGate(this->_mps,1,2,dt);
    g.squeeze(1e-14);
    _Gates.push_back(g);
  }
 
}

template<class KeyType,typename T>
void TEVEngine<KeyType,T>::doTimeSlice(std::pair<int,int>p,const uint t,const Real dt){
  this->computeGates(t,dt);
  if(p.second>=_N){
    cout<<"reached end of MPS-length"<<endl;
    abort();
  }
  //at this point, the gates and the swaps are already known;
  //alternate: gate(1,2),swap(1,2),gate(2,3),swap(2,3),gate(3,4),swap(3,4),...
  //the first time step is applied going from left to right; at the end, the impurity site is at the site L-2
  for(uint site=p.first;site<p.second;++site){
    applyGate(this->_mps[site],this->_mps[site+1],this->_Gates[site],this->_mps.getItoKey(site),this->_mps.getItoKey(site+1));
    if(site<_N-1)
      applyGate(this->_mps[site],this->_mps[site+1],_swap,this->_mps.getItoKey(site),this->_mps.getItoKey(site+1));
  }
  for(int site=p.second-1;site>=p.first;--site){//dont change int to uint!!!!
    applyGate(this->_mps[site],this->_mps[site+1],this->_Gates[site],this->_mps.getItoKey(site),this->_mps.getItoKey(site+1));
    if(site>0)
      applyGate(this->_mps[site-1],this->_mps[site],_swap,this->_mps.getItoKey(site-1),this->_mps.getItoKey(site));
    //    applyGate(this->_mps[site],this->_mps[site+1],_swap,this->_mps.getItoKey(site),this->_mps.getItoKey(site+1));
  }
}
//for now, oscillating bath parameters are not implemented
template<class KeyType,typename T>
void TEVEngine<KeyType,T>::sort(std::pair<int,int> &p,const uint t,const Real thresh){
  //check if the fist bond has dropped below threshhold
  if(_bathparams(t,p.first)<thresh){
    applyGate(this->_mps[p.first],this->_mps[p.first+1],_swap,this->_mps.getItoKey(p.first),this->_mps.getItoKey(p.first+1));
    p.first++;
  }
  if(_bathparams(t,p.second+1)>thresh){
    p.second++;
    if(p.second>=_N){
      cout<<"in sort(std::pair<int,int> p,const uint t,const Real thresh): cannot increase p.second further; reached end of MPS"<<endl;
    }
  }
  
}

template<class KeyType,typename T>
void TEVEngine<KeyType,T>::simulate(list<measurementObject<T> >&mList,const int start, const int stop,const Real dt, const int chimax,
				    const Real tr_weight,const int ncmax,const Real thresh){
  //the size of the gates is hard coded, this works only for spinless fermions
  std::pair<int,int> inds(0,0);
  cout<<"startin TEBD"<<endl;
  _swap.resize(4,4);
  std::vector<uint> in(2),out(2);

  in[0]=0;in[1]=0;out[0]=0;out[1]=0;
  OpKeyType<2>k1(out,in);
  _swap[k1]=1.0;

  in[0]=1;in[1]=0;out[0]=0;out[1]=1;
  OpKeyType<2>k2(out,in);
  _swap[k2]=1.0;

  in[0]=0;in[1]=1;out[0]=1;out[1]=0;
  OpKeyType<2>k3(out,in);
  _swap[k3]=1.0;

  in[0]=1;in[1]=1;out[0]=1;out[1]=1;
  OpKeyType<2>k4(out,in);
  _swap[k4]=-1.0;

  ofstream tw_outfile;
  std::string fname = "truncatedWeight" + this->_simName;
  std::vector<Real> trw;
  
  for(uint step=start;step!=stop;++step){
    //check if the any of the bonds is zero; if so, shift it to the left end of the mps
    sort(inds,step,thresh);
    doTimeSlice(inds,step,dt);
    //now compress; has to be checked!
    Complex conv=compress(this->_mps,chimax,ncmax,tr_weight,false,true);

    while(!checkMPS(this->_mps));

    if(step%_twMeasurementStep==0){	//append tw -lists
      trw.push_back(conv.real());
      trw.push_back(conv.imag());
    }
    if(step%this->_checkpointStep==0){	//access hdd only at those steps where it is anyway done by measurement
      tw_outfile.open(fname.c_str(),ios_base::app);
      for(uint s=0;s<trw.size()/2;++s)
	tw_outfile<<trw[2*s]<<trw[2*s+1]<<endl;
      tw_outfile.close();
      trw.clear();
    }
    //step%2 ensures that state is always measure in the left-to-right ordering situation (remember that we swap that state)
    if(step%2==1){
      doMPSmeasurement(mList);
    }
    //save the MPS for checkpointing
    if (step%this->_checkpointStep==0){
      std::string name=this->_simName+"_CP";
      this->save();
      this->write(name);
      writeMeasurements(mList);
      cout << endl;
    }

  }
  writeMeasurements(mList);
}



//THIS IS STILL UNDER CONSTRUCTION

//template<class KeyType,typename T>
//class TDVPengine:public Container<KeyType,T>{
//public:
//  TDVPengine():Container<KeyType,T>(){_blockthreads=1;_tebdthreads=1;_twMeasurementStep=1;_nameit=1;_storeallcp=false;_TDVPstep=1;}
//  TDVPengine(std::string simname):Container<KeyType,T>(simname)
//  {_blockthreads=1;_tebdthreads=1;_twMeasurementStep=1;_nameit=1;_storeallcp=false;_TDVPstep=1;}
//  TDVPengine(std::string simname,const bool storeallcp):Container<KeyType,T>(simname)
//  {_blockthreads=1;_tebdthreads=1;_twMeasurementStep=1;_nameit=1;_storeallcp=storeallcp;_TDVPstep=1;}
//
//  TDVPengine(const Container<KeyType,T>&c):Container<KeyType,T>(c){_blockthreads=1;_tebdthreads=1;
//    _twMeasurementStep=1;_nameit=1;_storeallcp=false;_TDVPstep=1;}
//  TDVPengine(const MPS<KeyType,T> &mps ,MPO<T> & mpo,const std::string simName,const bool storeallcp=false):Container<KeyType,T>(mps,mpo,simName){
//    this->_checkpointStep=50;
//    _blockthreads=1;
//    _tebdthreads=1;
//    _twMeasurementStep=1;
//    _nameit=1;
//    _storeallcp=storeallcp;
//    _TDVPstep=1;
//  }
//  
//  void simulateSingleSite(list<measurementObject<T> > &mList,const std::string order,const uint MAXSTEPS=100,const T dt=T(0.05),
//		const uint chimax=500, const bool doTrunc=false,const Real tr_weight=1e-12,const bool normalize=true,const bool recanonize=false,
//		const int blockthreads=1,const int tebdthreads=1);
//  void simulateTwoSite(list<measurementObject<T> > &mList,const std::string order,const uint MAXSTEPS=100,const T dt=T(0.05),
//		const uint chimax=500, const bool doTrunc=false,const Real tr_weight=1e-12,const bool normalize=true,const bool recanonize=false,
//		const int blockthreads=1,const int tebdthreads=1);
//
//  void measureCMPS(list<measurementObject<T> > &m){Container<KeyType,T>::measureCMPS(m);}
//  virtual void clear(){Container<KeyType,T>::clear();_NNGates.clear();_NNGatesHalf.clear();}
//  void setTwMeasurementStep(const uint gTwMeasurementStep){this->_twMeasurementStep=gTwMeasurementStep;}
//  uint getTwMeasurementStep(){return this->_twMeasurementStep;}
//
//  void storecheckpoints(const bool store){_storeallcp=store;}
//
//protected:
//  virtual void toarchive(boost::archive::binary_oarchive &oar){
//    Container<KeyType,T>::toarchive(oar);
//    oar<<_blockthreads;
//    oar<<_tebdthreads;
//    oar<<_twMeasurementStep;
//    oar<<_nameit;
//    oar<<_storeallcp;
//    oar<<_TDVPstep;
//    
//  }
//  virtual void fromarchive(boost::archive::binary_iarchive &iar){
//    Container<KeyType,T>::fromarchive(iar);
//    iar>>_blockthreads;
//    iar>>_tebdthreads;
//    iar>>_twMeasurementStep;
//    iar>>_nameit;
//    iar>>_storeallcp;
//    iar>>_TDVPstep;
//  }
//  void doCMPSmeasurement(list<measurementObject<T> > &m)const{Container<KeyType,T>::doCMPSmeasurement(m);}
//  bool doCMPSmeasurement(list<measurementObject<T> > &m,const uint chimax,const bool doTrunc,const Real tr_weight, bool saveIteration,VectorType &tw,
//			 const bool normalize=true,const bool recanonize=false,const bool skipsecondhalfstep=false);
//  //for Next NN interactions
//  bool doCMPSmeasurementNNN(list<measurementObject<T> > &m,const uint chimax,const bool doTrunc,const Real tr_weight, bool saveIteration,VectorType &twe,
//			    VectorType &two,const bool normalize=true,const bool recanonize=false,const bool skipsecondhalfstep=false);
//
//private:
//  int _blockthreads;
//  int _tebdthreads;
//  uint _twMeasurementStep;
//  uint _nameit;
//  bool _storeallcp;
//  uint _TDVPstep;//stores the current TDVP step of the simulate routine; when later restarted from a read() command, the routine starts from 
//                 //the corresponding position
//};
//
//template<class KeyType,typename T>
//class FineGrainer:
//public:
//FineGrainer(const MPS<KeyType,T>mps) {}
//private:
//MPS<KeyType,T>_mps;
//};
//
#endif
  
