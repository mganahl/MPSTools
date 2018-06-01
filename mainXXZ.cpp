#include "MPSTypes.hpp"
#include "MPSfunctions.hpp"
#include "helperfunctions.hpp"
#include "Hamiltonians.hpp"
#include "BlockTensor.hpp"
#include "tensormap.hpp"
#include "container.hpp"
#include "utilities.hpp"
#include "parparser.hpp"
#include <iostream>

int main(const int argc,const char **argv){
  parparser parser;
  parser.addoption("chiTEBD", Parameter_Integer,"100");
  parser.addoption("chiDMRG", Parameter_Integer,"100");
  parser.addoption("NmaxTEBD", Parameter_Integer,"1000");
  parser.addoption("NmaxDMRG", Parameter_Integer,"5");
  parser.addoption("deltaDMRG", Parameter_Float,"1e-12");
  parser.addoption("deltaTEBD", Parameter_Float,"1e-12");
  parser.addoption("dt", Parameter_Float,"0.05");
  parser.addoption("N", Parameter_Integer, "100");   
  parser.addoption("Jz", Parameter_Float, "1.0");      
  parser.addoption("Jxy", Parameter_Float, "1.0");
  parser.addoption("doTrunc", Parameter_Option);
  parser.addoption("simname", Parameter_String, "SpinlessFermions");
  parser.addoption("measure", Parameter_Integer,"10");
  parser.addoption("Nup", Parameter_Integer,"10");
  parser.addoption("load", Parameter_String,"n");
  parser.addoption("save", Parameter_Option);
  parser.addoption("density", Parameter_Option);
  parser.addoption("entropy", Parameter_Option);
  parser.addoption("overlap", Parameter_Option);
  parser.addoption("energy", Parameter_Option);
  parser.addoption("hole", Parameter_Option);
  parser.addoption("gauss", Parameter_Option);
  parser.addoption("flipsite", Parameter_Integer,"-1");
  parser.addoption("x0", Parameter_Float,"-1.0");
  parser.addoption("k0", Parameter_Float,"0.0");
  parser.addoption("sigma", Parameter_Float,"0.01");
  parser.addoption("target", Parameter_Integer,"0");

  if (!parser.load(argc, argv))
    return 0;
  std::cout << "# Options:\n" << parser << std::endl;
  const uint chiTEBD = parser.getint("chiTEBD");
  const uint chiDMRG = parser.getint("chiDMRG");
  const uint maxstepsTEBD = parser.getint("NmaxTEBD");
  const uint maxstepsDMRG = parser.getint("NmaxDMRG");
  const double deltaDMRG = parser.getfloat("deltaDMRG");
  const double deltaTEBD = parser.getfloat("deltaTEBD");
  const uint N = parser.getint("N");
  const double jz = parser.getfloat("Jz");
  const double jxy = parser.getfloat("Jxy");
  const double dt = parser.getfloat("dt");
  const bool doTrunc = parser.isgiven("doTrunc");
  const string simname = parser.getstring("simname");
  const uint measStep = parser.getint("measure");
  const uint QN = parser.getint("Nup");
  const string load = parser.getstring("load");
  const bool save = parser.isgiven("save");
  const bool dens = parser.isgiven("density");
  const bool ent = parser.isgiven("entropy");
  const bool overlap = parser.isgiven("overlap");
  const bool ham = parser.isgiven("energy");
  const bool part = !parser.isgiven("hole");
  const double x0 = parser.getfloat("x0");
  const double k0 = parser.getfloat("k0");
  const double sigma = parser.getfloat("sigma");
  const int Osite=parser.getint("flipsite");
  const bool gauss=parser.isgiven("gauss");
  const int target=parser.getint("target");
  //==========always save the parameters=============================
  std::ofstream pfile;
  string pfilename="params_"+simname;
  pfile.open(pfilename.c_str(),ios_base::app);
  pfile<<"====================="+simname+"================="<<endl;
  pfile<<parser<<endl;
  pfile.close();

  VectorType Jz(N-1),Jxy(N-1),B(N);
  boost_setTo(Jz,jz);  boost_setTo(Jxy,jxy);  boost_setTo(B,0.0);
  HeisenbergMPO<Real> mpo(Jz,Jxy,B);
  HeisenbergMPO<Complex> cmpo(Jz,Jxy,B);
  

  //============state initialization========
  UintVectorType dimP(N),localState(N),localState2(N);
  boost_setTo(dimP,uint(2));
  boost_setTo(localState,uint(0));
  boost_setTo(localState2,uint(0));
  int mod=static_cast<int>(floor(1.0*N/QN));
  std::vector<i_to_key<AbKeyType > > itokey(N);
  uint qn = 0;
  for(uint site = 0;site<N;site++){
    AbKeyType k1(0,1);    
    AbKeyType k2(1,1);
    itokey[site][0] =k1; 
    itokey[site][1] =k2;
    if(!(site%mod)&&qn<QN){
      localState(site) = 1;
      ++qn;
    }
    else{
      localState(site) = 0;
    }
  }
  if(qn!=QN){
    cout<<"could not generate filling"<<endl;
    return 0;
  }


  MPS<AbKeyType,Real>mps(dimP,localState,itokey);
  DMRGengine<AbKeyType,Real> HSim(mps,mpo,simname);
  if(maxstepsDMRG>0){
    HSim.simulate(maxstepsDMRG,1e-13,chiDMRG,doTrunc,deltaDMRG,1e-12,true,false,target);
  }
  if(save){
    cout<<"#saving SimContainer "<<simname<<endl;
    HSim.write();
  }
  if(load!="n"){
    cout<<"#loading SimContainer "<<load<<endl;
    HSim.setName(load);
    HSim.read();
  }
  
  MPS<AbKeyType,Complex>cmps = cast<AbKeyType,Complex,Real>(HSim.getMPS());
  MPS<AbKeyType,Complex>tmps;
  if(gauss&&(Osite<0)){
    gaussMPOXXZ gaussMPO(N,x0,k0,sigma);
    tmps = applyMPO(gaussMPO,cmps);
    cmps=tmps;
    while(!checkMPS(cmps));
  }
  if(Osite>=0&&!gauss){
    MPO<Complex> Op;  
    for(uint site=0;site<N;++site){
      MPOMat<1,Complex>mat;
      SparseLocalOperator<Complex> O_11(2,2),Sp(2,2),Sm(2,2);
      O_11[OpKeyType<1>(0,0)]=1.0;O_11[OpKeyType<1>(1,1)]=1.0;
      Sp[OpKeyType<1>(1,0)] = 1.0;
      if(site==Osite){
	mat.clear();
	cout<<"applying up particle creation"<<endl;
	mat[ukey<2>(0,0)]=Sp;
      }
      if(site<Osite)
	mat[ukey<2>(0,0)]=O_11;
      if(site>Osite)
	mat[ukey<2>(0,0)]=O_11;
      mat.Dl=1;mat.Dr=1; 
      Op[site]=mat;
    }
    tmps = applyMPO(Op,cmps);
    cmps=tmps;
    while(!checkMPS(cmps));
  }

  //==========================measurements================================
  std::vector<Operator<Complex>*> DENS,ENTROPY,OV,HAM;

  HAM.push_back(new MPO<Complex>(cmpo));
  for(int i = 0;i<N;++i){
    SparseLocalOperator<Complex> n(2,2);
    n[OpKeyType<1>(0,0)]=-0.5;n[OpKeyType<1>(1,1)]=0.5;
    n.setSite(i);
    DENS.push_back(new SparseLocalOperator<Complex>(n));
    if(i<(N-1)){
      //==============bipartite entanglement============
      Entropy<Complex> ent(i);
      ENTROPY.push_back(new Entropy<Complex>(ent));
    }
  }
  cout<<"#norm of Spr|0> before normalization: "<<sqrt(computeOverlap(cmps,cmps))<<endl;
  cout<<"#total quantum number of the state: "<<cmps.getQN()<<endl;

  orthonormalize(cmps,"rl");
  orthonormalize(cmps,"rl");
  orthonormalize(cmps,"lr");
  orthonormalize(cmps,"rl");
  orthonormalize(cmps,"lr");
  for(uint site=0;site<N;++site){
    DiagBlockMatrix<AbKeyType,Real > lam=cmps.lambda(site);
    lam.normalize();
    cmps.lambda(site)=lam;    
    cout<<"norm of lambda at site "<<site<<"  :"<<cmps.lambda(site).norm()<<endl;
  }
  while(!checkMPS(cmps));
  cout<<"#norm of Sp|0> after normalization: "<<sqrt(computeOverlap(cmps,cmps))<<endl;
  cout<<"#is state OK? "<<checkMPS(cmps)<<endl;
 
  TEBDengine<AbKeyType,Complex> HSimC(cmps,cmpo,simname);
  Overlap<Complex> ov(cmps);
  OV.push_back(new Overlap<Complex>(ov));
  string filename = simname;

  measurementObject<Complex> mDens(DENS,"n_"+simname,measStep);
  measurementObject<Complex> mEntropy(ENTROPY,"entropy_"+simname,measStep);
  measurementObject<Complex> mOverlap(OV,"overlap_"+simname,1);
  measurementObject<Complex> mHam(HAM,"Ham_"+simname,measStep);
  list<measurementObject<Complex> > mList;
  if(dens)
    mList.push_back(mDens);
  if(ent)
    mList.push_back(mEntropy);
  if(overlap)
    mList.push_back(mOverlap);
  if(ham)
      mList.push_back(mHam);

  orthonormalize(HSimC.getMPS(),"rl");
  orthonormalize(HSimC.getMPS(),"lr");
  HSimC.canonizeMPS();

  cout<<"#is state OK? "<<checkMPS(HSimC.getCMPS())<<endl;
  HSimC.measureCMPS(mList);
  HSimC.writeMeasurements(mList);
  HSimC.simulate(mList,"second",maxstepsTEBD,dt,chiTEBD,doTrunc,deltaTEBD);  
  cout<<"#finished successfully"<<endl;
  return 0;
}
