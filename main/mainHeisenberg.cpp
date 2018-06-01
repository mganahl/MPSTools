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
  parser.addoption("B0", Parameter_Float, "0.0");      
  parser.addoption("curvature", Parameter_Float, "0.0");      
  parser.addoption("n0", Parameter_Float, "0.0");      
  parser.addoption("doTrunc", Parameter_Option);
  parser.addoption("simname", Parameter_String, "HeisenbergXXZ");
  parser.addoption("measure", Parameter_Integer,"10");
  parser.addoption("2strsite", Parameter_Integer,"-1");
  parser.addoption("corrsite", Parameter_Integer,"-1");
  parser.addoption("Nup", Parameter_Integer,"10");
  parser.addoption("load", Parameter_Option);

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
  const double B0  = parser.getfloat("B0");
  const double curv  = parser.getfloat("curvature");
  const double n0  = parser.getfloat("n0");
  const bool doTrunc = parser.isgiven("doTrunc");
  const string simname = parser.getstring("simname");
  const uint measStep = parser.getint("measure");
  const int TwoStrSite = parser.getint("2strsite");
  const int corrSite = parser.getint("corrsite");
  const uint QN = parser.getint("Nup");
  const bool load = parser.isgiven("load");
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
  for(int n=0;n<B.size();++n)
    B(n) = curv*B0*(n-n0)*(n-n0);
  HeisenbergMPO<Real> mpoTrap(Jz,Jxy,B);
  HeisenbergMPO<Complex> mpoTrapComplex(Jz,Jxy,B);
  if(TwoStrSite>=0){
    B(TwoStrSite)+=10;
    B(TwoStrSite+1)+=10;
  }
  HeisenbergMPO<Real> mpoTrap2Str(Jz,Jxy,B);

//==========================measurements================================0
  std::vector<Operator<Complex>*> SZ,CURRENT,ENTROPY,UUPROJECTOR,UDPROJECTOR,DUPROJECTOR;
  std::vector<Operator<Real>*> CORR;

  for(int i=0;i<N;++i){
    SparseLocalOperator<Complex> Sz(2,2);
    SparseLocalOperator<Real> Szr(2,2);
    Sz[OpKeyType<1>(0,0)]=-0.5;Sz[OpKeyType<1>(1,1)]=0.5;
    Szr[OpKeyType<1>(0,0)]=-0.5;Szr[OpKeyType<1>(1,1)]=0.5;
    Sz.setSite(i);
    Szr.setSite(i);
    SZ.push_back(new SparseLocalOperator<Complex>(Sz));
    if(corrSite>=0){
      MPO<Real> corr;
      MPOMat<1,Real> mat;
      mat[ukey<2>(0,0)]=Szr;
      mat.Dl=1;mat.Dr=1;
      corr[corrSite]=mat;
      if(i!=corrSite){
	MPOMat<1,Real> mat;
	mat[ukey<2>(0,0)]=Szr;
	mat.Dl=1;mat.Dr=1;
	corr[i]=mat;
      }else{
	MPOMat<1,Real> mat;
	mat[ukey<2>(0,0)]=Szr*Szr;
	mat.Dl=1;mat.Dr=1;
	corr[i]=mat;
      }
      CORR.push_back(new MPO<Real>(corr));
    }
    if(i<(N-1)){
      //==============bipartite entanglement============
      Entropy<Complex> ent(i);
      ENTROPY.push_back(new Entropy<Complex>(ent));
      //==============current for Heisenenberg ===========================
      //the sign may be wrong!!!!!
      Sparse2SiteOperator<Complex> Current(4,4);
      std::vector<uint> in(2),out(2);

      out[0] = 1;out[1] = 0;
      in[0]  = 0;in[1]  = 1;
      Current[OpKeyType<2>(out,in)]=Complex(0,1);
      Current.setSites(i,i+1);

      out[0] = 0;out[1] = 1;
      in[0] = 1;in[1] = 0;
      Current[OpKeyType<2>(out,in)]=Complex(0,-1);
      
      CURRENT.push_back(new Sparse2SiteOperator<Complex>(Current));
      //===========================projectors============================
      Sparse2SiteOperator<Complex> UUProj(4,4),UDProj(4,4),DUProj(4,4);

      out[0] = 1;out[1] = 1;
      in[0]  = 1;in[1]  = 1;
      UUProj[OpKeyType<2>(out,in)]= Complex(1.0);
      UUProj.setSites(i,i+1);

      out[0] = 1;out[1] = 0;
      in[0]  = 1;in[1]  = 0;
      UDProj[OpKeyType<2>(out,in)]= Complex(1.0);
      UDProj.setSites(i,i+1);

      out[0] = 0;out[1] = 1;
      in[0]  = 0;in[1]  = 1;
      DUProj[OpKeyType<2>(out,in)]= Complex(1.0);
      DUProj.setSites(i,i+1);

      UUPROJECTOR.push_back(new Sparse2SiteOperator<Complex>(UUProj));
      UDPROJECTOR.push_back(new Sparse2SiteOperator<Complex>(UDProj));
      DUPROJECTOR.push_back(new Sparse2SiteOperator<Complex>(DUProj));
    }
  }
  string filename = simname;
  measurementObject<Complex> mSz(SZ,"Sz_"+simname,measStep);
  measurementObject<Complex> mCurrent(CURRENT,"Current_"+simname,measStep);
  measurementObject<Complex> mEntropy(ENTROPY,"entropy_"+simname,measStep);
  measurementObject<Complex> mUUPr(UUPROJECTOR,"UUProj_"+simname,measStep);
  measurementObject<Complex> mUDPr(UDPROJECTOR,"UDProj_"+simname,measStep);
  measurementObject<Complex> mDUPr(DUPROJECTOR,"DUProj_"+simname,measStep);
  measurementObject<Real> mCorr(CORR,"SzSzCorr_"+simname,measStep);
  
  list<measurementObject<Complex> > mList;
  list<measurementObject<Real> > mRList;
  mList.push_back(mSz);
  mList.push_back(mCurrent);
  mList.push_back(mEntropy);
  mList.push_back(mUUPr);
  mList.push_back(mUDPr);
  mList.push_back(mDUPr);
  //  mRList.push_back(mCorr);
  //============state initialization========
  UintVectorType dimP(N),localState(N);
  boost_setTo(dimP,uint(2));
  boost_setTo(localState,uint(0));
  int mod=static_cast<int>(floor(1.0*N/QN));
  std::vector<i_to_key<AbKeyType > > itokey(N);
  uint qn = 0;
  for(uint site = 0;site<N;site++){
    AbKeyType k1(0,1);    
    AbKeyType k2(1,1);
    itokey[site][0] =k1; 
    itokey[site][1] =k2;
    // itokey[site].setSign(-1);
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
  //  SimContainer<AbKeyType,Real> HSim(mps,mpo,simname);
  DMRGengine<AbKeyType,Real> HSim(mps,mpo,simname);
  if(maxstepsDMRG>0){
    HSim.simulate(5,1e-13,60,doTrunc,deltaDMRG);
    HSim.setMPO(mpoTrap);
    HSim.setName(simname+"withTrap");
    HSim.simulate(5,1e-13,60,doTrunc,deltaDMRG);
    HSim.setName(simname+"TwoStringwithTrap");
    HSim.setMPO(mpoTrap2Str);
    HSim.simulate(maxstepsDMRG,1e-13,chiDMRG,doTrunc,deltaDMRG);
  }
  HSim.canonizeMPS();
  mCorr.measure(HSim.getCMPS());
  mCorr.writeMeasurements();
  MPS<AbKeyType,Complex>cmps=cast<AbKeyType,Complex,Real>(HSim.getMPS());
  TEBDengine<AbKeyType,Complex> HSimComplex(cmps,mpoTrapComplex,simname);

  cout<<"total quantum number: "<<HSimComplex.getMPS().getQN()<<endl;
  orthonormalize(HSim.getMPS(),"lr");
  orthonormalize(HSimComplex.getMPS(),"lr");
  HSimComplex.canonizeMPS();
  if(load)
    HSimComplex.load();
  HSimComplex.setMPO(mpoTrapComplex);
  cout<<"is state OK? "<<checkMPS(HSimComplex.getMPS())<<endl;
  HSimComplex.measureCMPS(mList);
  HSimComplex.simulate(mList,"second",maxstepsTEBD,dt,chiTEBD,doTrunc,deltaTEBD);
  return 0;
}
