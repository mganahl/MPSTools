#include "MPSTypes.hpp"
#include "MPSfunctions.hpp"
#include "Hamiltonians.hpp"
#include "BlockTensor.hpp"
#include "tensormap.hpp"
#include "container.hpp"
#include "utilities.hpp"
#include "parparser.hpp"
#include <iostream>




int main(const int argc,const char **argv){
  parparser parser;

  parser.addoption("Nmax", Parameter_Integer,"100");
  parser.addoption("chiDMRG", Parameter_Integer,"100");
  parser.addoption("NmaxDMRG", Parameter_Integer,"5");
  parser.addoption("rescalingfactor", Parameter_Float, "2.0");
  parser.addoption("label", Parameter_String,"default_name");
  parser.addoption("ncomp", Parameter_Integer,"3");
  parser.addoption("write", Parameter_Integer,"5");
  parser.addoption("Esweeps", Parameter_Integer,"5");
  parser.addoption("Dmax", Parameter_Integer,"25");
  parser.addoption("Etruncstep", Parameter_Integer,"1");
  parser.addoption("chicomp", Parameter_Integer,"200");
  parser.addoption("compdelta", Parameter_Float,"1e-8");
  parser.addoption("landelta", Parameter_Float,"1e-10");
  parser.addoption("spin", Parameter_Integer,"1");//for the down-spin spectral function
  parser.addoption("branch", Parameter_Integer,"1");//for the hole part of the spectral function
  parser.addoption("verbosity", Parameter_Integer,"0");//for the hole part of the spectral function
  parser.addoption("cleanup", Parameter_Option);//for the hole part of the spectral function
  
  // parser.addoption("up",Parameter_Integer,"0");
  // parser.addoption("down",Parameter_Integer,"0");

  parser.addoption("testRes", Parameter_Option);
  parser.addoption("deltaDMRG", Parameter_Float,"1e-12");
  parser.addoption("conv", Parameter_Float,"1e-10");  

  if (!parser.load(argc, argv))
    return 0;
  std::cout << "# Options:\n" << parser << std::endl;


  const uint Nmax = parser.getint("Nmax");
  const uint chiDMRG = parser.getint("chiDMRG");
  const uint maxstepsDMRG = parser.getint("NmaxDMRG");
  const Real a = parser.getfloat("rescalingfactor");//rescaling parameter, scales the bandwidth by a*D
  const bool doTrunc = true;//parser.isgiven("doTrunc");
  string simname = parser.getstring("label");
  const uint ncomp = parser.getint("ncomp");
  const uint step = parser.getint("write");
  const uint etrstep = parser.getint("Etruncstep");
  const uint Esweeps = parser.getint("Esweeps");
  const uint Dmax = parser.getint("Dmax");
  const uint chicomp = parser.getint("chicomp");
  const double compdelta = parser.getfloat("compdelta");
  const double landelta = parser.getfloat("landelta");
  const int spin=parser.getint("spin");
  const int branch=parser.getint("branch");
  const int verbosity = parser.getint("verbosity");
  const uint verboseCheb=verbosity>=1?1:0;
  const uint verboseDMRG=verbosity>=2?1:0;
  const bool cleanup=parser.isgiven("cleanup");
  const Real W = 0.0;
  const bool recons = true;
  // const int changeup = parser.getint("up");
  // const int changedown = parser.getint("down");


  const bool testenergy=parser.isgiven("testRes")?false:true;
  const double deltaDMRG = parser.getfloat("deltaDMRG");
  const double conv = parser.getfloat("conv");


  //consistency check:


  
  string bath = "bath_params_"+simname;
  string H0 = "H_local_"+simname;


  bool part,up;
  if (spin==1)
    up=true;
  else if(spin==-1)
    up=false;
  else if((spin!=1)&&(spin!=-1)){
    cout<<"#unknown value of input variable spin; use 1 or -1"<<endl;
    abort();
  }
  if(branch==1)
    part=true;
  else if(branch==-1)
    part=false;
  else if((branch!=1)&&(branch!=-1)){
    cout<<"#unknown value of input variable branch; use 1 or -1"<<endl;
    abort();
  }

//  if(part&&up)
//    simname+="pu";
//  else if(part&&(!up))
//    simname+="pd";
//  else if(!part&&up)
//    simname+="hu";
//  else if((!part)&&(!up))
//    simname+="hd";
//
  //check if there's a simname.container file to be resumed ....
  string name=simname+".container";
  cout<<"#starting solver using "<<omp_get_max_threads()<<" parallel threads"<<endl;
  if(fexist(name)){
    ChebychevEngine<AbKeyType,Real> Cheb(simname,verboseCheb);
    Cheb.read();
    std::string oname=Cheb.getName();
    if(oname!=simname){
      cout<<"############################### .container file has been renamed. ########################## "<<endl;
      Cheb.setName(simname);
    }
    int N=Cheb.getMPS().size();
    if(N%2!=0){
      cout<<"#warning: using an odd chain length!"<<endl;
    }

    std::vector<Operator<Real>*> DEN,ENTROPY,HAM,NORM;
    for(int i = 0;i<N;++i){
      SparseLocalOperator<Real> n(2,2);
      n[OpKeyType<1>(1,1)]=(1.0);
      n.setSite(i);
      DEN.push_back(new SparseLocalOperator<Real>(n));
      if(i<(N-1)){
	//==============bipartite entanglement============
	Entropy<Real> ent(i);
	ENTROPY.push_back(new Entropy<Real>(ent));
      }
    }
    HAM.push_back(new MPO<Real>(Cheb.getMPO()));
    NORM.push_back(new Norm<Real>());

    measurementObject<Real> mN(DEN,"N"+simname,1);
    measurementObject<Real> mEntropy(ENTROPY,"entropy_"+simname,1);
    measurementObject<Real> mHam(HAM,"Ham_"+simname,10);
    measurementObject<Real> mNorm(NORM,"Norm_"+simname,10);

    list<measurementObject<Real> > mList;

    mList.push_back(mN);  
    mList.push_back(mEntropy);
    mList.push_back(mNorm);
    mList.push_back(mHam);


    
    cout<<"#resuming simulation ..."<<endl;
    Cheb.simulate(Nmax,chicomp,Esweeps,Dmax,mList,ncomp,compdelta,recons,true,etrstep,landelta,cleanup);
    cout<<"#finished simulation ..."<<endl;
  }else{
    //==========always save the parameters=============================
    std::ofstream pfile;
    string pfilename="params_"+simname;
    pfile.open(pfilename.c_str(),ios_base::app);
    pfile<<"====================="+simname+"================="<<endl;
    pfile<<parser<<endl;
    pfile.close();
    //================read the hopping parameters from infile========================
    std::fstream infile;
    if(!fexist(bath)){
      cout<<"#error reading "<<bath<<": no such file or directory"<<endl;
      return 1;
    }
    infile.open(bath.c_str(),ios_base::in);
    if(!infile.good()){
      cout<<"#error reading "<<bath<<": file is broken"<<endl;
      return 2;
    }
    cout<<"#reading data from "<<bath<<endl;
    int lbath;
    infile>>lbath;
    int N=lbath+1;
    VectorType hop(2*N-1),U(2*N),mu(2*N);
    boost_setTo(U,0.0);boost_setTo(mu,0.0);boost_setTo(hop,0.0);
    //================read hoppings from file ==========================
    Real hyb,D;  
    infile>>D;
    infile>>hyb;
    hop(N)= hyb;
    hop(N-2)= hyb;
    for(uint site=1;site<(N-1);++site){
      Real tn,en;
      infile>>tn;
      infile>>en;
      hop(N+site)=tn;
      hop(N-2-site)=tn;

      //because of an inconsistency the chemical potential for spinless fermions used for the unfolded case has to have the opposite sign to 
      //be equivalent to the SIAM case;
      mu(N+site)=-en;
      mu(N-2-site+1)=-en;
    }
    Real en;
    infile>>en;
    mu(mu.size()-1)=-en;
    mu(0)=-en;
    infile.close();
    cout<<"#done..."<<endl;
    //========read the local hamiltonian H0; order: first down-electron chemical potential, then the up-electron chemical potential, then hubbard U======


    if(!fexist(H0)){
      cout<<"#error reading "<<H0<<": no such file or directory"<<endl;
      return 1;
    }
    infile.open(H0.c_str(),ios_base::in);
    if(!infile.good()){
      cout<<"#error reading "<<H0<<": file is broken"<<endl;
      return 2;
    }
    cout<<"#reading data from "<<H0<<endl;
    Real e0d,e0u,u;
    infile>>e0d;
    infile>>e0u;
    infile>>u;
    infile.close();
    cout<<"#done..."<<endl;
    U(N-1)=u;
    //because of an inconsistency the chemical potential for spinless fermions used for the unfolded case has to have the opposite sign to 
    //be equivalent to the SIAM case;
    mu(N-1)=-e0d;
    mu(N)=-e0u;
    cout<<"#following parameters will be used in the solver:"<<endl;
    cout<<"#t: "<<hop<<endl;
    cout<<"#u: "<<U<<endl;
    cout<<"#mu: "<<mu<<endl;

    SpinlessFermionsMPO<Real> mpo(U,hop,mu);

    //============state initialization========
    UintVectorType dimP,localState;
    std::vector<i_to_key<AbKeyType > > itokey;
    dimP.resize(2*N);
    localState.resize(2*N);
    boost_setTo(dimP,uint(2));
    boost_setTo(localState,uint(0));
    itokey.resize(2*N);
    for(uint site = 0;site<2*N;site++){
      AbKeyType k1(0,0,2);    
      AbKeyType k2(1,0,2);
      AbKeyType k3(0,1,2);
      if(site<N){
	itokey[site][0] =k1; 
	itokey[site][1] =k2;
      }else if(site>=N){
	itokey[site][0] =k1; 
	itokey[site][1] =k2;
      }
    }
    for(uint site = 0;site<N;site++){
      localState(site)=site%2;
      localState(N+site)=site%2;//this is half filling
    }

    //localState(site)=site%2;
    //localState(N+site)=site%2;//this is half filling
//    switch (initial){
//    case 0:
//      for(uint site = 0;site<N;site++){
//	localState(site)=site%2;
//	localState(N+site)=site%2;//this is half filling
//      }
//      break;
//    case 1:
//      for(uint site = 0;site<N;site++){
//	localState(site)=site%2;
//	localState(N+site)=(site+1)%2;//this is half filling
//      }
//      break;
//    case 2:
//      for(uint site = 0;site<N;site++){
//	localState(site)=(site+1)%2;
//	localState(N+site)=(site+1)%2;//this is half filling
//      }
//      break;
//    case 3:
//      for(uint site = 0;site<N;site++){
//	localState(site)=(site+1)%2;
//	localState(N+site)=site%2;//this is half filling
//      }
//      break;
//    }
//
    cout<<"#the initial state for dmrg:"<<endl;
    cout<<localState<<endl;

    MPS<AbKeyType,Real>mps(dimP,localState,itokey),tmps;
    DMRGengine<AbKeyType,Real> dmrg(mps,mpo,simname,verboseDMRG);
    Real energy=0.0;
    if(maxstepsDMRG>0){
      energy=dmrg.simulate(maxstepsDMRG,conv,chiDMRG,doTrunc,deltaDMRG,landelta,testenergy);
    }
    cout<<"#GS-energy after two-site DMRG: "<<energy<<endl;
    std::string sssname=simname+"_singleSite";
    SSDMRGengine<AbKeyType,Real> SSdmrg(dmrg.getMPS(),mpo,sssname,verboseDMRG);
    if(maxstepsDMRG>0){
      energy=SSdmrg.simulate(maxstepsDMRG*2,1e-10,1e-12,true);
    }
    cout<<"#GS-energy after single-site DMRG: "<<energy<<endl;
    mps=SSdmrg.getMPS();
    energy=mpo.measure(mps);
    cout<<"#number of particles: "<<mps.getQN()<<endl;
    Real res=a*D;
    // cout<<"#got groundstate energy E0 of "<<energy<<endl;
    // cout<<"#additional groundstate shift W by "<<W<<endl;
    VectorType hop_shifted(2*N-1),U_shifted(2*N),mu_shifted(2*N);
    boost_setTo(U_shifted,0.0);boost_setTo(mu_shifted,0.0);boost_setTo(hop_shifted,0.0);
    for(uint s=0;s<U_shifted.size();++s){
      U_shifted(s)=U(s)/res;
      mu_shifted(s)=mu(s)/res;
    }
    for(uint s=0;s<hop_shifted.size();++s)
      hop_shifted(s)=hop(s)/res;
    SpinlessFermionsMPO<Real> mposhifted(U_shifted,hop_shifted,mu_shifted);
    mposhifted.shiftandrescale(energy/res,0.0,0.0);
    mpo=mposhifted;
    std::ofstream outfile;
    std::string fname="moments_"+simname+".dat";
    outfile.open(fname.c_str(),ios_base::out);

    outfile<<res<<" "<<W<<endl;

    outfile.close();

    fname="rmoments_"+simname+".dat";  
    outfile.open(fname.c_str(),ios_base::out);
    outfile<<res<<" "<<W<<endl;
    outfile.close();

    energy=mpo.measure(mps);
    //cout<<"#got rescaled energy of "<<energy<<endl;
    while(!checkMPS(mps));
    // cout<<"#norm of mps without normalization: "<<computeOverlap(mps,mps)<<endl;
    // cout<<"#total quantum number mps: "<<mps.getQN()<<endl;
    // cout<<"#is mps OK? "<<checkMPS(mps)<<endl;
    orthonormalize(mps,"lr");
    orthonormalize(mps,"rl");
    while(!checkMPS(mps));
    // cout<<"#norm of mps after normalization:"<<sqrt(computeOverlap(mps,mps))<<endl;

    MPO<Real> Op;
    for(uint site=0;site<2*N;++site){
      MPOMat<1,Real>mat;
      SparseLocalOperator<Real> O_11(2,2),cD(2,2),c(2,2),p(2,2);
      O_11[OpKeyType<1>(0,0)]=1.0;O_11[OpKeyType<1>(1,1)]=1.0;
      p[OpKeyType<1>(0,0)]=1.0;p[OpKeyType<1>(1,1)]=-1.0;
      cD[OpKeyType<1>(1,0)] = 1.0;
      c[OpKeyType<1>(0,1)] = 1.0;
      if(site<(N-1)){
	mat.clear();
	mat[ukey<2>(0,0)]=p;
	mat.Dl=1;mat.Dr=1;
	Op[site]=mat;
	continue;
      }
      if(site==N-1){
	mat.clear();
	if(part&&up){
	  cout<<"#applying up particle creation"<<endl;
	  mat[ukey<2>(0,0)]=cD;
	}
	else if(!part&&up){
	  cout<<"#applying up hole creation"<<endl;
	  mat[ukey<2>(0,0)]=c;
	}
	else if(!up)
	  mat[ukey<2>(0,0)]=p;
      }
      if(site==N){
	mat.clear();
	if(part&&!up){
	  cout<<"#applying down particle creation"<<endl;
	  mat[ukey<2>(0,0)]=cD;
	}
	else if(!part&&!up){
	  cout<<"#applying down hole creation"<<endl;
	  mat[ukey<2>(0,0)]=c;
	}
	else if(up)
	  mat[ukey<2>(0,0)]=O_11;
      }
      if(site>N){
	mat.clear();
	mat[ukey<2>(0,0)]=O_11;
      }
      mat.Dl=1;mat.Dr=1; 
      Op[site]=mat;
    }
  

    std::vector<Operator<Real>*> DEN,ENTROPY,HAM,NORM;
    for(int i = 0;i<2*N;++i){
      SparseLocalOperator<Real> n(2,2);
      n[OpKeyType<1>(1,1)]=(1.0);
      n.setSite(i);
      DEN.push_back(new SparseLocalOperator<Real>(n));
      if(i<(2*N-1)){
	//==============bipartite entanglement============
	Entropy<Real> ent(i);
	ENTROPY.push_back(new Entropy<Real>(ent));
      }
    }
    HAM.push_back(new MPO<Real>(mpo));
    NORM.push_back(new Norm<Real>());

    measurementObject<Real> mN(DEN,"N"+simname,10);
    measurementObject<Real> mEntropy(ENTROPY,"entropy_"+simname,10);
    measurementObject<Real> mHam(HAM,"Ham_"+simname,10);
    measurementObject<Real> mNorm(NORM,"Norm_"+simname,10);

    list<measurementObject<Real> > mList;

    mList.push_back(mN);  
    mList.push_back(mEntropy);
    mList.push_back(mNorm);
    mList.push_back(mHam);


    ChebychevEngine<AbKeyType,Real> Cheb(mps,mpo,Op,Op,simname,verboseCheb);
    Cheb.canonizeMPS();
    Cheb.measureCMPS(mList);
    Cheb.writeMeasurements(mList);
    Cheb.setCheckpointStep(step);
    Cheb.simulate(Nmax,chicomp,Esweeps,Dmax,mList,ncomp,compdelta,recons,false,etrstep,landelta,cleanup);
    cout<<"#finished simulation ..."<<endl;
  }
  cout<<"#exiting qis.exe ..."<<endl;
  return 0;
}
