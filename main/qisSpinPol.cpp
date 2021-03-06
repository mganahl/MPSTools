#include "MPSTypes.hpp"
#include "MPSfunctions.hpp"
#include "Hamiltonians.hpp"
#include "BlockTensor.hpp"
#include "tensormap.hpp"
#include "container.hpp"
#include "utilities.hpp"
#include "parparser.hpp"
#include <iostream>



//left chain is down, right chain is up
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
  const Real res = parser.getfloat("rescalingfactor");//rescaling parameter, scales the bandwidth by a*D
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


  
  string bath_up = "bath_params_up"+simname;
  string bath_down = "bath_params_down"+simname;
  string init_up = "filling_up"+simname;
  string init_down = "filling_down"+simname;
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
    std::fstream infile_up,infile_down,infile;
    if(!fexist(bath_up)){
      cout<<"#error reading "<<bath_up<<": no such file or directory"<<endl;
      return 1;
    }
    if(!fexist(bath_down)){
      cout<<"#error reading "<<bath_down<<": no such file or directory"<<endl;
      return 1;
    }

    infile_up.open(bath_up.c_str(),ios_base::in);
    if(!infile_up.good()){
      cout<<"#error reading "<<bath_up<<": file is broken"<<endl;
      return 2;
    }
    infile_down.open(bath_down.c_str(),ios_base::in);
    if(!infile_down.good()){
      cout<<"#error reading "<<bath_down<<": file is broken"<<endl;
      return 2;
    }

    cout<<"#reading data from "<<bath_up<<endl;
    int lbath_up,lbath_down;
    infile_up>>lbath_up;
    infile_down>>lbath_down;
    int N_up=lbath_up+1;
    int N_down=lbath_down+1;
    //int N=N_up+Ndown
    VectorType hop(N_up+N_down-1),U(N_up+N_down),mu(N_up+N_down);
    boost_setTo(U,0.0);boost_setTo(mu,0.0);boost_setTo(hop,0.0);
    //================read hoppings from file ==========================
    Real hyb_up,D_up,hyb_down,D_down;  
    infile_up>>D_up;
    infile_up>>hyb_up;

    infile_down>>D_down;
    infile_down>>hyb_down;

    hop(N_down)= hyb_up;
    hop(N_down-2)= hyb_down;
    for(uint site=1;site<(N_down-1);++site){
      Real tn_down,en_down;
      infile_down>>tn_down;
      infile_down>>en_down;
      hop(N_down-2-site)=tn_down;
      //because of an inconsistency the chemical potential for spinless fermions used for the unfolded case has to have the opposite sign to 
      //be equivalent to the SIAM case;
      mu(N_down-2-site+1)=-en_down;
    }

    for(uint site=1;site<(N_up-1);++site){
      Real tn_up,en_up;
      infile_up>>tn_up;
      infile_up>>en_up;
      hop(N_down+site)=tn_up;

      //because of an inconsistency the chemical potential for spinless fermions used for the unfolded case has to have the opposite sign to 
      //be equivalent to the SIAM case;
      mu(N_down+site)=-en_up;
    }

    Real en_up,en_down;
    infile_up>>en_up;
    infile_down>>en_down;
    mu(mu.size()-1)=-en_up;
    mu(0)=-en_down;
    infile_up.close();
    infile_down.close();
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
    U(N_down-1)=u;
    //because of an inconsistency the chemical potential for spinless fermions used for the unfolded case has to have the opposite sign to 
    //be equivalent to the SIAM case;
    mu(N_down-1)=-e0d;
    mu(N_down)=-e0u;
    cout<<"#following parameters will be used in the solver:"<<endl;
    cout<<"#t: "<<hop<<endl;
    cout<<"#u: "<<U<<endl;
    cout<<"#mu: "<<mu<<endl;

    SpinlessFermionsMPO<Real> mpo(U,hop,mu);

    //============state initialization========
    if(!fexist(init_up)){
      cout<<"#error reading "<<init_up<<": no such file or directory"<<endl;
      return 1;
    }
    infile_up.open(init_up.c_str(),ios_base::in);
    if(!infile_up.good()){
      cout<<"#error reading "<<init_up<<": file is broken"<<endl;
      return 2;
    }

    if(!fexist(init_down)){
      cout<<"#error reading "<<init_down<<": no such file or directory"<<endl;
      return 1;
    }
    infile_down.open(init_down.c_str(),ios_base::in);
    if(!infile_down.good()){
      cout<<"#error reading "<<init_down<<": file is broken"<<endl;
      return 2;
    }


    UintVectorType dimP,localState;
    std::vector<i_to_key<AbKeyType > > itokey;
    dimP.resize(N_up+N_down);
    localState.resize(N_up+N_down);
    boost_setTo(dimP,uint(2));
    boost_setTo(localState,uint(0));
    itokey.resize(N_up+N_down);
    for(uint site = 0;site<N_up+N_down;site++){
      AbKeyType k1(0,1);    
      AbKeyType k2(1,1);
      itokey[site][0] =k1; 
      itokey[site][1] =k2;
    }
    int sitedown=0;
    while (!infile_down.eof()){
      int occ;
      infile_down>>occ;
      localState(N_down-1-sitedown)=occ;
      sitedown++;
    }
    int siteup=0;
    while (!infile_up.eof()){
      int occ;
      infile_up>>occ;
      localState(N_down+siteup)=occ;
      siteup++;
    }

    if (N_down!=sitedown){
      cout<<"length of "<<init_down<<"="<<sitedown<<" is different from length given in "<<bath_down<<"="<<N_down<<endl;
      return 7;
    }

    if (N_up!=siteup){
      cout<<"length of "<<init_up<<"="<<siteup<<" is different from length given in "<<bath_up<<"="<<N_up<<endl;
      return 7;
    }

    cout<<"#the initial state for dmrg:"<<endl;
    cout<<localState<<endl;
    cout<<"lengths of chains: left chain (down): "<<N_down<<endl;
    cout<<"                    right chain (up): "<<N_up<<endl;
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

    // cout<<"#got groundstate energy E0 of "<<energy<<endl;
    // cout<<"#additional groundstate shift W by "<<W<<endl;
    VectorType hop_shifted(N_up+N_down-1),U_shifted(N_up+N_down),mu_shifted(N_up+N_down);
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
    for(uint site=0;site<N_up+N_down;++site){
      MPOMat<1,Real>mat;
      SparseLocalOperator<Real> O_11(2,2),cD(2,2),c(2,2),p(2,2);
      O_11[OpKeyType<1>(0,0)]=1.0;O_11[OpKeyType<1>(1,1)]=1.0;
      p[OpKeyType<1>(0,0)]=1.0;p[OpKeyType<1>(1,1)]=-1.0;
      cD[OpKeyType<1>(1,0)] = 1.0;
      c[OpKeyType<1>(0,1)] = 1.0;
      if(site<(N_down-1)){
	mat.clear();
	mat[ukey<2>(0,0)]=p;
	mat.Dl=1;mat.Dr=1;
	Op[site]=mat;
	continue;
      }
      if(site==N_down-1){
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
      if(site==N_down){
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
      if(site>N_down){
	mat.clear();
	mat[ukey<2>(0,0)]=O_11;
      }
      mat.Dl=1;mat.Dr=1; 
      Op[site]=mat;
    }
  

    std::vector<Operator<Real>*> DEN,ENTROPY,HAM,NORM;
    for(int i = 0;i<N_up+N_down;++i){
      SparseLocalOperator<Real> n(2,2);
      n[OpKeyType<1>(1,1)]=(1.0);
      n.setSite(i);
      DEN.push_back(new SparseLocalOperator<Real>(n));
      if(i<(N_up+N_down-1)){
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
