#include "MPSTypes.hpp"
#include "MPSfunctions.hpp"
#include "Hamiltonians.hpp"
#include "BlockTensor.hpp"
#include "tensormap.hpp"
#include "container.hpp"
//#include "containertypes.hpp"
#include "utilities.hpp"
#include "parparser.hpp"
#include <iostream>




int main(const int argc,const char **argv){
  cout<<"##########################################################################################################################"<<endl;
  cout<<"#######                     TWO-ORBITAL DMFT SOLVER FOR THE HUBBARD MODEL ON THE BETHE LATTICE                         ###"<<endl;
  cout<<"#######            IMPORTANT NOTE: FOR SPIN-Z=0 AND HALFFILLING, EACH CHAIN HAS TO BE OF EVEN LENGTH                   ###"<<endl;
  cout<<"##########################################################################################################################"<<endl;
  cout<<"#######  Usage: --label=name: label of the input files; bath_params_orbital_0_'name' and bath_params_orbital_1_'name'  ###"<<endl;
  cout<<"#######                       contain the bathparameters                                                               ###"<<endl;
  cout<<"#######                       H_local_'name' contains local Hamiltonian; structure of the file:                        ###"<<endl;
  cout<<"#######                       e0d,e0u,e1d,e1u,U,J                                                                      ###"<<endl;
  cout<<"##########################################################################################################################"<<endl;
  cout<<"##########################################################################################################################"<<endl;
  cout<<"##########################################################################################################################"<<endl;
  cout<<"##########################################################################################################################"<<endl;
  cout<<"##########################################################################################################################"<<endl;
  cout<<"##########################################################################################################################"<<endl;
  cout<<"##########################################################################################################################"<<endl;
  cout<<"##########################################################################################################################"<<endl;


  parparser parser;

  parser.addoption("NmaxTEBD", Parameter_Integer,"500");
  parser.addoption("chiDMRG", Parameter_Integer,"100");
  parser.addoption("NmaxDMRG", Parameter_Integer,"5");
  parser.addoption("NmaxSSDMRG", Parameter_Integer,"10");
  parser.addoption("label", Parameter_String,"default_name");
  parser.addoption("order", Parameter_String,"second");
  parser.addoption("measure", Parameter_Integer,"1");
  parser.addoption("checkpoint", Parameter_Integer,"50");
  parser.addoption("chiTEBD", Parameter_Integer,"200");
  parser.addoption("landelta", Parameter_Float,"1e-10");
  parser.addoption("dt", Parameter_Float,"0.05");
  parser.addoption("trweight", Parameter_Float,"1e-8");
  parser.addoption("spin", Parameter_Integer,"1");//for the down-spin spectral function
  parser.addoption("branch", Parameter_Integer,"1");//for the hole part of the spectral function
  parser.addoption("orbital", Parameter_Integer,"0");//choose an orbital
  parser.addoption("verbose", Parameter_Option);
  parser.addoption("testRes", Parameter_Option);
  parser.addoption("deltaDMRG", Parameter_Float,"1e-12");
  parser.addoption("conv", Parameter_Float,"1e-10");  
  parser.addoption("nonormalize", Parameter_Option);

  if (!parser.load(argc, argv))
    return 0;
  std::cout << "# Options:\n" << parser << std::endl;


  const uint NmaxTEBD = parser.getint("NmaxTEBD");
  const uint chiDMRG = parser.getint("chiDMRG");
  const uint maxstepsDMRG = parser.getint("NmaxDMRG");
  const uint maxstepsSSDMRG = parser.getint("NmaxSSDMRG");
  string simname = parser.getstring("label");
  string order = parser.getstring("order");
  const uint step = parser.getint("measure");
  const uint cpstep = parser.getint("checkpoint");
  const uint chimax = parser.getint("chiTEBD");
  const double landelta = parser.getfloat("landelta");
  const double dt = parser.getfloat("dt");
  const double tr_weight = parser.getfloat("trweight");
  const int spin=parser.getint("spin");
  const int branch=parser.getint("branch");
  const int orbital=parser.getint("orbital");
  const bool verbose=parser.isgiven("verbose");
  const bool testenergy=parser.isgiven("testRes")?false:true;
  const double deltaDMRG = parser.getfloat("deltaDMRG");
  const double conv = parser.getfloat("conv");
  const bool normalize=parser.isgiven("nonormalize")?false:true;

  //consistency check:
  string bath_0 = "bath_params_orbital_0_"+simname;
  string bath_1 = "bath_params_orbital_1_"+simname;
  string H0 = "H_local_"+simname;
  std::string temp;
  if(orbital==0)
    temp="_orbital_0_"+simname;
  if(orbital==1)
    temp="_orbital_1_"+simname;
  simname=temp;
  bool part,up,orb;
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

  //check if there's a simname.container file to be resumed ....
  string name=simname+"_CP.container";
  cout<<"#starting two-orbital solver using "<<omp_get_max_threads()<<" parallel threads"<<endl;
  if(fexist(name)){
    string loadname=simname+"_CP";
    TEBDengine<AbKeyType,Complex> SIAM(loadname);
    SIAM.read(loadname);
    int N=SIAM.getMPO().size();
    if(N%2!=0){
      cout<<"#warning: using an odd chain length!"<<endl;
    }



    SparseLocalOperator<Complex> nu(4,4),nd(4,4);   //elementary,identity and fermi sign operators
    //elementary,identity and fermi sign operators
    nu[OpKeyType<1>(2,2)]=(1.0);nu[OpKeyType<1>(3,3)]=(1.0);
    nd[OpKeyType<1>(1,1)]=(1.0);nd[OpKeyType<1>(3,3)]=(1.0);
    std::vector<Operator<Complex>*> NU,ND,ENTROPY,OV,NUND;
    for(int i = 0;i<N;++i){
      nu.setSite(i);
      nd.setSite(i);
      NU.push_back(new SparseLocalOperator<Complex>(nu));
      ND.push_back(new SparseLocalOperator<Complex>(nd));
      NUND.push_back(new SparseLocalOperator<Complex>(nu*nd));
      if(i<(N-1)){
	//==============bipartite entanglement============
	Entropy<Complex> ent(i);
	ENTROPY.push_back(new Entropy<Complex>(ent));
      }
    }
    MPS<AbKeyType,Real>mps;
    std::string name=simname+"_initial_state";
    mps.load(name);
    MPS<AbKeyType,Complex>cmps=cast<AbKeyType,Complex,Real>(mps);
    mps.clear();

    Overlap<Complex> ov(cmps);
    OV.push_back(new Overlap<Complex>(ov));
    measurementObject<Complex> mOverlap(OV,"overlap"+simname,step,true,true);
    measurementObject<Complex> mNu(NU,"Nu"+simname,step);
    measurementObject<Complex> mNuNd(NUND,"NuNd"+simname,step);
    measurementObject<Complex> mNd(ND,"Nd"+simname,step);
    measurementObject<Complex> mEntropy(ENTROPY,"entropy"+simname,step);

    list<measurementObject<Complex> > mList;
    mList.push_back(mOverlap);  
    mList.push_back(mNu);  
    mList.push_back(mNuNd);  
    mList.push_back(mNd);  
    mList.push_back(mEntropy);
    //bool normalize=true;
    bool recanonize=false;
    cout<<"#resuming simulation in "<<name<<endl;
    SIAM.simulate(mList,order,NmaxTEBD,dt,chimax,true,tr_weight,normalize,recanonize);
    cout<<"#finished simulation ..."<<endl;
  }else{
    //==========always save the parameters=============================
    std::ofstream pfile;
    string pfilename="params_"+simname;
    pfile.open(pfilename.c_str(),ios_base::app);
    pfile<<"====================="+simname+"================="<<endl;
    pfile<<parser<<endl;
    pfile.close();
    twoOrbitalSIAMMPO<Real> mpo(bath_0,bath_1,H0);
    twoOrbitalSIAMMPO<Complex> cmpo(bath_0,bath_1,H0);
    uint N=mpo.size();
    uint N0=mpo.N0();
    uint N1=mpo.N1();
    // uint N=8,N0=4,N1=4;
    // VectorType hop(N-1),U(N),chUp(N),chDown(N);

    // boost_setTo(hop,1.0);
    // boost_setTo(U,0.0);
    // boost_setTo(chUp,0.0);
    // boost_setTo(chDown,0.0);
    // FermiHubbardMPO<Real>mpo(hop,U,chUp,chDown);
    // FermiHubbardMPO<Complex>cmpo(hop,U,chUp,chDown);

    //============state initialization========
    UintVectorType dimP,localState;
    std::vector<i_to_key<AbKeyType > > itokey;
    dimP.resize(N);
    localState.resize(N);
    boost_setTo(dimP,uint(4));
    boost_setTo(localState,uint(0));
    itokey.resize(N);
    for(uint site = 0;site<N;site++){
      AbKeyType k1(0,0,2);
      AbKeyType k2(0,1,2);
      AbKeyType k3(1,0,2);
      AbKeyType k4(1,1,2);
      itokey[site][0] = k1;
      itokey[site][1] = k2;
      itokey[site][2] = k3;
      itokey[site][3] = k4;
    }
    for(uint site = 0;site<N;site++){
      localState(site)=site%2+1;
    }

    cout<<"#the initial state for dmrg:"<<endl;
    cout<<localState<<endl;

    MPS<AbKeyType,Real>mps(dimP,localState,itokey),tmps;
    DMRGengine<AbKeyType,Real> dmrg(mps,mpo,simname,verbose);
    Real energy=0.0;
    if(maxstepsDMRG>0){
      bool doTrunc=true;
      energy=dmrg.simulate(maxstepsDMRG,conv,chiDMRG,doTrunc,deltaDMRG,landelta,testenergy);
    }
    cout<<"#GS-energy after two-site DMRG: "<<energy<<endl;
    std::string sssname=simname+"_singleSite";
    SSDMRGengine<AbKeyType,Real> SSdmrg(dmrg.getMPS(),mpo,sssname,verbose);
    if(maxstepsDMRG>0){
      energy=SSdmrg.simulate(maxstepsSSDMRG,1e-10,1e-12,true);
    }
    mpo.shiftandrescale(energy,0.0,0.0);
    cmpo.shiftandrescale(Complex(energy,0.0),Complex(0.0,0.0),Complex(0.0,0.0));	
    mps=SSdmrg.getMPS();

    Real shifted_energy=mpo.measure(mps);
    cout<<"#GS energy has been shifted from "<<energy<<" to "<<shifted_energy<<"!"<<endl;
    cout<<"#number of particles: "<<mps.getQN()<<endl;
    std::ofstream outfile;
    while(!checkMPS(mps));
    orthonormalize(mps,"lr");
    orthonormalize(mps,"rl");
    while(!checkMPS(mps));


    MPO<Real> Op;
    for(uint site=0;site<N;++site){
      SparseLocalOperator<Real>cuD(4,4),cdD(4,4),cu(4,4),cd(4,4),p(4,4),O_11(4,4);
      cuD[OpKeyType<1>(2,0)] = (1.0);cuD[OpKeyType<1>(3,1)] = (1.0); 
      cu[OpKeyType<1>(0,2)] = (1.0);cu[OpKeyType<1>(1,3)] = (1.0);   
								   
      cdD[OpKeyType<1>(1,0)] = (1.0);cdD[OpKeyType<1>(3,2)] = (-1.0);
      cd[OpKeyType<1>(0,1)] = (1.0);cd[OpKeyType<1>(2,3)] = (-1.0);  

      p[OpKeyType<1>(0,0)]=(1.0);p[OpKeyType<1>(1,1)]=(-1.0);p[OpKeyType<1>(2,2)]=(-1.0);p[OpKeyType<1>(3,3)]=(1.0);

      O_11[OpKeyType<1>(0,0)]=1.0;O_11[OpKeyType<1>(1,1)]=1.0;O_11[OpKeyType<1>(2,2)]=1.0;O_11[OpKeyType<1>(3,3)]=1.0;

      MPOMat<1,Real>mat;
      if(site<(N0-1)){
	mat.clear();
	p.setSite(site);
	mat[ukey<2>(0,0)]=p;
      }
      if(site==(N0-1)){
	mat.clear();
	if(part&&up&&(orbital==0)){
	  cout<<"#applying up particle creation in orbital "<<orbital<<endl;
	  cuD.setSite(site);
	  mat[ukey<2>(0,0)]=cuD;
	}else if((!part)&&up&&(orbital==0)){
	  cout<<"#applying up hole creation in orbital "<<orbital<<endl;
	  cu.setSite(site);
	  mat[ukey<2>(0,0)]=cu;
	}else if(part&&(!up)&&(orbital==0)){
	  cout<<"#applying down particle creation in orbital "<<orbital<<endl;
	  cdD.setSite(site);
	  mat[ukey<2>(0,0)]=cdD;
	}else if((!part)&&(!up)&&(orbital==0)){
	  cout<<"#applying down hole creation in orbital "<<orbital<<endl;
	  cd.setSite(site);
	  mat[ukey<2>(0,0)]=cd;
	}else if(orbital==1){
	  p.setSite(site);
	  mat[ukey<2>(0,0)]=p;
	}
      }
      if(site==N0){
	mat.clear();
	if(part&&up&&(orbital==1)){
	  cout<<"#applying up particle creation in orbital "<<orbital<<endl;
	  cuD.setSite(site);
	  mat[ukey<2>(0,0)]=cuD;
	}else if((!part)&&up&&(orbital==1)){
	  cout<<"#applying up hole creation in orbital "<<orbital<<endl;
	  cu.setSite(site);
	  mat[ukey<2>(0,0)]=cu;
	}else if(part&&(!up)&&(orbital==1)){
	  cout<<"#applying down particle creation in orbital "<<orbital<<endl;
	  cdD.setSite(site);
	  mat[ukey<2>(0,0)]=cdD;
	}else if((!part)&&(!up)&&(orbital==1)){
	  cout<<"#applying down hole creation in orbital "<<orbital<<endl;
	  cd.setSite(site);
	  mat[ukey<2>(0,0)]=cd;
	}else if(orbital==0){
	  O_11.setSite(site);
	  mat[ukey<2>(0,0)]=O_11;
	}
      }
      if(site>N0){
	mat.clear();
	O_11.setSite(site);
	mat[ukey<2>(0,0)]=O_11;
      }
      mat.Dl=1;mat.Dr=1; 
      Op[site]=mat;
    }
    
    tmps = applyMPO(Op,mps);
    while(!checkMPS(tmps));
    cout<<"before normalization: <0|c* c|0>= "<<computeOverlap(tmps,tmps)<<endl;
    orthonormalize(tmps,"lr");
    orthonormalize(tmps,"rl");
    cout<<"after normalization: <0|c* c|0> = "<<computeOverlap(tmps,tmps)<<endl;
    std::string name=simname+"_initial_state";
    MPS<AbKeyType,Complex>cmps=cast<AbKeyType,Complex,Real>(tmps);
    tmps.clear();
    mps.clear();

    SparseLocalOperator<Complex> nu(4,4),nd(4,4);

    nu[OpKeyType<1>(2,2)]=(1.0);nu[OpKeyType<1>(3,3)]=(1.0);
    nd[OpKeyType<1>(1,1)]=(1.0);nd[OpKeyType<1>(3,3)]=(1.0);

    std::vector<Operator<Complex>*> NU,ND,ENTROPY,OV,NUND,NU0NU1,NU0ND1,ND0ND1,ND0NU1;
    MPO<Complex> nu0nu1,nd0nd1,nu0nd1,nd0nu1;
    MPOMat<1,Complex> mat;

    nu.setSite(N/2-1);
    mat[ukey<2>(0,0)]=nu;
    mat.Dl=1;mat.Dr=1;
    nu0nu1[N/2-1]=mat;
    nu0nd1[N/2-1]=mat;
    mat.clear();

    nd.setSite(N/2-1);
    mat[ukey<2>(0,0)]=nd;
    mat.Dl=1;mat.Dr=1;
    nd0nu1[N/2-1]=mat;
    nd0nd1[N/2-1]=mat;
    mat.clear();

    nu.setSite(N/2);
    mat[ukey<2>(0,0)]=nu;
    mat.Dl=1;mat.Dr=1;
    nu0nu1[N/2]=mat;
    nd0nu1[N/2]=mat;
    mat.clear();

    nd.setSite(N/2);
    mat[ukey<2>(0,0)]=nd;
    mat.Dl=1;mat.Dr=1;
    nu0nd1[N/2]=mat;
    nd0nd1[N/2]=mat;
    mat.clear();
    
    NU0NU1.push_back(new MPO<Complex>(nu0nu1));
    NU0ND1.push_back(new MPO<Complex>(nu0nd1));
    ND0NU1.push_back(new MPO<Complex>(nd0nu1));
    ND0ND1.push_back(new MPO<Complex>(nd0nd1));

    for(int i = 0;i<N;++i){
      nu.setSite(i);
      nd.setSite(i);
      NU.push_back(new SparseLocalOperator<Complex>(nu));
      NUND.push_back(new SparseLocalOperator<Complex>((nu*nd)));
      ND.push_back(new SparseLocalOperator<Complex>(nd));
      if(i<(N-1)){
	//==============bipartite entanglement============
	Entropy<Complex> ent(i);
	ENTROPY.push_back(new Entropy<Complex>(ent));
      }
    }
    Overlap<Complex> ov(cmps);
    OV.push_back(new Overlap<Complex>(ov));

    measurementObject<Complex> mOverlap(OV,"overlap"+simname,step,true,true);
    measurementObject<Complex> mNu(NU,"Nu"+simname,step);
    measurementObject<Complex> mNd(ND,"Nd"+simname,step);
    measurementObject<Complex> mNuNd(NUND,"NuNd"+simname,step);
    measurementObject<Complex> mNu0Nd1(NU0ND1,"Nu0Nd1"+simname,step);
    measurementObject<Complex> mNu0Nu1(NU0NU1,"Nu0Nu1"+simname,step);
    measurementObject<Complex> mNd0Nd1(ND0ND1,"Nd0Nd1"+simname,step);
    measurementObject<Complex> mNd0Nu1(ND0NU1,"Nd0Nu1"+simname,step);
    measurementObject<Complex> mEntropy(ENTROPY,"entropy"+simname,step);

    list<measurementObject<Complex> > mList;
    //mList.push_back(mOverlap);  
    mList.push_back(mNu);  
    mList.push_back(mNd);  
    mList.push_back(mNuNd);  
    mList.push_back(mNu0Nd1);  
    mList.push_back(mNu0Nu1);  
    mList.push_back(mNd0Nd1);  
    mList.push_back(mNd0Nu1);  
    mList.push_back(mEntropy);
    cout<<cmpo._hop<<endl;
    TEBDengine<AbKeyType,Complex> SIAM(cmps,cmpo,simname);
    SIAM.canonizeMPS();
    SIAM.measureCMPS(mList);
    SIAM.writeMeasurements(mList);
    SIAM.setCheckpointStep(cpstep);
    //bool normalize=true;
    bool recanonize=false;
    
    SIAM.simulate(mList,order,NmaxTEBD,dt,chimax,true,tr_weight,normalize,recanonize);
    cout<<"#finished simulation ..."<<endl;
  }
  cout<<"#exiting qis.exe ..."<<endl;
  return 0;
}
