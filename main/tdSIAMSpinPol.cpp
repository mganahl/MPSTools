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
  cout<<"####################################################################################################################################"<<endl;
  cout<<"#######                                                                                                                          ###"<<endl;
  cout<<"#######                               SINGLE-ORBITAL DMFT SOLVER FOR THE HUBBARD MODEL ON THE BETHE LATTICE                      ###"<<endl;
  cout<<"#######                                                                                                                          ###"<<endl;
  cout<<"        ########################################################################################################################"<<endl;
  cout<<"#######                                                                                                                          ###"<<endl;
  cout<<"#######    Usage: --label='name': name of the bath_params_'name' input file and the local Hamiltonian H_local_'name'             ###"<<endl;
  cout<<"#######                           bath_params_'name' is a column vector; first entry is the number of bath orbitals              ###"<<endl;
  cout<<"#######                           second is the bandwidth, third is the hybridization hopping (sqrt{xi_0/pi} in                  ###"<<endl;
  cout<<"#######                           Rev. Mod. Phys 80, 395); rest alternating bath-hopping and bath-onsite potential               ###"<<endl;
  cout<<"#######                           (the last two entries are then two bath-onsite potentials).                                    ###"<<endl;
  cout<<"#######                           H_local_'name' is a column vector with three entries: first is the down-potential              ###"<<endl;
  cout<<"#######                           at the impurity, second the up-potential at the impurity, third is U                           ###"<<endl;
  cout<<"#######                                                                                                                          ###"<<endl;
  cout<<"####################################################################################################################################"<<endl;
  cout<<"####################################################################################################################################"<<endl;

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
  parser.addoption("verbose", Parameter_Option);
  parser.addoption("SingleParticleGreensfunction", Parameter_Option);
  parser.addoption("SpinStructureFactor", Parameter_Option);
  parser.addoption("ChargeStructureFactor", Parameter_Option);
  parser.addoption("DoubleOccupancyStructureFactor", Parameter_Option);
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
  const bool spg=parser.isgiven("SingleParticleGreensfunction");
  const bool ssf=parser.isgiven("SpinStructureFactor");
  const bool csf=parser.isgiven("ChargeStructureFactor");
  const bool dosf=parser.isgiven("DoubleOccupancyStructureFactor");
  const bool verbose=parser.isgiven("verbose");
  const bool testenergy=parser.isgiven("testRes")?false:true;
  const double deltaDMRG = parser.getfloat("deltaDMRG");
  const double conv = parser.getfloat("conv");
  const bool normalize=parser.isgiven("nonormalize")?false:true;


  //consistency check:
    
  string bath_up = "bath_params_up"+simname;
  string bath_down = "bath_params_down"+simname;
  string init_up = "filling_up"+simname;
  string init_down = "filling_down"+simname;
  string H0 = "H_local_"+simname;


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
  cout<<"#starting 1-orbital tdSIAM.exe solver using "<<omp_get_max_threads()<<" parallel threads"<<endl;
  if(fexist(name)){
    string loadname=simname+"_CP";
    TEBDengine<AbKeyType,Complex> SIAM(loadname);
    SIAM.read(loadname);
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


    std::vector<Operator<Complex>*> DEN,ENTROPY,HAM,NORM,OV,REDRHO;
    MPO<Complex> e_e,d_d,ud_ud,u_u;;
    for(int i = 0;i<N_up+N_down;++i){
      SparseLocalOperator<Complex> n(2,2),e(2,2);
      e[OpKeyType<1>(0,0)]=(1.0);
      n[OpKeyType<1>(1,1)]=(1.0);
      n.setSite(i);
      e.setSite(i);
      DEN.push_back(new SparseLocalOperator<Complex>(n));
      if(i<(N_up+N_down-1)){
	//==============bipartite entanglement============
	Entropy<Complex> ent(i);
	ENTROPY.push_back(new Entropy<Complex>(ent));
      }
      if(i==N_down-1){

	MPOMat<1,Complex> m;
	m[ukey<2>(0,0)]=n;
	m.Dl=1;m.Dr=1;
	ud_ud[i]=m;

	m.clear();
	m[ukey<2>(0,0)]=e;
	m.Dl=1;m.Dr=1;
	u_u[i]=m;

	m[ukey<2>(0,0)]=n;
	m.Dl=1;m.Dr=1;
	d_d[i]=m;

	m.clear();
	m[ukey<2>(0,0)]=e;
	m.Dl=1;m.Dr=1;
	e_e[i]=m;
      }
      if(i==N_down){
	MPOMat<1,Complex> m;
       	n.setSite(i);
	e.setSite(i);
	m[ukey<2>(0,0)]=n;
	m.Dl=1;m.Dr=1;
	ud_ud[i]=m;

	m[ukey<2>(0,0)]=n;
	m.Dl=1;m.Dr=1;
	u_u[i]=m;

	m.clear();
	m[ukey<2>(0,0)]=e;
	m.Dl=1;m.Dr=1;
	d_d[i]=m;

	m.clear();
	m[ukey<2>(0,0)]=e;
	m.Dl=1;m.Dr=1;
	e_e[i]=m;
      }
    }
    // HAM.push_back(new MPO<Complex>(SIAM.getMPO()));
    // NORM.push_back(new Norm<Complex>());
    REDRHO.push_back(new MPO<Complex>(ud_ud));
    REDRHO.push_back(new MPO<Complex>(u_u));
    REDRHO.push_back(new MPO<Complex>(d_d));
    REDRHO.push_back(new MPO<Complex>(e_e));

    MPS<AbKeyType,Real>mps;
    std::string name=simname+"_initial_state";
    mps.load(name);
    MPS<AbKeyType,Complex>cmps=cast<AbKeyType,Complex,Real>(mps);
    mps.clear();

    //Overlap<Complex> ov(cmps);
    TRIOverlap<Complex> ov;
    OV.push_back(new TRIOverlap<Complex>(ov));
    measurementObject<Complex> mOverlap(OV,"overlap_"+simname,step,true,true);
    measurementObject<Complex> mN(DEN,"N"+simname,step);
    measurementObject<Complex> mEntropy(ENTROPY,"entropy_"+simname,step);
    measurementObject<Complex> mRedrho(REDRHO,"ReducedDensityMatrix_"+simname,step,true,true);
    // measurementObject<Complex> mHam(HAM,"Ham_"+simname,step);
    // measurementObject<Complex> mNorm(NORM,"Norm_"+simname,step);
    list<measurementObject<Complex> > mList;
    mList.push_back(mOverlap);  
    mList.push_back(mN);  
    mList.push_back(mEntropy);
    mList.push_back(mRedrho);
    //mList.push_back(mNorm);
    //mList.push_back(mHam);
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

    SpinlessFermionsMPO<Real> mpo(U,hop,mu);
    SpinlessFermionsMPO<Complex> cmpo(U,hop,mu);
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
    MPO<Real> SpinOp;
    for(uint site=0;site<N_up+N_down;++site){
      MPOMat<1,Real>mat;
      SparseLocalOperator<Real> O_11(2,2),cD(2,2),c(2,2),p(2,2),n(2,2);
      O_11[OpKeyType<1>(0,0)]=1.0;O_11[OpKeyType<1>(1,1)]=1.0;
      p[OpKeyType<1>(0,0)]=1.0;p[OpKeyType<1>(1,1)]=-1.0;
      cD[OpKeyType<1>(1,0)] = 1.0;
      c[OpKeyType<1>(0,1)] = 1.0;
      n[OpKeyType<1>(1,1)]=1.0;
      if(spg){
	if(site<(N_down-1)){
	  mat.clear();
	  mat[ukey<2>(0,0)]=p;
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
      }else if(ssf){
	if(site<(N_down-1)){
	  mat.clear();
	  mat[ukey<2>(0,0)]=O_11;
	  mat.Dl=1;mat.Dr=1;
	}
	if(site==N_down-1){
	  cout<<"#calculating spin-structure factor"<<endl;
	  mat.clear();
	  mat[ukey<2>(0,0)]=n;
	  mat[ukey<2>(0,1)]=O_11;
	  mat.Dl=1;mat.Dr=2;
	}
	if(site==N_down){
	  mat.clear();
	  mat[ukey<2>(0,0)]=O_11;
	  mat[ukey<2>(1,0)]=(-1.0)*n;
	  mat.Dl=2;mat.Dr=1;
	}
	if(site>N_down){
	  mat.clear();
	  mat[ukey<2>(0,0)]=O_11;
	  mat.Dl=1;mat.Dr=1; 
	}
	Op[site]=mat;
      }else if(csf){

	if(spin==-1){
	  if(site<(N_down-1)){
	    mat.clear();
	    mat[ukey<2>(0,0)]=O_11;
	    mat.Dl=1;mat.Dr=1;
	  }
	  if(site==(N_down-1)){
	    cout<<"#calculating charge-structure factor for "<<spin<<" spins"<<endl;
	    mat.clear();
	    mat[ukey<2>(0,0)]=n;
	    mat[ukey<2>(0,1)]=O_11;
	    mat.Dl=1;mat.Dr=2;
	  }
	  if(site>(N_down-1)){
	    mat.clear();
	    mat[ukey<2>(0,0)]=O_11;
	    mat.Dl=1;mat.Dr=1; 
	  }
	}else if(spin==1){
	  if(site<N_down){
	    mat.clear();
	    mat[ukey<2>(0,0)]=O_11;
	    mat.Dl=1;mat.Dr=1;
	  }
	  if(site==N_down){
	    mat.clear();
	    mat[ukey<2>(0,0)]=n;
	    mat[ukey<2>(0,1)]=O_11;
	    mat.Dl=1;mat.Dr=2;
	  }
	  if(site>N_down){
	    mat.clear();
	    mat[ukey<2>(0,0)]=O_11;
	    mat.Dl=1;mat.Dr=1; 
	  }
	}
	Op[site]=mat;
      }else if(dosf){

	if(site<(N_down-1)){
	  mat.clear();
	  mat[ukey<2>(0,0)]=O_11;
	  mat.Dl=1;mat.Dr=1;
	}
	if(site==(N_down-1)){
	  cout<<"#calculating doubleoccupancy-structure factor"<<endl;
	  mat.clear();
	  mat[ukey<2>(0,0)]=n;
	  mat.Dl=1;mat.Dr=1;
	}
	if(site==N_down){
	  mat.clear();
	  mat[ukey<2>(0,0)]=n;
	  mat.Dl=1;mat.Dr=1;
	}
	if(site>N_down){
	  mat.clear();
	  mat[ukey<2>(0,0)]=O_11;
	  mat.Dl=1;mat.Dr=1; 
	}
	Op[site]=mat;
      }
    }
    tmps = applyMPO(Op,mps);
    while(!checkMPS(tmps));
    cout<<"before normalization: <0|c* c|0>= "<<computeOverlap(tmps,tmps)<<endl;
    orthonormalize(tmps,"lr");
    orthonormalize(tmps,"rl");
    cout<<"after normalization: <0|c* c|0> = "<<computeOverlap(tmps,tmps)<<endl;
    std::string name=simname+"_initial_state";
    tmps.save(name);
    MPS<AbKeyType,Complex>cmps=cast<AbKeyType,Complex,Real>(tmps);
    MPS<AbKeyType,Complex>mps_GS=cast<AbKeyType,Complex,Real>(mps);
    orthonormalize(mps_GS,"lr");
    orthonormalize(mps_GS,"rl");
    CanonizedMPS<AbKeyType,Complex>canmps_GS=canonize(mps_GS);
    tmps.clear();
    mps.clear();

    std::vector<Operator<Complex>*> DEN,ENTROPY,HAM,NORM,OV,REDRHO;
    MPO<Complex> e_e,d_d,u_u,ud_ud,ud_e,ud_u,ud_d,u_d,u_e,d_e;
    for(int i = 0;i<N_up+N_down;++i){
      SparseLocalOperator<Complex> n(2,2),e(2,2);
      e[OpKeyType<1>(0,0)]=(1.0);
      n[OpKeyType<1>(1,1)]=(1.0);
      n.setSite(i);
      e.setSite(i);
      DEN.push_back(new SparseLocalOperator<Complex>(n));
      if(i<(N_up+N_down-1)){
	//==============bipartite entanglement============
	Entropy<Complex> ent(i);
	ENTROPY.push_back(new Entropy<Complex>(ent));
      }
      if(i==N_down-1){
	MPOMat<1,Complex> m;
	m[ukey<2>(0,0)]=n;
	m.Dl=1;m.Dr=1;
	ud_ud[i]=m;

	m.clear();
	m[ukey<2>(0,0)]=e;
	m.Dl=1;m.Dr=1;
	u_u[i]=m;

	m[ukey<2>(0,0)]=n;
	m.Dl=1;m.Dr=1;
	d_d[i]=m;

	m.clear();
	m[ukey<2>(0,0)]=e;
	m.Dl=1;m.Dr=1;
	e_e[i]=m;
      }
      if(i==N_down){
	MPOMat<1,Complex> m;
       	n.setSite(i);
	e.setSite(i);
	m[ukey<2>(0,0)]=n;
	m.Dl=1;m.Dr=1;
	ud_ud[i]=m;

	m[ukey<2>(0,0)]=n;
	m.Dl=1;m.Dr=1;
	u_u[i]=m;

	m.clear();
	m[ukey<2>(0,0)]=e;
	m.Dl=1;m.Dr=1;
	d_d[i]=m;

	m.clear();
	m[ukey<2>(0,0)]=e;
	m.Dl=1;m.Dr=1;
	e_e[i]=m;
      }

    }
    // HAM.push_back(new MPO<Complex>(cmpo));
    // NORM.push_back(new Norm<Complex>());
    REDRHO.push_back(new MPO<Complex>(ud_ud));
    REDRHO.push_back(new MPO<Complex>(u_u));
    REDRHO.push_back(new MPO<Complex>(d_d));
    REDRHO.push_back(new MPO<Complex>(e_e));
    Overlap<Complex> ov(cmps);
    OV.push_back(new Overlap<Complex>(ov));

    measurementObject<Complex> mOverlap(OV,"overlap_"+simname,step,true,true);
    measurementObject<Complex> mN(DEN,"N"+simname,step);
    measurementObject<Complex> mEntropy(ENTROPY,"entropy_"+simname,step);
    //measurementObject<Complex> mHam(HAM,"Ham_"+simname,step);
    //measurementObject<Complex> mNorm(NORM,"Norm_"+simname,step);
    measurementObject<Complex> mRedrho(REDRHO,"ReducedDensityMatrix_"+simname,step,true,true);

    list<measurementObject<Complex> > mList;
    mList.push_back(mN);  
    mList.push_back(mEntropy);
    //mList.push_back(mNorm);
    //mList.push_back(mHam);
    mList.push_back(mRedrho);

    //GS measurements=============================
    
    typename list<measurementObject<Complex> >::iterator it;
    for(it = mList.begin();it!=mList.end();++it){
      it->measure(canmps_GS);
    }
    

    mps_GS.clear();    
    canmps_GS.clear();

    mList.push_back(mOverlap);  

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
