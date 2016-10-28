// #include <boost/python/module.hpp>
// #include <boost/python/def.hpp>
#include "MPSTypes.hpp"

#include "MPSfunctions.hpp"
#include "Hamiltonians.hpp"
#include "BlockTensor.hpp"
#include "tensormap.hpp"
#include "container.hpp"
#include "utilities.hpp"
#include <iostream>




int main(const int argc,const char **argv){
  cout<<"##########################################################################################################################"<<endl;
  cout<<"#######                     TWO-ORBITAL DMFT SOLVER FOR THE HUBBARD MODEL ON THE BETHE LATTICE                         ###"<<endl;
  cout<<"#######            IMPORTANT NOTE: FOR SPIN-Z=0 AND HALFFILLING, EACH CHAIN HAS TO BE OF EVEN LENGTH                   ###"<<endl;
  cout<<"##########################################################################################################################"<<endl;
  cout<<"#######  Usage: --label=name: label of the input files; bath_params_orbital_0_'label' and bath_params_orbital_1_'label'###"<<endl;
  cout<<"#######                       contain the bathparameters                                                               ###"<<endl;
  cout<<"#######                       H_local_'label' contains local Hamiltonian; structure of the file:                       ###"<<endl;
  cout<<"#######                       e0d                                                                                      ###"<<endl;
  cout<<"#######                       e0u                                                                                      ###"<<endl;
  cout<<"#######                       e1d                                                                                      ###"<<endl;
  cout<<"#######                       e1u                                                                                      ###"<<endl;
  cout<<"#######                       U                                                                                        ###"<<endl;
  cout<<"#######                       J                                                                                        ###"<<endl;
  cout<<"##########################################################################################################################"<<endl;
  cout<<"#######                                                                                                                ###"<<endl;
  cout<<"#######  Parameters: NmaxTEBD: number of time steps to be done                                                         ###"<<endl;
  cout<<"#######              chiDMRG: matrix dimension of DMRG                                                                 ###"<<endl;
  cout<<"#######              NmaxDMRG: number of two-side DMRG step                                                            ###"<<endl;
  cout<<"#######              NmaxSSDMRG: number of single site DMRG step                                                       ###"<<endl;
  cout<<"#######              label: name of the simulation (see above)                                                         ###"<<endl;
  cout<<"#######              measure: measurement interval in TEBD-timesteps                                                   ###"<<endl;
  cout<<"#######              chiTEBD: matrix dimension during TEBD                                                             ###"<<endl;
  cout<<"#######              dt:      time step for trotterization; it can be positive or negative; a negative dt run          ###"<<endl;
  cout<<"#######                       creates 'Numfiles' MPS files at increasint times and stores them on the hard disc        ###"<<endl;
  cout<<"#######                       (can be very large!). If dt>0, a collection run is started which does a TEBD time        ###"<<endl;
  cout<<"#######                       evolution and which collects the formerly produces MPS files to obtain the Greens        ###"<<endl;
  cout<<"#######                       function at the corresponding time difference. If it can not find an MPS file, it        ###"<<endl;
  cout<<"#######                       waits until it will be created (-> make sure the production run is running).             ###"<<endl;
  cout<<"#######                                                                                                                ###"<<endl;
  cout<<"#######              trweight: the truncated weight that should be obtained (soft bound, will be                       ###"<<endl;
  cout<<"#######                        broken when chiTEBD is reached                                                          ###"<<endl;
  cout<<"#######              spin,brach,orbital: determines which Greens function to be calculated                             ###"<<endl;
  cout<<"#######              verbose: output flag                                                                              ###"<<endl;
  cout<<"#######              Numfiles: the number of of MPS-files that the collection run can expect to be there               ###"<<endl;
  cout<<"#######                                                                                                                ###"<<endl;
  cout<<"#######                                                                                                                ###"<<endl;
  cout<<"#######  Output: the collection run with dt>0 (see above) if run correctly, produces files with names                  ###"<<endl;
  cout<<"#######          overlap_'label'dt'dt' of the Greensfunction at the corresponding times. These can                     ###"<<endl;
  cout<<"#######          then be further used in subsequent data analysis.                                                     ###"<<endl;
  cout<<"#######                                                                                                                ###"<<endl;
  cout<<"##########################################################################################################################"<<endl;

  
  parser.addoption("NmaxTEBD", Parameter_Integer,"500");
  parser.addoption("chiDMRG", Parameter_Integer,"100");
  parser.addoption("NmaxDMRG", Parameter_Integer,"5");
  parser.addoption("NmaxSSDMRG", Parameter_Integer,"10");
  parser.addoption("label", Parameter_String,"default_name");
  parser.addoption("measure", Parameter_Integer,"1");
  parser.addoption("chiTEBD", Parameter_Integer,"200");
  parser.addoption("dt", Parameter_Float,"0.05");
  parser.addoption("trweight", Parameter_Float,"1e-8");
  parser.addoption("spin", Parameter_Integer,"1");//for the down-spin spectral function
  parser.addoption("branch", Parameter_Integer,"1");//for the hole part of the spectral function
  parser.addoption("orbital", Parameter_Integer,"0");//choose an orbital
  parser.addoption("verbose", Parameter_Option);
  parser.addoption("deltaDMRG", Parameter_Float,"1e-12");
  parser.addoption("Numfiles", Parameter_Integer,"10");

  if (!parser.load(argc, argv))
    return 0;
  std::cout << "# Options:\n" << parser << std::endl;


  const uint NmaxTEBD = parser.getint("NmaxTEBD");
  const uint chiDMRG = parser.getint("chiDMRG");
  const uint maxstepsDMRG = parser.getint("NmaxDMRG");
  const uint maxstepsSSDMRG = parser.getint("NmaxSSDMRG");
  string simname = parser.getstring("label");
  const uint step = parser.getint("measure");
  const uint chimax = parser.getint("chiTEBD");
  const double dt = parser.getfloat("dt");
  const double tr_weight = parser.getfloat("trweight");
  const int spin=parser.getint("spin");
  const int branch=parser.getint("branch");
  const int orbital=parser.getint("orbital");
  const bool verbose=parser.isgiven("verbose");
  const double deltaDMRG = parser.getfloat("deltaDMRG");
  const int Numfiles=parser.getint("Numfiles");

  
  
  const string order="second";
  const uint cpstep = step;
  const double landelta = 1e-10;
  const bool testenergy=true;
  const double conv = 1e-10;
  const bool normalize=true;


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

  //check if there's a simname_CP.container checkpoint file of a time evolution which is to be resumed
  std::ostringstream convert;
  convert<<dt;
  string name=simname+"dt"+convert.str()+"_CP.container";
  cout<<"#starting two-orbital solver using "<<omp_get_max_threads()<<" parallel threads"<<endl;
  if(!fexist(name)){
    cout<<"starting fresh run"<<endl;
    //bath parameters are read in by the mpo class
    twoOrbitalSIAMMPO<Real> mpo(bath_0,bath_1,H0);
    twoOrbitalSIAMMPO<Complex> cmpo(bath_0,bath_1,H0);
    uint N=mpo.size();
    uint N0=mpo.N0();
    uint N1=mpo.N1();

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
    //now check if a dmrg groundstate calculation already exists. if so, then load it, if not then do a
    //dmrg and save all the relevant data


    string sssname=simname+"_singleSite";
    string dmrgGS=simname+"_singleSite.container";
    if(!fexist(dmrgGS)){
      cout<<"#Starting dmrg run"<<endl;
      if(maxstepsDMRG>0){
	bool doTrunc=true;
	dmrg.simulate(maxstepsDMRG,conv,chiDMRG,doTrunc,deltaDMRG,landelta,testenergy);
      }
      SSDMRGengine<AbKeyType,Real> SSdmrg(dmrg.getMPS(),mpo,sssname,verbose);
      if(maxstepsDMRG>0){
	SSdmrg.simulate(maxstepsSSDMRG,1e-10,1e-12,true);
      }
      SSdmrg.write();
      mps=SSdmrg.getMPS();
      Real energy=mpo.measure(mps);
      cout<<"#GS-energy after two-site and single-site DMRG: "<<energy<<endl;
      mpo.shiftandrescale(energy,0.0,0.0);
      cmpo.shiftandrescale(Complex(energy,0.0),Complex(0.0,0.0),Complex(0.0,0.0));
      Real shifted_energy=mpo.measure(mps);
      cout<<"#GS energy has been shifted from "<<energy<<" to "<<shifted_energy<<"!"<<endl;
      cout<<"#number of particles: "<<mps.getQN()<<endl;
      std::ofstream outfile;
      while(!checkMPS(mps));
      orthonormalize(mps,"lr");
      orthonormalize(mps,"rl");
      while(!checkMPS(mps));
    }else if(fexist(dmrgGS)){
      cout<<"#found "<<dmrgGS<<" file. Reading in data:"<<endl;
      SSDMRGengine<AbKeyType,Real> SSdmrg(sssname);
      SSdmrg.read();
      mps=SSdmrg.getMPS();
      Real energy=mpo.measure(mps);
      mpo.shiftandrescale(energy,0.0,0.0);
      cmpo.shiftandrescale(Complex(energy,0.0),Complex(0.0,0.0),Complex(0.0,0.0));
      Real shifted_energy=mpo.measure(mps);
      cout<<"#GS energy has been shifted from "<<energy<<" to "<<shifted_energy<<"!"<<endl;
      cout<<"#number of particles: "<<mps.getQN()<<endl;
      std::ofstream outfile;
      while(!checkMPS(mps));
      orthonormalize(mps,"lr");
      orthonormalize(mps,"rl");
      while(!checkMPS(mps));
    }

    //creation and annihilation operators
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
    //now create the initial state; this state is common to the two tebd runs 
    tmps = applyMPO(Op,mps);

    while(!checkMPS(tmps));
    cout<<"before normalization: <0|c* c|0>= "<<computeOverlap(tmps,tmps)<<endl;
    orthonormalize(tmps,"lr");
    orthonormalize(tmps,"rl");
    cout<<"after normalization: <0|c* c|0> = "<<computeOverlap(tmps,tmps)<<endl;
    MPS<AbKeyType,Complex>cmps=cast<AbKeyType,Complex,Real>(tmps);
    MPS<AbKeyType,Complex>mps_GS=cast<AbKeyType,Complex,Real>(mps);
    tmps.clear();
    mps.clear();
    orthonormalize(mps_GS,"lr");
    orthonormalize(mps_GS,"rl");
    CanonizedMPS<AbKeyType,Complex>canmps_GS=canonize(mps_GS);//for groundstate measurements
    mps_GS.clear();    

    SparseLocalOperator<Complex> nu(4,4),nd(4,4);

    nu[OpKeyType<1>(2,2)]=(1.0);nu[OpKeyType<1>(3,3)]=(1.0);
    nd[OpKeyType<1>(1,1)]=(1.0);nd[OpKeyType<1>(3,3)]=(1.0);

    std::vector<Operator<Complex>*> NU,ND,ENTROPY,OV;
    for(int i = 0;i<N;++i){
      nu.setSite(i);
      nd.setSite(i);
      NU.push_back(new SparseLocalOperator<Complex>(nu));
      ND.push_back(new SparseLocalOperator<Complex>(nd));
      if(i<(N-1)){
	//==============bipartite entanglement============
	Entropy<Complex> ent(i);
	ENTROPY.push_back(new Entropy<Complex>(ent));
      }
    }

    //========================================================================================================================================================================================================
    //now decide wether to do a storage-run or a measurement run; if dt>0, a measurement run is done which collects .mps files from a 
    //previously started dt<0 run.
    if(dt<0.0){
      //this run produces .mps files with a running index; they can be collected
      //note that only a single .container file is written
      const bool store=true;
      convert.str("");
      convert<<dt;

      //this run only stores the .container files at different breakpoints, given by checkpointstep
      measurementObject<Complex> mNd(NU,"Nup"+simname+"dt"+convert.str(),step,true,false);
      measurementObject<Complex> mNu(ND,"Ndown"+simname+"dt"+convert.str(),step,true,false);
      measurementObject<Complex> mEntropy(ENTROPY,"entropy_"+simname+"dt"+convert.str(),step,true,false);

      list<measurementObject<Complex> > mList;
      mList.push_back(mNu);  
      mList.push_back(mNd);  
      mList.push_back(mEntropy);

      //GS measurements=============================
      typename list<measurementObject<Complex> >::iterator it;
      for(it = mList.begin();it!=mList.end();++it)
	it->measure(canmps_GS);

      canmps_GS.clear();
      if(NmaxTEBD!=0){
	TEBDengine<AbKeyType,Complex> SIAM(cmps,cmpo,simname+"dt"+convert.str(),store);
	SIAM.canonizeMPS();
	SIAM.measureCMPS(mList);
	SIAM.writeMeasurements(mList);
	SIAM.setCheckpointStep(cpstep);

	bool recanonize=false;
	//the first measurement at t=0;
	SIAM.save(simname+"dt"+convert.str()+"_CP0");
	SIAM.simulate(mList,order,NmaxTEBD,dt,chimax,true,tr_weight,normalize,recanonize);
      }
      cout<<"#finished simulation ..."<<endl;
    }else if(dt>0.0){
      //this run collects .mps files from a previous run and computes the overlap with the corresponding mps from a 
      //run with dt>0; 
      convert.str("");
      convert<<(-dt);
      std::string filelabel=simname+"dt"+convert.str();
      convert.str("");
      convert<<(dt);
      OverlapWithStored<Complex> ov(filelabel,Numfiles);
      OV.push_back(new OverlapWithStored<Complex>(ov));

      measurementObject<Complex> mOverlap(OV,"overlap_"+simname+"dt"+convert.str(),cpstep,true,true);
      measurementObject<Complex> mNd(NU,"Nup"+simname+"dt"+convert.str(),step,true,false);
      measurementObject<Complex> mNu(ND,"Ndown"+simname+"dt"+convert.str(),step,true,false);
      measurementObject<Complex> mEntropy(ENTROPY,"entropy_"+simname+"dt"+convert.str(),step,true,false);
      
      list<measurementObject<Complex> > mList;
      mList.push_back(mNu);  
      mList.push_back(mNd);  
      mList.push_back(mEntropy);


      //GS measurements=============================
      typename list<measurementObject<Complex> >::iterator it;
      for(it = mList.begin();it!=mList.end();++it)
	it->measure(canmps_GS);

      canmps_GS.clear();
      if(NmaxTEBD!=0){
	mList.push_back(mOverlap);  
	const bool store=false;
	TEBDengine<AbKeyType,Complex> SIAM(cmps,cmpo,simname+"dt"+convert.str());    
	SIAM.canonizeMPS();
	SIAM.measureCMPS(mList);
	SIAM.writeMeasurements(mList);
	SIAM.setCheckpointStep(cpstep);
	bool recanonize=false;
	SIAM.simulate(mList,order,NmaxTEBD,dt,chimax,true,tr_weight,normalize,recanonize);
      }
      cout<<"#finished simulation ..."<<endl;
    }

  }else if(fexist(name)){
    //resumes the .mps file production/collection run, depending on wether dt>0 or dt<0
    const bool store=true;//can be omitted; the loaded container already contains the flag
    string loadname=simname+"dt"+convert.str()+"_CP";
    TEBDengine<AbKeyType,Complex> SIAM(loadname);
    SIAM.read();
    int N=SIAM.getMPO().size();
    if(N%2!=0){
      cout<<"#warning: using an odd chain length!"<<endl;
    }

    SparseLocalOperator<Complex> nu(4,4),nd(4,4);   //elementary,identity and fermi sign operators
    //elementary,identity and fermi sign operators
    nu[OpKeyType<1>(2,2)]=(1.0);nu[OpKeyType<1>(3,3)]=(1.0);
    nd[OpKeyType<1>(1,1)]=(1.0);nd[OpKeyType<1>(3,3)]=(1.0);
    std::vector<Operator<Complex>*> NU,ND,ENTROPY,OV;
    for(int i = 0;i<N;++i){
      nu.setSite(i);
      nd.setSite(i);
      NU.push_back(new SparseLocalOperator<Complex>(nu));
      ND.push_back(new SparseLocalOperator<Complex>(nd));
      if(i<(N-1)){
	//==============bipartite entanglement============
	Entropy<Complex> ent(i);
	ENTROPY.push_back(new Entropy<Complex>(ent));
      }
    }

    convert.str("");
    if(dt>0.0)
      convert<<(-dt);
    std::string filelabel=simname+"dt"+convert.str();

    OverlapWithStored<Complex> ov(filelabel,Numfiles);
    OV.push_back(new OverlapWithStored<Complex>(ov));


    convert.str("");
    convert<<(dt);
    measurementObject<Complex> mOverlap(OV,"overlap_"+simname+"dt"+convert.str(),cpstep,true,true);
    measurementObject<Complex> mNu(NU,"Nu"+simname,step);
    measurementObject<Complex> mNd(ND,"Nd"+simname,step);
    measurementObject<Complex> mEntropy(ENTROPY,"entropy_"+simname+"dt"+convert.str(),step);


    list<measurementObject<Complex> > mList;
    mList.push_back(mNu);  
    mList.push_back(mNd);  
    mList.push_back(mEntropy);
    if(dt>0.0)
      mList.push_back(mOverlap);  
    bool recanonize=false;
    cout<<"#resuming simulation in "<<name<<endl;
    SIAM.simulate(mList,order,NmaxTEBD,dt,chimax,true,tr_weight,normalize,recanonize);
    cout<<"#finished simulation ..."<<endl;
  }
  cout<<"#exiting qis.exe ..."<<endl;
  return 0;
}

int main(const int argc,const char **argv){
  TOSolver(200,400,10, 10,"test_label",1,50,-0.05,1e-10,1,1,1,true,1e-10,10);
  TOSolver(200,400,10, 10,"test_label",1,50,0.05,1e-10,1,1,1,true,1e-10,200);
}

// BOOST_PYTHON_MODULE(TO_solver_ext)
// {
//     using namespace boost::python;
//     def("TOSolver", TOSolver);
// }

