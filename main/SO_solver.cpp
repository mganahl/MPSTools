#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include "MPSTypes.hpp"
#include "MPSfunctions.hpp"
#include "Hamiltonians.hpp"
#include "BlockTensor.hpp"
#include "tensormap.hpp"
#include "container.hpp"
#include "utilities.hpp"
#include <iostream>




int SOSolver(uint NmaxTEBD,uint chiDMRG,uint maxstepsDMRG, uint maxstepsSSDMRG, string simname, string order, uint step, uint chimax, Real dt, Real tr_weight, int spin, int branch, bool verbose, 
	     Real deltaDMRG, int Numfiles){
  cout<<"####################################################################################################################################"<<endl;
  cout<<"#######                                                                                                                          ###"<<endl;
  cout<<"#######                         IMPROVED SINGLE-ORBITAL DMFT SOLVER FOR THE HUBBARD MODEL ON THE BETHE LATTICE                   ###"<<endl;
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
  cout<<"#######                                                                                                                          ###"<<endl;
  cout<<"#######                                                                                                                          ###"<<endl;
  cout<<"#######                                                                                                                          ###"<<endl;
  cout<<"#######                                                                                                                          ###"<<endl;
  cout<<"####################################################################################################################################"<<endl;
  cout<<"####################################################################################################################################"<<endl;

  std::ofstream file;
  std::string paramfile="Parameters_"+simname;
  file.open(paramfile.c_str());
  file<<"NmaxTEBD:          "<<NmaxTEBD<<endl;
  file<<"chiDMRG:           "<<chiDMRG<<endl;
  file<<"maxstepsDMRG:      "<<maxstepsDMRG<<endl;
  file<<"maxstepsSSDMRG:    "<<maxstepsSSDMRG<<endl;
  file<<"simname:           "<<simname<<endl;
  file<<"order:             "<<order<<endl;
  file<<"measurement:       "<<step<<endl;
  file<<"chiTEBD:           "<<chimax<<endl;
  file<<"dt:                "<<dt<<endl;
  file<<"truncated weight:  "<<tr_weight<<endl;
  file<<"spin:              "<<spin<<endl;
  file<<"branch:            "<<branch<<endl;
  file<<"verbose:           "<<verbose<<endl;
  file<<"deltaDMRG:         "<<deltaDMRG<<endl;
  file<<"Numfiles:          "<<Numfiles<<endl;
  file.close();
    

  const uint cpstep = step;
  const bool spg=true;
  const bool ssf=false;
  const bool csf=false;
  const bool dosf=false;
  const bool testenergy=true;
  const double conv = 1e-10;
  const bool normalize=true;
  const double landelta = 1e-10;
  //consistency check:
  string bath = "bath_params_"+simname;
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

  //check if there's a simnamedt..._CP.container file to be resumed ....
  std::ostringstream convert;
  convert<<dt;
  string name=simname+"dt"+convert.str()+"_CP.container";
  cout<<"#starting 1-orbital tdSIAM.exe solver using "<<omp_get_max_threads()<<" parallel threads"<<endl;
  if(!fexist(name)){
    // ========================================================   Fresh run ===============================================================
    cout<<"starting fresh run"<<endl;
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
    SpinlessFermionsMPO<Complex> cmpo(U,hop,mu);

    //============state initialization========
    UintVectorType dimP,localState;
    std::vector<i_to_key<AbKeyType > > itokey;
    dimP.resize(2*N);
    localState.resize(2*N);
    boost_setTo(dimP,uint(2));
    boost_setTo(localState,uint(0));
    itokey.resize(2*N);
    for(uint site = 0;site<2*N;site++){
      AbKeyType k1(0,1);    
      AbKeyType k2(1,1);
      itokey[site][0] =k1; 
      itokey[site][1] =k2;
    }
    for(uint site = 0;site<N;site++){
      localState(site)=site%2;
      localState(N+site)=site%2;//this is half filling
    }

    cout<<"#the initial state for dmrg:"<<endl;
    cout<<localState<<endl;

    MPS<AbKeyType,Real>mps(dimP,localState,itokey),tmps;
    DMRGengine<AbKeyType,Real> dmrg(mps,mpo,simname,verbose);
    string dmrgGS=simname+"_singleSite.container";
    string sssname=simname+"_singleSite";
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

    //define the particle creation/annihilation operators to create the initial state 
    MPO<Real> Op;
    MPO<Real> SpinOp;
    for(uint site=0;site<2*N;++site){
      MPOMat<1,Real>mat;
      SparseLocalOperator<Real> O_11(2,2),cD(2,2),c(2,2),p(2,2),n(2,2);
      O_11[OpKeyType<1>(0,0)]=1.0;O_11[OpKeyType<1>(1,1)]=1.0;
      p[OpKeyType<1>(0,0)]=1.0;p[OpKeyType<1>(1,1)]=-1.0;
      cD[OpKeyType<1>(1,0)] = 1.0;
      c[OpKeyType<1>(0,1)] = 1.0;
      n[OpKeyType<1>(1,1)]=1.0;
      if(spg){
	if(site<(N-1)){
	  mat.clear();
	  mat[ukey<2>(0,0)]=p;
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
      }else if(ssf){
	if(site<(N-1)){
	  mat.clear();
	  mat[ukey<2>(0,0)]=O_11;
	  mat.Dl=1;mat.Dr=1;
	}
	if(site==N-1){
	  cout<<"#calculating spin-structure factor"<<endl;
	  mat.clear();
	  mat[ukey<2>(0,0)]=n;
	  mat[ukey<2>(0,1)]=O_11;
	  mat.Dl=1;mat.Dr=2;
	}
	if(site==N){
	  mat.clear();
	  mat[ukey<2>(0,0)]=O_11;
	  mat[ukey<2>(1,0)]=(-1.0)*n;
	  mat.Dl=2;mat.Dr=1;
	}
	if(site>N){
	  mat.clear();
	  mat[ukey<2>(0,0)]=O_11;
	  mat.Dl=1;mat.Dr=1; 
	}
	Op[site]=mat;
      }else if(csf){

	if(spin==-1){
	  if(site<(N-1)){
	    mat.clear();
	    mat[ukey<2>(0,0)]=O_11;
	    mat.Dl=1;mat.Dr=1;
	  }
	  if(site==(N-1)){
	    cout<<"#calculating charge-structure factor for "<<spin<<" spins"<<endl;
	    mat.clear();
	    mat[ukey<2>(0,0)]=n;
	    mat[ukey<2>(0,1)]=O_11;
	    mat.Dl=1;mat.Dr=2;
	  }
	  if(site>(N-1)){
	    mat.clear();
	    mat[ukey<2>(0,0)]=O_11;
	    mat.Dl=1;mat.Dr=1; 
	  }
	}else if(spin==1){
	  if(site<N){
	    mat.clear();
	    mat[ukey<2>(0,0)]=O_11;
	    mat.Dl=1;mat.Dr=1;
	  }
	  if(site==N){
	    mat.clear();
	    mat[ukey<2>(0,0)]=n;
	    mat[ukey<2>(0,1)]=O_11;
	    mat.Dl=1;mat.Dr=2;
	  }
	  if(site>N){
	    mat.clear();
	    mat[ukey<2>(0,0)]=O_11;
	    mat.Dl=1;mat.Dr=1; 
	  }
	}
	Op[site]=mat;
      }else if(dosf){

	if(site<(N-1)){
	  mat.clear();
	  mat[ukey<2>(0,0)]=O_11;
	  mat.Dl=1;mat.Dr=1;
	}
	if(site==(N-1)){
	  cout<<"#calculating doubleoccupancy-structure factor"<<endl;
	  mat.clear();
	  mat[ukey<2>(0,0)]=n;
	  mat.Dl=1;mat.Dr=1;
	}
	if(site==N){
	  mat.clear();
	  mat[ukey<2>(0,0)]=n;
	  mat.Dl=1;mat.Dr=1;
	}
	if(site>N){
	  mat.clear();
	  mat[ukey<2>(0,0)]=O_11;
	  mat.Dl=1;mat.Dr=1; 
	}
	Op[site]=mat;
      }
    }
    //now create the initial state
    tmps = applyMPO(Op,mps);

    //some checks:
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
    

    //Measurements
    std::vector<Operator<Complex>*> DEN,ENTROPY,HAM,NORM,OV;
    for(int i = 0;i<2*N;++i){
      SparseLocalOperator<Complex> n(2,2),e(2,2);
      e[OpKeyType<1>(0,0)]=(1.0);
      n[OpKeyType<1>(1,1)]=(1.0);
      n.setSite(i);
      e.setSite(i);
      DEN.push_back(new SparseLocalOperator<Complex>(n));
      if(i<(2*N-1)){
	//==============bipartite entanglement============
	Entropy<Complex> ent(i);
	ENTROPY.push_back(new Entropy<Complex>(ent));
      }
    }
    //now decide wether to do a storage-run or a measurement run; if dt>0, a measurement run is done which collects .mps files from a 
    //previously started dt<0 run.
    if(dt<0.0){
      //this run produces .mps files with a running index; they can be collected
      //note that only a single .container file is written
      const bool store=true;
      convert.str("");
      convert<<dt;

      //this run only stores the .container files at different breakpoints, given by checkpointstep
      measurementObject<Complex> mN(DEN,"N"+simname+"dt"+convert.str(),step,true,false);
      measurementObject<Complex> mEntropy(ENTROPY,"entropy_"+simname+"dt"+convert.str(),step,true,false);

      list<measurementObject<Complex> > mList;
      mList.push_back(mN);  
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
      measurementObject<Complex> mN(DEN,"N"+simname+"dt"+convert.str(),step,true,false);
      measurementObject<Complex> mEntropy(ENTROPY,"entropy_"+simname+"dt"+convert.str(),step,true,false);
      
      list<measurementObject<Complex> > mList;
      mList.push_back(mN);  
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

    cout<<"#exiting qis.exe ..."<<endl;
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
    //===============================measurements =======================================================
    std::vector<Operator<Complex>*> DEN,ENTROPY,OV;
    for(int i = 0;i<N;++i){
      SparseLocalOperator<Complex> n(2,2),e(2,2);
      e[OpKeyType<1>(0,0)]=(1.0);
      n[OpKeyType<1>(1,1)]=(1.0);
      n.setSite(i);
      e.setSite(i);
      DEN.push_back(new SparseLocalOperator<Complex>(n));
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
    measurementObject<Complex> mN(DEN,"N"+simname+"dt"+convert.str(),step);
    measurementObject<Complex> mEntropy(ENTROPY,"entropy_"+simname+"dt"+convert.str(),step);

    list<measurementObject<Complex> > mList;
    mList.push_back(mN);  
    mList.push_back(mEntropy);
    if(dt>0.0)
      mList.push_back(mOverlap);  
    bool recanonize=false;
    SIAM.simulate(mList,order,NmaxTEBD,dt,chimax,true,tr_weight,normalize,recanonize);
    cout<<"#finished simulation ..."<<endl;
  }


  return 0;
}

BOOST_PYTHON_MODULE(SO_solver_ext)
{
    using namespace boost::python;
    def("SOSolver", SOSolver);
}


