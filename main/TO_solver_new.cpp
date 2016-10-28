#include "MPSTypes.hpp"
#include "MPSfunctions.hpp"
#include "Hamiltonians.hpp"
#include "BlockTensor.hpp"
#include "tensormap.hpp"
#include "container.hpp"
#include "utilities.hpp"
#include <iostream>
#include "parparser.hpp"


std::pair<double,std::vector<int> > findenergyminimum(double energy,int N0down,int N0up,int N0,int N1down,int N1up,int N1,MPO<Real> mpo,int dir,uint which,UintVectorType dimP ,std::vector<i_to_key<AbKeyType > > itokey,uint maxstepsDMRG,int chi,double conv,double deltaDMRG,double landelta){
  bool converged=false;
  std::vector<int> vec(4);
  double energynew;
  if(which==0){
    int N0downtemp=N0down;
    while(converged==false){
      if(dir>0){
	N0downtemp++;
      }else if(dir<0){
	N0downtemp--;
      }
      if ((N0downtemp<0)||(N0downtemp>(N0+1))){
	cout<<"can't change N0downtemp any more: aborting() at "<<N0downtemp<<" at direction "<<dir<<endl; 
      }
      UintVectorType localState=createTwoorbitalfilling(N0downtemp,N0up,N0,N1down,N1up,N1);
      MPS<AbKeyType,Real>mpstemp(dimP,localState,itokey);
      DMRGengine<AbKeyType,Real> dmrgtemp(mpstemp,mpo,"minimumsearch",false);

      energynew=dmrgtemp.simulate(maxstepsDMRG,conv,chi,true,deltaDMRG,landelta,true);
      cout<<"finished simulating dmrg at N0downtemp= "<<N0downtemp<<" at direction "<<dir<<" and energy "<<energynew<<endl;
      if (energy<energynew){
	converged=true;
	break;
      }
      energy=energynew;
    }
    if(dir>0)
      N0downtemp--;
    else if(dir<0)
      N0downtemp++;
    
    vec[0]=N0downtemp,    vec[1]=N0up,    vec[2]=N1down,    vec[3]=N1up;
  }
  if(which==1){
    int N0uptemp=N0up;
    while(converged==false){
      if(dir>0){
	N0uptemp++;
      }else if(dir<0){
	N0uptemp--;
      }
      if ((N0uptemp<0)||(N0uptemp>(N0+1))){
	cout<<"can't change N0uptemp any more: aborting() at "<<N0uptemp<<" at direction "<<dir<<endl; 
      }
      UintVectorType localState=createTwoorbitalfilling(N0down,N0uptemp,N0,N1down,N1up,N1);
      MPS<AbKeyType,Real>mpstemp(dimP,localState,itokey);
      DMRGengine<AbKeyType,Real> dmrgtemp(mpstemp,mpo,"minimumsearch",false);
      energynew=dmrgtemp.simulate(maxstepsDMRG,conv,30,true,deltaDMRG,landelta,true);
      cout<<"finished simulating dmrg at N0uptemp= "<<N0uptemp<<" at direction "<<dir<<" and energy "<<energynew<<endl;
      if (energy<energynew){
	converged=true;
	break;
      }
      energy=energynew;
    }
    if(dir>0)
      N0uptemp--;
    else if(dir<0)
      N0uptemp++;

    vec[0]=N0down,    vec[1]=N0uptemp,    vec[2]=N1down,    vec[3]=N1up;
  }
  if(which==2){
    int N1downtemp=N1down;
    while(converged==false){
      if(dir>0){
	N1downtemp++;
      }else if(dir<0){
	N1downtemp--;
      }
      if ((N1downtemp<0)||(N1downtemp>(N0+1))){
	cout<<"can't change N1downtemp any more: aborting() at "<<N1downtemp<<" at direction "<<dir<<endl; 
      }
      UintVectorType localState=createTwoorbitalfilling(N0down,N0up,N0,N1downtemp,N1up,N1);
      MPS<AbKeyType,Real>mpstemp(dimP,localState,itokey);
      DMRGengine<AbKeyType,Real> dmrgtemp(mpstemp,mpo,"minimumsearch",false);
      energynew=dmrgtemp.simulate(maxstepsDMRG,conv,30,true,deltaDMRG,landelta,true);
      cout<<"finished simulating dmrg at N1downtemp= "<<N1downtemp<<" at direction "<<dir<<" and energy "<<energynew<<endl;
      if (energy<energynew){
	converged=true;
	break;
      }
      energy=energynew;
    }
    if(dir>0)
      N1downtemp--;
    else if(dir<0)
      N1downtemp++;

    vec[0]=N0down,    vec[1]=N0up,    vec[2]=N1downtemp,    vec[3]=N1up;
  }

  if(which==3){
    int N1uptemp=N1up;
    while(converged==false){
      if(dir>0){
	N1uptemp++;
      }else if(dir<0){
	N1uptemp--;
      }
      if ((N1uptemp<0)||(N1uptemp>(N0+1))){
	cout<<"can't change N1uptemp any more: aborting() at "<<N1uptemp<<" at direction "<<dir<<endl; 
      }
      UintVectorType localState=createTwoorbitalfilling(N0down,N0up,N0,N1down,N1uptemp,N1);
      MPS<AbKeyType,Real>mpstemp(dimP,localState,itokey);
      DMRGengine<AbKeyType,Real> dmrgtemp(mpstemp,mpo,"minimumsearch",false);
      energynew=dmrgtemp.simulate(maxstepsDMRG,conv,30,true,deltaDMRG,landelta,true);
      cout<<"finished simulating dmrg at N1uptemp= "<<N1uptemp<<" at direction "<<dir<<" and energy "<<energynew<<endl;
      if (energy<energynew){
	break;
	converged=true;
      }
      energy=energynew;
    }
    if(dir>0)
      N1uptemp--;
    else if(dir<0)
      N1uptemp++;

    vec[0]=N0down,    vec[1]=N0up,    vec[2]=N1down,    vec[3]=N1uptemp;
  }
  std::pair<double,std::vector<int> > out;
  out.first=energy;
  out.second=vec;

  
  return out;
}



//reads in bath parameters in file "bath" and local Hamiltonian from file "Hlocal" and diagonalizes the non-interacting 
//bath + local site Hamiltonian for each bath.
//structure of bath:
//bath parameters for orbital (s,o) are arranged columnwise; s is the spin index, (-1 is down, 1 is up), and o is the orbital 
//index (ranging from 1 to Norb)
//
//  (-1,1)     (-1,1)  ... (-1,Norb)    (1,1)      (1,1)     ...  (1,Norb)
//  t    eps   t    eps    t     eps    t    eps   t    eps       t    eps
//
//t's are a column vectors of hoppings, eps are column vectors of local potentials in the corresponding baths

std::vector<int> bathoccupations(std::string bath,std::string Hlocal){
  MatrixType bathparams=readMatrix(bath);
  MatrixType H=readMatrix(Hlocal);

  uint x;
  int Nspbath=bathparams.size2()/2;
  int Norb=Nspbath/2;
   if(Nspbath==0){
    cout<<"number of baths is 0, or file "<<bath<<"has an odd number of columns"<<endl;
    abort();
  }
  std::vector<int>occs(Nspbath),lengths(Nspbath);
  for(uint b=0;b<Nspbath;b++){
    int x=0;
    while((fabs(bathparams(x,2*b))>1E-10)&&(x<bathparams.size1())){
      x++;
    }
    lengths[b]=x;
  }
  for(int b=0;b<Norb;++b){
    int N=std::max(lengths[b],lengths[b+Norb]);
    lengths[b]=N;
    lengths[b+Norb]=N;
  }
  for(uint b=0;b<Nspbath;++b){
    //note that the first hopping the bath is the hybridization with the bath
    MatrixTypeColMaj mat(lengths[b]+1,lengths[b]+1);
    boost_setTo(mat,0.0);
    VectorType E(lengths[b]);
    mat(0,0)=H(b,0);
    for(uint n=0;n<lengths[b];n++){
      mat(n+1,n+1)=bathparams(n,2*b+1);
    }
    for(uint n=0;n<lengths[b];n++){
      mat(n+1,n)=bathparams(n,2*b);
      mat(n,n+1)=bathparams(n,2*b);
    }
    BoostDiag(mat,E);
    uint occ=0;
    while(E(occ)<0.0){
      occ++;
    }
    occs[b]=occ;
  }
  return occs;
}
//int TOSolver(uint NmaxTEBD,uint chiDMRG,uint maxstepsDMRG, uint maxstepsSSDMRG, string sname, uint step, uint chimax, Real dt, Real tr_weight, int spin, int branch, int orbital, bool verbose,Real deltaDMRG,int Numfiles);

//#########################################################################################################################################
//sign of (branch*dt) determines wether the started run is a production (dt*branch<0) or collection (dt*branch>0) run
//sign of dt determines wether run is done forward (dt>0) or backward (dt<0) in time.
//
  //int TOSolver(uint NmaxTEBD,uint chiDMRG,uint maxstepsDMRG, uint maxstepsSSDMRG, string sname, uint step, uint chimax, Real dt, Real tr_weight, int spin, int branch, int orbital, bool verbose,Real deltaDMRG,int Numfiles){

int main(const int argc,const char **argv){

  parparser parser;
  parser.addoption("NmaxTEBD", Parameter_Integer,"500");
  parser.addoption("chiDMRG", Parameter_Integer,"100");
  parser.addoption("chiDMRG_scanning", Parameter_Integer,"10");
  parser.addoption("NmaxDMRG", Parameter_Integer,"5");
  parser.addoption("scan_particlenumber", Parameter_Option);
  parser.addoption("NmaxDMRGscanning", Parameter_Integer,"5");
  parser.addoption("NmaxSSDMRG", Parameter_Integer,"10");
  parser.addoption("label", Parameter_String,"default_name");
  parser.addoption("measure", Parameter_Integer,"1");
  parser.addoption("chiTEBD", Parameter_Integer,"200");
  parser.addoption("landelta", Parameter_Float,"1e-5");
  parser.addoption("dt", Parameter_Float,"0.05");
  parser.addoption("trweight", Parameter_Float,"1e-8");
  parser.addoption("spin", Parameter_Integer,"1");//for the down-spin spectral function
  parser.addoption("branch", Parameter_Integer,"1");//for the hole part of the spectral function
  parser.addoption("orbital", Parameter_Integer,"1");//choose an orbital
  parser.addoption("verbose", Parameter_Option);
  parser.addoption("info", Parameter_Option);
  parser.addoption("deltaDMRG", Parameter_Float,"1e-12");
  parser.addoption("Numfiles", Parameter_Integer,"10");

  //int TOSolver(uint NmaxTEBD,uint chiDMRG,uint maxstepsDMRG, uint maxstepsSSDMRG, string sname, uint step, uint chimax, Real dt, Real tr_weight, int spin, int branch, int orbital, bool verbose,Real deltaDMRG,int Numfiles){
  if (!parser.load(argc, argv))
    return 0;
  std::cout << "# Options:\n" << parser << std::endl;
  const uint NmaxTEBD = parser.getint("NmaxTEBD");
  const uint chiDMRG = parser.getint("chiDMRG");
  const uint chi = parser.getint("chiDMRG_scanning");
  const uint maxstepsDMRGinit = parser.getint("NmaxDMRG");
  const uint maxstepsDMRG = parser.getint("NmaxDMRGscanning");
  const uint maxstepsSSDMRG = parser.getint("NmaxSSDMRG");
  string simname = parser.getstring("label");
  const uint step = parser.getint("measure");
  const uint chimax = parser.getint("chiTEBD");
  const double landelta = parser.getfloat("landelta");
  const double dt = parser.getfloat("dt");
  const double tr_weight = parser.getfloat("trweight");
  const int spin=parser.getint("spin");
  const int branch=parser.getint("branch");
  const int orbital=parser.getint("orbital");
  const bool verbose=parser.isgiven("verbose");
  const bool info=parser.isgiven("info");
  const bool scan=parser.isgiven("scan_particlenumber");
  const double deltaDMRG = parser.getfloat("deltaDMRG");
  const int Numfiles=parser.getint("Numfiles");

  if(info==true){
    cout<<"##########################################################################################################################"<<endl;
    cout<<"#######                                            TWO-ORBITAL DMFT SOLVER                                             ###"<<endl;
    cout<<"#######            IMPORTANT NOTE: FOR SPIN-Z=0 AND HALFFILLING, EACH CHAIN HAS TO BE OF EVEN LENGTH                   ###"<<endl;
    cout<<"##########################################################################################################################"<<endl;
    cout<<"#######  Usage: --orbital: takes values 1,2; the orbital index of the Greens function                                  ###"<<endl;
    cout<<"#######                                                                                                                ###"<<endl;
    cout<<"#######         --spin: takes values -1,1;   the spin index of the Greens function                                     ###"<<endl;
    cout<<"#######                                                                                                                ###"<<endl;
    cout<<"#######         --branch: takes values -1,1;   the branch index of the Greens function (particle of hole part)         ###"<<endl;
    cout<<"#######                                                                                                                ###"<<endl;
    cout<<"#######         --scan: if given, the solver scans the particle number sectors to find the absolute minimum of energy  ###"<<endl;
    cout<<"#######                 an initial guess for the particle numbers is obtained by diagonalizing the non-interacting     ###"<<endl;
    cout<<"#######                 system. starting from that guess, it keeps on adding or removing particles of type alpha until ###"<<endl;
    cout<<"#######                 energy does not decrease anymore; it does so for all flavours alpha                            ###"<<endl;
    cout<<"#######                                                                                                                ###"<<endl;
    cout<<"#######         --label=name: label of the input files; bath_params_'name' (contains the bathparameters)               ###"<<endl;
    cout<<"#######                       structure:                                                                               ###"<<endl;
    cout<<"#######                                                                                                                ###"<<endl;
    cout<<"#######                        parameters for orbital (s,o) are arranged columnwise;                                   ###"<<endl;
    cout<<"#######                        s is the spin index, (-1 is down, 1 is up), and o is the orbital index                  ###"<<endl;
    cout<<"#######                        (ranging from 1 to Norb)                                                                ###"<<endl;
    cout<<"#######                        (-1,1)     (-1,1)  ... (-1,Norb)    (1,1)      (1,1)     ...  (1,Norb)                  ###"<<endl;
    cout<<"#######                        t    eps   t    eps    t     eps    t    eps   t    eps       t    eps                  ###"<<endl;
    cout<<"#######                                                                                                                ###"<<endl;
    cout<<"#######                        t's are a column vectors of hoppings, eps are column vectors of local potentials        ###"<<endl;
    cout<<"#######                        in the corresponding baths; the first element of each t is the hybridization of the     ###"<<endl;
    cout<<"#######                        of the local orbital o with the bath                                                    ###"<<endl;
    cout<<"#######                        H_local_'name' contains local Hamiltonian; structure of the file:                       ###"<<endl;
    cout<<"#######                        e0d                                                                                     ###"<<endl;
    cout<<"#######                        e1d                                                                                     ###"<<endl;
    cout<<"#######                        e0d                                                                                     ###"<<endl;
    cout<<"#######                        e1u                                                                                     ###"<<endl;
    cout<<"#######                        U                                                                                       ###"<<endl;
    cout<<"#######                        J                                                                                       ###"<<endl;
    cout<<"#######                                                                                                                ###"<<endl;
    cout<<"#######                                                                                                                ###"<<endl;
    cout<<"#######  GENERAL:                                                                                                      ###"<<endl;
    cout<<"#######  All output files are appended by an identifier branch*(orbital+(spin+1))                                      ###"<<endl;
    cout<<"#######  calculations of one lesser greens function should split up into a positive and negative time evolution:       ###"<<endl;
    cout<<"#######  if (dt<0 and branch=1) or (dt>0 and branch=-1), the solver produces .mps files that have to be read by a      ###"<<endl;
    cout<<"#######  subsequent run which the sign of the corresponding dt reversed. this second run produces the lesser Greens    ###"<<endl;
    cout<<"#######  function. the second run also produces file lasterased. DO NOT DELETE THIS FILE BEFORE THE RUN IS FINISHED    ###"<<endl;
    cout<<"#######  the file bath_params_'name' must not contain empty lines at the bottom; this will confuse the read-in         ###"<<endl;
    cout<<"#######  function!                                                                                                     ###"<<endl;
    cout<<"#######  If bath_params_'name' has hoppings that are 0, the baths are cut at the smallest possible position such       ###"<<endl;
    cout<<"#######  that baths belonging to the same orbital and opposite spin have the same length                               ###"<<endl;
    cout<<"##########################################################################################################################"<<endl;
    cout<<"##########################################################################################################################"<<endl;
    return 0;
  }



  //string simname(sname);
  std::ofstream file;
  std::string paramfile="Parameters_"+simname;
  file.open(paramfile.c_str());
  file<<"NmaxTEBD:          "<<NmaxTEBD<<endl;
  file<<"chiDMRG:           "<<chiDMRG<<endl;
  file<<"maxstepsDMRG:      "<<maxstepsDMRG<<endl;
  file<<"maxstepsSSDMRG:    "<<maxstepsSSDMRG<<endl;
  file<<"simname:           "<<simname<<endl;
  file<<"measurement:       "<<step<<endl;
  file<<"chiTEBD:           "<<chimax<<endl;
  file<<"dt:                "<<dt<<endl;
  file<<"truncated weight:  "<<tr_weight<<endl;
  file<<"spin:              "<<spin<<endl;
  file<<"branch:            "<<branch<<endl;
  file<<"orbital:           "<<orbital<<endl;
  file<<"verbose:           "<<verbose<<endl;
  file<<"deltaDMRG:         "<<deltaDMRG<<endl;
  file<<"Numfiles:          "<<Numfiles<<endl;
  file.close();
  
  const string order="second";
  const uint cpstep = step;
  //const double landelta = 1e-4;
  const bool testenergy=true;
  const double conv = 1e-10;
  const bool normalize=true;
  const bool doTrunc=true;


  string bath = "bath_params_"+simname;
  string H0 = "H_local_"+simname;
  std::ostringstream convert;
  convert<<branch*(orbital+(spin+1));
  std::string temp;
  temp="_spinorbital"+convert.str()+simname;
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
  convert.str("");
  convert<<dt;
  string name=simname+"dt"+convert.str()+"_CP.container";
  MPS<AbKeyType,Real>mpsfinal,tmps;
  cout<<"#starting two-orbital solver using "<<omp_get_max_threads()<<" parallel threads"<<endl;
  if(!fexist(name)){
    cout<<"starting ground state run"<<endl;
    twoOrbitalSIAMMPO2<Real> mpo(bath,H0);
    uint N=mpo.size();
    uint N0=mpo.N0();
    uint N1=mpo.N1();
    
    twoOrbitalSIAMMPO2<Complex> cmpo(bath,H0);

    std::vector<int> occs=bathoccupations(bath,H0);
    int N0down=occs[0];
    int N1down=occs[1];
    int N0up=occs[2];
    int N1up=occs[3];
    cout<<"following bath occupations have been determined from diagonalization of the bath"<<endl;
    cout<< occs<<endl;


    cout<<"#following parameters will be used in the solver:"<<endl;
    cout<<"#t-up   : "<<mpo._hopu<<endl;
    cout<<"#t-down : "<<mpo._hopd<<endl;
    cout<<"#mu-up  : "<<mpo._muu<<endl;
    cout<<"#mu-down: "<<mpo._mud<<endl;
 

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

    string sssname=simname+"_singleSite";
    string dmrgGS=simname+"_singleSite.container";
    double energy;
    if(!fexist(dmrgGS)){
      localState=createTwoorbitalfilling(N0down,N0up,N0,N1down,N1up,N1);
      //now check if a dmrg groundstate calculation already exists. if so, then load it, if not then do a
      //dmrg and save all the relevant data

      if(scan==true){
	cout<<"determining particle numbers .... "<<endl;
	MPS<AbKeyType,Real>mpstemp(dimP,localState,itokey);
	DMRGengine<AbKeyType,Real> dmrgtemp(mpstemp,mpo,simname,false);
	energy=dmrgtemp.simulate(maxstepsDMRGinit,conv,chi,doTrunc,deltaDMRG,landelta,testenergy);
	dmrgtemp.getMPS().clear();
	mpstemp.clear();
	cout<<"initial fillings: N0u "<<N0up<<", N0d "<<N0down<<", N1u "<<N1up<<", N1d "<<N1down<<endl;

	std::pair<double,std::vector<int> >out= findenergyminimum(energy,N0down,N0up,N0,N1down,N1up,N1,mpo,-1,0,dimP,itokey,maxstepsDMRGinit,chi,conv,deltaDMRG,landelta);
	N0down=out.second[0];
	energy=out.first;
	cout<<"current lowest energy "<<energy<<"and current fillings: N0u "<<N0up<<", N0d "<<N0down<<", N1u "<<N1up<<", N1d "<<N1down<<endl;
	cout<<endl;

	out= findenergyminimum(energy,N0down,N0up,N0,N1down,N1up,N1,mpo,1,0,dimP,itokey,maxstepsDMRGinit,chi,conv,deltaDMRG,landelta);
	N0down=out.second[0];
	energy=out.first;
	cout<<"current lowest energy "<<energy<<"and current fillings: N0u "<<N0up<<", N0d "<<N0down<<", N1u "<<N1up<<", N1d "<<N1down<<endl;
	cout<<endl;

	out= findenergyminimum(energy,N0down,N0up,N0,N1down,N1up,N1,mpo,-1,1,dimP,itokey,maxstepsDMRGinit,chi,conv,deltaDMRG,landelta);
	N0up=out.second[1];
	energy=out.first;
	cout<<"current lowest energy "<<energy<<"and current fillings: N0u "<<N0up<<", N0d "<<N0down<<", N1u "<<N1up<<", N1d "<<N1down<<endl;
	cout<<endl;

	out= findenergyminimum(energy,N0down,N0up,N0,N1down,N1up,N1,mpo,1,1,dimP,itokey,maxstepsDMRGinit,chi,conv,deltaDMRG,landelta);
	N0up=out.second[1];
	energy=out.first;
	cout<<"current lowest energy "<<energy<<"and current fillings: N0u "<<N0up<<", N0d "<<N0down<<", N1u "<<N1up<<", N1d "<<N1down<<endl;
	cout<<endl;

	out= findenergyminimum(energy,N0down,N0up,N0,N1down,N1up,N1,mpo,-1,2,dimP,itokey,maxstepsDMRGinit,chi,conv,deltaDMRG,landelta);
	N1down=out.second[2];
	energy=out.first;
	cout<<"current lowest energy "<<energy<<"and current fillings: N0u "<<N0up<<", N0d "<<N0down<<", N1u "<<N1up<<", N1d "<<N1down<<endl;
	cout<<endl;

	out= findenergyminimum(energy,N0down,N0up,N0,N1down,N1up,N1,mpo,1,2,dimP,itokey,maxstepsDMRGinit,chi,conv,deltaDMRG,landelta);
	N1down=out.second[2];
	energy=out.first;
	cout<<"current lowest energy "<<energy<<"and current fillings: N0u "<<N0up<<", N0d "<<N0down<<", N1u "<<N1up<<", N1d "<<N1down<<endl;
	cout<<endl;

	out= findenergyminimum(energy,N0down,N0up,N0,N1down,N1up,N1,mpo,-1,3,dimP,itokey,maxstepsDMRGinit,chi,conv,deltaDMRG,landelta);
	N1up=out.second[3];
	energy=out.first;
	cout<<"current lowest energy "<<energy<<"and current fillings: N0u "<<N0up<<", N0d "<<N0down<<", N1u "<<N1up<<", N1d "<<N1down<<endl;
	cout<<endl;

	out= findenergyminimum(energy,N0down,N0up,N0,N1down,N1up,N1,mpo,1,3,dimP,itokey,maxstepsDMRGinit,chi,conv,deltaDMRG,landelta);
	N1up=out.second[3];
	energy=out.first;
	cout<<"final lowest energy "<<energy<<" and final fillings: N0u "<<N0up<<", N0d "<<N0down<<", N1u "<<N1up<<", N1d "<<N1down<<endl;
	cout<<endl;
      }

      cout<<"starting procution runs with fillings N0u "<<N0up<<", N0d "<<N0down<<", N1u "<<N1up<<", N1d "<<N1down<<endl;
      localState=createTwoorbitalfilling(N0down,N0up,N0,N1down,N1up,N1);
      cout<<"#the initial state for dmrg:"<<endl;
      cout<<localState<<endl;
      cout<<"fillings: N0u "<<N0up<<", N0d "<<N0down<<", N1u "<<N1up<<", N1d "<<N1down<<endl;
      MPS<AbKeyType,Real>mps(dimP,localState,itokey);
      DMRGengine<AbKeyType,Real> dmrg(mps,mpo,simname,verbose);

      cout<<"#Starting dmrg run"<<endl;
      if(maxstepsDMRG>0){
	dmrg.simulate(maxstepsDMRG,conv,chiDMRG,doTrunc,deltaDMRG,landelta,testenergy);
      }
      SSDMRGengine<AbKeyType,Real> SSdmrg(dmrg.getMPS(),mpo,sssname,verbose);
      if(maxstepsDMRG>0){
	SSdmrg.simulate(maxstepsSSDMRG,1e-10,1e-12,true);
      }
      SSdmrg.write();

      mpsfinal=SSdmrg.getMPS();
      Real energy=mpo.measure(mpsfinal);
      cout<<"#GS-energy after two-site and single-site DMRG: "<<energy<<endl;
      mpo.shiftandrescale(energy,0.0,0.0);
      cmpo.shiftandrescale(Complex(energy,0.0),Complex(0.0,0.0),Complex(0.0,0.0));
      Real shifted_energy=mpo.measure(mpsfinal);
      cout<<"#GS energy has been shifted from "<<energy<<" to "<<shifted_energy<<"!"<<endl;
      cout<<"#number of particles: "<<mpsfinal.getQN()<<endl;
      std::ofstream outfile;
      while(!checkMPS(mpsfinal));
      orthonormalize(mpsfinal,"lr");
      orthonormalize(mpsfinal,"rl");
      while(!checkMPS(mpsfinal));


    }else if(fexist(dmrgGS)){
      cout<<"#found "<<dmrgGS<<" file. Reading in data:"<<endl;
      SSDMRGengine<AbKeyType,Real> SSdmrg(sssname);
      SSdmrg.read();

      mpsfinal=SSdmrg.getMPS();
      Real energy=mpo.measure(mpsfinal);
      mpo.shiftandrescale(energy,0.0,0.0);
      cmpo.shiftandrescale(Complex(energy,0.0),Complex(0.0,0.0),Complex(0.0,0.0));
      Real shifted_energy=mpo.measure(mpsfinal);
      cout<<"#GS energy has been shifted from "<<energy<<" to "<<shifted_energy<<"!"<<endl;
      cout<<"#number of particles: "<<mpsfinal.getQN()<<endl;
      std::ofstream outfile;
      while(!checkMPS(mpsfinal));
      orthonormalize(mpsfinal,"lr");
      orthonormalize(mpsfinal,"rl");
      while(!checkMPS(mpsfinal));
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
      if(site<(N0)){
	mat.clear();
	p.setSite(site);
	mat[ukey<2>(0,0)]=p;
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
	}else if(orbital==2){
	  p.setSite(site);
	  mat[ukey<2>(0,0)]=p;
	}
      }
      if(site==(N0+1)){
	mat.clear();
	if(part&&up&&(orbital==2)){
	  cout<<"#applying up particle creation in orbital "<<orbital<<endl;
	  cuD.setSite(site);
	  mat[ukey<2>(0,0)]=cuD;
	}else if((!part)&&up&&(orbital==2)){
	  cout<<"#applying up hole creation in orbital "<<orbital<<endl;
	  cu.setSite(site);
	  mat[ukey<2>(0,0)]=cu;
	}else if(part&&(!up)&&(orbital==2)){
	  cout<<"#applying down particle creation in orbital "<<orbital<<endl;
	  cdD.setSite(site);
	  mat[ukey<2>(0,0)]=cdD;
	}else if((!part)&&(!up)&&(orbital==2)){
	  cout<<"#applying down hole creation in orbital "<<orbital<<endl;
	  cd.setSite(site);
	  mat[ukey<2>(0,0)]=cd;
	}else if(orbital==1){
	  O_11.setSite(site);
	  mat[ukey<2>(0,0)]=O_11;
	}
      }
      if(site>(N0+1)){
	mat.clear();
	O_11.setSite(site);
	mat[ukey<2>(0,0)]=O_11;
      }
      mat.Dl=1;mat.Dr=1; 
      Op[site]=mat;
    }
    //now create the initial state; this state is common to the two tebd runs 
    tmps = applyMPO(Op,mpsfinal);

    while(!checkMPS(tmps));
    cout<<"before normalization: <0|c* c|0>= "<<computeOverlap(tmps,tmps)<<endl;
    orthonormalize(tmps,"lr");
    orthonormalize(tmps,"rl");
    cout<<"after normalization: <0|c* c|0> = "<<computeOverlap(tmps,tmps)<<endl;
    MPS<AbKeyType,Complex>cmps=cast<AbKeyType,Complex,Real>(tmps);
    MPS<AbKeyType,Complex>mps_GS=cast<AbKeyType,Complex,Real>(mpsfinal);
    tmps.clear();
    mpsfinal.clear();
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
    if(((dt<0.0)&&(branch==1))||((dt>0.0)&&(branch==-1))){
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
	cout<<"starting simulation"<<endl;
	SIAM.save(simname+"dt"+convert.str()+"_CP0");
	SIAM.simulate(mList,order,NmaxTEBD,branch*dt,chimax,true,tr_weight,normalize,recanonize);
      }
      cout<<"#finished simulation ..."<<endl;
    }else if(((dt<0.0)&&(branch==-1))||((dt>0.0)&&(branch==1))){
      //this run collects .mps files from a previous run and computes the overlap with the corresponding mps from a 
      //run with dt>0; 
      convert.str("");
      convert<<-dt;
      std::string filelabel=simname+"dt"+convert.str();
      convert.str("");
      convert<<(dt);
      OverlapWithStored<Complex> ov(filelabel,Numfiles);
      
      OV.push_back(new OverlapWithStored<Complex>(ov));

      measurementObject<Complex> mOverlap(OV,"overlap"+simname,cpstep,true,true);
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
	SIAM.simulate(mList,order,NmaxTEBD,branch*dt,chimax,true,tr_weight,normalize,recanonize);
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
    measurementObject<Complex> mOverlap(OV,"overlap"+simname,cpstep,true,true);
    measurementObject<Complex> mNd(NU,"Nup"+simname+"dt"+convert.str(),step,true,false);
    measurementObject<Complex> mNu(ND,"Ndown"+simname+"dt"+convert.str(),step,true,false);
    measurementObject<Complex> mEntropy(ENTROPY,"entropy_"+simname+"dt"+convert.str(),step);


    list<measurementObject<Complex> > mList;
    mList.push_back(mNu);  
    mList.push_back(mNd);  
    mList.push_back(mEntropy);
    if(dt>0.0)
      mList.push_back(mOverlap);  
    bool recanonize=false;
    cout<<"#resuming simulation in "<<name<<endl;
    SIAM.simulate(mList,order,NmaxTEBD,branch*dt,chimax,true,tr_weight,normalize,recanonize);
    cout<<"#finished simulation ..."<<endl;
  }
  cout<<"#exiting qis.exe ..."<<endl;
  return 0;
}



//int main(const int argc,const char **argv){
//  TOSolver(200,400,10, 10,"test_label",1,50,-0.05,1e-10,1,1,1,true,1e-10,10);
//  TOSolver(200,400,10, 10,"test_label",1,50,0.05,1e-10,1,1,1,true,1e-10,200);
//}

// BOOST_PYTHON_MODULE(TO_solver_ext)
// {
//     using namespace boost::python;
//     def("TOSolver", TOSolver);
// }

