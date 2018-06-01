/*
this is a example file for how setting up a dmrg followed by tebd tine evolution
if anything is unclear or bugs are found, contact me (martin.ganahl@gmail.com) 
any time.
I extensively use std and boost libraries. for the latter, i made a lot of typedefs. these are all in
boost_typedefs.hpp file. you may want to check this file.

General things: 
I use a lot of different classes. A crude structure of the program is:

MPS related: 
There is a class called MPSBase, which is the base class for all Matrix Product states. It's a map which maps some key to boost matrices
the key is related to the particle number conservation. MPSBase is templated class: MPSBase<KeyType,T>, where KeyType is template paremeter
specifying the type of the key. The idea was to be able to include abelian and non-abelain quantum numbers. However, the latter ones would
require a major rewriting of the whole code. The template parameter T can be either Real or Complex.
Derived from MPSBase class are MPS<KeyType,T> and CanonizedMPS<KeyType,T>. 
MPSBase also contains a set of lambdas (schmidt-eigenvalues).
For a full list of feature of these classes, check out the file MPSTypes.hpp.


there are a number of operations one can perform on MPS, listed in the following:

orthonormalize(MPSBase<KeyType,T>mps,string direction,bool storenorm=false):
Arguments: 
mps: is the mps; note that both CanonizedMPS and MPS can be passed. For the first, orthonormalization takes place on the Gamma matrices, which means
that the state is actually changed!!!
dir: direction of orthonormalization: either "lr" or "rl" (left-right of right-left); IMPORTANT: DO NOT ORTHONORMALIZE TWICE IN THE SAME DIRECTION!!!! this will result in a phase 
inconsistency of the lambdas and the Gamma matrices to the left and to the right of each lambda. This inconsistency can be removed by applying a single orthonormalization into the other 
direction after any number of orthonormalizations in the other direction have been applied. 
storenorm: if true, the norm of the mps is stored in the mps-variable _norm, which can be assessed by mps.norm().

The function sets a mps-intern flag so that the mps knows in which state it is (left or right orthonormal). The final state is either in A...AAA ((dir="lr") or B...BBB (dir="rl") form.
lambdas are always stored in the mps; note that if mps is in no particular orthonormal form prior to orthonormalize, the lambdas after orthonormalize cannot be interpreted as schmidtvalues.
the norm of mps after this routine is 1. The old norm can be assessed by mps.norm(). The function prepareSite is defined in MPSTypes.hpp.

CanonizedMPWS<KeyType,T> = canonize(MPS<KeyType,T> mps):
this brings mps brings it into a so called canonical form (Gamma Lambda Gamma Lamba ... Lambda Gamma). 
canonize doesn't change the state but gives back an object of class CanonizedMPS<KeyType,T>, where T is either Real or Complex. This can then be used for
time evolution.

All Operator related things are in MPOTypes.hpp. There is class MPO<T> with T either real or complex, which is the base class for a matrix Product Operators
there are no operations defined for this class. Hamiltonians are then defined in Hamiltonians.hpp. every hamiltonian has to be derived from MPO<T>

for efficiency, i use sparse operators. In the end, every operator is a map from some key to a real or complex number. I know that the labeling of all the keys 
is confusing, sorry for that.... However, besides that, it's pretty straightforward: you can insert elements into all kinds of operators and MPO-matrices by
using the corresponding key and a matrix-like indexing. see the Hamiltonians.hpp file, it's all there ...

the third big junk is defined in container.hpp.
Every kind of simulation, like DMRG or TEBD got it's own engine. So to run a simulation, you have to create an object of the corresponding class and initialize it
(which is in most cases done by the constructor). the you have to invoke the memberfunction .simulate(some parameters) which will start the simulation.
note that all engines are derived from a class called Container<KeyType,T>, where KeyType refers to certain type of key used (at the moment we can only use AbKeyType)
and T is either Real or Complex. Container contains objects of type MPS,MPO and also other stuff. simulation will in general change these objects

You can read and write Containers, and there is a checkpointing for the time evolution. Note that the members save and load only saves and loads the MPS object within 
the Container.

Stuff about DMRG and lanzcos:
I wrote a class called LanzcosEngine, and from that derived classes TWLanzcosEngine and SSLanzcosEngine which are used in the two site and single site dmrg approach, 
respectively. both can be found in lanzcos.hpp. they are initialized with some single or two site mps and mpo matrices, and by calling simulate(...) the lanzcos 
run starts.
The dmrg is then actually be done in container.hpp within the DMRGengine. 
All function needed for contractions and such stuff are in MPSfunctions.hpp, also the applyMPO routine. 

Although i use the boost library, i actually coded the matrix mutiplications myself. I use BLAS and LAPACK for all vectors-space algebraic stuff like
diagonalizations or singular value decompositions (SVD). These are defined in blas_routines.hpp/.cpp and lapack_routines.hpp/.cpp.

I also have designed a class BlockTensor, defined in BlockTensor.hpp, which is kind of block diagonal matrix, and is used as the mps-matrix in class
MPS<KeyType,T>.

adding two MPS is also possible, there is an operator + for that. use mps3=mps1+mps2; for that. afterwards you may call orthonormalize(mps3,"lr");
orthonormalize(mps3,"rl");

applying an MPO with bond dimension D to an mps with bond dimension chi increases the bond dimension of the resulting mps to 
chi*D. multiple appications are thus infeasible. you can however compress the mps again by using the routine 
compress(mps,chimax,maxstep,delta); it can only handle MPS, not CanonizedMPS!. it reduces the dimension of an mps to chimax by an iterative procedure very
much like a dmrg run. the number of sweeps is specified by maxstep. the last parameter delta is a truncation parameter. to initialize the compression,
an initial guess for the compressed state is needed. it's smart to take the original state and on every bond reduce the size of the lambdas, which means
throwing away states 

I will perhaps write a manual at some stage, but at the moment, you will just have to ask me if anything is unclear.



 */

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

  //first, parse all the parameters needed; I'll not explain anything about model parameters in this file
  parparser parser;
  
  //DMRG parameters:
  parser.addoption("chiDMRG", Parameter_Integer,"100");
  //this is the maximum matrix dimension of the mps for dmrg runs

  parser.addoption("chiTEBD", Parameter_Integer,"100");
  //chiTEBD is the maximum matrix dimension of the mps during time evolution. large chiTEBD gives accurate results but also slows the simulation down (by a factor chi^3).



  parser.addoption("NmaxTEBD", Parameter_Integer,"1000");//the number of time steps you wich to do 
  parser.addoption("NmaxDMRG", Parameter_Integer,"5");   //the number of dmrg sweeps
  parser.addoption("doTrunc", Parameter_Option);//see one below
  parser.addoption("deltaDMRG", Parameter_Float,"1e-12");//if doTrunc is given, then the program tries to keep the truncated weight below this value, 
                                                         //until chiDMRG is reached.
  parser.addoption("doTrunc", Parameter_Option);//see one below
  parser.addoption("deltaTEBD", Parameter_Float,"1e-12");//that's like deltaDMRG, just for the TEBD runs with chiTEBD
  parser.addoption("dt", Parameter_Float,"0.05");//time step of the trotter expansion
  parser.addoption("N", Parameter_Integer, "100");   //length of the chain
  parser.addoption("V", Parameter_Float, "1.0");      //interaction parameter for the spinless fermions hamiltonian
  parser.addoption("t", Parameter_Float, "-1.0");    //hopping parameter
  parser.addoption("kd", Parameter_Float, "0.0");//this is for creation of an excitation. kd=a_k^dagger e.g. craetes a block electron in state k
  parser.addoption("k", Parameter_Float, "0.0"); //as kd but annihilation of a bloch electron
  parser.addoption("Osite", Parameter_Integer, "-1"); //that's for application of a localized excitation at site=Ositexb
  parser.addoption("msite", Parameter_Integer, "0"); //that's for something else

  parser.addoption("simname", Parameter_String, "SpinlessFermions");//that's obvious
  parser.addoption("measure", Parameter_Integer,"10");//this sets the measurement step; measurement will be down every measure-th step (10 is default) 
                                                      //during the time evolution
  parser.addoption("Nfermions", Parameter_Integer,"10");//that's the number of fermions in the system. will be used in localState-initialization
  parser.addoption("load", Parameter_String,"n");//you can save and load simulations
  parser.addoption("save", Parameter_Option);
  parser.addoption("density", Parameter_Option);//following flags specify if you want to measure the corresponding quantity during time evolution
  parser.addoption("entropy", Parameter_Option);
  parser.addoption("overlap", Parameter_Option);
  parser.addoption("energy", Parameter_Option);
  parser.addoption("hole", Parameter_Option);//adds a hole instead of a particle at Osite

  if (!parser.load(argc, argv))
    return 0;
  std::cout << "# Options:\n" << parser << std::endl;

  //now read in the parameters
  const uint chiTEBD = parser.getint("chiTEBD");
  const uint chiDMRG = parser.getint("chiDMRG");
  const uint maxstepsTEBD = parser.getint("NmaxTEBD");
  const uint maxstepsDMRG = parser.getint("NmaxDMRG");
  const double deltaDMRG = parser.getfloat("deltaDMRG");
  const double deltaTEBD = parser.getfloat("deltaTEBD");
  const uint N = parser.getint("N");
  const double jz = parser.getfloat("V");
  const double jxy = parser.getfloat("t");
  const double dt = parser.getfloat("dt");
  const double kd  = parser.getfloat("kd");
  const int Osite = parser.getint("Osite");
  const double k  = parser.getfloat("k");
  const bool doTrunc = parser.isgiven("doTrunc");
  const string simname = parser.getstring("simname");
  const uint measStep = parser.getint("measure");
  const uint QN = parser.getint("Nfermions");
  const uint msite = parser.getint("msite");
  const string load = parser.getstring("load");
  const bool save = parser.isgiven("save");
  const bool dens = parser.isgiven("density");
  const bool ent = parser.isgiven("entropy");
  const bool overlap = parser.isgiven("overlap");
  const bool ham = parser.isgiven("energy");
  const bool part = !parser.isgiven("hole");


  //==========always save the parameters in a parameter file=============================
  std::ofstream pfile;
  string pfilename="params_"+simname;
  pfile.open(pfilename.c_str(),ios_base::app);
  pfile<<"====================="+simname+"================="<<endl;
  pfile<<parser<<endl;
  pfile.close();

  VectorType Jz(N-1),Jxy(N-1),B(N);
  boost_setTo(Jz,jz);  boost_setTo(Jxy,jxy);  boost_setTo(B,0.0);
  // for(uint n=0;n<N-1;++n)
  //   Jz(n)+=c_rand(-0.1,0.1);
  // cout<<Jz<<endl;
  //this is the matrix product operator. there's a Real=double one, and Complex=std::complex one. the first is for the dmrg, which is twice as fast if you use real 
  //rather than complex numbers. the second one is for real time evolution, which has to be done with complex numbers.
  SpinlessFermionsMPO<Real> mpo(Jz,Jxy,B);
  SpinlessFermionsMPO<Complex> cmpo(Jz,Jxy,B);

  //============state initialization========
  //OK, so now let's define the parameters needed for initialization of an MPS:
  //dimP is a vector which specifies the local hilbert space dimension living at each bond; it can be any distribution
  //localState is a vector of length N of uints and is passed to the constructor of the mps. it specifies in which state 
  //a certain site is. To do that, I use a number coding (starting at 0) which maps every local state to a number: for the case of spinless fermios for
  //example, we have as local states |0>, |1> which means zero and one particle respectively. the apparent mapping is
  //|0>->0 and |1>->1. Thus localState contais a distribution of 0's and 1's, defining a product state where every site n (running from 0 to N-1) 
  //is in state localState(n). This initial distribution defines the number of conserved particles.

  //Conservation of particles is a symmetry of the Hamilonian. In general, symmetries can be exploited to improve efficiency of the program.
  //There are many symmetries that could be incorporated, e.g. particle-hole, parity, total momentum, z-komponent of angular momentum (total S_z), particle number ...
  //I only have implemented abelian symmetries: symmetries in general are groups. if you think of space rotations, they form a group. this group is however not
  //abelian which means that different elements of the group in general do not commute. In many cases such groups contain subgroups which do form an abelian 
  //group (for example rotations about a certain axis). Conservation of z-component of spin is an example of the conservation law accompanying such a symmetry.
  //conservation of number of electrons (or charge conservation) is another one. Non abelian symmetries are harder to implement (but will improve much on speed).

  //during the simulation, the program needs to know how the local states of the system (e.g. the states |0> and |1>) are connected with the conserved quantities
  //so you have to define a mapping which maps from the local states (more precisely from the numbers representing the local states) to the local symmetry of the 
  //local state. This is what the std::vector<i_to_key<AbKeyType> > itokey does. now for all abelian types of keys, there is class called AbKeyType with which this
  //can be done. AbKeyType contains an array which stores the number of particles in a certain local state (this can be many if there are many different particles
  //present). in our case, it just contains one number which is either 0 or 1, meaning that there can be zero or one particle on a certain site.
  //the AbKeyType object is created by AbKeyType k1(0,1) for the case of the empty state, the first number gives the number of particles in the state, and the 
  //second gives the number of different flavours (or particles) which are actually there in the system. if there were both down and up spin electrons (which were 
  //both conserved) then it should be two.
 
  UintVectorType dimP(N),localState(N),localState2(N);
  boost_setTo(dimP,uint(2));
  boost_setTo(localState,uint(0));
  boost_setTo(localState2,uint(0));
  int mod=static_cast<int>(floor(1.0*N/QN));
  std::vector<i_to_key<AbKeyType > > itokey(N);
  uint qn = 0;
  for(uint site = 0;site<N;site++){
    //since sites can in general differ from each other, we have to define this mapping for every site of the system
    AbKeyType k1(0,1);    
    AbKeyType k2(1,1);
    itokey[site][0] =k1; //here the mapping at site "site" is defined: state "0" (=|0>) has an AbKeyType k1, which means that there is no particle there
    itokey[site][1] =k2;//state "1" (=|1>) on the other hand has an AbKeyType k2, thus there is particle
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
  //wo now we can define the mps.
  MPS<AbKeyType,Real>mps(dimP,localState,itokey);
  //and initialize a dmrgengine with it.
  DMRGengine<AbKeyType,Real> HSim(mps,mpo,simname);

  //and finally, we just have to run it.
  if(maxstepsDMRG>0){
    HSim.simulate(maxstepsDMRG,1e-8,chiDMRG,doTrunc,deltaDMRG);
  }
  //some notes on reading and writing the container: it works! you may however get segmentation faults when reading data from a previous run IF the system size N
  //given by input differs from the N of the previous run, because N goes into the operator definition of for example the gaussian creation operator. If its size
  //differs from that of the mps, it can of course not be applied to it.
  if(save){
    cout<<"#saving Container "<<simname<<endl;
    HSim.write();
  }
  if(load!="n"){
    cout<<"#loading Container "<<load<<endl;
    HSim.setName(load);
    HSim.read();
  }
  // to do a time evolution with tebd we have to cast the mps from Real to Complex: that is done here  ...
  MPS<AbKeyType,Complex>cmps = cast<AbKeyType,Complex,Real>(HSim.getMPS());
  MPS<AbKeyType,Complex>tmps;
  //the next two guys are creation and annihilation operators of "bloch" states (on open boundary systems). They are defined in Hamiltonians.hpp
  //both are derived from class MPO
  ck_dagger ck_d(N,kd);
  ck c_k(N,k);
  //in the following, i construct a localized excitation which either creates or annihilates (depending on input parameters) a particle at site Osite
  //there is a function called applyMPO (defined in MPSfunctions.hpp) which can apply an MPO to an MPS. It can ONLY apply an MPO, so i first have to construct
  //the MPO which creates or annihilates a particle
  MPO<Complex> Op,Op2,Op3;  
  if(Osite>=0){
    for(uint site=0;site<N;++site){
      MPOMat<1,Complex>mat,mat2,mat3;
      SparseLocalOperator<Complex> O_11(2,2),cD(2,2),c(2,2),p(2,2),bla(2,2);
      O_11[OpKeyType<1>(0,0)]=1.0;O_11[OpKeyType<1>(1,1)]=1.0;
      p[OpKeyType<1>(0,0)]=1.0;p[OpKeyType<1>(1,1)]=-1.0;
      cD[OpKeyType<1>(1,0)] = 1.0;
      c[OpKeyType<1>(0,1)] = 1.0;
      bla[OpKeyType<1>(1,0)] = 1.0;bla[OpKeyType<1>(0,1)] = 1.0;
      if(site==Osite&&part){
	mat.clear();
	cout<<"applying particle creation"<<endl;
	mat[ukey<2>(0,0)]=cD;
	mat2[ukey<2>(0,0)]=c;
	mat3[ukey<2>(0,0)]=bla;
      }
      if(site==Osite&&(!part)){
	mat.clear();
	cout<<"applying particle destruction"<<endl;
	mat[ukey<2>(0,0)]=c;
	mat2[ukey<2>(0,0)]=cD;
	mat3[ukey<2>(0,0)]=bla;
      }
      if(site<Osite){
	mat[ukey<2>(0,0)]=p;
	mat2[ukey<2>(0,0)]=p;
	mat3[ukey<2>(0,0)]=p;
      }
      if(site>Osite){
	mat[ukey<2>(0,0)]=O_11;
	mat2[ukey<2>(0,0)]=O_11;
	mat3[ukey<2>(0,0)]=O_11;
      }
      mat.Dl=1;mat.Dr=1; 
      mat2.Dl=1;mat2.Dr=1; 
      mat3.Dl=1;mat3.Dr=1; 
      Op[site]=mat;
      Op2[site]=mat2;
      Op3[site]=mat3;
    }
  }
  //the actual applycation of all the above defined MPO's is done here
  if(kd!=0.0&&k==0.0){
    tmps = applyMPO(ck_d,cmps);
    while(!checkMPS(tmps));
    orthonormalize(tmps,"lr");
    orthonormalize(tmps,"rl");
    orthonormalize(tmps,"lr");
    cmps=tmps;
  }else if(kd==0.0&&k!=0.0){
    tmps = applyMPO(c_k,cmps);
    while(!checkMPS(tmps));
    orthonormalize(tmps,"lr");
    orthonormalize(tmps,"rl");
    orthonormalize(tmps,"lr");
    cmps=tmps;
  }else if(kd!=0.0&&k!=0.0){
    cout<<"#kd and k are both !=0; set one of them to zero"<<endl;
  }
  if(kd==0.0&&k==0.0){
    if(Osite>=0){
      tmps = applyMPO(Op,cmps);
      while(!checkMPS(tmps));
      orthonormalize(tmps,"lr");
      orthonormalize(tmps,"rl");
      orthonormalize(tmps,"lr");
      cmps=tmps;
    }
  }
  //the function checkMPS takes an mps and fixes certain things that can emerge after application of operators. 
  //run it whenever you think that the state may be messed up for any reason. it also prints some information of 
  //what exactly is wrong with it. if the state is OK, it returns true, if not false is returned
  while(!checkMPS(cmps));
  cout<<"#norm of ck_dagger|0> before normalization: "<<sqrt(computeOverlap(cmps,cmps))<<endl;
  cout<<"#total quantum number of the state: "<<cmps.getQN()<<endl;

  //==========================measurements================================
  //This is where all measurements are defined. I tried to to implement measurments to be as close as possible to the way one does it in 
  //quantum mechanics. So measurements are connected with operators. if you want to measure anything, you first have to define the according
  //operator. this has already been done for the case of SparseLocalOperator's and MPO's (they are defined in MPOTypes.hpp). 
  //It's very unlikely that you will have to define a new operator
  //since most of the measurements can be done with SparseLocalOperator's and MPO's. Both these classes are derived from Operators. suppose you have constructed
  //a sparse local operator,let's say a density at site k: SparseLocalOperator<Real> n_k(2,2); n_k.setSite(k);n_k(OpKeyType<1>(1,1)]=1.0; suppose you have a state
  //|psi> stored as an MPS in a variable "mps" and want to measure <psi|n_k|psi>. The class SparseLocalOperator (as well as MPO) has a public member named 
  //measure(CanonizedMPS<KeyType,T>) which takes a canonized mps and computes <psi|n_k|psi>. so you just have to call Real expectation_value=n_k.measure(mps). 
  //that's all. a very important point is the fact that SparseLocalOperator
  //ONLY can handle CanonizedMPS at the moment. please note that the function measure(CanonizedMPS) is overloaded with measure(MPS), so it will compile if you pass
  //an MPS instead of an CanonizedMPS, but it will give back the value -1000 (or +1000, i don't recall), indipendent of the operator details.
  //measuring an MPO is done in EXACTLY the same way. The only difference here is that MPO can handle BOTH MPS and CanonizedMPS, so you can actually pass either of
  //the two.

  //to make life easier (although at first glance this might not seem to be the case), you can pack different measurements together into larger objects.
  //consider the follwing scenario: you want to measure local density at every site, every n_d-th step, then you want to measure spin at every second site 
  //at every n_s-th step. If you want to measure something, you have to take the time evolved state and for every single operator have to call .measure(canonizedmps)
  //at the desired time step. this can become pretty unhandy if you have large systems and many different measurements. so i designed a class 
  //measurementObject<T>, with T either Real or Complex, which collects such measurements. every measurementObject has a certain name which you can specify (for
  //example "density"). It has been designed to collect measurments which are of the same type into a single object. by "same type" i mean for example 
  //all local density measurements: they acutally differ only in the site at which they measure the density. so i pack all density measurements into a single 
  //measurementObject with name density, and this object then measures all the different density-operators and writes them into a single file, one after the other.
  //this is done by invoking the measurementObject-member measure(CanonizedMPS) or measure(MPS) (the two functions are overloaded), just like in the case
  //of Operators (forgive the extensive use of the name "measure" as function name).
  //the way this is achieved is, for example for the density operators, by storing all these operators in a vector, put this vector together with an appropriate
  //name for the measurements and the desired time steps at which you want to measure into a measurementObject and pass it to the simulation container
  //which will call the memberfunction measure(mps) as indicated by the desired time step. But because i want to use this measurementObject not only for 
  //SparseLocalOperators but for ALL kinds of operators, i have to store the POINTERS to the operators instead of the operators itself. so i use a   
  //std::vector<Operator<Complex>*> instead of   std::vector<Operator<Complex> >  (note the missing *). This in turn forces me to use the C++ operator "new" below, 
  //in case you may wonder about all this ... And because one usually has many different kinds of such measurements (like for example local density AND local spin)
  //I create a measurementObject for each with a certain filename and a certain mesaurement step, put them all into a std::list (called mList below) and pass this
  //monstrum to the TEBDcontainer which then does the measurements as described above.

  //by calling the Container-memberfunction measureCMPS(std::list<measurementObject<T> >&), the container does the measurement and stores it in a buffer.
  //by calling Container member writeMeasurements(std::list<measurementObject<T> >&), it writes the buffers to files.
  //remark:every measurmentObject has an internal incremental variable "inc" which can be incremented. This can be used to remember how often on tried to measure 
  //the object, and e.g. measure it only at every n_d-th or n_s-th step.
  
  std::vector<Operator<Complex>*> DENS,ENTROPY,OV,HAM,MNO,MPODENS;
  std::vector<Operator<Real>*> RDENS,RMPODENS;

  HAM.push_back(new MPO<Complex>(cmpo));
  MPO<Real>mpodens;
  for(int i = 0;i<N;++i){
    MPOMat<1,Real>mat;
    SparseLocalOperator<Complex> n(2,2);
    SparseLocalOperator<Real> nr(2,2),id(2,2);
    n[OpKeyType<1>(0,0)]=0.0;n[OpKeyType<1>(1,1)]=1.0;
    nr[OpKeyType<1>(0,0)]=0.0;nr[OpKeyType<1>(1,1)]=1.0;
    id[OpKeyType<1>(0,0)]=1.0;id[OpKeyType<1>(1,1)]=1.0;
    n.setSite(i);//this is important because the measurement-function lookes for the mps-matrices at site i to compute the density at site i
    nr.setSite(i);

    DENS.push_back(new SparseLocalOperator<Complex>(n));
    RDENS.push_back(new SparseLocalOperator<Real>(nr));
    for(uint s=0;s<N;++s){
      id.setSite(i);
      nr.setSite(i);
      if(s!=i){
	mat.clear();
	mat[ukey<2>(0,0)]=id;
      }
      if(s==i){
	mat.clear();
	mat[ukey<2>(0,0)]=nr;
      }
      mat.Dl=1;mat.Dr=1;
      mpodens[s]=mat;
    }
    RMPODENS.push_back(new MPO<Real>(mpodens));
    if(i<(N-1)){
      //==============bipartite entanglement============
      Entropy<Complex> ent(i);
      ENTROPY.push_back(new Entropy<Complex>(ent));
    }
  }
  orthonormalize(cmps,"rl");
  orthonormalize(cmps,"lr");
  while(!checkMPS(cmps));
  cout<<"#norm of ck_dagger|0> after normalization: "<<sqrt(computeOverlap(cmps,cmps))<<endl;
  cout<<"#is state OK? "<<checkMPS(cmps)<<endl;

  //now we can take the state, obtained by application of an operator to the dmrg-groundstate, and put it into a TEBDengine to do time evolution with it.
  //note that the constructor of TEBDengine takes an MPS<KeyType,T> and NOT a CanonizedMPS<KeyType,T>, although to do time evolution, a canonized MPS
  //would really be needed. It stores this mps in an internal private object _mps of type MPS<KeyType,T>.
  //TEBDengine however has also a private member object _cmps of type CanonizedMPS<KeyType,T>, and by calling TEBDengine<KeyType,T>.canonizeMPS(), the state
  //_mps is canonized and stored in _cmps, which is then time evolved.
  TEBDengine<AbKeyType,Complex> HSimC(cmps,cmpo,simname);
  //this is only another measurement; Overlap is also derived from Operator
  //The overlap object contains an MPS<KeyType,T>_mps, and by calling Overlap.measure(mps), it returns <_mps|mps>
  Overlap<Complex> ov(cmps);
  OV.push_back(new Overlap<Complex>(ov));
  string filename = simname;
  measurementObject<Complex> mDens(DENS,"n_"+simname,measStep);
  measurementObject<Real> mRDens(RDENS,"nreal_"+simname,measStep);
  measurementObject<Real> mRMPODens(RMPODENS,"nmpo_"+simname,measStep);
  measurementObject<Complex> mEntropy(ENTROPY,"entropy_"+simname,measStep);
  measurementObject<Complex> mOverlap(OV,"overlap_"+simname,1);
  measurementObject<Complex> mHam(HAM,"Ham_"+simname,measStep);
  list<measurementObject<Complex> > mList;
  list<measurementObject<Real> > mRList;

  mRList.push_back(mRMPODens);
  if(dens)
    mList.push_back(mDens);
  if(ent)
    mList.push_back(mEntropy);
  if(overlap)
    mList.push_back(mOverlap);
  if(ham)
      mList.push_back(mHam);

  //that's the orthonormalization. "rl" means orthonormalizing form right to left, thus leaving the state in BBBBB form, where as 
  //"lr" goes from left to right and leaves the state in AAAAA form.
  orthonormalize(HSimC.getMPS(),"rl");
  orthonormalize(HSimC.getMPS(),"lr");
  //here the function canonize is invoked. it creates a canonized MPS from the usual MPS (stored in _mps in the TEBDengine)
  HSimC.canonizeMPS();
  maxNorm<Complex> mnorm (HSimC.getCMPS()[msite]*HSimC.getCMPS().lambda(msite),msite);
  MNO.push_back(new maxNorm<Complex>(mnorm));
  measurementObject<Complex> mMaxNorm(MNO,"maxNorm_"+simname,measStep);
  if(msite>0)
    mList.push_back(mMaxNorm);
  cout<<"#is state OK? "<<checkMPS(HSimC.getCMPS())<<endl;
  HSimC.measureCMPS(mList);
  HSimC.writeMeasurements(mList);
  //arguments of simulate:
  //the mList is the list of measurementObjects desribed above
  //maxstepsTEBD is the number of time steps to be done
  //third argument HAS TO BE 1! At the moment, ONLY NEAREST NEIGHBOR INTERACTIONS can be time evolved. Someday in future it might be extended to 
  //                arbitrary interactions. The third argument tells the program how longranged the interactions are.
  //dt is the time step (usually around .01 -.05)
  //chiTEBD is the maximum chi at which the program starts to truncate the bond dimension. take it as large as possible to get good results.
  //   but note that computation time scales like chi^3.
  //doTrunc is a bool; if it's true then the program tries truncate in such a way that the truncated weight is smaller than deltaTEBD, under the constraint
  //   that the resulting chi is <= chiTEBD.
  HSimC.simulate(mList,"second",maxstepsTEBD,dt,chiTEBD,doTrunc,deltaTEBD,true,false);


  cout<<"#finished successfully"<<endl;
  return 0;
}
