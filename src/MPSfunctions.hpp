#ifndef MPSFUNCTIONS_H__
#define MPSFUNCTIONS_H__
#include "BlockTensor.hpp"
#include <set>
#include "MPSTypes.hpp"
#include "MPOTypes.hpp"
#include "blasroutines.hpp"
#include "lapack_routines.hpp"
#include "utilities.hpp"
#include "keyTypes.hpp"
#include <random>
#include <chrono>
//#include "mkl_boost_ublas_matrix_prod.hpp"

using std::set;
//using std::real;
using std::cout;using std::endl;
using boost::numeric::ublas::identity_matrix;


template<class KeyType,typename T>
bool prepareBlock(const set<TensorKeyType<2,1,KeyType> >out,BlockMatrix<1,KeyType,T>& BL1,BlockMatrix<1,KeyType,T>&BL2,BlockMatrix<1,KeyType,T>&BLp,DiagBlockMatrix<KeyType,Real>&lambda,const string dir,
		  const bool multiply,const bool truncate=false,const Real delta=0.0){  
  if (dir != "r" && dir != "l") {
    cout << "prepareBlock: wrong input parameter dir";
    abort();
  }
  bool erase = true;
  if (dir == "l") {
    KeyType keyout = (*out.begin())[1];
    typename BlockMatrix<1,KeyType,T>::iterator itBL;
    typename set<TensorKeyType<2,1,KeyType> >::const_iterator it=out.begin(),sit;;
    uint SIZE2=BL1[*it].size2(),SIZE1=0,start=0;
    range_map<TensorKeyType<2,1,KeyType> > rmap;
    for(it = out.begin();it!=out.end();++it){
      assert(BL1[*it].size2()==SIZE2);
      assert((*it)[1]==keyout);
      rmap[(*it)]=Range(SIZE1,SIZE1+BL1[*it].size1());
      SIZE1+=BL1[*it].size1();
    }
    matrix<T,column_major> tempmat(SIZE1,SIZE2);
    boost_setTo(tempmat,T(0.0));
    for(it = out.begin();it!=out.end();++it){
      matrix_range<matrix<T, column_major> >(tempmat,rmap[(*it)], Range(0,SIZE2))= BL1[*it];
    }
    matrix<T, column_major> U,Vd,Utmp,Vdtmp;
    VectorType D;
    //VectorType D;
    //TODO replace svd by QR decomposition
    //there is no truncation in prepareBlock!!!! this function ONLY normalizes a block!!!!!! truncation has to be done in another place!
    BoostSVD(tempmat, U, D, Vd);
    //discard all values smaller than delta
    //adapt size of matrices
    if(truncate){
      VectorType Dtmp;
      Dtmp=cut(D,delta);D = Dtmp;
    }
    SIZE1 = U.size1(), SIZE2 = D.size();
    if(truncate){
      cut(U, Utmp, SIZE1, SIZE2);U = Utmp;
      SIZE1 = D.size();SIZE2 =Vd.size2();
      cut(Vd, Vdtmp, SIZE1, SIZE2);Vd = Vdtmp;
    }
    SIZE2 = D.size();
    //this is kind of stupid, but i have no real tensor classes
    if(SIZE2){
      erase=false;
      lambda[keyout]=D;
      //    cout<<"the norm of the lambda"<<dot_prod(D.size(),D.data().begin(),D2.data().begin())<<endl;
      //the new U is inserted into the map by adressing each iterator to the right blocks
      //if D.size()==0, then the block has to be discarded at all.
      for(it = out.begin();it!=out.end();++it){
	matrix<T,column_major>tmp(matrix_range<matrix<T, column_major> >(U,rmap[(*it)], Range(0,SIZE2)));
	//matrix_range<matrix<T, column_major> >(U,rmap[itBL->first[0]], Range(0,SIZE2))
	BLp[*it]=tmp;
      } 
    }
    if(multiply&&SIZE2){
      for(itBL = BL2.begin();itBL!=BL2.end();++itBL){
       	if(itBL->first[0]==keyout){
       	  assert(itBL->second.size1() == Vd.size2());
       	  MatMatProd(T(1.0),Vd, itBL->second, tempmat, 'n', 'n');
	  MatMatProd(T(1.0),D, tempmat, itBL->second, 'n', 'n');
       	}
      }
    }
  }//end of dir==l
  else if (dir == "r") {
    KeyType keyout = (*out.begin())[0];
    typename set<TensorKeyType<2,1,KeyType> >::const_iterator it=out.begin();
    typename BlockMatrix<1,KeyType,T>::iterator itBL;
    uint SIZE1=BL2[*it].size1(),SIZE2=0,start=0;
    range_map<TensorKeyType<2,1,KeyType> > rmap;

    for(it = out.begin();it!=out.end();++it){
      assert(BL2[*it].size1()==SIZE1);
      assert(keyout==(*it)[0]);
      rmap[(*it)]=Range(SIZE2,SIZE2+BL2[*it].size2());
      SIZE2+=BL2[*it].size2();
    }
    //TODO: this loop can be parallelized
    matrix<T,column_major> tempmat(SIZE1,SIZE2);
    boost_setTo(tempmat,T(0.0));
    for(it = out.begin();it!=out.end();++it){
      matrix_range<matrix<T, column_major> >(tempmat, Range(0,SIZE1),rmap[(*it)])= BL2[*it];
    }
    matrix<T, column_major> U,Vd,Utmp,Vdtmp;
    VectorType D;
    //TODO replace svd by QR decomposition
    //there is no truncation in prepareBlock!!!! this function ONLY normalizes a block!!!!!! truncation has to be done in another place!
    BoostSVD(tempmat, U, D, Vd);
    //discard all values smaller than delta
    //adapt size of matrices
    if(truncate){
      VectorType Dtmp;
      Dtmp=cut(D,delta);D = Dtmp;
      SIZE1 = U.size1(), SIZE2 = D.size();
      cut(U, Utmp, SIZE1, SIZE2);U = Utmp;
      SIZE1 = D.size();SIZE2 = Vd.size2();
      cut(Vd, Vdtmp, SIZE1, SIZE2);Vd = Vdtmp;
    }
    SIZE1 = D.size();
    if(SIZE1){
      erase=false;
      lambda[keyout]=D;
      //the new U is inserted into the map by adressing each iterator to the right blocks
      //if D.size()==0, then the block has to be discarded at all.
      for(it =out.begin();it!=out.end();++it){
	matrix<T,column_major> tmp(matrix_range<matrix<T, column_major> >(Vd, Range(0,SIZE1),rmap[(*it)]));
	BLp[*it] = tmp;
      } 
      for(itBL = BL1.begin();itBL!=BL1.end();++itBL){
       	if(keyout==itBL->first[1]){
       	  assert(itBL->second.size2() == U.size1());
       	  MatMatProd(T(1.0),itBL->second, U, tempmat, 'n', 'n');
	  MatMatProd(T(1.0),tempmat, D, itBL->second, 'n', 'n');
       	}
      }
    }
  }//end of dir==r
  return erase;
}


//This routine prepares an MPS-matrix in either left (dir="l") or right (dir="r") form.
//Arguments: BL1, BL2 are two mps-matrices; if dir="l", then BL1 is left-orthonormalized, if dir="r", then BL2 is right-orthonormalized; if multiply=true, then the 
//remaining part of the svd of the MPS-matrix is multiplied to the other matrix (so as not to change the full mps-state). if truncate=true, then D-values of the SVD below delta are discarded
//(note that they not necessarily need to be schmidt-values)
template<class KeyType,typename T>
bool prepareMatrix(BlockMatrix<1,KeyType,T>& BL1, DiagBlockMatrix<KeyType,Real>&lambda,BlockMatrix<1,KeyType,T>&BL2,const string dir, const bool multiply,const bool truncate,const Real delta,const int maxthreads=omp_get_max_threads()){

  if(dir!="l"&&dir!="r"){cout<<"prepareMatrix: wrong input string dir = "<<dir<<"; use \"l\" or \"r\" instead"<<endl;abort();}
  bool erase = false;
  if(dir=="l"){
    typename BlockMatrix<1,KeyType,T>::iterator itBL;
    //uint _maxthreads=omp_get_max_threads();
    uint _maxthreads=maxthreads;
    const int oldmaxthreads=omp_get_max_threads();
    omp_set_num_threads(_maxthreads);
    DiagBlockMatrix<KeyType,Real> Lamp[_maxthreads];
    BlockMatrix<1,KeyType,T> BLp[_maxthreads];
    std::vector<KeyType> possiblekeys;
    std::set<KeyType>Er[_maxthreads];
    map_container<KeyType,set<TensorKeyType<2,1,KeyType> > >Tkeys;
    typename map_container<KeyType,set<TensorKeyType<2,1,KeyType> > >::iterator Tkit;
    for(itBL=BL1.begin();itBL!=BL1.end();++itBL){
      Tkeys[itBL->first[1]].insert(itBL->first);
    }
    for(Tkit=Tkeys.begin();Tkit!=Tkeys.end();++Tkit){
      possiblekeys.push_back(Tkit->first);
    }
#pragma omp parallel for
    for(uint i = 0;i<possiblekeys.size();i++){
      const uint jthread = omp_get_thread_num();
	std::set<TensorKeyType<2,1,KeyType> >keyset=Tkeys[possiblekeys[i]];
	std::set<KeyType>&er = Er[jthread];
	BlockMatrix<1,KeyType,T>& bl_p=BLp[jthread];
	DiagBlockMatrix<KeyType,Real>& lam_p = Lamp[jthread];
	erase = prepareBlock(keyset,BL1,BL2,bl_p,lam_p,dir,multiply,truncate,delta);
	if(erase)er.insert(possiblekeys[i]);
      }
    omp_set_num_threads(oldmaxthreads);
    //erase all blocks with zero size:
    typename std::set<KeyType>::iterator erit;
    std::vector<typename BlockMatrix<1,KeyType,T>::iterator> eraseBl;
    typename std::vector<typename BlockMatrix<1,KeyType,T>::iterator>::iterator eBl;
    for(uint i = 0;i<_maxthreads;++i){
      for(erit=Er[i].begin();erit!=Er[i].end();++erit)
	for(itBL=BL2.begin();itBL!=BL2.end();++itBL){
	  if (itBL->first[0]==*erit)eraseBl.push_back(itBL);
	}
    }
    for(eBl = eraseBl.begin();eBl!=eraseBl.end();++eBl){
      BL2.erase(*eBl);
    }
    BL1.clear();
    lambda.clear();
    //merge the blocks from each thread into one map
    for(uint i=0;i<_maxthreads;++i){
      merge(BL1,BLp[i]);
      merge(lambda,Lamp[i]);
      BLp[i].clear();
    }
    // BlockMatrix<1,KeyType,T> BLt;
    // if(multiply){
    //   BLt = lambda*BL2;
    // }
    // BL2=BLt;
  }//end of dir=="l"
  else if(dir=="r"){
    typename BlockMatrix<1,KeyType,T>::iterator itBL;
    uint _maxthreads=maxthreads;
    const int oldmaxthreads=omp_get_max_threads();
    omp_set_num_threads(_maxthreads);
    DiagBlockMatrix<KeyType,Real> Lamp[_maxthreads];
    BlockMatrix<1,KeyType,T> BLp[_maxthreads];
    std::vector<KeyType> possiblekeys;
    std::set<KeyType>Er[_maxthreads];
    map_container<KeyType,set<TensorKeyType<2,1,KeyType> > >Tkeys;
    typename map_container<KeyType,set<TensorKeyType<2,1,KeyType> > >::iterator Tkit;
    for(itBL=BL2.begin();itBL!=BL2.end();++itBL){
      KeyType key=itBL->first[0];
      Tkeys[key].insert(itBL->first);
    }
    for(Tkit=Tkeys.begin();Tkit!=Tkeys.end();++Tkit){
      possiblekeys.push_back(Tkit->first);
    }
#pragma omp parallel for
    for(uint i = 0;i<possiblekeys.size();i++){
	std::set<TensorKeyType<2,1,KeyType> >keyset=Tkeys[possiblekeys[i]];
	const uint jthread = omp_get_thread_num();
	std::set<KeyType>&er = Er[jthread];
	BlockMatrix<1,KeyType,T>& bl_p=BLp[jthread];
	DiagBlockMatrix<KeyType,Real>& lam_p = Lamp[jthread];
	erase = prepareBlock(keyset,BL1,BL2,bl_p,lam_p,dir,multiply,truncate,delta);
	if(erase)er.insert(possiblekeys[i]);
      }
    omp_set_num_threads(oldmaxthreads);
    //erase all blocks with zero size:
    typename std::set<KeyType>::iterator erit;
    std::vector<typename BlockMatrix<1,KeyType,T>::iterator> eraseBl;
    typename std::vector<typename BlockMatrix<1,KeyType,T>::iterator>::iterator eBl;
    for(uint i = 0;i<_maxthreads;++i){
      for(erit=Er[i].begin();erit!=Er[i].end();++erit)
	for(itBL=BL1.begin();itBL!=BL1.end();++itBL){
	  if (itBL->first[1]==*erit)eraseBl.push_back(itBL);
	}
    }
    for(eBl = eraseBl.begin();eBl!=eraseBl.end();++eBl){
      BL1.erase(*eBl);
    }

    BL2.clear();
    lambda.clear();
    //merge the blocks from each thread into one map
    for(uint i=0;i<_maxthreads;++i){
      merge(BL2,BLp[i]);
      merge(lambda,Lamp[i]);
      BLp[i].clear();
      Lamp[i].clear();
    }
    // BlockMatrix<1,KeyType,T> BLt;
    // if(multiply){
    //   BLt = BL1*lambda;
    // }
    // BL1=BLt;
  }//end of dir=="r"
  return erase;
}



//this is a subroutine used by prepareMatrix; it should not really matter to the user what it does
//WARNING: pass BLits NOT as reference; the #pragma parallel for which calls add_out_to_sub somehow causes the thing to crash
//when all threads operate on the same BLits; these are only iterators anyways, so copying them is cheap
template<class KeyType,typename T>
BlockMatrix<1,KeyType,T> add_out_sub_2(std::map<uint, std::vector<typename BlockMatrix<1,KeyType,T>::const_iterator> > BLits,
				       const BlockMatrix<1,KeyType,T>& A,const MPOMat<1,T>& mpo, const KeyType p,i_to_key<KeyType> itokey,
				       const int direction){
  assert(itokey.size());
  const int idir = (direction > 0 ? 1 : 0);
  BlockMatrix<1,KeyType,T> B;
  const int out_a = (direction > 0 ? 1 : 0);
  typename BlockMatrix<1,KeyType,T>::const_iterator a;
  //sort all matrices with outgoing p in B
  for(a=A.begin();a!=A.end();++a){
    if (a->first[out_a]==p)
      B[a->first] = a->second;
  }
  BlockMatrix<1,KeyType,T> Bout;
  //choose an mpo entry:
  for (typename MPOMat<1,T> ::const_iterator o=mpo.begin();o!=mpo.end();++o) {
    const ukey<2>& key_o = o->first;
    const SparseOperator<1,T>& lo = o->second;
    const uint o_out = key_o[idir], o_in = key_o[1-idir];
    //for that entry, run through the physical operator entries; this defines you the a lower tensorkeytype
    for(typename SparseOperator<1,T>::const_iterator oi=lo.begin();oi!=lo.end();++oi) {
      TensorKeyType<2,1,KeyType> key_b;
      //here the tensor keytype is defined; it's not sure yet if it can be contracted
      KeyType kyout=itokey[oi->first(0)];
      if (direction > 0){
	key_b = TensorKeyType<2,1,KeyType> (p-kyout,p,oi->first.getOut());
      }else {
	key_b = TensorKeyType<2,1,KeyType> (p, p + kyout,oi->first.getOut());
      }
      typename BlockMatrix<1,KeyType,T>::iterator posb = B.find(key_b);//now see if there really is a matrix with that key
      if (posb==B.end()) continue;//of not, then go to the next Operator entry
      // now run through all keys of BLits: so itBL is an iterator of BlockMatrix with outgoing auxmpo-ind o_in
      // and see if the key_b also matches any of the BL's
      typename std::vector<typename BlockMatrix<1,KeyType,T>::const_iterator>::iterator itBL;
      //            cout<<BLits.size()<<endl;
      for(itBL=BLits[o_in].begin();itBL!=BLits[o_in].end();++itBL) {
	const TensorKeyType<2,1,KeyType>& key_itBL = (*itBL)->first;
	if (direction >0 && (key_itBL(0)!=key_o[0]) ) continue;//compares if mpo-aux-ind of BL-block matches the mpo-aux-ind of the mpo
	if (direction <=0 && (key_itBL(0)!=key_o[1]) ) continue;//same as above, for the other direction
	if (direction >0 && (key_itBL[1]!=key_b[0]) ) continue;//look if the lower KeyType of the BL matches the chosen left-keytype of key_b of lower mat
	if (direction <=0 && (key_itBL[1]!=key_b[1]) ) continue;//same as above for other direction
	// finally, it's the right one.
	//now find an upper matrix that matches into the construction
	//the key is now defined by the other indices:
	TensorKeyType<2,1,KeyType> key_a;
	KeyType kyin=itokey[oi->first[0]];
	if (direction > 0)
	  key_a = TensorKeyType<2,1,KeyType> (key_itBL[0], key_itBL[0] + kyin,oi->first.getIn());
	else
	  key_a = TensorKeyType<2,1,KeyType> (key_itBL[0] - kyin, key_itBL[0],oi->first.getIn());
	typename BlockMatrix<1,KeyType,T>::const_iterator posa = A.find(key_a);
	if (posa==A.end()) continue;

	// ok. we've got everything.
	const matrix<T,column_major>&B  = (*itBL)->second;
	const matrix<T,column_major>&au = posa->second;
	const matrix<T,column_major>&ad = posb->second;
	matrix<T,column_major>Tmp1,Tmp2;
	if (direction > 0) {
	  MatMatProd(T(1.0),B, conj(au),Tmp1, 't', 'n');
	  MatMatProd(T(1.0),Tmp1,ad, Tmp2,'t', 'n');
	} else {
	  MatMatProd(T(1.0),conj(au),B,Tmp1, 'n', 'n');
	  MatMatProd(T(1.0),Tmp1,ad, Tmp2,'n', 't');
	}
	Tmp2 *= (oi->second);
	TensorKeyType<2,1,KeyType> ky(key_a[idir], key_b[idir], key_o[idir]);
	typename BlockMatrix<1,KeyType,T>::iterator ypos = Bout.find(ky);
	if (ypos != Bout.end() )
	  Tmp2 += ypos->second;
	Bout[ky] = Tmp2;
      }
    }
  }
  return Bout;    
}




//WARNING: pass BLits NOT as reference; the #pragma parallel for which calls add_out_to_sub somehow causes the thing to crash
//when all threads operate on the same BLits; these are only iterators anyways, so copying them is cheap
template<class KeyType,typename T>
BlockMatrix<1,KeyType,T> add_out_sub(std::map<uint, std::vector<typename BlockMatrix<1,KeyType,T>::const_iterator> > BLits,
				     const BlockMatrix<1,KeyType,T>& A,const MPOMat<1,T>& mpo,const BlockMatrix<1,KeyType,T>& Ap,
				     const KeyType p,i_to_key<KeyType> itokey,const int direction){
  assert(itokey.size());
  const int idir = (direction > 0 ? 1 : 0);
  BlockMatrix<1,KeyType,T> B;
  const int out_a = (direction > 0 ? 1 : 0);
  typename BlockMatrix<1,KeyType,T>::const_iterator a;
  //sort all matrices with outgoing p in B
  for(a=A.begin();a!=A.end();++a){
    if (a->first[out_a]==p)
      B[a->first] = a->second;
  }
  BlockMatrix<1,KeyType,T> Bout;
  //choose an mpo entry:
  for (typename MPOMat<1,T> ::const_iterator o=mpo.begin();o!=mpo.end();++o) {
    const ukey<2>& key_o = o->first;
    const SparseOperator<1,T>& lo = o->second;
    const uint o_out = key_o[idir], o_in = key_o[1-idir];
    //for that entry, run through the physical operator entries; this defines you the a lower tensorkeytype
    for(typename SparseOperator<1,T>::const_iterator oi=lo.begin();oi!=lo.end();++oi) {
      TensorKeyType<2,1,KeyType> key_b;
      //here the tensor keytype is defined; it's not sure yet if it can be contracted
      KeyType kyout=itokey[oi->first(0)];
      if (direction > 0){
	key_b = TensorKeyType<2,1,KeyType> (p-kyout,p,oi->first.getOut());
      }else {
	key_b = TensorKeyType<2,1,KeyType> (p, p + kyout,oi->first.getOut());
      }
      typename BlockMatrix<1,KeyType,T>::iterator posb = B.find(key_b);//now see if there really is a matrix with that key
      if (posb==B.end()) continue;//of not, then go to the next Operator entry
      // now run through all keys of BLits: so itBL is an iterator of BlockMatrix with outgoing auxmpo-ind o_in
      // and see if the key_b also matches any of the BL's
      typename std::vector<typename BlockMatrix<1,KeyType,T>::const_iterator>::iterator itBL;
      //            cout<<BLits.size()<<endl;
      for(itBL=BLits[o_in].begin();itBL!=BLits[o_in].end();++itBL) {
	const TensorKeyType<2,1,KeyType>& key_itBL = (*itBL)->first;
	if (direction >0 && (key_itBL(0)!=key_o[0]) ) continue;//compares if mpo-aux-ind of BL-block matches the mpo-aux-ind of the mpo
	if (direction <=0 && (key_itBL(0)!=key_o[1]) ) continue;//same as above, for the other direction
	if (direction >0 && (key_itBL[1]!=key_b[0]) ) continue;//look if the lower KeyType of the BL matches the chosen left-keytype of key_b of lower mat
	if (direction <=0 && (key_itBL[1]!=key_b[1]) ) continue;//same as above for other direction
	// finally, it's the right one.
	//now find an upper matrix that matches into the construction
	//the key is now defined by the other indices:
	TensorKeyType<2,1,KeyType> key_a;
	KeyType kyin=itokey[oi->first[0]];
	if (direction > 0)
	  key_a = TensorKeyType<2,1,KeyType> (key_itBL[0], key_itBL[0] + kyin,oi->first.getIn());
	else
	  key_a = TensorKeyType<2,1,KeyType> (key_itBL[0] - kyin, key_itBL[0],oi->first.getIn());
	typename BlockMatrix<1,KeyType,T>::const_iterator posa = Ap.find(key_a);
	if (posa==Ap.end()) continue;

	// ok. we've got everything.
	const matrix<T,column_major>&B  = (*itBL)->second;
	const matrix<T,column_major>&au = posa->second;//upper matrix
	const matrix<T,column_major>&ad = posb->second;//lower matrix
	matrix<T,column_major>Tmp1,Tmp2;
	if (direction > 0) {
	  MatMatProd(T(1.0),B, conj(au),Tmp1, 't', 'n');
	  MatMatProd(T(1.0),Tmp1,ad, Tmp2,'t', 'n');
	} else {
	  MatMatProd(T(1.0),conj(au),B,Tmp1, 'n', 'n');
	  MatMatProd(T(1.0),Tmp1,ad, Tmp2,'n', 't');
	}
	Tmp2 *= (oi->second);
	TensorKeyType<2,1,KeyType> ky(key_a[idir], key_b[idir], key_o[idir]);
	typename BlockMatrix<1,KeyType,T>::iterator ypos = Bout.find(ky);
	if (ypos != Bout.end() )
	  Tmp2 += ypos->second;
	Bout[ky] = Tmp2;
      }
    }
  }
  return Bout;    
}


//this function is for L and R contractions and should be used instead of addToR and addToL
//it contracts a single matrices into network; this is needed in dmrg to build the block-representation of the Hamiltonian.
//arguments: BL is an E-matrix, either of the left or the right type; A is an mps matrix, mpo is an mpo-matrix, itokey is the mapping of local
//holbertspace variables to local quantum numbers; if direction is >0, the routine assumes that BL is a left-block, otherwise it assumes it to be right block type.
template<typename T,typename KeyType>
BlockMatrix<1,KeyType,T> addContractionLayer(const BlockMatrix<1,KeyType,T>&BL,const BlockMatrix<1,KeyType,T>&A,const MPOMat<1,T>&mpo,
					     i_to_key<KeyType> itokey, const int direction,const int maxthreads=omp_get_max_threads()){
  if (!A.size() ) {
    throw mps_exception("Error in addContractionLayer: A is empty.");
  }
  typename BlockMatrix<1,KeyType,T>::const_iterator a;
  typename MPOMat<1,T>::const_iterator o;
  // get all outgoing legs of A.
  const int idir= (direction > 0 ? 1 : 0);
  std::set<KeyType> outgoing_set;//the outgoing keys after contraction
  std::set<uint> ingoing_o_set;//the incoming mpo-aux indices
  for(a=A.begin();a!=A.end();++a)
    outgoing_set.insert(a->first[idir]);
  for(o=mpo.begin();o!=mpo.end();++o)
    ingoing_o_set.insert(o->first[1-idir]);
  BlockMatrix<1,KeyType,T>Bout;    //this will be the resulting new BL's
      
  std::vector<KeyType> outgoing = set2array(outgoing_set);//array is needed for #pragma omp paralell for 
    
  const int n_out = outgoing.size();


  const int _maxthreads = maxthreads;
  const int oldmaxthreads=omp_get_max_threads();
  omp_set_num_threads(maxthreads);
  BlockMatrix<1,KeyType,T> bouts[_maxthreads];
    
  //sort all BL blocks according to their mpo-aux index;BLits maps the aux inds to 
  //a set of BL iterators, all with the same mpo aux ind
  std::map<uint,std::vector<typename BlockMatrix<1,KeyType,T>::const_iterator > >BLits;
  typename BlockMatrix<1,KeyType,T> ::const_iterator BLit = BL.begin();


  for(BLit=BL.begin();BLit!=BL.end();++BLit){
    const uint o_in = BLit->first(0);
    BLits[o_in].push_back(BLit);
  }
    
#pragma omp parallel for
  for(int i_out=0;i_out<n_out;++i_out) {
    const KeyType out = outgoing[i_out];
    const int jthr = omp_get_thread_num();
    BlockMatrix<1,KeyType,T> B_i = add_out_sub_2(BLits,A,mpo,out,itokey,direction);
    typename BlockMatrix<1,KeyType,T> ::iterator it;
    merge(bouts[jthr], B_i);
  }
  omp_set_num_threads(oldmaxthreads);
  for(int i=0;i<_maxthreads;++i)
    merge(Bout, bouts[i]);
  return Bout;
}


//this function is for L and R contractions and should be used instead of addToR and addToL
//it contracts two different matrices (A and Ap) into the network. Naming of the matrices is inconsistent with usual literature, sorry.
//else, it does everything like the above addContractionLayer routine.
template<typename T,typename KeyType>
BlockMatrix<1,KeyType,T> addContractionLayer(const BlockMatrix<1,KeyType,T>&BL,const BlockMatrix<1,KeyType,T>&A,const MPOMat<1,T>&mpo,
 					     const BlockMatrix<1,KeyType,T>&Ap,i_to_key<KeyType> itokey, const int direction,
 					     const int maxthreads=omp_get_max_threads()){
  if (!A.size() ) {
    throw mps_exception("Error in addContractionLayer: A is empty.");
  }
  if (!Ap.size() ) {
    throw mps_exception("Error in addContractionLayer: Ap is empty.");
  }

  typename BlockMatrix<1,KeyType,T>::const_iterator a;
  typename MPOMat<1,T>::const_iterator o;
  // get all outgoing legs of A.
  const int idir= (direction > 0 ? 1 : 0);
  std::set<KeyType> outgoing_set;//the outgoing keys after contraction
  std::set<uint> ingoing_o_set;//the incoming mpo-aux indices
  for(a=A.begin();a!=A.end();++a)
    outgoing_set.insert(a->first[idir]);
  for(o=mpo.begin();o!=mpo.end();++o)
    ingoing_o_set.insert(o->first[1-idir]);
  BlockMatrix<1,KeyType,T>Bout;    //this will be the resulting new BL's
    
  std::vector<KeyType> outgoing = set2array(outgoing_set);//array is needed for #pragma omp paralell for 
  
  const int n_out = outgoing.size();


  const int _maxthreads = maxthreads;
  const int oldmaxthreads=omp_get_max_threads();
  omp_set_num_threads(maxthreads);
  BlockMatrix<1,KeyType,T> bouts[_maxthreads];
  
  //sort all BL blocks according to their mpo-aux index;BLits maps the aux inds to 
  //a set of BL iterators, all with the same mpo aux ind
  std::map<uint,std::vector<typename BlockMatrix<1,KeyType,T>::const_iterator > >BLits;
  typename BlockMatrix<1,KeyType,T> ::const_iterator BLit = BL.begin();


  for(BLit=BL.begin();BLit!=BL.end();++BLit){
    const uint o_in = BLit->first(0);
    BLits[o_in].push_back(BLit);
  }
  
#pragma omp parallel for
  for(int i_out=0;i_out<n_out;++i_out) {
    const KeyType out = outgoing[i_out];
    const int jthr = omp_get_thread_num();
    BlockMatrix<1,KeyType,T> B_i = add_out_sub(BLits,A,mpo,Ap,out,itokey,direction);
    typename BlockMatrix<1,KeyType,T> ::iterator it;
    merge(bouts[jthr], B_i);
  }
  omp_set_num_threads(oldmaxthreads);
  for(int i=0;i<_maxthreads;++i)
    merge(Bout, bouts[i]);
  return Bout;
}



//it contracts different matrices into the network
 // template<typename T,typename KeyType>
 // BlockMatrix<1,KeyType,T> addContractionLayer(const BlockMatrix<1,KeyType,T>&BL,const BlockMatrix<1,KeyType,T>&A,const MPOMat<1,T>&mpo,
 // 					     const BlockMatrix<1,KeyType,T>&B,i_to_key<KeyType> itokey, const int direction,
 // 					     const int maxthreads=omp_get_max_threads()){
 //   if (!A.size() ) {
 //     throw mps_exception("Error in addContractionLayer: A is empty.");
 //   }
 //   if (!B.size() ) {
 //     throw mps_exception("Error in addContractionLayer: B is empty.");
 //   }

 //   std::set<KeyType> key_A;
 //   const uint idir=direction>0?1:0;
 //   std::map<unsigned int,std::set<ukey<2> > > leg_O;
 //   std::map<KeyType,std::vector<typename BlockMatrix<1,KeyType,T>::const_iterator> > leg_A,leg_BL;

 //   //first get all mpo-index-pairs, ordered after outgoing
 //   for(typename MPOMat<1,T>::const_iterator o=mpo.begin();o!=mpo.end();++o) {
 //     leg_O[ (unsigned int)(o->first[idir]) ].insert(o->first);
 //   }
 //   //get all A-matrix iterators belonging to a certain outgoing Key
 //   for(typename BlockMatrix<1,KeyType,T>::const_iterator a=A.begin();a!=A.end();++a) {
 //     key_A.insert( a->first[idir]);
 //     leg_A[a->first[idir]].push_back(a);        
 //   }
 //   //store BL-iterators according to their upper Key
 //   for(typename BlockMatrix<1,KeyType,T>::const_iterator bl=BL.begin();bl!=BL.end();++bl){
 //     leg_BL[bl->first[0]].push_back(bl);
 //   }

 //   std::vector<KeyType> keys_A = set2array(key_A);
 //   const int _maxthreads = maxthreads;
 //   const int oldmaxthreads=omp_get_max_threads();
 //   omp_set_num_threads(_maxthreads);
 //   BlockMatrix<1,KeyType,T> BOUT[_maxthreads];

 // #pragma omp parallel for
 //   for(int k=0;k<keys_A.size();++k) {
 //     const KeyType kyout=keys_A[k];
 //     const int jthr = omp_get_thread_num();
 //     BlockMatrix<1,KeyType,T>& Bout = BOUT[jthr];

 //     std::vector<typename BlockMatrix<1,KeyType,T>::const_iterator> va=leg_A[kyout],vbl;
 //     typename std::vector<typename BlockMatrix<1,KeyType,T>::const_iterator>::iterator ia,ibl;
 //     for(ia=va.begin();ia!=va.end();++ia){
 //       KeyType kyin=(*ia)->first[1-idir];
 //       if(leg_BL.find(kyin)==leg_BL.end())
 // 	continue;	//OK, we got a key which is an mps key which is not in the BL's ...
    
 //       vbl=leg_BL.find(kyin)->second;
 //       for(ibl=vbl.begin();ibl!=vbl.end();++ibl){
 // 	//now the mpo keys...
 // 	for(std::set<ukey<2> >::iterator o=leg_O[(*ibl)->first(0)].begin();o!=leg_O[(*ibl)->first(0)].end();++o){
 // 	  //and now the local operators ...
 // 	  const SparseOperator<1,T> &lo=mpo.find(*o)->second;
 // 	  typename SparseOperator<1,T>::const_iterator ilo;
 // 	  KeyType lkey,rkey;
 // 	  for(ilo=lo.begin();ilo!=lo.end();++ilo){
 // 	    if(ilo->first[0]!=(*ia)->first(0)){
 // 	      continue;
 // 	    }
 // 	    if(direction>0){
 // 	      lkey = (*ibl)->first[1];
 // 	      rkey = (*ibl)->first[1]+itokey[ilo->first(0)];
 // 	    }
 // 	    else if(direction<=0){
 // 	      rkey = (*ibl)->first[1];
 // 	      lkey = (*ibl)->first[1]-itokey[ilo->first(0)];
 // 	    }
 // 	    //now we got everything ...	    
 // 	    TensorKeyType<2,1,KeyType> kyB(lkey,rkey,ilo->first(0));
 // 	    if(B.find(kyB)==B.end())
 // 	      continue;

 // 	    TensorKeyType<2,1,KeyType> kyBout(kyout,direction>0?rkey:lkey,(*o)[idir]);
 // 	    matrix<T,column_major> mA=(*ia)->second;
 // 	    matrix<T,column_major> mB=B.find(kyB)->second;
 // 	    matrix<T,column_major> mBL=(*ibl)->second;
 // 	    matrix<T,column_major> temp;
 // 	    mB*=ilo->second;
 // 	    if(direction>0){
 // 	      MatMatProd(T(1.0),mBL,conj(mA),temp,'t','n');
 // 	      MatMatProd(T(1.0),temp,mB,mA,'t','n');
 // 	    }else if(direction<=0){
 // 	      MatMatProd(T(1.0),conj(mA),mBL,temp,'n','n');
 // 	      MatMatProd(T(1.0),temp,mB,mA,'n','t');
 // 	    }
 // 	    typename BlockMatrix<1,KeyType,T>::iterator ypos=Bout.find(kyBout);
 // 	    if(ypos!=Bout.end() )
 // 	      mA += ypos->second;
 // 	    Bout[kyBout] = mA;
 // 	  }
 // 	}
 //       }
 //     }
  
 //   }//end of pragma
 //   omp_set_num_threads(oldmaxthreads);
 //   BlockMatrix<1,KeyType,T> out;
 //   for(int i=0;i<_maxthreads;++i)
 //     merge(out, BOUT[i]);
 //   return out;
 // }


//    this does contraction with a two site mpo, but it's wrong for the fermi hubbard
// template<class KeyType,typename T>
// void doContraction(BlockMatrix<2,KeyType,T>&Aout,const std::vector<typename BlockMatrix<1,KeyType,T>::const_iterator> &BLits,
//       		   const BlockMatrix<1,KeyType,T>&BR,const BlockMatrix<2,KeyType,T>&Ain,const MPOMat<1,T>&mpol,const MPOMat<1,T>&mpor
//       		   ,std::vector<i_to_key<KeyType> >itokey)
// {
//   const MPOMat<2,T>mpo(mpol,mpor,itokey[0],itokey[1]);
//   typename std::vector<typename BlockMatrix<1,KeyType,T>::const_iterator>::const_iterator itBL;
//   typename MPOMat<2,T> ::const_iterator o;
//   for(itBL=BLits.begin();itBL!=BLits.end();++itBL){
//     for (o=mpo.begin();o!=mpo.end();++o) {
//       const ukey<2>& key_o = o->first;
//       const SparseOperator<2,T>& lo = o->second;
//       const uint o_out = key_o[1], o_in = key_o[0];
//       if(o_in!=(*itBL)->first(0))continue;
//       //for that entry, run through the physical operator entries; this defines you the a lower tensorkeytype
//       //Key of sparseLocalOp is OpKeyType;[] gives back the incoming, () the outgoing index
//       for(typename SparseOperator<2,T>::const_iterator oi=lo.begin();oi!=lo.end();++oi) {
//       	TensorKeyType<2,2,KeyType> key_mat;
//       	//here the tensorkeytype is defined; it's not sure yet if it can be contracted
//    	//get the lower key:
//    	KeyType kylo=itokey[0][oi->first(0)];
//    	for(uint u = 1;u<itokey.size();++u)kylo+=itokey[u][oi->first(u)];

//    	key_mat = TensorKeyType<2,2,KeyType> ((*itBL)->first[1],(*itBL)->first[1]+kylo,oi->first.getOut());
            
//    	typename BlockMatrix<2,KeyType,T>::const_iterator pos = Ain.find(key_mat);//now see if there really is a matrix with that key
//    	if (pos==Ain.end()) continue;//of not, then go to the next Operator entry
//    	// now run through all keys of BLits: so itBL is an iterator of BlockMatrix with outgoing auxmpo-ind o_in
//    	// and see if the key_b also matches any of the BL's
//    	KeyType kyup=itokey[0][oi->first[0]];
//    	for(uint u = 1;u<itokey.size();++u)kyup+=itokey[u][oi->first[u]];

//    	TensorKeyType<2,1,KeyType> kyR((*itBL)->first[0]+kyup,(*itBL)->first[1]+kylo,o_out);
//    	typename  BlockMatrix<1,KeyType,T>::const_iterator posR = BR.find(kyR);//now see if there really is a matrix with that key
//    	if(posR==BR.end())continue;
//    	//now we got all
//    	const matrix<T,column_major>&Lm  = (*itBL)->second;
//    	matrix<T,column_major>a = pos->second;
//    	const matrix<T,column_major>&Rm = posR->second;
//    	matrix<T,column_major>Tmp1,Tmp2;
//    	assert(lo.getSign());
//    	// KeyType kyleft= (*itBL)->first[1];
//    	// for(uint ind = 0;ind<kyleft.size();++ind)
//    	//   ql += kyleft(ind);
//    	a*=((oi->second));
//    	TensorKeyType<2,2,KeyType> key_out((*itBL)->first[0],posR->first[0],oi->first.getIn());
//    	MatMatProd(Lm,a,Tmp1,'n','n');
//    	MatMatProd(Tmp1,Rm,Tmp2,'n','t');
//    	typename  BlockMatrix<2,KeyType,T>::iterator ypos = Aout.find(key_out);
//    	if (ypos != Aout.end() )
//    	  Tmp2 += ypos->second;
//    	Aout[key_out] = Tmp2;
//       }
//     } 
//   }
// }




template<class KeyType,typename T>
void doContraction(BlockMatrix<2,KeyType,T>&Aout,const std::vector<typename BlockMatrix<1,KeyType,T>::const_iterator> &BLits,
		   const BlockMatrix<1,KeyType,T>&BR,const BlockMatrix<2,KeyType,T>&Ain,const MPOMat<1,T> &mpol,const MPOMat<1,T> &mpor,
		   std::vector<i_to_key<KeyType> >itokey,const bool debug=false){
	    
  std::vector<uint> out(2);
  std::vector<uint> in(2);
  KeyType kylo,kyup;
  typename BlockMatrix<2,KeyType,T>::const_iterator pos;
  typename  BlockMatrix<2,KeyType,T>::iterator ypos;
  typename BlockMatrix<1,KeyType,T>::const_iterator posR;
  typename std::vector<typename BlockMatrix<1,KeyType,T>::const_iterator>::const_iterator itBL;
  typename MPOMat<1,T> ::const_iterator o1,o2;
  for(itBL=BLits.begin();itBL!=BLits.end();++itBL){
    for (o1=mpol.begin();o1!=mpol.end();++o1){
      for (o2=mpor.begin();o2!=mpor.end();++o2) {
	if(o1->first[0]!=(*itBL)->first(0))continue;
	if(o2->first[0]!=o1->first[1])continue;
	const ukey<2>& key_r = o2->first,key_l = o1->first;
	//const SparseOperator<1,T>& lol = o1->second,lor = o2->second;
	const uint o_out = key_r[1], o_in = key_l[0];
	//for that entry, run through the physical operator entries; this defines you the a lower tensorkeytype
	//Key of sparseLocalOp is OpKeyType;[] gives back the incoming, () the outgoing index
	for(typename SparseOperator<1,T>::const_iterator oil=o1->second.begin();oil!=o1->second.end();++oil) {
	  for(typename SparseOperator<1,T>::const_iterator oir=o2->second.begin();oir!=o2->second.end();++oir) {
	    TensorKeyType<2,2,KeyType> key_mat;
	    //here the tensorkeytype is defined; it's not sure yet if it can be contracted
	    //get the lower key:
	    //std::vector<uint> out(2);
	    //std::vector<uint> in(2);
	    //KeyType kylo=itokey[0][oil->first(0)];
	    kylo.reset();
	    kylo=itokey[0][oil->first(0)];
	    kylo+=itokey[1][oir->first(0)];
	    out[0]= oil->first(0);out[1] = oir->first(0);
	    in[0]= oil->first[0];in[1] = oir->first[0];
	    //bool indexcheck=(in[0]==11&&in[1]==8);
	    key_mat = TensorKeyType<2,2,KeyType> ((*itBL)->first[1],(*itBL)->first[1]+kylo,out);
	    pos = Ain.find(key_mat);//now see if there really is a matrix with that key
	    if (pos==Ain.end()){
	      continue;//of not, then go to the next Operator entry
	    }

	    // now run through all keys of BLits: so itBL is an iterator of BlockMatrix with outgoing auxmpo-ind o_in
	    // and see if the key_b also matches any of the BL's
	    //KeyType kyup=itokey[0][oil->first[0]];
	    kyup.reset();
	    kyup=itokey[0][oil->first[0]];
	    kyup+=itokey[1][oir->first[0]];
	    //	    for(uint u = 1;u<itokey.size();++u)

	    TensorKeyType<2,1,KeyType> kyR((*itBL)->first[0]+kyup,(*itBL)->first[1]+kylo,o_out);
	    posR = BR.find(kyR);//now see if there really is a matrix with that key
	    if(posR==BR.end()){
	      continue;
	    }
	    
	    //now we got all
	    //const matrix<T,column_major>&Lm  = (*itBL)->second;
	    //matrix<T,column_major>&a = pos->second;
	    //const matrix<T,column_major>&Rm = posR->second;
	    matrix<T,column_major>Tmp1,Tmp2;
	    T num =((oil->second)*(oir->second));
	    // for(uint k=0;k<a.size1()*a.size2();++k)
	    //   a.data()[k]*=num;
	    //a*=((oil->second)*(oir->second));
	    TensorKeyType<2,2,KeyType> key_out((*itBL)->first[0],posR->first[0],in);
	    MatMatProd(num,(*itBL)->second,pos->second,Tmp1,'n','n');
	    MatMatProd(T(1.0),Tmp1,posR->second,Tmp2,'n','t');
	    ypos = Aout.find(key_out);
	    if (ypos != Aout.end() )
	      Tmp2 += ypos->second;
	    Aout[key_out] = Tmp2;
	  }
	}
      } 
    }
  }
}


//does the single site contraction
template<class KeyType,typename T>
void doContraction(BlockMatrix<1,KeyType,T>&Aout,const std::vector<typename BlockMatrix<1,KeyType,T>::const_iterator>  &BLits,
		   const BlockMatrix<1,KeyType,T>&BR,const BlockMatrix<1,KeyType,T>&Ain,const MPOMat<1,T> &mpo,i_to_key<KeyType>itokey)
{
  typename std::vector<typename BlockMatrix<1,KeyType,T>::const_iterator>::const_iterator itBL;
  typename MPOMat<1,T> ::const_iterator o;
  for(itBL=BLits.begin();itBL!=BLits.end();++itBL){
    for (o=mpo.begin();o!=mpo.end();++o){
      const ukey<2>& key_o = o->first;
      const SparseOperator<1,T>& lo = o->second;
      const uint o_out = key_o[1], o_in = key_o[0];
      if(o_in!=(*itBL)->first(0))continue;
      //for that entry, run through the physical operator entries; this defines you the a lower tensorkeytype
      //Key of sparseLocalOp is OpKeyType;[] gives back the incoming, () the outgoing index
      for(typename SparseOperator<1,T>::const_iterator oi=lo.begin();oi!=lo.end();++oi){
	TensorKeyType<2,1,KeyType> key_mat;
	//here the tensorkeytype is defined; it's not sure yet if it can be contracted
	//get the lower key:
	KeyType kylo=itokey[oi->first(0)];
	key_mat = TensorKeyType<2,1,KeyType> ((*itBL)->first[1],(*itBL)->first[1]+kylo,oi->first(0));
            
	typename BlockMatrix<1,KeyType,T>::const_iterator pos = Ain.find(key_mat);//now see if there really is a matrix with that key
	if (pos==Ain.end()) continue;//if not, then go to the next Operator entry
	// now run through all keys of BLits: so itBL is an iterator of BlockMatrix with outgoing auxmpo-ind o_in
	// and see if the key_b also matches any of the BL's
	KeyType kyup=itokey[oi->first[0]];
	TensorKeyType<2,1,KeyType> kyR((*itBL)->first[0]+kyup,(*itBL)->first[1]+kylo,o_out);
	typename  BlockMatrix<1,KeyType,T>::const_iterator posR = BR.find(kyR);//now see if there really is a matrix with that key
	if(posR==BR.end())continue;
	//now we got all
	const matrix<T,column_major>&Lm  = (*itBL)->second;
	matrix<T,column_major>a = pos->second;
	const matrix<T,column_major>&Rm = posR->second;
	matrix<T,column_major>Tmp1,Tmp2;
	T num=(oi->second);
	//a*=(oi->second);
	TensorKeyType<2,1,KeyType> key_out((*itBL)->first[0],posR->first[0],oi->first[0]);
	MatMatProd(num,Lm,a,Tmp1,'n','n');
	MatMatProd(T(1.),Tmp1,Rm,Tmp2,'n','t');
	typename  BlockMatrix<1,KeyType,T>::iterator ypos = Aout.find(key_out);
	if (ypos != Aout.end() )
	  Tmp2 += ypos->second;
	Aout[key_out] = Tmp2;
      }
    }
  }
}

//same as contractTwosite, but with a single-site mps wave-function.
template<class KeyType,typename T>
void contractSingle(BlockMatrix<1,KeyType,T>&Aout ,const BlockMatrix<1,KeyType,T>&BL,const BlockMatrix<1,KeyType,T> &A,const MPOMat<1,T>&mpo,
					const BlockMatrix<1,KeyType,T> &BR,i_to_key<KeyType> itokey,const int maxthreads=omp_get_max_threads()){
  Aout.clear();
  std::set<KeyType>pleft;
  std::map<KeyType,std::vector<typename BlockMatrix<1,KeyType,T>::const_iterator> >BLits;
  typename BlockMatrix<1,KeyType,T>::const_iterator l;
  //Store the left blocks according to the left key
  for(l=BL.begin();l!=BL.end();++l){
    BLits[l->first[0]].push_back(l);
    pleft.insert(l->first[0]);
  }
  std::vector<KeyType>lefts =set2array(pleft);
  const uint n = lefts.size();
  const uint max_threads=maxthreads;
  const int oldmaxthreads=omp_get_max_threads();
  omp_set_num_threads(max_threads);
  BlockMatrix<1,KeyType,T> aout[max_threads];
#pragma omp parallel for
  for(uint i = 0;i<n;++i){
    KeyType p = lefts[i];
    const uint jthread=omp_get_thread_num();
    BlockMatrix<1,KeyType,T>&_aout_=aout[jthread];
    doContraction(_aout_,BLits[p],BR,A,mpo,itokey);
  }
  omp_set_num_threads(oldmaxthreads);
  for(uint i=0;i<max_threads;++i)
    merge(Aout,aout[i]);
}


//this is the primitive version which ...
template<class KeyType,typename T>
void contractTwosite2(BlockMatrix<2,KeyType,T>&Aout,const BlockMatrix<1,KeyType,T>&L,const BlockMatrix<1,KeyType,T>&R,const BlockMatrix<2,KeyType,T>&Ain,
		const MPOMat<2,T>&mpo,std::vector<i_to_key<KeyType> >itokey){
  Aout.clear();
  typename MPOMat<2,T>::const_iterator itM;
  typename BlockMatrix<1,KeyType,T>::const_iterator itL,itR;
  typename BlockMatrix<2,KeyType,T>::const_iterator itAin;
  typename BlockMatrix<2,KeyType,T>::iterator itAout;
  UintVectorType ql,qr,ql_prime,qr_prime;
  boost::numeric::ublas::matrix<T,boost::numeric::ublas::column_major>tmp,tmp1,tmp2;
  for (itL = L.begin(); itL!= L.end(); itL++)
    {
      KeyType ql = itL->first[0];
      KeyType ql_prime = itL->first[1];
      for(itM = mpo.begin();itM!=mpo.end();itM++)
	{
	  if(itM->first(0) == itL->first(0))
	    {
	      const SparseOperator<2,T>& lo = itM->second;	      
	      typename SparseOperator<2,T>::const_iterator li;
	      for(li = lo.begin();li!=lo.end();++li){
		KeyType qr = ql+itokey[0][li->first[0]]+itokey[1][li->first[1]];
		KeyType qr_prime =ql_prime+itokey[0][li->first(0)]+itokey[1][li->first(1)];
		TensorKeyType<2,1,KeyType>kr(qr,qr_prime,itM->first(1));
		itR = R.find(kr);
		TensorKeyType<2,2,KeyType>kA(ql_prime,qr_prime,li->first.getOut());
		itAin = Ain.find(kA);
		if(itAin!=Ain.end()&&itR!=R.end())
		  {
		    MatMatProd(T(1.),itL->second,(itAin->second),tmp1,'n','n');
		    MatMatProd(T(1.),tmp1,itR->second,tmp2,'n','t');
		    tmp2=tmp2*li->second;
		    TensorKeyType<2,2,KeyType>kAout(ql,qr,li->first.getIn());
		    itAout = Aout.find(kAout);
		    //if no entry is found then insert it
		    if(itAout==Aout.end())
		      {
			Aout.insert(pair<TensorKeyType<2,2,KeyType>,boost::numeric::ublas::matrix<T,column_major> >(kAout,tmp2));
		      }
		    //if there exist an entry for this key, then add tmp2 to it
		    else
		      {
			tmp1 = itAout->second;
			itAout->second = tmp1+tmp2;
		      }
		  }
	      }
	    }
	}
    }
}

//this should work for single and two side product, but is not
//templated yet. don't forget that itokey maps a std::vector<uint>
//(=local indicces) to KeyType; there is no conjugation of the vector
//involved!  does a two-site contraction of a left E-matrix BL, a
//right E-matrix BR, a two-site mps wave function A and two mpo
//matrices mpol and mpor and stores the resulting two-site mps
//wavefunction in Aout this routine is needed in the lanzcos algorithm
//in lanzcos.hpp
template<class KeyType,typename T>
void contractTwosite(BlockMatrix<2,KeyType,T>&Aout,const BlockMatrix<1,KeyType,T>&BL,const BlockMatrix<1,KeyType,T>&BR,
		     const BlockMatrix<2,KeyType,T>&A,const MPOMat<1,T>&mpol,const MPOMat<1,T>&mpor,std::vector<i_to_key<KeyType> >itokey,
		     const int maxthreads=omp_get_max_threads(),const bool debug=false){
  Aout.clear();
  std::set<KeyType>pleft;
  std::map<KeyType,std::vector<typename BlockMatrix<1,KeyType,T>::const_iterator> >BLits;
  typename BlockMatrix<1,KeyType,T>::const_iterator l;
  //Store the left blocks according to the left key
  for(l=BL.begin();l!=BL.end();++l){
    BLits[l->first[0]].push_back(l);
    pleft.insert(l->first[0]);
  }

  std::vector<KeyType>lefts =set2array(pleft);
  const uint n = lefts.size();
  const uint max_threads = maxthreads;
  const int oldmaxthreads=omp_get_max_threads();
  omp_set_num_threads(max_threads);
  BlockMatrix<2,KeyType,T> aout[max_threads];
#pragma omp parallel for
  for(uint i = 0;i<n;++i){
    KeyType p = lefts[i];
    const uint jthread=omp_get_thread_num();
    BlockMatrix<2,KeyType,T>&_aout_=aout[jthread];

    doContraction(_aout_,BLits[p],BR,A,mpol,mpor,itokey,false);
  }
  omp_set_num_threads(oldmaxthreads);
  for(uint i=0;i<max_threads;++i)
    merge(Aout,aout[i]);

}

//this should work for single and two side product if doContraction was adapted ...
// template<int Ns,class KeyType,typename T>
// void contract(BlockMatrix<Ns,KeyType,T>&Aout,const BlockMatrix<1,KeyType,T>&BL,const BlockMatrix<1,KeyType,T>&BR,const BlockMatrix<Ns,KeyType,T>&A,
// 	      const std::vector<MPOMat<1,T> >&mpos,std::vector<i_to_key<KeyType> >itokey)
// {
//   Aout.clear();
//   std::set<KeyType>pleft;
//   std::map<KeyType,std::vector<typename BlockMatrix<1,KeyType,T>::const_iterator> >BLits;
//   typename BlockMatrix<1,KeyType,T>::const_iterator l;
//   //Store the left blocks according to the left key
//   for(l=BL.begin();l!=BL.end();++l){
//     BLits[l->first[0]].push_back(l);
//     pleft.insert(l->first[0]);
//   }
//   std::vector<KeyType>lefts =set2array(pleft);
//   const uint n = lefts.size();
//   const uint max_threads = omp_get_max_threads();
//   BlockMatrix<Ns,KeyType,T> aout[max_threads];
// #pragma omp parallel for
//   for(uint i = 0;i<n;++i){
//     KeyType p = lefts[i];
//     const uint jthread=omp_get_thread_num();
//     BlockMatrix<Ns,KeyType,T>&_aout_=aout[jthread];
//     doContraction(_aout_,BLits[p],BR,A,mpos,itokey);
//   }
//   for(uint i=0;i<max_threads;++i)
//     merge(Aout,aout[i]);
// }




template<class KeyType,typename T>
triple<T,Real,uint> doTwoSiteLanzcos(const BlockMatrix<1,KeyType,T>&BL,BlockMatrix<1,KeyType,T> &mpsl,BlockMatrix<1,KeyType,T> &mpsr,
				     const BlockMatrix<1,KeyType,T>&BR,const MPOMat<1,T>&mpol,const MPOMat<1,T>&mpor,DiagBlockMatrix<KeyType,Real>&lambda,
				     i_to_key<KeyType> itokeyl,i_to_key<KeyType> itokeyr,const int direction,Real delta=1e-12,bool testenergy=true,
				     const uint chimax=100,const bool doTrunc=false,const Real tr_weight = 1e-12)
{
  BlockMatrix<2,KeyType,T> Ainit(mpsl,mpsr);
  mpsl.clear();mpsr.clear();
  std::vector<BlockMatrix<1,KeyType,T> > mps(2);
  std::vector<DiagBlockMatrix<KeyType,Real> > lam(1);
  std::vector<i_to_key<KeyType> >itokey;
  itokey.push_back(itokeyl);
  itokey.push_back(itokeyr);
  triple<T,Real,uint>out=doLanzcos(BL,Ainit,BR,mpol,mpor,mps,lam,itokey,direction,delta,testenergy,500,chimax,doTrunc,tr_weight);
  mpsl = mps[0];
  mpsr=mps[1];
  lambda = lam[0];
  return out;
  
}


//this is a two-site lanzcos procedure for dmrg ...
template<class KeyType,typename T>
triple<T,Real,uint> doLanzcos(const BlockMatrix<1,KeyType,T>&BL,BlockMatrix<2,KeyType,T>&A1,const BlockMatrix<1,KeyType,T>&BR,
			      const MPOMat<1,T>&mpol,const MPOMat<1,T>&mpor,std::vector<BlockMatrix<1,KeyType,T> >&mps,
			      std::vector<DiagBlockMatrix<KeyType,Real> >&lambda,std::vector<i_to_key<KeyType> >itokey,const int direction,
			      const Real delta,const bool testenergy,const int itmax=500,const uint chimax=100,const bool doTrunc=false,
			      const Real tr_weight=1e-12,const bool normalize=true){
  mps.clear();
  lambda.clear();
  //bring everything into right form:
  int it = 0,Ndiag = 5;
  //clock_t P0,P1;
  T k = 0,e,eold = T(300.0);
  double *eig = NULL;
  bool converged = false, first = true;
  std::vector<BlockMatrix<2,KeyType,T> >As;
  vector<T> es,ks,estmp,kstmp;
  BlockMatrix<2,KeyType,T>A0,A2,A3;
  es.clear();
  estmp.clear();
  ks.clear();
  kstmp.clear();
  Real resi=1e5;
  A0.clear();A2.clear();A3.clear();
  while(!converged){
    k=(sqrt(A1*A1));
    if(abs(k)<delta){converged = true;break;}
    if(!first){ks.push_back(k);}
    A1/=k;
    contractTwosite(A2,BL,BR,A1,mpol,mpor,itokey);
    As.push_back(A1);
    e = A2*A1;
    es.push_back(e);
    if(it%Ndiag==0&&it!=0){
      if(testenergy){
	kstmp = ks;
	estmp = es;
	tridiag(&estmp[0],estmp.size(),&kstmp[0],eig,'N');
	cout.precision(18);
	if(abs(estmp[0]-eold)<delta){cout<<"E="<<estmp[0]<<" after "<<it<<" iterations; ";converged = true;break;}
	eold = estmp[0];
      }
      else{
	kstmp = ks;
	estmp = es;
	BlockMatrix<2,KeyType,T>A4,A5,A6;
	double*eig2 = new double[estmp.size()*estmp.size()];
	tridiag(&estmp[0],estmp.size(),&kstmp[0],eig2,'V');
	A4 = As[0]*T(eig2[0]);
	for(uint i = 1;i<estmp.size();i++){
	  A5 = (A4+(As[i]*T(eig2[i])));
	  A4 = A5;
	}
	delete[] eig2;
	//now test for convergence;
	contractTwosite(A5,BL,BR,A4,mpol,mpor,itokey);
	A6=A5+A4*(-estmp[0]);
	resi = abs(sqrt(A6*A6));
	if(abs(resi)<delta){cout<<"state converged at "<<estmp[0]<<" after "<<it<<" iterations at a residuum of "<<resi<<"; ";converged = true;break;}
      }
    }
    if(!first){
      A3 = (((A1*(-e))+=A2)+=(A0*=(-k)));
      A0 = A1;
      A1=A3;
    }
    if(first){
      A0 = A1;
      A3= ((A1*(-e))+A2);
      A1=A3;
      first = false;
    }
    it++;
    if(it>=itmax){
      if(testenergy){cout<<"doTwoSiteLanzcos:: no convergence after "<<it<<" iteration steps; exit Lanzcos loop at energy "<<estmp[0]<<endl;
	converged = true;break;}
      else{cout<<"doTwoSiteLanzcos:: no convergence after "<<it<<" iteration steps; exit Lanzcos loop at energy "<<estmp[0]
	       <<" and a residuum of "<<resi<<endl;converged = true;break;}
    }
  }
  if(converged&&es.size()!=0){
    double*eig2 = new double[es.size()*es.size()];
    tridiag(&es[0],es.size(),&ks[0],eig2,'V');
    A3 = As[0]*T(eig2[0]);
    for(uint i = 1;i<es.size();i++){
      A2.clear();
      A2 = (A3+(As[i]*T(eig2[i])));
      A3 = A2;
    }
    delete[] eig2;
  }
  if(abs(k)<delta&&es.size()!=0){cout<<"k < "<<delta<<": E="<<es[0]<<"; ";}
  A1 = A3;
  //A1.squeeze(1e-16);//split has problems when there are too many empty blocks
  pair<Real,uint> sret=split(A1,mps,lambda,itokey,direction,chimax,doTrunc,normalize,tr_weight);
  return triple<T,Real,uint>(es[0],sret.first,sret.second);
}



//that's to initialize the recursive template;
//splits up a two-site matrix A and returns a vector containing two single-site blockmatrices, and a vector containing the
//lambda belonging to the splitting. Note that "direction" is there because in the recursive template call, you need the same
//arguments for the call. "direction" has no effect when calling this function. Formerly, it multiplied "lambda" onto the left (direction <0)
//or the right (direction >0) matrix. This is NOT done anymore!!!!! The effect is now the same as direction=0 was earlier.
template<class KeyType,typename T>
pair<std::vector<Real>,std::vector<uint> > split(BlockMatrix<2,KeyType,T>&A,std::vector<BlockMatrix<1,KeyType,T> >&mps,
						 std::vector<DiagBlockMatrix<KeyType,Real> >&lambda,std::vector<i_to_key<KeyType> > itokey,
						 const int direction,uint chimax,bool doTrunc,const Real tr_weight,const bool normalize=true,
						 const int maxthreads=omp_get_max_threads()){
  assert(itokey.size()==2);
  typename BlockMatrix<2,KeyType,T>::iterator a;
  //if it's called recursively, then clearing is bad ...
  //lambda.clear();
  //mps.clear();
  //store all blocks according to their qc number
  std::map<KeyType, std::vector<typename BlockMatrix<2,KeyType,T>::iterator> >blocks;
  std::set<KeyType> QC;
  //the ranges: given a qc, all blocks are sorted into a boost matrix
  //range_map maps the template parameter to a boost::range
  std::map<KeyType,range_map<std::vector<uint> > >RY;
  std::map<KeyType,range_map<uint> >RX;
  typename std::map<KeyType,range_map<std::vector<uint> > >::iterator iry;
  typename std::map<KeyType,range_map<uint > >::iterator irx;
  for(a=A.begin();a!=A.end();++a){
    const uint xp=a->first(0);
    std::vector<uint> yp(2-1);
    for(uint s = 1;s<(2);s++)yp[s-1] = a->first(s);
    KeyType qc = a->first[0]+itokey[0][xp];
    blocks[qc].push_back(a);
    RX[qc][xp]=Range(0,a->second.size1());
    RY[qc][yp]=Range(0,a->second.size2());
    QC.insert(qc);
  }
  //every qc matrix got its own partition:
  map_container<KeyType,uint>sizex,sizey;
  for(irx=RX.begin();irx!=RX.end();++irx){
    uint start = 0,size=0;
    range_map<uint>::iterator i1;
    for(i1 = irx->second.begin();i1!=irx->second.end();++i1){
      size = i1->second.size();
      i1->second = Range(start,start+size);
      start +=size;
    }
    sizex[irx->first] = start;
  }
  for(iry=RY.begin();iry!=RY.end();++iry){
    uint start = 0,size=0;
    range_map<std::vector<uint> >::iterator i2;
    for(i2 = iry->second.begin();i2!=iry->second.end();++i2){
      size = i2->second.size();
      i2->second = Range(start,start+size);
      start +=size;
    }
    sizey[iry->first] = start;
  }
  const uint max_threads = maxthreads;
  const int oldmaxthreads=omp_get_max_threads();
  omp_set_num_threads(max_threads);
  std::vector<KeyType> qcv = set2array(QC);
  uint N = qcv.size();
  BlockMatrix<1,KeyType,T> PMPS[max_threads];
  BlockMatrix<2-1,KeyType,T> RMPS[max_threads];
  DiagBlockMatrix<KeyType,Real>LAMS[max_threads];      

#pragma omp parallel for
  for(uint i = 0;i<N;++i){
    const uint jthread = omp_get_thread_num();
    BlockMatrix<1,KeyType,T>&pmps = PMPS[jthread];
    BlockMatrix<2-1,KeyType,T>&rmps = RMPS[jthread];
    DiagBlockMatrix<KeyType,Real>&lam = LAMS[jthread];      
    std::vector<typename BlockMatrix<2,KeyType,T>::iterator> qcB;
    qcB=blocks[qcv[i]];
    matrix<T,column_major>tmp(sizex[qcv[i]],sizey[qcv[i]]),U,Vd;
    VectorType D;
    boost_setTo(tmp,T(0.0));
    typename std::vector<typename BlockMatrix<2,KeyType,T>::iterator>::iterator b;
    for(b=qcB.begin();b!=qcB.end();++b){
      const uint xp=(*b)->first(0);
      std::vector<uint> yp(2-1);
      for(uint s = 1;s<2;++s)yp[s-1] = (*b)->first(s);
      Range rx = RX[qcv[i]][xp];
      Range ry = RY[qcv[i]][yp];
      matrix_range<matrix<T,column_major> >(tmp,rx,ry) = (*b)->second;
    }
    BoostSVD(tmp, U, D, Vd);
    uint chi = D.size();
    range_map<uint>::iterator i1;
    for(i1 = RX[qcv[i]].begin();i1 != RX[qcv[i]].end();++i1){
      TensorKeyType<2,1,KeyType> ky(qcv[i]-itokey[0][i1->first],qcv[i],i1->first);
      pmps[ky]=matrix_range<matrix<T,column_major> >(U,i1->second,Range(0,chi));
    }
    lam[qcv[i]]=D;
    range_map<std::vector<uint> >::iterator i2;
    assert(2-1>0);
    for(i2 = RY[qcv[i]].begin();i2 != RY[qcv[i]].end();++i2){
      std::vector<uint> yp(2-1);//without the first site; it has been split off
      KeyType kyte = itokey[1][i2->first[0]];
      for(uint s = 2;s<(2-1);++s)kyte+= itokey[s][i2->first[s-1]];
      TensorKeyType<2,2-1,KeyType> ky(qcv[i],qcv[i]+kyte,i2->first);
      rmps[ky]=matrix_range<matrix<T,column_major> >(Vd,Range(0,chi),i2->second);
    }
  }
  omp_set_num_threads(oldmaxthreads);
  //now clear A; this is save (it's been split up)
  A.clear();
  BlockMatrix<2-1,KeyType,T>A1;
  BlockMatrix<1,KeyType,T>mps0;
  DiagBlockMatrix<KeyType,Real>lambda0;
  for(uint k =0;k<max_threads;++k){
    merge(A1,RMPS[k]);
    merge(mps0,PMPS[k]);
    merge(lambda0,LAMS[k]);
  }

  pair<std::vector<Real>,std::vector<uint> > out;
  if(!chimax&&doTrunc){cout<<"no chimax specified: setting chimax to 100"<<endl;chimax=100;}
  if((!chimax)&&(!doTrunc)){cout<<"no truncation made at all! this will result in a memory crash!!!; setting chimax = 200;"<<endl;chimax =200;}

  if(chimax&&(!doTrunc)){
    pair<Real,uint> t=lambda0.truncate(chimax,normalize);
    out.first.push_back(t.first);
    out.second.push_back(t.second);
  }

  if(chimax&&doTrunc){
    pair<Real,uint> t=lambda0.truncate(tr_weight,chimax,normalize);
    out.first.push_back(t.first);
    out.second.push_back(t.second);
  }


  //cout<<out.first[0]<<endl;
  
  dropDangling(mps0,lambda0,A1);
  lambda.push_back(lambda0);
  mps.push_back(mps0);
  mps.push_back(A1);
  

  return out;
}



//An Ns-site BM is split-up into Ns single-site BM's, stored in "mps". the lambdas are normalized, everything is truncated (if specified)
//an important remark: if you want this function go give back true schmidt-coefficients, then the passed Ns-site BlockMatrix HAS to be of the 
//form Lambda Gamma ... Gamma Lambda!, e.g the boundary gammas have to be included! Furthermore, direction has to be >0, so that after splitting
//off a single site, the lambda at the corresponding bond is multiplied to the other (yet-to-decompose) part! In the end, to get the Gamma's back 
//you HAVE to divide the lambdas back out for EVERY single entry in the vector "mps", depending on direction. If direction > 0, then "mps" 
//contains matrices of type lambda*Gamma, except for the last entry in "mps", that's a lambda*Gamma*lambda type.
//if direction<=0, nonsense will be returned.
//the routine has been tested and seems to work

template<int Ns,class KeyType,typename T>
pair<std::vector<Real>,std::vector<uint> >  split(BlockMatrix<Ns,KeyType,T>&A,std::vector<BlockMatrix<1,KeyType,T> >&mps,
						  std::vector<DiagBlockMatrix<KeyType,Real> >&lambda,std::vector<i_to_key<KeyType> > itokey,
						  const int direction,uint chimax,bool doTrunc,const Real tr_weight,const bool normalize=true,
						  const int maxthreads=omp_get_max_threads()){

  assert(itokey.size()==Ns);
  typename BlockMatrix<Ns,KeyType,T>::iterator a;
  //if it's called recurseively, then clearing is bad ...
  //lambda.clear();
  //mps.clear();
  //store all blocks according to their qc number
  std::map<KeyType, std::vector<typename BlockMatrix<Ns,KeyType,T>::iterator> >blocks;
  std::set<KeyType> QC;
  //the ranges: given a qc, all blocks are sorted into a boost matrix
  //range_map maps the template parameter to a boost::range
  std::map<KeyType,range_map<std::vector<uint> > >RY;
  std::map<KeyType,range_map<uint> >RX;
  typename std::map<KeyType,range_map<std::vector<uint> > >::iterator iry;
  typename std::map<KeyType,range_map<uint > >::iterator irx;
  for(a=A.begin();a!=A.end();++a){
    const uint xp=a->first(0);
    std::vector<uint> yp(Ns-1);
    for(uint s = 1;s<(Ns);s++)yp[s-1] = a->first(s);
    KeyType qc = a->first[0]+itokey[0][xp];
    blocks[qc].push_back(a);
    RX[qc][xp]=Range(0,a->second.size1());
    RY[qc][yp]=Range(0,a->second.size2());
    QC.insert(qc);
  }
  //every qc matrix got its own partition:
  map_container<KeyType,uint>sizex,sizey;
  for(irx=RX.begin();irx!=RX.end();++irx){
    uint start = 0,size=0;
    range_map<uint>::iterator i1;
    for(i1 = irx->second.begin();i1!=irx->second.end();++i1){
      size = i1->second.size();
      i1->second = Range(start,start+size);
      start +=size;
    }
    sizex[irx->first] = start;
  }
  for(iry=RY.begin();iry!=RY.end();++iry){
    uint start = 0,size=0;
    range_map<std::vector<uint> >::iterator i2;
    for(i2 = iry->second.begin();i2!=iry->second.end();++i2){
      size = i2->second.size();
      i2->second = Range(start,start+size);
      start +=size;
    }
    sizey[iry->first] = start;
  }
  const uint max_threads = maxthreads;
  const int oldmaxthreads=omp_get_max_threads();
  omp_set_num_threads(max_threads);
  std::vector<KeyType> qcv = set2array(QC);
  uint N = qcv.size();
  BlockMatrix<1,KeyType,T> PMPS[max_threads];
  BlockMatrix<Ns-1,KeyType,T> RMPS[max_threads];
  DiagBlockMatrix<KeyType,Real>LAMS[max_threads];      

#pragma omp parallel for
  for(uint i = 0;i<N;++i){
    const uint jthread = omp_get_thread_num();
    BlockMatrix<1,KeyType,T>&pmps = PMPS[jthread];
    BlockMatrix<Ns-1,KeyType,T>&rmps = RMPS[jthread];
    DiagBlockMatrix<KeyType,Real>&lam = LAMS[jthread];      
    std::vector<typename BlockMatrix<Ns,KeyType,T>::iterator> qcB;
    qcB=blocks[qcv[i]];
    matrix<T,column_major>tmp(sizex[qcv[i]],sizey[qcv[i]]),U,Vd;
    VectorType D;
    boost_setTo(tmp,T(0.0));
    typename std::vector<typename BlockMatrix<Ns,KeyType,T>::iterator>::iterator b;
    for(b=qcB.begin();b!=qcB.end();++b){
      const uint xp=(*b)->first(0);
      std::vector<uint> yp(Ns-1);
      for(uint s = 1;s<Ns;++s)yp[s-1] = (*b)->first(s);
      Range rx = RX[qcv[i]][xp];
      Range ry = RY[qcv[i]][yp];
      matrix_range<matrix<T,column_major> >(tmp,rx,ry) = (*b)->second;
    }
    BoostSVD(tmp, U, D, Vd);
    uint chi = D.size();
    range_map<uint>::iterator i1;
    for(i1 = RX[qcv[i]].begin();i1 != RX[qcv[i]].end();++i1){
      TensorKeyType<2,1,KeyType> ky(qcv[i]-itokey[0][i1->first],qcv[i],i1->first);
      pmps[ky]=matrix_range<matrix<T,column_major> >(U,i1->second,Range(0,chi));
    }
    lam[qcv[i]]=D;
    range_map<std::vector<uint> >::iterator i2;
    assert(Ns-1>0);
    for(i2 = RY[qcv[i]].begin();i2 != RY[qcv[i]].end();++i2){
      std::vector<uint> yp(Ns-1);//without the first site; it has been split off
      KeyType kyte = itokey[1][i2->first[0]];
      for(uint s = 2;s<(Ns-1);++s)kyte+= itokey[s][i2->first[s-1]];
      TensorKeyType<2,Ns-1,KeyType> ky(qcv[i],qcv[i]+kyte,i2->first);
      rmps[ky]=matrix_range<matrix<T,column_major> >(Vd,Range(0,chi),i2->second);
    }
  }
  omp_set_num_threads(oldmaxthreads);
  //now clear A; this is save (it's been split up)
  A.clear();
  BlockMatrix<Ns-1,KeyType,T>A1;
  BlockMatrix<1,KeyType,T>mps0;
  DiagBlockMatrix<KeyType,Real>lambda0;
  for(uint k =0;k<max_threads;++k){
    merge(A1,RMPS[k]);
    merge(mps0,PMPS[k]);
    merge(lambda0,LAMS[k]);
  }


  pair<std::vector<Real>,std::vector<uint> > out,tout;
  if(!chimax&&doTrunc){cout<<"no chimax specified: setting chimax to 100"<<endl;chimax=100;}
  if((!chimax)&&(!doTrunc)){cout<<"no truncation made at all! this will result in a memory crash!!!; setting chimax = 200;"<<endl;chimax =200;}
  if(chimax&&(!doTrunc)){
    pair<Real,uint> t=lambda0.truncate(chimax,normalize);
    out.first.push_back(t.first);
    out.second.push_back(t.second);
  }
  if(chimax&&doTrunc){
    pair<Real,uint> t=lambda0.truncate(tr_weight,chimax,normalize);
    out.first.push_back(t.first);
    out.second.push_back(t.second);
  }
  dropDangling(mps0,lambda0,A1);
  lambda.push_back(lambda0);
  if(direction>0){
    mps.push_back(mps0);
    A1*=lambda0;
  }
  if(direction<0){
    mps.push_back(mps0*lambda0);
  }
  if(direction==0){
    mps.push_back(mps0);
  }

  if(Ns>2){
    //we need to remove the first entry in itokey, since we've split it up already ...
    itokey.erase(itokey.begin());
    tout=split(A1,mps,lambda,itokey,direction,chimax,doTrunc,tr_weight,normalize,maxthreads);
    for(uint u=0;u<tout.first.size();++u)
      out.first.push_back(tout.first[u]);
    for(uint u=0;u<tout.second.size();++u)
      out.second.push_back(tout.second[u]);
  }
  return out;
}



//adapts A_L and A_R to match the sizes of lambda, and erases all blocks in A_L and A_R which do not have a corresponding lambda block
template<int Ns1,int Ns2,class KeyType,typename T>
bool dropDangling(BlockMatrix<Ns1,KeyType,T>&L,DiagBlockMatrix<KeyType,Real>&lambda,BlockMatrix<Ns2,KeyType,T> &R){
  bool out = false;
  typename DiagBlockMatrix<KeyType,Real>::iterator e = lambda.end(),lam;
  boost::numeric::ublas::matrix<T,column_major> tmp;
  typename BlockMatrix<Ns1,KeyType,T>::iterator l;  
  typename BlockMatrix<Ns2,KeyType,T>::iterator r;  
  std::vector<typename BlockMatrix<Ns1,KeyType,T>::iterator> el;
  std::vector<typename BlockMatrix<Ns2,KeyType,T>::iterator> er;
  for(l=L.begin();l!=L.end();++l){
    lam = lambda.find(l->first[1]);
    if(lam==e){out=true;el.push_back(l);continue;}
    if(l->second.size2()>lam->second.size())
      l->second.resize(l->second.size1(),lam->second.size(),true);
    //this is now done somewhere else
    // if(direction<0&&lam!=e){
    //   MatMatProd(T(1.),l->second,lam->second,tmp,'n','n');
    //   l->second = tmp;
    // }
  }
  for(uint i = 0;i<el.size();++i)
    L.erase(el[i]);

  el.clear();

  for(r=R.begin();r!=R.end();++r){
    lam = lambda.find(r->first[0]);
    if(lam==e){out=true;er.push_back(r);continue;}
    if(r->second.size1()>lam->second.size())
      r->second.resize(lam->second.size(),r->second.size2(),true);
    //this is now done somewhere else
    // if(direction>0&&lam!=e){
    //   MatMatProd(T(1.),lam->second,r->second,tmp,'n','n');
    //   r->second = tmp;
    // }

  }
  for(uint i = 0;i<er.size();++i)
    R.erase(er[i]);
  return out;
}






//application of a two site gate
template<class KeyType,typename T>
pair<std::vector<Real>,std::vector<uint> > applyGate(DiagBlockMatrix<KeyType,Real> &lambdal,BlockMatrix<1,KeyType,T> &Gammal,
						     DiagBlockMatrix<KeyType,Real> &lambdac, BlockMatrix<1,KeyType,T> &Gammar,
						     DiagBlockMatrix<KeyType,Real>  &lambdar,const SparseOperator<2,T> &Gate,
						     i_to_key<KeyType>itokeyl,i_to_key<KeyType>itokeyr,const uint chimax,const bool doTrunc,
						     const Real tr_weight,const bool normalize=true,const int maxthreads=omp_get_max_threads())
{
  BlockMatrix<1,KeyType,T>MPSl,MPSr,MPST;
  typename  BlockMatrix<2,KeyType,T>::iterator itAin,ittheta;
  typename SparseOperator<2,T>::const_iterator itG;
  boost::numeric::ublas::matrix<T,boost::numeric::ublas::column_major>tmp1;
  //first create the A and B matrices:

  MPSl = lambdal*Gammal*lambdac;
  MPSr = Gammar*lambdar;

  //now create TwoSiteBlockmatrix
  BlockMatrix<2,KeyType,T> Ain(MPSl,MPSr),theta;
  while(Ain.size()!=0){
    itAin = Ain.begin();
    for(itG = Gate.begin();itG!=Gate.end();itG++){
      if(itG->first[0]==itAin->first(0)&&itG->first[1]==itAin->first(1)){
	tmp1 = itG->second*itAin->second;
	std::vector<uint>v;v.push_back(itG->first(0));v.push_back(itG->first(1));
	TensorKeyType<2,2,KeyType> ky(itAin->first[0],itAin->first[1],v);
	//asserts conservation of quantum numbers by the gate:
	assert(itAin->first[1]==itAin->first[0]+itokeyl[v[0]]+itokeyr[v[1]]);
	
	ittheta=theta.find(ky);
	if(ittheta!=theta.end())
	  tmp1+=ittheta->second;
	theta[ky] = tmp1;
      }
    }
    Ain.erase(itAin);
  }
  theta.squeeze(1e-15);//to be sure ....
  std::vector<BlockMatrix<1,KeyType,T> >mps;
  std::vector<DiagBlockMatrix<KeyType,Real> >lambda;
  std::vector<i_to_key<KeyType> >itokey(2);
  itokey[0] = itokeyl;  itokey[1] = itokeyr;
  pair<std::vector<Real>,std::vector<uint> > out =split(theta,mps,lambda,itokey,0,chimax,doTrunc,tr_weight,normalize,maxthreads);
  //and now the lambdal and lambdar have to be reintroduced

  mps[0].squeeze(1e-15);  //get rid of empty entries:
  mps[1].squeeze(1e-15);  //get rid of empty entries:

  Gammal = lambdal/mps[0];
  Gammar = mps[1]/lambdar;
  lambdac=lambda[0];
  return out;
  //thats it, now the state is canonized again
}





//application of a two site gate, but it does not assume anything about the matrices; they don't have to be orthonormalized, and they will not be after
//application of applyGate. it still uses SVD, which is an overkill. qr will be implemented later.
//it also DOES NOT TRUNCATE THE STATE!!! successive application without compression is infeasible!!!!!
template<class KeyType,typename T>
pair<std::vector<Real>,std::vector<uint> > applyGate(BlockMatrix<1,KeyType,T> &MPSl,
						     BlockMatrix<1,KeyType,T> &MPSr,const SparseOperator<2,T> &Gate,
						     i_to_key<KeyType>itokeyl,i_to_key<KeyType>itokeyr,const int maxthreads=omp_get_max_threads())
{
  BlockMatrix<1,KeyType,T>MPST;
  typename  BlockMatrix<2,KeyType,T>::iterator itAin,ittheta;
  typename SparseOperator<2,T>::const_iterator itG;
  boost::numeric::ublas::matrix<T,boost::numeric::ublas::column_major>tmp1;
  
  //now create TwoSiteBlockmatrix
  BlockMatrix<2,KeyType,T> Ain(MPSl,MPSr),theta;
  while(Ain.size()!=0){
    itAin = Ain.begin();
    for(itG = Gate.begin();itG!=Gate.end();itG++){
      if(itG->first[0]==itAin->first(0)&&itG->first[1]==itAin->first(1)){
	tmp1 = itG->second*itAin->second;
	std::vector<uint>v;v.push_back(itG->first(0));v.push_back(itG->first(1));
	TensorKeyType<2,2,KeyType> ky(itAin->first[0],itAin->first[1],v);
	assert(itAin->first[1]==itAin->first[0]+itokeyl[v[0]]+itokeyr[v[1]]);
	ittheta=theta.find(ky);
	if(ittheta!=theta.end())
	  tmp1+=ittheta->second;
	theta[ky] = tmp1;
      }
    }
    Ain.erase(itAin);
  }
  theta.squeeze(1e-14);//that erases very small (unphysical) entries of theta and IS crucial!! 
  std::vector<BlockMatrix<1,KeyType,T> >mps;
  std::vector<DiagBlockMatrix<KeyType,Real> >lambda;
  std::vector<i_to_key<KeyType> >itokey(2);
  itokey[0] = itokeyl;  itokey[1] = itokeyr;
  uint chimax=itokeyl.size()*MPSl.size2();
  //there is no compression done in here
  pair<std::vector<Real>,std::vector<uint> > out =split(theta,mps,lambda,itokey,0,chimax,true,1e-16,false,maxthreads);
  mps[0].squeeze(1e-15);  //get rid of empty entries:
  mps[1].squeeze(1e-15);  //get rid of empty entries:
  MPSl = mps[0];
  MPSr = lambda[0]*mps[1];
  return out;
}


//this is under construction: the prolongation (like bela bauer didi it)

//template<class KeyType,typename T>
//pair<std::vector<Real>,std::vector<uint> > prolongation(BlockMatrix<1,KeyType,T> &MPS,i_to_key<KeyType>itokey,const int maxthreads=omp_get_max_threads())
//{
//  BlockMatrix<1,KeyType,T>Ain(MPS);
//  BlockMatrix<2,KeyType,T> theta;
//  typename  BlockMatrix<2,KeyType,T>::iterator itAin,ittheta;
//  boost::numeric::ublas::matrix<T,boost::numeric::ublas::column_major>tmp1;
//
//  while(Ain.size()!=0){
//    itAin = Ain.begin();
//    for(it = Isometry.begin();it!=Isometry.end();it++){
//      if(it->first[0]==itAin->first(0)){
//	tmp1 = it->second*itAin->second;
//	std::vector<uint>v;v.push_back(it->first[1]);v.push_back(it->first[2]);
//	TensorKeyType<2,2,KeyType> ky(itAin->first[0],itAin->first[1],v);
//	assert(itAin->first[1]==itAin->first[0]+itokeyl[v[0]]+itokeyr[v[1]]);
//	ittheta=theta.find(ky);
//	if(ittheta!=theta.end())
//	  tmp1+=ittheta->second;
//	theta[ky] = tmp1;
//      }
//    }
//    Ain.erase(itAin);
//  }
//  theta.squeeze(1e-14);//that erases very small (unphysical) entries of theta and IS crucial!! 
//  std::vector<BlockMatrix<1,KeyType,T> >mps;
//  std::vector<DiagBlockMatrix<KeyType,Real> >lambda;
//  std::vector<i_to_key<KeyType> >itokey(2);
//  itokey[0] = itokeyl;  itokey[1] = itokeyr;
//  uint chimax=itokeyl.size()*MPSl.size2();
//  //there is no compression done in here
//  pair<std::vector<Real>,std::vector<uint> > out =split(theta,mps,lambda,itokey,0,chimax,true,1e-16,false,maxthreads);
//  mps[0].squeeze(1e-15);  //get rid of empty entries:
//  mps[1].squeeze(1e-15);  //get rid of empty entries:
//  MPSl = mps[0];
//  MPSr = lambda[0]*mps[1];
//  return out;
//}


//that's application of a three-site gate...
template<class KeyType,typename T>
pair<std::vector<Real>,std::vector<uint> >
applyGate(DiagBlockMatrix<KeyType,Real> &lambdal,BlockMatrix<1,KeyType,T> &Gammal,DiagBlockMatrix<KeyType,Real> &lambdacl,
	  BlockMatrix<1,KeyType,T> &Gammac,DiagBlockMatrix<KeyType,Real> &lambdacr,BlockMatrix<1,KeyType,T> &Gammar,
	  DiagBlockMatrix<KeyType,Real> &lambdar,
	  const SparseOperator<3,T> &Gate,i_to_key<KeyType>itokeyl,i_to_key<KeyType>itokeyc,i_to_key<KeyType>itokeyr,
	  const uint chimax,const bool doTrunc,const Real tr_weight,const bool normalize=true,const int maxthreads=omp_get_max_threads()){

  BlockMatrix<1,KeyType,T> MPSl,MPSc,MPSr;
  typename BlockMatrix<1,KeyType,T>::iterator itl,itr;
  typename SparseOperator<3,T>::const_iterator itG;
  boost::numeric::ublas::matrix<T,boost::numeric::ublas::column_major>tmp1;

  MPSl = lambdal*Gammal;
  MPSc = lambdacl*Gammac*lambdacr;
  MPSr = Gammar*lambdar;

  BlockMatrix<3,KeyType,T> Ain(MPSl,MPSc,MPSr),theta;
  typename BlockMatrix<3,KeyType,T> ::iterator itAin,ittheta;
  while(Ain.size()!=0){
    itAin = Ain.begin();
    for(itG = Gate.begin();itG!=Gate.end();itG++){
      if(itG->first[0]==itAin->first(0)&&itG->first[1]==itAin->first(1)&&itG->first[2]==itAin->first(2)){
	tmp1 = itG->second*itAin->second;
	std::vector<uint>v;v.push_back(itG->first(0));v.push_back(itG->first(1));v.push_back(itG->first(2));
	TensorKeyType<2,3,KeyType> ky(itAin->first[0],itAin->first[1],v);
	assert(itAin->first[1]==itAin->first[0]+itokeyl[v[0]]+itokeyc[v[1]]+itokeyr[v[2]]);
	ittheta=theta.find(ky);
	if(ittheta!=theta.end())
	  tmp1+=ittheta->second;
	theta[ky] = tmp1;
      }
    }
    Ain.erase(itAin);
  }

  theta.squeeze(1e-15);//to be sure ....
  std::vector<BlockMatrix<1,KeyType,T> >mps;
  std::vector<DiagBlockMatrix<KeyType,Real> >lambda;

  std::vector<i_to_key<KeyType> >itokey(3);
  itokey[0] = itokeyl;  itokey[1] = itokeyc; itokey[2] = itokeyr;
  pair<std::vector<Real>,std::vector<uint> >  out =split(theta,mps,lambda,itokey,1,chimax,doTrunc,tr_weight,normalize,maxthreads);

  mps[0].squeeze(1e-15);  //get rid of empty entries:
  mps[1].squeeze(1e-15);  //get rid of empty entries:
  mps[2].squeeze(1e-15);  //get rid of empty entries:
  lambdacl=lambda[0];
  lambdacr=lambda[1];

  Gammal = lambdal/mps[0];
  Gammac = lambdacl/mps[1];
  Gammar =mps[2]/lambdar;
  
  return out;
}



//truncates the mps; make sure it is orthonormalized!!!!!
//chimax is the hard limit; if delta!=0, then its truncated down to delta, but NEVER higher than chimax
//Be aware that the state maybe not normalized (I'm not completely sure since i normalize the lambdas...)
//also note that in a recursive sceme like the Chebychev stuff, the state may get a very very very small norm, thus 
//it could happen that on some bonds the state is reduced to zero size by truncate
template<typename T,class KeyType>
void truncate(MPSBase<KeyType,T>&mps,const uint chimax,const Real delta=0.0){
  if(!chimax){
    cout<<"no chimax given  ..."<<endl;
    return;
  }
  //now truncate ...
  if((chimax)&&(delta==0.0)){
    for(uint site=0;site<(mps.size()-1);++site){
      // if(mps.chi(site)<=chimax)
      // continue;
      mps.lambda(site).truncate(delta,chimax,false);
      mps.lambda(site).normalize();
      dropDangling(mps[site],mps.lambda(site),mps[site+1]);
    }
  }
  else if((chimax)&&(delta!=0.0)){
    for(uint site=0;site<(mps.size()-1);++site){
      //if(mps.chi(site)<=chimax)
      //continue;
      mps.lambda(site).truncate(delta,chimax,false);
      mps.lambda(site).normalize();
      dropDangling(mps[site],mps.lambda(site),mps[site+1]);
    }
  }
  //NO!!! DON'T ORTHONORMALIZE IT; this has to be done after this routine ....
  //and then orthonormalize it again; that's because the norm of the state is wrong now ...
  //maybe two sweeps is too much; for critical applications, one may only need one of the two.
  // if(mps.getROrtho()){
  //   orthonormalize(mps,"lr");
  //   orthonormalize(mps,"rl");
  //   orthonormalize(mps,"lr");
  // }
  // else{
  //   orthonormalize(mps,"rl");
  //   orthonormalize(mps,"lr");
  // }
}

//this function does an orthonormalization of an mps. 
//Arguments: 
//mps: is the mps; note that both CanonizedMPS and MPS can be passed. For the first, orthonormalization takes place on the Gamma matrices, which means
//that the state is actually changed!!!
//dir: direction of orthonormalization: either "lr" or "rl" (left-right of right-left); IMPORTANT: DO NOT ORTHONORMALIZE TWICE IN THE SAME DIRECTION!!!! this will result in a phase 
//inconsistency of the lambdas and the Gamma matrices to the left and to the right of each lambda. This inconsistency can be removed by applying a single orthonormalization into the other 
//direction after any number of orthonormalizations in the other direction have been applied. 
//storenorm: if true, the norm of the mps is stored in the mps-variable _norm, which can be assessed by mps.norm().

//The function sets a mps-intern flag so that the mps knows in which state it is (left or right orthonormal). The final state is either in A...AAA ((dir="lr") or B...BBB (dir="rl") form.
//lambdas are always stored in the mps; note that if mps is in no particular orthonormal form prior to orthonormalize, the lambdas after orthonormalize cannot be interpreted as schmidtvalues.
//the norm of mps after this routine is 1. The old norm can be assessed by mps.norm(). The function prepareSite is defined in MPSTypes.hpp.

template<class KeyType,typename T>
void orthonormalize(MPSBase<KeyType,T> &mps,const string dir,const bool storenorm=false)
{
  //  cout<<"orthonormalizing ...."<<endl;
  const uint N = mps.size();
  if(dir=="lr"){
    if((!mps.getLOrtho()&&mps.getROrtho())||(!mps.getLOrtho()&&!mps.getROrtho())){
      for(int site = 0;site<N;++site){
	prepareSite(mps,site,"l",true,storenorm,true,1e-16);
      }
      //      cout<<"the norm after ortho: "<<mps.lambda(N-1).norm()<<endl;
      mps.setROrtho(false);
      mps.setLOrtho(true);
    }else if(mps.getLOrtho()&&!mps.getROrtho()){
      cout<<"orthonormalize(mps,\"lr\"): detected left orthonormalized state; skipping orthonormalization"<<endl;return;
    }else if(mps.getLOrtho()&&mps.getROrtho()){
      cout<<"orthonormalize(mps,\"lr\"): detected LEFT AND RIGHT orthonormalized state; strange ... skipping"<<endl;return;
    }
  }
  else if(dir=="rl"){
    if((mps.getLOrtho()&&!mps.getROrtho())||(!mps.getLOrtho()&&!mps.getROrtho())){
      for(int site = N-1;site>=0;--site){
	if(site<0)cerr<<"site < 0!"<<endl;
	prepareSite(mps,site,"r",true,storenorm,true,1e-16);
      }
      //      cout<<"the norm after ortho: "<<mps.lambda(-1).norm()<<endl;
      mps.setROrtho(true);
      mps.setLOrtho(false);
    }else if(!mps.getLOrtho()&&mps.getROrtho()){
      cout<<"orthonormalize(mps,\"rl\"): detected right orthonormalized state; skipping orthonormalization"<<endl;return;
    }else if(mps.getLOrtho()&&mps.getROrtho()){
      cout<<"orthonormalize(mps,\"rl\"): detected LEFT AND RIGHT orthonormalized state; strange ... skipping sicherheitshalber"<<endl;return;
    }
  }else if(dir!="rl"&&dir!="lr"){
    cerr<<"orthonormalize: wrong input string dir; aborting ... "<<dir<<endl;
    abort();
  }
}

//raw orthonormalization, keeps norm of the mps (e.g. at the boundary, it does not throw away anything!). The reslulting state should be in the form MA...A or B...BM, where
//M is a non-orthonormali matrix.
template<class KeyType,typename T>
void rawOrthonormalize(MPS<KeyType,T> &mps,const int dir){

  if(dir>=0){
    for(uint site = 0;site<mps.size()-1;++site){
      prepareSite(mps,site,"l",true);
    }
    BlockMatrix<1,KeyType,T>a=mps[mps.size()-1];
    prepareMatrix(a,mps.lambda(mps.size()-1),mps._rightDummy,"l",true,false,1e-12);
    for(int site = mps.size()-1;site>0;--site){
      prepareSite(mps,site,"r",true);
    }
    a=mps[0];
    prepareMatrix(mps._leftDummy,mps.lambda(-1),a,"r",true,false,1e-12);

  }
  if(dir<0){
    for(int site = mps.size()-1;site>0;--site){
      prepareSite(mps,site,"r",true);
    }
    BlockMatrix<1,KeyType,T>a=mps[0];
    prepareMatrix(mps._leftDummy,mps.lambda(-1),a,"r",true,false,1e-12);
    for(uint site = 0;site<mps.size()-1;++site){
      prepareSite(mps,site,"l",true);
    }
    a=mps[mps.size()-1];
    prepareMatrix(a,mps.lambda(mps.size()-1),mps._rightDummy,"l",true,false,1e-12);
  }
}

//canonizes an mps. Assumes nothing about "mps" 
template<class KeyType,typename T>
CanonizedMPS<KeyType,T> canonize(MPS<KeyType,T>&mps){
  //  cout<<"canonizing"<<endl;
  CanonizedMPS<KeyType,T> cmps;
  if(mps.getROrtho()){
    //copy the lambdas to cmps; the leftmost lambda is a bit stupid to obtain if it's not there from initialization of the mps;
    orthonormalize(mps,"lr");
    cmps=mps;
    cmps.clearMatrices();
    cmps[0] = mps[0];
    for(int site = 1;site<mps.size();++site){
      cmps[site]=mps.lambda(site-1)/mps[site];
    }
  }
  if(mps.getLOrtho()){
    orthonormalize(mps,"rl");
    //copy the lambdas to cmps; the leftmost lambda is a bit stupid to obtain if it's not there from initialization of the mps;
    cmps=mps;
    cmps.clearMatrices();
    cmps[mps.size()-1] = mps[mps.size()-1];
    for(int site = (mps.size()-2);site>=0;--site){
      cmps[site]=mps[site]/mps.lambda(site);
    }
  }
  if(!mps.getLOrtho()&&!mps.getROrtho()){
    cout<<"in canonize: doing orthonormalization"<<endl;
    orthonormalize(mps,"lr");
    canonize(mps);
  }
  return cmps;
}


//canonizes an "mps" according to direction "dir". Assumes that "mps" has been orthonormalized! 
template<class KeyType,typename T>
CanonizedMPS<KeyType,T> rawCanonize(MPS<KeyType,T> &mps,const int dir){
  cout<<"#rawCanonizing ..."<<endl;
  CanonizedMPS<KeyType,T> cmps;
  cmps.clearMatrices();
  cmps.resize(mps.size());
  if(dir<0){
    //copy the lambdas to cmps; the leftmost lambda is a bit stupid to obtain if it's not there from initialization of the mps;
    cmps=mps;
    cmps.clearMatrices();
    cmps[mps.size()-1] = mps[mps.size()-1];
    for(int site = (mps.size()-2);site>=0;--site){
      if(mps.lambda(site).size()==0){
	cout<<"#rawCanonize: input mps contains no lambdas at site "<<site<<endl;
	abort();
      }
      cmps[site]=mps[site]/mps.lambda(site);
      cmps.lambda(site)=mps.lambda(site);
    }
    // cmps.lambda(mps.size()-1)=mps.lambda(mps.size()-1);
    // cmps.lambda(-1)=mps.lambda(-1);
  }  
  if(dir>=0){
    //copy the lambdas to cmps; the leftmost lambda is a bit stupid to obtain if it's not there from initialization of the mps;
    cmps=mps;
    cmps.clearMatrices();
    cmps[0] = mps[0];
    for(int site = 1;site<mps.size();++site){
      if(mps.lambda(site-1).size()==0){
	cout<<"#in rawCanonize: input mps contains no lambdas at site "<<site<<endl;
	abort();
      }
      cmps[site]=mps.lambda(site-1)/mps[site];
      cmps.lambda(site-1)=mps.lambda(site-1);
    }
    // cmps.lambda(mps.size()-1)=mps.lambda(mps.size()-1);
    // cmps.lambda(-1)=mps.lambda(-1);
  }
  cout<<"#done ... "<<endl;
  return cmps;
}
  
template<class KeyType,typename T>
BlockMatrix<1,KeyType,T> getBoundaryBL(const BlockMatrix<1,KeyType,T> &A,const uint dir){
  typename BlockMatrix<1,KeyType,T>::const_iterator a;
  BlockMatrix<1,KeyType,T> BL;
  for(a=A.begin();a!=A.end();++a){
    uint size=(dir<=0)?a->second.size1():a->second.size2();
    TensorKeyType<2,1,KeyType>ky(a->first[dir],a->first[dir],0);
    identity_matrix<T,column_major> mat(size);
    BL[ky] = mat;
  }
  return BL;
}

template<class KeyType,typename T>
BlockMatrix<1,KeyType,T> getLMin(const MPSBase<KeyType,T>&mps){
  BlockMatrix<1,KeyType,T> LMin;
  KeyType key = mps.getItoKey(0).begin()->second;
  reset(key);
  TensorKeyType<2,1,KeyType> Tkey(key,key,0);
  boost::numeric::ublas::matrix<T,column_major> mat(1,1);mat(0,0)=1;
  LMin[Tkey] = mat;
  return LMin;
}


template<class KeyType,typename T>
BlockMatrix<1,KeyType,T> fixEnd(const MPSBase<KeyType,T>&mps){
  const int n=mps.size();
  BlockMatrix<1,KeyType,T> A=mps[n-1];
  typename BlockMatrix<1,KeyType,T>::iterator a;
  int m=0;
  for(a=A.begin();a!=A.end();++a){
    m=m>a->second.size2()?m:a->second.size2();
  }
  for(a=A.begin();a!=A.end();++a){
    if(a->second.size2()==m)
      continue;
    matrix<T,column_major> t(a->second.size1(),m);
    boost_setTo(t,T(0.0));
    matrix_range<matrix<T,column_major> >(t,Range(0,a->second.size1()),Range(0,a->second.size2()))=a->second;
    a->second=t;
  }
  return A;
}

//returns the right boundary R; this may have many different qn's
// template<class KeyType,typename T>
// BlockMatrix<1,KeyType,T> getRMax(const MPSBase<KeyType,T>&mps){
//   // get all right keys.
//   const int N=mps.size();
//   BlockMatrix<1,KeyType,T> A=mps[N-1]; //better not here ...
//   std::set<KeyType> ps;
//   typename BlockMatrix<1,KeyType,T>::iterator a;
//   for(a=A.begin();a!=A.end();++a)
//     ps.insert(a->first[1]);
//   BlockMatrix<1,KeyType,T> RMax;
//   typename std::set<KeyType>::iterator i,j;
//   for(i=ps.begin();i!=ps.end();++i) {
//     const KeyType & p1 = *i;
//     for(j=ps.begin();j!=ps.end();++j) {
//       const KeyType & p2 = *j;
//       // now get shape.
//       unsigned int m1=0,m2=0;
//       bool found1=false,found2=false;
//       for(a=A.begin();a!=A.end();++a) {
// 	if (a->first[1]==p1) {
// 	  m1=a->second.size2();
// 	  found1=true;
// 	}
// 	if (a->first[1]==p2) {
// 	  m2=a->second.size2();
// 	  found2=true;
// 	}
// 	if(found1&&found2)
// 	  break;
//       }
//       if(m1!=m2){
// 	cout<<"in getRMax: the left most mps-matrix has differnt blocks of different size2()"<<endl;
// 	abort();
//       }
//       assert( m1 );
//       assert( m2 );
//       assert( m1==m2 );
//       TensorKeyType<2,1,KeyType>ky(p1,p2,uint(0));
//       if(!RMax.count(ky))
// 	RMax[ky] = identity_matrix<T>(m1);
//     }
//   }
//   return RMax;
// }


template<class KeyType,typename T>
BlockMatrix<1,KeyType,T> getRMax(const MPSBase<KeyType,T>&mps){
  // get all right keys.
  const int N= mps.size();
  std::set<KeyType> ps;
  typename BlockMatrix<1,KeyType,T>::const_iterator a;
  for(a=mps[N-1].begin();a!=mps[N-1].end();++a)
    ps.insert(a->first[1]);
    
  BlockMatrix<1,KeyType,T> RMax;
  BlockMatrix<1,KeyType,T> A=mps[mps.size()-1];
  typename std::set<KeyType>::iterator i;
  for(i=ps.begin();i!=ps.end();++i) {
    const KeyType & p = *i;
    // now get shape.
    unsigned int m =0;
    for(a=A.begin();a!=A.end();++a) {
      if (a->first[1]==p) {
	m = a->second.size2();
	break;
      }
    }
    assert(m);
    RMax[TensorKeyType<2,1,KeyType>(p,p,uint(0))] = identity_matrix<T>(m);
  }
  return RMax;
}

//computes all R's and stores them in R, site included
template<class KeyType,typename T>
void computeR(const MPS<KeyType,T>&mps,const MPO<T> &mpo,const int site,map_container<int,BlockMatrix<1,KeyType,T> >& R){
  const int N=mps.size();
  MPO<T>_mpo=mpo;
  assert(N==_mpo.size());
  //if(!mps.getROrtho()){cerr<<"computeR: mps is not r-orthonormal"<<endl;}
  R.clear();
  R[N]= getRMax(mps);
  for(int s = N-1;s>=site;--s){
    assert(s>=0);
    R[s] = addContractionLayer(R[s+1],mps[s],_mpo[s],mps.getItoKey(s),0);
  }
}

//computes all L's and stores them in L; site has to be N-1 at most; site is included
template<class KeyType,typename T>
void computeL(const MPS<KeyType,T>&mps,const MPO<T> &mpo,const int site,map_container<int,BlockMatrix<1,KeyType,T> >& L){
  const int N=mps.size();
  MPO<T>_mpo=mpo;
  assert(N==_mpo.size());
  //if(!mps.getLOrtho()){cerr<<"computeL: mps is not l-orthonormal"<<endl;}
  L.clear();
  L[-1]= getLMin(mps);
  for(int s = 0;s<=site;++s){
    assert(s!=N);
    L[s] = addContractionLayer(L[s-1],mps[s],_mpo[s],mps.getItoKey(s),1);
  }
}

//checks MPS and also removes unconnected matrices
template<class KeyType,typename T>
bool checkMPS(MPSBase<KeyType,T>&mps){
  bool out=true;
  for(int site = 0;site<(mps.size()-1);++site){
    if(!checkBMs(mps[site],mps[site+1]))
      out = false;
  }
  return out;
}

//note that functions computes the complex cojugate ...
template<typename T,class KeyType>
T computeOverlap(const MPSBase<KeyType,T>&a,const MPSBase<KeyType,T>&b){
  BlockMatrix<1,KeyType,T>L=getLMin(a);
  typename BlockMatrix<1,KeyType,T>::iterator l;
  assert(a.size()==b.size());
  //assert(a.getPSpace()==b.getPSpace()); this is now working properly ...
  for(uint site=0;site<a.size();++site){
    MPOMat<1,T> ID;
    SparseOperator<1,T> id(a.getPSpace()(site),a.getPSpace()(site));
    for(uint ind=0;ind<a.getPSpace()(site);++ind)
      id[OpKeyType<1>(ind,ind)]= T(1.);
    ID[ukey<2>(0,0)]=id;
    ID.Dl=1;ID.Dr=1;
    L=addContractionLayer(L,a[site],ID,b[site],a.getItoKey(site),1);
  }
  //cout<<"L.size() = "<<L.size()<<endl;
  //cout<<L<<endl;
  return L.size()>0?L.begin()->second(0,0):T(0.0);
}


//this is a monster ... but nice one :)
template<class KeyType,typename T>
MPS<KeyType,T> applyMPO(MPO<T> &mpo,MPS<KeyType,T> &mps)
{
  //  cout<<"#applying MPO ..."<<endl;
  using boost::numeric::ublas::matrix;  
  using boost::numeric::ublas::matrix_range;
  using boost::numeric::ublas::column_major;
  MPS<doubleKey<KeyType,uint>,T>FiveLegs(mps.size());
  MPS<KeyType,T>omps(mps.getPSpace(),mps.getItoKey());
  typename BlockMatrix<1,KeyType,T>::iterator itA;
  typename MPOMat<1,T>::const_iterator itM;
  //this is a vector containing maps for MPS-MPO matrices; not the matrices are mappde, but the iterators to the individual matrices and MPO entries
  //this has to be a map cause a index-pair (q,beta) can ONLY have ONE quantum number!!!
  std::map<doubleKey<KeyType,uint>,KeyType> newQNLeft,newQNRight;//this is the new number_map for the enlarged MPS
  //this on the other hand has to be a multimap cause ONE q can have many index pairs (q_l,beta)
  std::multimap<KeyType,doubleKey<KeyType,uint> >invNewQNLeft,invNewQNRight;
  std::map<KeyType,UintVectorType> leftSizes,rightSizes;
  typename std::map<KeyType,UintVectorType> ::iterator itrightSize;

  typename std::vector<typename BlockMatrix<1,doubleKey<KeyType,uint>,T>::iterator> Matsvec;
  typename std::vector<typename BlockMatrix<1,doubleKey<KeyType,uint>,T>::iterator>::iterator Matsvecit;
  i_to_key<KeyType> lmap;
  typename i_to_key<KeyType>::iterator lmapit;
  UintVectorType largestsize1,largestsize2,startx,starty;
  //this will be the new quantum numbers ...
  KeyType ql=mps.getItoKey(0).begin()->second,qr;
  reset(ql);
  typename BlockMatrix<1,doubleKey<KeyType,uint>,T> ::iterator BMit;
  matrix<T,column_major> mat,mat2;
  //first compute all MPO-MPS matrices; this is not really efficient, but it is save ...
  //this can be parallelized
  bool first = true;
  for(uint site = 0;site<mps.size();site++){
    BlockMatrix<1,doubleKey<KeyType,uint>,T> &A = FiveLegs[site];
    for(itM = mpo[site].begin();itM != mpo[site].end();++itM){
      SparseOperator<1,T> lo = itM->second;
      for(typename SparseOperator<1,T>::iterator li=lo.begin();li!=lo.end();++li){
	for(itA = mps[site].begin();itA!=mps[site].end();++itA){
	  //is the p-ind of the mps the same as the one of the Local operator entry?
	  if(itA->first(0)==li->first[0]){
	    mat=itA->second*li->second;
	    doubleKey<KeyType,uint> keyLeft,keyRight;
	    keyLeft.first = itA->first[0];
	    keyLeft.second =itM->first[0];
	    keyRight.first = itA->first[1];
	    keyRight.second =itM->first[1];
	    TensorKeyType<2,1,doubleKey<KeyType,uint> >key(keyLeft,keyRight,li->first(0));
	    if(site==0&&first){
	      if(newQNLeft.find(keyLeft)!=newQNLeft.end()){
		cout<<"bad ...at site 0, there should only be one bond"<<endl;
		abort();
	      }
	      newQNLeft[keyLeft]=ql;
	      invNewQNLeft.insert(pair<KeyType,doubleKey<KeyType,uint> >(ql,keyLeft));
	      first = false;
	    }

	    BMit = A.find(key);
	    if(BMit==A.end()){
	      A[key]=mat;
	    }else{
	      // cout<<"being here means that the input MPO is probably NOT HERMITIAN. I really advice you to stop here and rethink what you are doing, unless you REALLY know what you're doing. i dont know what will happen now, be prepared for anything :)"<<endl;
	      // my_pause();
	      const uint x = BMit->second.size1()>mat.size1()?BMit->second.size1():mat.size1();
	      const uint y = BMit->second.size2()>mat.size2()?BMit->second.size2():mat.size2();
	      matrix<T,column_major> m1(x,y),m2(x,y);
	      boost_setTo(m1,T(0.0));
	      boost_setTo(m2,T(0.0));
	      matrix_range<matrix<T,column_major> >(m1,Range(0,BMit->second.size1()),Range(0,BMit->second.size2()))=BMit->second;
	      matrix_range<matrix<T,column_major> >(m2,Range(0,BMit->second.size1()),Range(0,BMit->second.size2()))=mat;
	      A[key]=m1+m2;
	    }
	  }
	}
      }
    }
  }
  //FiveLegs containes the new block matrices with the enlarged keys of a MPS-MPO matrix
  typename BlockMatrix<1,doubleKey<KeyType,uint>,T> ::iterator flit;
  typename std::map<doubleKey<KeyType,uint>,KeyType>::iterator it;
  leftSizes.clear();
  reset(ql);
  largestsize1.resize(mpo[0].Dl);
  boost_setTo(largestsize1,uint(1));
  leftSizes.insert(std::pair<KeyType,UintVectorType>(ql,largestsize1));

  for(uint site = 0;site<mps.size();++site){
    //    cout<<"at site "<<site<<endl;
    BlockMatrix<1,doubleKey<KeyType,uint>,T > &A = FiveLegs[site];
    //cout<<"entries in FiveLeg at site "<<site<<"before inserting them: "<<A.size()<<endl;
    // cout<<"the left quantum numbers"<<endl;
    // for(it=newQNLeft.begin();it!=newQNLeft.end();++it)
    //   cout<<it->first<<"->"<<it->second<<endl;

    // cout<<"that's before doing anything, at site "<<site<<endl;
    // for(flit=FiveLegs[site].begin();flit!=FiveLegs[site].end();++flit){
    //   cout<<flit->first<<endl;
    //   //      cout<<flit->second<<endl;
    // }

    largestsize1.resize(mpo[site].Dl);
    startx.resize(mpo[site].Dl);
    largestsize2.resize(mpo[site].Dr);
    starty.resize(mpo[site].Dr);
    lmap = mps.getItoKey(site);
    Matsvec.clear();
    for(BMit = A.begin();BMit!=A.end();++BMit){
      //Look up if the new matrix has a corresponding entry in the left new numbertable; if not, maybe it should be erased ....
      if(newQNLeft.find(BMit->first[0])!=newQNLeft.end()){
	//now store for every matrix its new qr value together with its pair-indices, this will later be used 
	qr = newQNLeft[BMit->first[0]]+lmap[BMit->first(0)];
	if(!newQNRight.count(BMit->first[1])){
	  newQNRight[BMit->first[1]]=qr;
	  invNewQNRight.insert(pair<KeyType,doubleKey<KeyType,uint> >(qr,BMit->first[1]));
	}
      }else if(newQNLeft.find(BMit->first[0])==newQNLeft.end()){
	Matsvec.push_back(BMit);

      }
    }//end of for(BMit = ...)
    for(Matsvecit=Matsvec.begin();Matsvecit!=Matsvec.end();++Matsvecit)
      A.erase((*Matsvecit));
    typename std::multimap<KeyType,doubleKey<KeyType,uint> >::iterator inql,inqll,inqlr;
    typename std::map<doubleKey<KeyType,uint>,KeyType>::iterator nql,nqr;
    pair<typename std::multimap<KeyType,doubleKey<KeyType,uint> >::iterator,typename std::multimap<KeyType,doubleKey<KeyType,uint> >::iterator>
      itrangeLeft,itrangeRight;  
    inql=invNewQNLeft.begin();
    //we now have to define the coarse grid into which the FiveLegs will be inserted ...
    //the point is, that into one new MPS matrix, many different FiveLegs will be inserted, according to the grid defined by the mpo-aux indices.
    //if the mpo is e.g. 2x2, then the new MPS matrix is a 2x2 matrix of FiveLegs. if the five leg has mpo_aux indices 0,1, then it is inserted at 
    //position (0,1) into the new coarse grid. These FiveLegs however may have different old quantum numbers, e.g.
    //not all ingoing old quantum numbers of the FiveLegs are equal (the same goes for the old outgoing quantum numbers). Thus they may have different sizes,
    //and we need to now them to adapt the coarse grid, e.g. if a FiveLeg has a dimension 1 which is smaller than all others, it has to be resized and 
    //filled with 0.0 at the resized indices.
    //note that the x-sizes of the coarse grid are already known from the last iteration, so only the y-sizes have to be determined
    while(inql!=invNewQNLeft.end()){
      //get all left double indices with the same quantum numbers ...
      itrangeLeft = invNewQNLeft.equal_range(inql->first);
      ql = inql->first;
      for(lmapit = lmap.begin();lmapit!=lmap.end();++lmapit){
	qr = ql+lmapit->second;
	itrangeRight = invNewQNRight.equal_range(qr);
	for(inqll = itrangeLeft.first;inqll!=itrangeLeft.second;++inqll){
	  for(inqlr = itrangeRight.first;inqlr!=itrangeRight.second;inqlr++){
	    TensorKeyType<2,1,doubleKey<KeyType,uint> > key(inqll->second,inqlr->second,lmapit->first);
	    BMit = A.find(key);
	    if(BMit!=A.end()){
	      itrightSize = rightSizes.find(qr);
	      if(itrightSize==rightSizes.end()){
		boost_setTo(largestsize2,uint(0));
		largestsize2(BMit->first[1].second) = BMit->second.size2()>largestsize2(BMit->first[1].second)?BMit->second.size2():
		  largestsize2(BMit->first[1].second);
		rightSizes[qr]=largestsize2;
	      }
	      else{
		largestsize2 = itrightSize->second;
		largestsize2(BMit->first[1].second) = BMit->second.size2()>largestsize2(BMit->first[1].second)?BMit->second.size2():
		  largestsize2(BMit->first[1].second);
		itrightSize->second = largestsize2;
	      }
	    }
	  }
	}
      }
      inql = itrangeLeft.second;
    }
    //now rightSizes contains for every qr the largest right sizes, taken over all FiveLegs having a right QN qr
    while(invNewQNLeft.size()!=0){
      //this is the new ql (happens to be the same as the old one)
      ql = invNewQNLeft.begin()->first;
      //get all left-keys of the new MPS with the same ql
      itrangeLeft = invNewQNLeft.equal_range(invNewQNLeft.begin()->first);
      for(lmapit = lmap.begin();lmapit!=lmap.end();++lmapit){
	qr = ql+lmapit->second;
	TensorKeyType<2,1,KeyType> newkey(ql,qr,lmapit->first);
	bool insert=false;//are we going to insert anything?
	itrangeRight = invNewQNRight.equal_range(qr);
	Matsvec.clear();
	for(inqll = itrangeLeft.first;inqll!=itrangeLeft.second;++inqll){
	  for(inqlr = itrangeRight.first;inqlr!=itrangeRight.second;++inqlr){
	    TensorKeyType<2,1,doubleKey<KeyType,uint> > key(inqll->second,inqlr->second,lmapit->first);
	    BMit=A.find(key);
	    if(BMit!=A.end()){
	      insert=true;//bingo!!!
	      Matsvec.push_back(BMit);
	    }
	  }
	}
	if(leftSizes.find(ql)!=leftSizes.end())
	  largestsize1 = leftSizes.find(ql)->second;
   	
	if(rightSizes.find(qr)!=rightSizes.end())
	  largestsize2 = rightSizes.find(qr)->second;
   	
	if(insert==true){
	  uint SIZEX=0,SIZEY=0;
	  uint cnt = 0;
	  boost_setTo(startx,uint(0));
	  boost_setTo(starty,uint(0));
	  for(UintVectorType::iterator it = largestsize1.begin();it!=largestsize1.end();it++){
	    SIZEX+=*it;
	    if(cnt>0)
	      startx(cnt) = startx(cnt-1)+largestsize1(cnt-1);
	    cnt++;
	  }
	  cnt = 0;
	  for(UintVectorType::iterator it = largestsize2.begin();it!=largestsize2.end();it++){
	    SIZEY+=*it;
	    if(cnt>0)
	      starty(cnt) = starty(cnt-1)+largestsize2(cnt-1);
	    cnt++;
	  }
	  mat2.resize(SIZEX,SIZEY);
	  boost_setTo(mat2,T(0.0));
	  for(Matsvecit=Matsvec.begin();Matsvecit!=Matsvec.end();Matsvecit++){
	    //first adapt the FiveLeg to the coarse-grid sizes ...
	    uint size1=(*Matsvecit)->second.size1();
	    uint size2 = (*Matsvecit)->second.size2();
	    mat.resize(largestsize1((*Matsvecit)->first[0].second),largestsize2((*Matsvecit)->first[1].second));
	    boost_setTo(mat,T(0.0));
	    matrix_range<matrix<T,column_major> >(mat,Range(0,size1), Range(0,size2))= (*Matsvecit)->second;
	    (*Matsvecit)->second = mat;
	    //and then insert it into the coarse-grid ...
	    uint start1 = startx((*Matsvecit)->first[0].second);
	    uint stop1 = start1+largestsize1((*Matsvecit)->first[0].second);
	    uint start2 = starty((*Matsvecit)->first[1].second);
	    uint stop2 = start2+largestsize2((*Matsvecit)->first[1].second);
	    matrix_range<matrix<T,column_major> > (mat2,Range(start1,stop1), Range(start2,stop2))= (*Matsvecit)->second;
	  }
	  //and now the new block is finished!
	  omps[site][newkey]=mat2;
	  for(Matsvecit=Matsvec.begin();Matsvecit!=Matsvec.end();Matsvecit++){
	    if((*Matsvecit)==A.end())cout<<"what the hell ... "<<endl;
	    A.erase((*Matsvecit));
	  }
	}//end of if(insert==true)
      }//end of lmapit
      invNewQNLeft.erase(itrangeLeft.first,itrangeLeft.second);
    }//End of while(invNewQNLeft.size()!=0)-loop
    // cout<<"the right quantum numbers"<<endl;
    // for(it=newQNRight.begin();it!=newQNRight.end();++it)
    //   cout<<it->first<<"->"<<it->second<<endl;
    // cin.get();
    leftSizes = rightSizes;
    rightSizes.clear();
    newQNLeft=newQNRight;
    invNewQNLeft=invNewQNRight;
    newQNRight.clear();
    invNewQNRight.clear();
    //    cout<<"entries left in FiveLeg at site "<<site<<": "<<A.size()<<endl;
    // cout<<A<<endl;
    // cin.get();
    omps.setItoKey(site,mps.getItoKey(site));
  }
  //FINALLY, WE'RE DONE!!!
  //  cout<<"#done ..."<<endl;
  return omps;
}



//this does a variational single site compression; it conserves the norm of mps;
//mps has structure MBB...B after compress; M will not be normalized to identity, but it will
//be proportional to identity
template<typename T,class KeyType>
T compress(MPSBase<KeyType,T> &mps,const uint chimax=100,const uint maxstep=5,const Real delta=1e-10,const bool verbose=true,
	   const bool normalize=false){
  
  if(verbose)
    cout<<"preparing input mps...."<<endl;
  int Nmin=0,Nmax=0;
  if(normalize==false){
    Nmax=mps.size()-1;
    Nmin=1;
  }else if(normalize==true){
    Nmax=mps.size();
    Nmin=0;
  }

  for(int site = 0;site<Nmax;++site){
    prepareSite(mps,site,"l",true);
  }
  for(int site = Nmax-1;site>=Nmin;--site){
    prepareSite(mps,site,"r",true);
  }

  if(verbose)
    cout<<"done"<<endl;
  //  cout<<"#finished orthonormalizing"<<endl;
  MPSBase<KeyType,T> compmps=mps;
  if(verbose)
    cout<<"initializing mps ...."<<endl;
  // cout<<"                ||compmps|| before tuncate: "<<sqrt(computeOverlap(compmps,compmps))<<endl;
  // cout<<"                ||mps|| before tuncate:     "<<sqrt(computeOverlap(mps,mps))<<endl;
  truncate(compmps,chimax,delta);
  //cout<<"                ||compmps|| AFTER tuncate:  "<<sqrt(computeOverlap(compmps,compmps))<<endl;
  for(uint site = 0;site<mps.size()-1;++site){
    prepareSite(compmps,site,"l",true);
  }
  for(uint site = mps.size()-1;site>0;--site){
    prepareSite(compmps,site,"r",true);
  }
  if(verbose)
    cout<<"done"<<endl;
  const T n1=sqrt(computeOverlap(mps,mps));
  //cout<<"                ||compmps|| AFTER tuncate and AFTER prepareSite (lr and rl sweep): "<<sqrt(computeOverlap(compmps,compmps))<<endl;
  T n2=0;
  map_container<int,BlockMatrix<1,KeyType,T> >BL,BR;
  BL[-1] = getLMin(mps);
  BR[mps.size()]=getRMax(mps);
  bool converged = maxstep==0?true:false;
  if(verbose)
    cout<<"initializing reduced blocks ...."<<endl;
  if(!converged){
    for(uint site=compmps.size()-1;site>0;--site){
      MPOMat<1,T> id;
      SparseOperator<1,T> O_11(mps.getPSpace()(site),mps.getPSpace()(site));
      for(uint ind=0;ind<mps.getPSpace()(site);++ind)
	O_11[OpKeyType<1>(ind,ind)]=T(1.0);
      id[ukey<2>(0,0)]=O_11;
      id.Dl=1;id.Dr=1;
      //the right order of contracting the matrices into the network is not always straight forward. I put both orders here, for fast switching
      BR[site]=addContractionLayer(BR[site+1],mps[site],id,compmps[site],mps.getItoKey(site),-1);
      //BR[site]=addContractionLayer(BR[site+1],compmps[site],id,mps[site],mps.getItoKey(site),-1);
    }
  }
  if(verbose)
    cout<<"done"<<endl;
  //got all BRs ...now start compressing

  uint step=1;
  Real conv,conv_old=1.0;
  while(!converged){
    if(verbose){
      cout<<".";
      cout.flush();
    }
    for(int site=0;site<(compmps.size()-1);++site){
      MPOMat<1,T> id;
      SparseOperator<1,T> O_11(mps.getPSpace()(site),mps.getPSpace()(site));
      for(uint ind=0;ind<mps.getPSpace()(site);++ind)
	O_11[OpKeyType<1>(ind,ind)]=T(1.0);
      id[ukey<2>(0,0)]=O_11;
      id.Dl=1;id.Dr=1;
      contractSingle(compmps[site],BL[site-1],mps[site],id,BR[site+1],mps.getItoKey(site));
      //conv=1.0-computeOverlap(mps,compmps)/(n1*n1);
      //prepareSite(mps,site,"l",true);
      prepareSite(compmps,site,"l",true);
      BL[site]=addContractionLayer(BL[site-1],mps[site],id,compmps[site],mps.getItoKey(site),1);
      //BL[site]=addContractionLayer(BL[site-1],compmps[site],id,mps[site],mps.getItoKey(site),1);
    }
    for(int site=(compmps.size()-1);site>0;--site){
      MPOMat<1,T> id;
      SparseOperator<1,T> O_11(mps.getPSpace()(site),mps.getPSpace()(site));
      for(uint ind=0;ind<mps.getPSpace()(site);++ind)
	O_11[OpKeyType<1>(ind,ind)]=T(1.0);
      id[ukey<2>(0,0)]=O_11;
      id.Dl=1;id.Dr=1;
      contractSingle(compmps[site],BL[site-1],mps[site],id,BR[site+1],mps.getItoKey(site));
      //      cout<<"#overlap: "<<1.0-computeOverlap(mps,compmps)/(n1*n1)<<endl;
      //prepareSite(mps,site,"r",true);
      prepareSite(compmps,site,"r",true);
      BR[site]=addContractionLayer(BR[site+1],mps[site],id,compmps[site],mps.getItoKey(site),-1);
      //BR[site]=addContractionLayer(BR[site+1],compmps[site],id,mps[site],mps.getItoKey(site),-1);
    }

    //    conv=T(1.0)-computeOverlap(compmps,mps)/(sqrt(computeOverlap(mps,mps))*sqrt(computeOverlap(compmps,compmps)));
    n2=sqrt(computeOverlap(compmps,compmps));
    conv=(std::real(n1)*std::real(n1)+std::real(n2)*std::real(n2)-2*abs(computeOverlap(compmps,mps)))/(std::real(n1)*std::real(n1));
    //    conv=1.0-abs(computeOverlap(compmps,mps))/abs((n1*sqrt(computeOverlap(compmps,compmps))));
    if(conv_old<conv){
      if(verbose)
	cout<<"reached limit of convergence; stopping compression;"<<endl;
      break;
    }
    conv_old=conv;

    if(verbose)
      cout<<conv<<"  ";
    if(step>=maxstep||abs(conv)<delta){
      converged=true;
    }
    ++step;
  }
  //THE LAST SITE!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  {
    uint site =0;
    MPOMat<1,T> id;
    SparseOperator<1,T> O_11(mps.getPSpace()(site),mps.getPSpace()(site));
    for(uint ind=0;ind<mps.getPSpace()(site);++ind)
      O_11[OpKeyType<1>(ind,ind)]=T(1.0);
    id[ukey<2>(0,0)]=O_11;
    id.Dl=1;id.Dr=1;
    contractSingle(compmps[site],BL[site-1],mps[site],id,BR[site+1],mps.getItoKey(site));
  }

  //this is a test: it should not be needed!
  for(uint site = 0;site<compmps.size()-1;++site){
    prepareSite(compmps,site,"l",true);
  }
  for(uint site = compmps.size()-1;site>0;--site){
    prepareSite(compmps,site,"r",true);
  }

  uint chim=0;
  if(verbose)
    cout<<endl;
  for(uint s=0;s<compmps.size();++s)
    chim=chim>compmps.chi(s)?chim:compmps.chi(s);
  cout.precision(10);
  n2=sqrt(computeOverlap(compmps,compmps));
  conv=(std::real(n1)*std::real(n1)+std::real(n2)*std::real(n2)-2*abs(computeOverlap(compmps,mps)))/(std::real(n1)*std::real(n1));
  //conv=1.0-abs(computeOverlap(compmps,mps))/abs(n1*sqrt(computeOverlap(compmps,compmps)));
  // if(verbose)
  //   cout<<"#rel. trun. weight = "<<conv;
  mps.clear();
  mps=compmps;
  compmps.setLOrtho(false);
  compmps.setROrtho(true);
  return conv;
  //  cout<<"#done ..."<<endl;
}



//this does a variational single site compression; it conserves the norm of mps;
//mps has structure MBB...B after compress; M will not be normalized to identity, but it will
//be proportional to identity
template<typename T,class KeyType>
T compressTest(MPSBase<KeyType,T> &mps,const uint chimax=100,const uint maxstep=5,const Real delta=1e-10,const bool verbose=true){
  // if(verbose)
  //   cout<<"#compressing to chimax "<<chimax<<": ";
  //OK, now you really have to be careful when using this function; the state has to be PROPERLY normalized ...
  //cout<<"entering compress: "<<endl;
  if(verbose)
    cout<<"preparing input mps...."<<endl;
  for(uint site = 0;site<mps.size()-1;++site){
    prepareSite(mps,site,"l",true);
  }
  for(uint site = mps.size()-1;site>0;--site){
    prepareSite(mps,site,"r",true);
  }
  if(verbose)
    cout<<"done"<<endl;
  //  cout<<"#finished orthonormalizing"<<endl;
  MPSBase<KeyType,T> compmps=mps;
  if(verbose)
    cout<<"initializing mps ...."<<endl;
  // cout<<"                ||compmps|| before tuncate: "<<sqrt(computeOverlap(compmps,compmps))<<endl;
  // cout<<"                ||mps|| before tuncate:     "<<sqrt(computeOverlap(mps,mps))<<endl;
  truncate(compmps,chimax,delta);
  // cout<<"                ||compmps|| AFTER tuncate:  "<<sqrt(computeOverlap(compmps,compmps))<<endl;
  for(uint site = 0;site<mps.size()-1;++site){
    prepareSite(compmps,site,"l",true);
  }
  for(uint site = mps.size()-1;site>0;--site){
    prepareSite(compmps,site,"r",true);
  }
  if(verbose)
    cout<<"done"<<endl;
  const T n1=sqrt(computeOverlap(mps,mps));
  //cout<<"                ||compmps|| AFTER tuncate and AFTER prepareSite (lr and rl sweep): "<<sqrt(computeOverlap(compmps,compmps))<<endl;
  T n2=0;
  map_container<int,BlockMatrix<1,KeyType,T> >BL,BR;
  BL[-1] = getLMin(mps);
  BR[mps.size()]=getRMax(mps);
  bool converged = maxstep==0?true:false;
  if(verbose)
    cout<<"initializing reduced blocks ...."<<endl;
  if(!converged){
    for(uint site=compmps.size()-1;site>0;--site){
      MPOMat<1,T> id;
      SparseOperator<1,T> O_11(mps.getPSpace()(site),mps.getPSpace()(site));
      for(uint ind=0;ind<mps.getPSpace()(site);++ind)
	O_11[OpKeyType<1>(ind,ind)]=T(1.0);
      id[ukey<2>(0,0)]=O_11;
      id.Dl=1;id.Dr=1;
      //the right order of contracting the matrices into the network is not always straight forward. I put both orders here, for fast switching
      BR[site]=addContractionLayer(BR[site+1],mps[site],id,compmps[site],mps.getItoKey(site),-1);
      //BR[site]=addContractionLayer(BR[site+1],compmps[site],id,mps[site],mps.getItoKey(site),-1);
    }
  }
  if(verbose)
    cout<<"done"<<endl;
  //got all BRs ...now start compressing

  uint step=1;
  Real conv,conv_old=1.0;
  while(!converged){
    if(verbose){
      cout<<".";
      cout.flush();
    }
    for(int site=0;site<(compmps.size()-1);++site){
      MPOMat<1,T> id;
      SparseOperator<1,T> O_11(mps.getPSpace()(site),mps.getPSpace()(site));
      for(uint ind=0;ind<mps.getPSpace()(site);++ind)
	O_11[OpKeyType<1>(ind,ind)]=T(1.0);
      id[ukey<2>(0,0)]=O_11;
      id.Dl=1;id.Dr=1;
      contractSingle(compmps[site],BL[site-1],mps[site],id,BR[site+1],mps.getItoKey(site));
      //conv=1.0-computeOverlap(mps,compmps)/(n1*n1);
      prepareSite(mps,site,"l",true);
      prepareSite(compmps,site,"l",true);
      BL[site]=addContractionLayer(BL[site-1],mps[site],id,compmps[site],mps.getItoKey(site),1);
      //BL[site]=addContractionLayer(BL[site-1],compmps[site],id,mps[site],mps.getItoKey(site),1);
    }
    for(int site=(compmps.size()-1);site>0;--site){
      MPOMat<1,T> id;
      SparseOperator<1,T> O_11(mps.getPSpace()(site),mps.getPSpace()(site));
      for(uint ind=0;ind<mps.getPSpace()(site);++ind)
	O_11[OpKeyType<1>(ind,ind)]=T(1.0);
      id[ukey<2>(0,0)]=O_11;
      id.Dl=1;id.Dr=1;
      contractSingle(compmps[site],BL[site-1],mps[site],id,BR[site+1],mps.getItoKey(site));
      //      cout<<"#overlap: "<<1.0-computeOverlap(mps,compmps)/(n1*n1)<<endl;
      prepareSite(mps,site,"r",true);
      prepareSite(compmps,site,"r",true);
      BR[site]=addContractionLayer(BR[site+1],mps[site],id,compmps[site],mps.getItoKey(site),-1);
      //BR[site]=addContractionLayer(BR[site+1],compmps[site],id,mps[site],mps.getItoKey(site),-1);
    }
    //    conv=T(1.0)-computeOverlap(compmps,mps)/(sqrt(computeOverlap(mps,mps))*sqrt(computeOverlap(compmps,compmps)));
    n2=sqrt(computeOverlap(compmps,compmps));
    conv=(n1*n1+n2*n2-2*abs(computeOverlap(compmps,mps)))/(n1*n1);
    //    conv=1.0-abs(computeOverlap(compmps,mps))/abs((n1*sqrt(computeOverlap(compmps,compmps))));
    if(conv_old<conv){
      if(verbose)
	cout<<"reached limit of convergence; stopping compression;"<<endl;
      break;
    }
    conv_old=conv;

    if(verbose)
      cout<<conv<<"  ";
    if(step>=maxstep||abs(conv)<delta){
      converged=true;
    }
    ++step;
  }
  //damn, i forgot the last site!!!!!!!!!!!!
  {
    uint site =0;
    MPOMat<1,T> id;
    SparseOperator<1,T> O_11(mps.getPSpace()(site),mps.getPSpace()(site));
    for(uint ind=0;ind<mps.getPSpace()(site);++ind)
      O_11[OpKeyType<1>(ind,ind)]=T(1.0);
    id[ukey<2>(0,0)]=O_11;
    id.Dl=1;id.Dr=1;
    contractSingle(compmps[site],BL[site-1],mps[site],id,BR[site+1],mps.getItoKey(site));
  }

  //this is a test: it should not be needed!
  for(uint site = 0;site<compmps.size()-1;++site){
    prepareSite(compmps,site,"l",true);
  }
  for(uint site = compmps.size()-1;site>0;--site){
    prepareSite(compmps,site,"r",true);
  }

  uint chim=0;
  if(verbose)
    cout<<endl;
  for(uint s=0;s<compmps.size();++s)
    chim=chim>compmps.chi(s)?chim:compmps.chi(s);
  cout.precision(10);
  n2=sqrt(computeOverlap(compmps,compmps));
  conv=(n1*n1+n2*n2-2*abs(computeOverlap(compmps,mps)))/(n1*n1);
  //conv=1.0-abs(computeOverlap(compmps,mps))/abs(n1*sqrt(computeOverlap(compmps,compmps)));
  // if(verbose)
  //   cout<<"#rel. trun. weight = "<<conv;
  mps.clear();
  mps=compmps;
  compmps.setLOrtho(false);
  compmps.setROrtho(true);
  return conv;
  //  cout<<"#done ..."<<endl;
}


//this does application of an mpo + variational single site compression; it conserves the norm of mps;
//mps has structure MBB...B after compress; M will not be normalized to identity, but it will
//be proportional to identity
template<typename T,class KeyType>
T compressSingleSite(MPS<KeyType,T> &mps,MPO<T> mpo,const uint chimax=100,const uint maxstep=5,const Real delta=1e-10,const bool verbose=true){

  //initial guess is state before mpo-application
  MPS<KeyType,T> tmps,compmps;

  //get the norm of the initial state (mpo|mps>)
  if(verbose){
    cout<<"computing norm of full mps ...."<<endl;
  }
  tmps=applyMPO(mpo,mps);
  while(!checkMPS(tmps));
  const T n1=sqrt(computeOverlap(tmps,tmps));  
  if(verbose)
    cout<<"done"<<endl;

  tmps.clear();
  if(verbose)
    cout<<"initializing mps ...."<<endl;
  tmps=mps;
  uint D=mpo[1].Dr;
  uint chi=floor(double(chimax)/double(D));
  for(uint site = 0;site<tmps.size()-1;++site){
    prepareSite(tmps,site,"l",true);
  }
  for(uint site = tmps.size()-1;site>0;--site){
    prepareSite(tmps,site,"r",true);
  }
  truncate(tmps,chi,delta);
  compmps=applyMPO(mpo,tmps);
  tmps.clear();
  for(uint site = 0;site<compmps.size()-1;++site){
    prepareSite(compmps,site,"l",true);
  }
  for(uint site = compmps.size()-1;site>0;--site){
    prepareSite(compmps,site,"r",true);
  }
  if(verbose)
    cout<<"done"<<endl;
  // orthonormalize(compmps,"lr"); 
  // orthonormalize(compmps,"rl");
  map_container<int,BlockMatrix<1,KeyType,T> >BL,BR;
  BL[-1] = getLMin(mps);
  BR[mps.size()]=getRMax(mps);
  bool converged = maxstep==0?true:false;
  if(verbose)
    cout<<"initializing reduced blocks ...."<<endl;
  if(!converged){
    for(uint site=compmps.size()-1;site>0;--site){
      //BR[site]=addContractionLayer(BR[site+1],compmps[site],mpo[site],mps[site],mps.getItoKey(site),-1);
      BR[site]=addContractionLayer(BR[site+1],mps[site],mpo[site],compmps[site],mps.getItoKey(site),-1);
    }
  }
  if(verbose)
    cout<<"done"<<endl;
  //got all BRs ...now start compressing

  BlockMatrix<1,KeyType,T> mat;
  uint step=1;
  Real conv,conv_old=1.0;
  T overlap;
  
  while(!converged){
    if(verbose){
      cout<<".";
      cout.flush();
    }
    for(int site=0;site<(compmps.size()-1);++site){
      contractSingle(compmps[site],BL[site-1],mps[site],mpo[site],BR[site+1],mps.getItoKey(site));
      prepareSite(compmps,site,"l",true);
      //BL[site]=addContractionLayer(BL[site-1],compmps[site],mpo[site],mps[site],mps.getItoKey(site),1);
      BL[site]=addContractionLayer(BL[site-1],mps[site],mpo[site],compmps[site],mps.getItoKey(site),1);
    }
    for(int site=(compmps.size()-1);site>0;--site){
      if(site==1)
	mat=compmps[site];
      contractSingle(compmps[site],BL[site-1],mps[site],mpo[site],BR[site+1],mps.getItoKey(site));
      if(site==1){
	overlap=mat*compmps[1];
	mat.clear();
      }
      prepareSite(compmps,site,"r",true);
      //BR[site]=addContractionLayer(BR[site+1],compmps[site],mpo[site],mps[site],mps.getItoKey(site),-1);
      BR[site]=addContractionLayer(BR[site+1],mps[site],mpo[site],compmps[site],mps.getItoKey(site),-1);
    }
    T n2=sqrt(computeOverlap(compmps,compmps));
    conv=(n1*n1+n2*n2-2*abs(computeOverlap(compmps,mps)))/(n1*n1);
    //conv=1.0-abs(overlap)/(n1*sqrt(computeOverlap(compmps,compmps)));
    if(conv_old<conv){
      if(verbose)
	cout<<"reached limit of convergence; stopping compression;"<<endl;
      break;
    }
    conv_old=conv;
    if(verbose)
      cout<<conv<<"  ";
  
    if(step>=maxstep||abs(conv)<delta){
      converged=true;
    }
    ++step;
  }
  uint chim=0;

  for(uint s=0;s<compmps.size();++s)
    chim=chim>compmps.chi(s)?chim:compmps.chi(s);
  cout.precision(10);
  T n2=sqrt(computeOverlap(compmps,compmps));
  conv=(n1*n1+n2*n2-2*abs(computeOverlap(compmps,mps)))/(n1*n1);

  //  conv=1.0-abs(overlap)/(n1*sqrt(computeOverlap(compmps,compmps)));
  if(verbose)
    cout<<endl;

  mps=compmps;
  compmps.setLOrtho(false);
  compmps.setROrtho(true);
  return conv;
  //  cout<<"#done ..."<<endl;
}




//this does application of an mpo + variational two site compression; it conserves the norm of mps;
//mps has structure MBB...B after compress; M will not be normalized to identity, but it will
//be proportional to identity
template<typename T,class KeyType>
T compressTwoSite(MPS<KeyType,T> &mps,MPO<T> mpo,const uint chimax=100,const uint maxstep=5,const Real delta=1e-10,const bool doTrunc=true,const Real tr_weight=1e-10,
		  const bool verbose=true){

  //initial guess is state before mpo-application
  MPS<KeyType,T> tmps,compmps;

  //get the norm of the initial state (mpo|mps>): first apply mpo to mps:
  if(verbose)
    cout<<"computing norm of full mps ...."<<endl;
  tmps=applyMPO(mpo,mps);
  while(!checkMPS(tmps));
  //then get norm of the resulting state
  const T n1=sqrt(computeOverlap(tmps,tmps));  
  T n2=0;
  if(verbose)
    cout<<"done"<<endl;
  //clear it!!!! it's huge!!!
  tmps.clear();
  //get a good initial state: first cut mps down to very small chi...
  if(verbose)
    cout<<"initializing mps ...."<<endl;
  tmps=mps;
  uint D=mpo[1].Dr;
  uint chi=floor(double(chimax)/double(D));
  for(uint site = 0;site<tmps.size()-1;++site){
    prepareSite(tmps,site,"l",true);
  }
  for(uint site = tmps.size()-1;site>0;--site){
    prepareSite(tmps,site,"r",true);
  }
  truncate(tmps,chi,delta);
  //then apply mpo to this truncated state. result is a (good) initial state with small bond dimension
  compmps=applyMPO(mpo,tmps);
  tmps.clear();

  for(uint site = 0;site<compmps.size()-1;++site){
    prepareSite(compmps,site,"l",true);
  }
  for(uint site = compmps.size()-1;site>0;--site){
    prepareSite(compmps,site,"r",true);
  }
  // orthonormalize(compmps,"lr"); 
  // orthonormalize(compmps,"rl");
  if(verbose)
    cout<<"done"<<endl;
  map_container<int,BlockMatrix<1,KeyType,T> >BL,BR;
  BL[-1] = getLMin(mps);
  BR[mps.size()]=getRMax(mps);
  bool converged = maxstep==0?true:false;
  if(verbose)
    cout<<"initializing reduced blocks ...."<<endl;
  if(!converged){
    for(uint site=compmps.size()-1;site>1;--site){
      //      BR[site]=addContractionLayer(BR[site+1],compmps[site],mpo[site],mps[site],mps.getItoKey(site),-1);
      BR[site]=addContractionLayer(BR[site+1],mps[site],mpo[site],compmps[site],mps.getItoKey(site),-1);
    }
  }
  if(verbose)
    cout<<"done"<<endl;
  //got all BRs ...now start compressing

  BlockMatrix<2,KeyType,T> mat,out,in;
  uint step=1;
  Real conv,conv_old=1.0;
  T overlap;
  std::vector<i_to_key<KeyType> >itokey(2);
  while(!converged){
    if(verbose){
      cout<<".";
      cout.flush();
    }
    for(int site=0;site<(mps.size()-2);++site){
      BlockMatrix<2,KeyType,T> in(mps[site],mps[site+1]);
      itokey[0]=mps.getItoKey(site);
      itokey[1]=mps.getItoKey(site+1);
      contractTwosite(out,BL[site-1],BR[site+2],in,mpo[site],mpo[site+1],itokey);
      std::vector<BlockMatrix<1,KeyType,T> >outsplit;
      std::vector<DiagBlockMatrix<KeyType,Real> > lambda;
      pair<std::vector<Real>,std::vector<uint> >sre=split(out,outsplit,lambda,itokey,0,chimax,doTrunc,tr_weight,false);
      compmps[site]=outsplit[0];
      compmps[site+1]=lambda[0]*outsplit[1];
      //      BL[site]=addContractionLayer(BL[site-1],compmps[site],mpo[site],mps[site],mps.getItoKey(site),1);      
      BL[site]=addContractionLayer(BL[site-1],mps[site],mpo[site],compmps[site],mps.getItoKey(site),1);
      out.clear();
      in.clear();
    }
    for(int site=(mps.size()-2);site>0;--site){
      BlockMatrix<2,KeyType,T> in(mps[site],mps[site+1]);
      itokey[0]=mps.getItoKey(site);
      itokey[1]=mps.getItoKey(site+1);
      contractTwosite(out,BL[site-1],BR[site+2],in,mpo[site],mpo[site+1],itokey);
      if(site==1){
	mat=BlockMatrix<2,KeyType,T>(compmps[site],compmps[site+1]);
	overlap=mat*out;
	mat.clear();
      }
      std::vector<BlockMatrix<1,KeyType,T> >outsplit;
      std::vector<DiagBlockMatrix<KeyType,Real> > lambda;
      pair<std::vector<Real>,std::vector<uint> >sre=split(out,outsplit,lambda,itokey,0,chimax,doTrunc,tr_weight,false);
      compmps[site]=outsplit[0]*lambda[0];
      compmps[site+1]=outsplit[1];
      //      BR[site+1]=addContractionLayer(BR[site+2],compmps[site+1],mpo[site+1],mps[site+1],mps.getItoKey(site+1),-1);
      BR[site+1]=addContractionLayer(BR[site+2],mps[site+1],mpo[site+1],compmps[site+1],mps.getItoKey(site+1),-1);
      out.clear();
      in.clear();
    }

    n2=sqrt(computeOverlap(compmps,compmps));
    conv=(n1*n1+n2*n2-2*abs(computeOverlap(compmps,mps)))/(n1*n1);
    //conv=1.0-abs(overlap)/(n1*sqrt(computeOverlap(compmps,compmps)));
    if(conv_old<conv){
      if(verbose)
	cout<<"reached limit of convergence; stopping compression;"<<endl;
      break;
    }
    conv_old=conv;

    if(verbose)
      cout<<conv<<"  ";
  
    if(step>=maxstep||abs(conv)<delta){
      converged=true;
    }
    ++step;
  }
  uint chim=0;

  for(uint s=0;s<compmps.size();++s)
    chim=chim>compmps.chi(s)?chim:compmps.chi(s);
  cout.precision(10);
  n2=sqrt(computeOverlap(compmps,compmps));
  conv=(n1*n1+n2*n2-2*abs(computeOverlap(compmps,mps)))/(n1*n1);
  //conv=1.0-abs(overlap)/(n1*sqrt(computeOverlap(compmps,compmps)));
  if(verbose)
    cout<<endl;

  mps=compmps;
  compmps.setLOrtho(false);
  compmps.setROrtho(true);
  return conv;
  //  cout<<"#done ..."<<endl;
}
//constructor for product states where each site is in a superposition of states given by the list entries in local-state; e.g. if at site n
//the list has entries 0,1,2,3, then the local state at site n is |0>+|1>+|2>+|3> (no normalization!!!).
template<class KeyType,typename T>
MPS<KeyType, T > FermionicSuperOperatorPerturbation(uint operatorsite,uint ostate,UintVectorType dimSiteSpace,UintVectorType localState1,UintVectorType localState2,std::vector<i_to_key<KeyType> > itokey)
{
  MPS<KeyType,T> mps;
  int N=localState1.size();
  mps.resize(N);
  mps.setItoKey(itokey);
  mps.setPSpace(dimSiteSpace);
  matrix<T,column_major> mat(1,1);
  mat(0,0) = T(1);


  KeyType keyl=itokey[0].begin()->second,keyr;
  reset(keyl);
  TensorKeyType<2,1,KeyType>TK(keyl,keyl,0);
  mps._leftDummy[TK] = mat;
  //this code works only if the local quantum numbers associated with localState1(site) and localState2(site) are the same!!!
  for(int site=0;site<N;++site){
    mps[site].clear();
    if(site<operatorsite){
      if(site==0){reset(keyl);
      }else if(site<N){keyl = keyr;}
      keyr=keyl+itokey[site][localState1[site]];
      TensorKeyType<2,1,KeyType>Tkey1(keyl,keyr,localState1(site));
      TensorKeyType<2,1,KeyType>Tkey2(keyl,keyr,localState2(site));
      mps[site][Tkey1]= mat;
      mps[site][Tkey2]= mat;
    }else if(site==operatorsite){
      if(site<N){keyl = keyr;}
      keyr=keyl+itokey[site][ostate];
      TensorKeyType<2,1,KeyType>Tkey(keyl,keyr,ostate);
      mps[site][Tkey]= mat;
    }else if(site>operatorsite){
      if(site<N){keyl = keyr;}
      keyr=keyl+itokey[site][localState1[site]];
      TensorKeyType<2,1,KeyType>Tkey1(keyl,keyr,localState1(site));
      TensorKeyType<2,1,KeyType>Tkey2(keyl,keyr,localState2(site));
      mps[site][Tkey1]= mat;
      mps[site][Tkey2]= mat;

    }
  }
  TensorKeyType<2,1,KeyType>TKm(keyr,keyr,0);
  mps._rightDummy[TKm] = mat;
  return mps;
}
template<class KeyType,typename T>
MPS<KeyType, T > FermionicSuperOperatorIdentity(UintVectorType dimSiteSpace,UintVectorType localState1,UintVectorType localState2,std::vector<i_to_key<KeyType> > itokey)
{
  MPS<KeyType,T> mps;
  int N=localState1.size();
  mps.resize(N);
  mps.setItoKey(itokey);
  mps.setPSpace(dimSiteSpace);
  matrix<T,column_major> mat(1,1);
  mat(0,0) = T(1);

  KeyType keyl=itokey[0].begin()->second,keyr;
  reset(keyl);
  TensorKeyType<2,1,KeyType>TK(keyl,keyl,0);
  mps._leftDummy[TK] = mat;

  for(int site=0;site<N;++site){
    mps[site].clear();
    if(site==0){reset(keyl);
    }else if(site<N){keyl = keyr;}
    //this code works only if the local quantum numbers associated with localState1(site) and localState2(site) are the same!!!
    keyr=keyl+itokey[site][localState1[site]];
    TensorKeyType<2,1,KeyType>Tkey1(keyl,keyr,localState1(site));
    TensorKeyType<2,1,KeyType>Tkey2(keyl,keyr,localState2(site));
    mps[site][Tkey1]= mat;
    mps[site][Tkey2]= mat;
  }
  TensorKeyType<2,1,KeyType>TKm(keyr,keyr,0);
  mps._rightDummy[TKm] = mat;
  return mps;
}


template<class KeyType,typename T>
MPS<KeyType, T > FermionicSuperOperator(UintVectorType dimSiteSpace,UintVectorType localState1,UintVectorType localState2,std::vector<i_to_key<KeyType> > itokey)
{
  MPS<KeyType,T> mps;

  int N=localState1.size();
  mps.resize(N);
  mps.setItoKey(itokey);
  mps.setPSpace(dimSiteSpace);
  matrix<T,column_major> mat(1,1);
  mat(0,0) = T(1);

  KeyType keyl=itokey[0].begin()->second,keyr;
  reset(keyl);
  TensorKeyType<2,1,KeyType>TK(keyl,keyl,0);
  mps._leftDummy[TK] = mat;

  for(int site=0;site<N;++site){
    mps[site].clear();
    if(site==0){reset(keyl);
    }else if(site<N){keyl = keyr;}
    //this code works only if the local quantum numbers associated with localState1(site) and localState2(site) are the same!!!
    keyr=keyl+itokey[site][localState1[site]];
    TensorKeyType<2,1,KeyType>Tkey1(keyl,keyr,localState1(site));
    TensorKeyType<2,1,KeyType>Tkey2(keyl,keyr,localState2(site));
    mps[site][Tkey1]= mat;
    mps[site][Tkey2]= mat;
  }
  TensorKeyType<2,1,KeyType>TKm(keyr,keyr,0);
  mps._rightDummy[TKm] = mat;
  return mps;
}


//currently only works with no quantum numbers!!!
template<class KeyType,typename T>
MPO<T> MPStoMPO(MPS<KeyType,T>&mps,UintVectorType opdim){
  MPO<T>mpo;
  uint N=mps.size();
  for(uint site=0;site<N;++site){
    MPOMat<1,T> M;
    typename BlockMatrix<1,KeyType,T>::iterator it;
    std::vector<uint> dims(2);
    dims[0]=opdim(site);
    dims[1]=opdim(site);
    uint chi1=site>0?mps.getAuxSpace()(site-1):1;
    uint chi2=site<(N-1)?mps.getAuxSpace()(site):1;
    for(uint c1=0;c1<chi1;++c1){
      for(uint c2=0;c2<chi2;++c2){
	SparseLocalOperator<T> op(dims[0],dims[1]);
	op.setSite(site);
	for(it=mps[site].begin();it!=mps[site].end();++it){
	  uint s=it->first(0);//that is the physical index of the Blockmatrix
	  // std::vector<uint> inds(2);
	  // if(s==0)
	  //   inds[0]=0;inds[1]=0;
	  // if(s==1)
	  //   inds[0]=1;inds[1]=0;
	  // if(s==2)
	  //   inds[0]=0;inds[1]=1;
	  // if(s==3)
	  //   inds[0]=1;inds[1]=1;

	  std::vector<uint> inds=SiToMuIndex(s,dims),i2(2);
	  if(s==1||s==2){
	    i2[0]=abs(1-inds[0]);
	    i2[1]=abs(1-inds[1]);
	    // cout<<"inds: mapping "<<s<<" to "<<inds<<endl;
	    // cout<<"i2:   mapping "<<s<<" to "<<i2<<endl;
	    // cin.get();

	    inds=i2;
	  }
	  op[OpKeyType<1>(inds[0],inds[1])]=it->second(c1,c2);
	}
	M[ukey<2>(c1,c2)]=op;
      }
    }
    M.squeeze(1e-15);
    M.Dl=chi1;M.Dr=chi2;
    mpo[site]=M;
  }
  return mpo;


  
}


//this is a dangerous routine!!! it works only for states with no quantum numbers;
template<class KeyType,typename T>
void grandomize(MPS<KeyType,T>&mps,Real mean,Real var){
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::normal_distribution<double> nd(mean,var);
  std::default_random_engine generator (seed);
  //std::default_random_engine generator;

  for(uint s=0;s<mps.size();++s){
    typename BlockMatrix<1,KeyType,T>::iterator it=mps[s].begin();
    //VectorType v(it->second.size1()*it->second.size2());
    //for(uint n=0;n<it->second.size1()*it->second.size2();++n)
    //v[n]=nd(generator);
    for(it=mps[s].begin();it!=mps[s].end();++it){
      for(uint n=0;n<it->second.size1()*it->second.size2();++n){
      //if(it->second.data()[n]>0.0)
      //it->second.data()[n]+=fabs(v[n]);
      //if(it->second.data()[n]<0.0)
      //it->second.data()[n]-=fabs(v[n]);
	it->second.data()[n]+=nd(generator);
      }
    }
  }
}
#endif
