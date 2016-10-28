/*
 * helperfunctions.cpp
 *
 *  Created on: Nov 3, 2009
 *      Author: martinjg
 */
#include "helperfunctions.hpp"

//using namespace BoostTools;
using namespace BoostTypedefs;


Real myFabs(const Real num){return fabs(num);}
Real myFabs(const Complex num){return sqrt(pow(num.real(),2)+pow(num.imag(),2));}
bool isequal(UintVectorType vec1,UintVectorType vec2)
{
	assert(vec1.size()==vec2.size());
	bool comp=true;
	UintVectorType::iterator ita = vec1.begin();
	UintVectorType::iterator itb = vec2.begin();
	while(ita!=vec1.end())
	{
		if(*ita!=*itb){comp = false;break;}
		else{ita++;itb++;}
	}
	return comp;
}
//returns TRUE only if vec1(i)<=vec2(i) for all i
bool VectorCompareLEQ(UintVectorType vec1,UintVectorType vec2)
{
	assert(vec1.size()==vec2.size());
	bool comp=true;
	UintVectorType::iterator ita = vec1.begin();
	UintVectorType::iterator itb = vec2.begin();
	while(ita!=vec1.end())
	{
		if(*ita>*itb){comp = false;break;}
		else{ita++;itb++;}
	}
	return comp;
}
//checks if vec1<vec2
bool VectorCompareLESS(UintVectorType vec1,UintVectorType vec2)
{
	assert(vec1.size()==vec2.size());
	bool comp=true;
	UintVectorType::iterator ita = vec1.begin();
	UintVectorType::iterator itb = vec2.begin();
	while(ita!=vec1.end())
	{
		if(*ita>=*itb){comp = false;break;}
		else{ita++;itb++;}
	}
	return comp;
}

bool VectorCompareLEQ( VectorType vec1,VectorType vec2)
{
	assert(vec1.size()==vec2.size());
	bool comp=true;
	VectorType::iterator ita = vec1.begin();
	VectorType::iterator itb = vec2.begin();
	while(ita!=vec1.end())
	{
		if(*ita>*itb){comp = false;break;}
		else{ita++;itb++;}
	}
	return comp;
}
void print2dboostCarray(const boost::multi_array<Complex,2> &ar)
{using std::endl;
using std::cout;
	uint size1 = ar.shape()[0],size2 = ar.shape()[1];
	cout<<"Ar(:,:)="<<endl<<endl;for(uint x = 0;x<size1;x++){for(uint y=0;y<size2;y++){std::cout<<ar[x][y].real()<<"+1i*"<<ar[x][y].imag()<<" ";}cout<<endl;}cout<<"]"<<endl;
}
void printboostCmatrix(const CMatrixType&ar)
{using std::endl;
using std::cout;
	uint size1 = ar.size1(),size2 = ar.size2();
	cout<<"Ar(:,:)=["<<endl<<endl;for(uint x = 0;x<size1;x++){for(uint y=0;y<size2;y++){std::cout<<ar(x,y).real()<<"+1i*"<<ar(x,y).imag()<<" ";}cout<<endl;}cout<<"]"<<endl;
}
void printboostCmatrix(const CMatrixTypeColMaj&ar)
{using std::endl;
using std::cout;
	uint size1 = ar.size1(),size2 = ar.size2();
	cout<<"Ar(:,:)=["<<endl<<endl;for(uint x = 0;x<size1;x++){for(uint y=0;y<size2;y++){std::cout<<ar(x,y).real()<<"+1i*"<<ar(x,y).imag()<<" ";}cout<<endl;}cout<<"]"<<endl;
}

Real Anorm(const boost::multi_array<Complex,3> &ar)
{
	Real val = 0;
	for(uint cnt = 0;cnt<ar.num_elements();cnt++)
	{
		val+=pow(ar.data()[cnt].real(),2)+pow(ar.data()[cnt].imag(),2);
	}
	return sqrt(val);
}

Real Anorm(const boost::multi_array<Complex,2> &ar)
{
	Real val = 0;
	for(uint cnt = 0;cnt<ar.num_elements();cnt++)
	{
		val+=pow(ar.data()[cnt].real(),2)+pow(ar.data()[cnt].imag(),2);
	}
	return sqrt(val);
}
Complex AScalarProd(const boost::multi_array<Complex,3>&ar1,boost::multi_array<Complex,3>&ar2)
{
	assert(ar1.num_elements()==ar2.num_elements());
	Complex val = 0;
		for(uint cnt = 0;cnt<ar1.num_elements();cnt++)
		{
			val=val+(ar1.data()[cnt])*ar2.data()[cnt];
		}
		return val;
}
Complex AScalarProd(const boost::multi_array<Complex,2>&ar1,boost::multi_array<Complex,2>&ar2)
{
	assert(ar1.num_elements()==ar2.num_elements());
	Complex val = 0;
		for(uint cnt = 0;cnt<ar1.num_elements();cnt++)
		{
			val=val+(ar1.data()[cnt])*ar2.data()[cnt];
		}
		return val;
}


//input: uint x, std::vector<uint> max_occ, pointer to an EMPTY (= no elements) std::vector<uint> &occ
//output: std::vector<uint> &occ
//the index-integer x is mapped onto a multi-index occ = (i_1,i_2,...,i_(N-1),i_N) of size N=max_occ.size().
//Mapping is specified by max_occ by the following equation:
//x = i_N*dim(i_1,...,i_(n-1))+i_(N-1)*dim(i_1,...,i_(N-2))+...+ i_1
//where dim(i_1,...,i_l),... denotes the number of different l-tupels where each  i_j (1<=j<=l) of the l-tupel
//runs from 0 to max_occ[j].
void compute_multi_index(uint x,std::vector<uint> max_occ,std::vector<uint> &occ)
{
	uint x_max=1;
	for(uint cnt=0;cnt<max_occ.size();cnt++)x_max*=max_occ[cnt];
	//cout<<"x_max="<<x_max<<endl;
	assert(x<x_max);
	std::vector<uint>::iterator it = max_occ.end();
	uint dim= 1;
	if(max_occ.size()>1){
	it--;
	it--;
	while(it!=max_occ.begin())
	{
		dim*=*it;
		it--;
	}
	dim*=*it;

	}
	max_occ.pop_back();
	if(max_occ.size()>0){
	compute_multi_index(x%dim,max_occ,occ);
	occ.push_back((x-occ[max_occ.size()-1])/dim);
	}
	else{ occ.push_back(x);}

}

//input: std::vector<uint> &occ, std::vector<uint> max_occ
//output: uint x

//multi index occ = (i_1,i_2,...,i_(N-1),i_N) is mapped onto a single index x according to the equation
//x = i_N*dim(i_1,...,i_(n-1))+i_(N-1)*dim(i_1,...,i_(N-2))+...+ i_1
//where dim(i_1,...,i_l),... denotes the number of different l-tupels where each  i_j (1<=j<=l) of the l-tupel
//runs from 0 to max_occ[j].
uint compute_single_index(std::vector<uint> occ,std::vector<uint> max_occ)
{
	uint x;
	std::vector<uint>::iterator it1,it2;
	if(occ.size()>1){
		it1 = occ.end();
		it2 = max_occ.end();
		it1--;it2--;it2--;
		uint dim=1;
		while(it2!=max_occ.begin()){dim*=(*it2);it2--;}
		dim*=*it2;
		x=*it1*dim;
		max_occ.pop_back();occ.pop_back();
		x+=compute_single_index(occ,max_occ);
	}
	else x=occ[0];
	return x;

}

void generate_next_occupation(std::vector<uint> &occ,std::vector<uint> &max_occ)
{
	bool finish = false,found_empty_site=false;
	uint index = 0;
	while(!finish)
	{
		if(occ[index]<max_occ[index]-1){ occ[index]+=1;finish =true;}
		else {index++; if(index==occ.size())return; if(occ[index]<max_occ[index]-1) found_empty_site = true;}
		if(found_empty_site)for(uint cnt = 0;cnt<index;cnt++)occ[cnt] = 0;

	}
}


void my_pause()
{
	std::string unused_var;
	std::cout<<"press key to continue";
	std::cin>>unused_var;

}

//template<class T>
	std::string itoa(int num,int base,int precision)
	{
		std::string a_string;
		uint i = 0 ;
		int temp = static_cast<int>(floor(num*pow(10.0,precision)));
		uint cnt = 0;
	if(temp!=0)
	{	while(temp!=0)
		{
			i = temp%base;
			//std::cout<<i<<std::endl;
			temp = (temp-i)/10;
			a_string.insert(a_string.begin(),static_cast<char>(i+48));
			cnt++;
			//std::cout<<a_string<<std::endl;
		}

		std::string::iterator it;
		it = a_string.end();
		while(precision>0)
		{
			precision--;
			if (it !=a_string.begin())
			{
				it--;
			}
			else
			{
				a_string.insert(it,'0');
				it--;
			}
		}


		if (it==a_string.begin())
		{

			a_string.insert(it,'.');

			it=a_string.begin();

			a_string.insert(it,'0');
		}
		else
		{
			a_string.insert(it,'.');
		}

	}
	else
	{
	a_string ="0";
	}

		return a_string;
}
	/// calculate Sum A[i]*A+[i] and print results for every site
	/// should give identity matrix in relevant subspace


	// -----------------------------------------
	/// calculate Sum B+[i]*B[i] and print results for every site
	/// should give identity matrix in relevant subspace
