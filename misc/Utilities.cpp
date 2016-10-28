/*
 * Utilities.cpp
 *
 *  Created on: Mar 18, 2009
 *      Author: michael
 */

#include "Utilities.h"
#include <math.h>

#include <iostream>
#include <fstream>
#include <strstream>

#include <boost/numeric/ublas/io.hpp>
#include <boost/lexical_cast.hpp>

using namespace boost::numeric::ublas;

Utilities::Utilities() {

}

Utilities::~Utilities() {

}

double Utilities::roundDouble(double n, int places){
	double t;
	t=n*pow(10,places)-floor(n*pow(10,places));
	if (t>=0.5) {
		  n*=pow(10,places);
		  n=ceil(n);
		  n/=pow(10,places);
	} else {
		  n*=pow(10,places);
		  n=floor(n);
		  n/=pow(10,places);
	}
	return n;
}

vector<double> Utilities::roundDoubleVector(vector<double> v, int places){
	for (uint i=0; i<v.size(); i++){
		v(i)=roundDouble(v(i),places);
	}
	return v;
}

std::vector<double> Utilities::roundDoubleVector(std::vector<double> v, int places){
	for (uint i=0; i<v.size(); i++){
		v[i]=roundDouble(v[i],places);
	}
	return v;
}


vector<double> Utilities::convertStringToVector(std::string mString){
	std::vector<std::string> data = divideStringAtBlanc(mString);
	vector<double> res(data.size()); res.clear();
	for (std::size_t i=0; i<data.size(); i++){
		res(i)=convertStringToNumber(data[i]);
	}
	return res;
}

std::vector<std::string> Utilities::divideStringAtBlanc(std::string mString){
	std::vector<std::string> data;

	while(mString.length()!=0){
		// equal to trim left
		mString.erase(0, mString.find_first_not_of(" "));
		if (mString.find(";")==std::string::npos){
			data.push_back(mString);
			break;
		}else{
			data.push_back(mString.substr(0, mString.find(";")));
			mString.erase(0,mString.find(";")+1);
		}
	}

	return data;
}

matrix<double> Utilities::convertStringToMatrix(std::string mString){
	matrix<double> m;
	std::vector<std::string> data;
	std::vector<std::vector<double> > dData;

	if (mString.length()>0){
		data=divideStringAtBlanc(mString);

		dData.resize(data.size());
		for (uint i=0; i<data.size(); i++){
			std::string line = data[i];
			while(line.length()!=0){
				// equal to trim left
				line.erase(0, line.find_first_not_of(" "));

				// extract string
				uint length=line.find(" ");
				if (length==std::string::npos){
					length=line.length();
				}
				std::string element=line.substr(0, length);

				//std::cout << element << ", " << convertStringToNumber(element) << std::endl;

				dData[i].push_back(convertStringToNumber(element));
				if (line.find(" ")==std::string::npos){
					break;
				}
				line.erase(0, line.find(" ")+1);
			}
		}

		m.resize(dData.size(), dData[0].size());
		for (uint i=0; i<dData.size(); i++){
			for (uint j=0; j<dData[i].size(); j++){
				m(i,j)=dData[i][j];
			}
		}
	}
	return m;
}

double Utilities::convertFullStringToNumber(std::string str){
	// find + or -
	uint indexP=str.find("+",1);
	uint indexM=str.find("-",1);
	uint index=indexP<indexM?indexP:indexM;

	std::vector<std::string> part;
	std::vector<bool> pm;
	pm.push_back(indexP<indexM);

	while (index!=std::string::npos){
		// extract part
		part.push_back(str.substr(0, index));
		//std::cout << *part.end() << std::endl;

		uint indexP=str.find("+",1);
		uint indexM=str.find("-",1);
		index=indexP<indexM?indexP:indexM;
		pm.push_back(indexP<indexM);
	}
}

double Utilities::convertStringToNumber(std::string element){
	double res=0.0;
	std::string e=element;

	try{
		res=boost::lexical_cast<double>(e);
		return res;
	}catch(...){   }

	uint index=e.find("sqrt(");
	std::vector<std::string> sqrtVector;
	std::string firstPart="";
	std::string lastPart="";
	if (index!=std::string::npos){
		firstPart=e.substr(0,index);
		//std::cout << "fp=" << firstPart << std::endl;
		while (index!=std::string::npos){
			uint index2=e.find(")");
			if (index2==std::string::npos){
				throw std::runtime_error("Utilities::convertStringToNumber: Wrong syntax!");
			}
			index+=5;
			sqrtVector.push_back(e.substr(index, index2-index));
			//std::cout << "sqrt=" << e.substr(index, index2-index) << std::endl;
			e=e.erase(0,index2+1);
			index=e.find("sqrt(");
		}
	}
	lastPart=e;
	//std::cout << "lp=" << lastPart << std::endl;

	// first part
	if (firstPart.length()>0){
		// check for /
		res = splitDivString(firstPart);
	}

	// middle part = sqrt part
	for (uint i=0; i<sqrtVector.size(); i++){
		double d=sqrt(splitDivString(sqrtVector[i]));
		if (res==0){
			res=d;
		}else{
			res/=d;
		}
	}

	// last part
	if (lastPart.length()>0){
		if (res==0){
			res=splitDivString(lastPart);
		}else{
			res/=splitDivString(lastPart);
		}
	}

	return res;
}

std::vector<double> Utilities::convertRangeToVector(std::string range){
	uint index=range.find(":");
	double from=boost::lexical_cast<double>(range.substr(0,index));

	range=range.erase(0,index+1);
	index=range.find(":");
	int nSteps=boost::lexical_cast<int>(range.substr(0,index));

	range=range.erase(0,index+1);
	double to=boost::lexical_cast<double>(range);

	double delta=(to-from)/(double)(nSteps-1.);

	std::vector<double> res;
	for (int i=0; i<nSteps; i++){
		res.push_back(from);
		from+=delta;
	}
	return res;
}

double Utilities::splitDivString(std::string element){
	uint index=element.find("/");

	// remove / if it is at first position
	if (index==0){
		element=element.erase(0,index+1);
		index=element.find("/");
	}

	std::vector<std::string> res;
	while (index!=std::string::npos){
		std::string part;
		part=element.substr(0,index);
		element.erase(0,index+1);
		res.push_back(part);
		index=element.find("/");
	}
	std::string part=element;
	if (part.size()!=0){
		res.push_back(part);
	}

	double d=0.0;
	if (res.size()>0){
		d=boost::lexical_cast<double>(res[0]);
		for (uint i=1; i<res.size(); i++){
			d/=boost::lexical_cast<double>(res[i]);
		}

	}
	return d;
}

int Utilities::pointInPolygon(matrix<double> vertices, vector<double> x){
	int nvert=vertices.size2();
	int i, j, c = 0;
	for (i = 0, j = nvert-1; i < nvert; j = i++) {
		if ( ((vertices(1,i)>x(1)) != (vertices(1,j)>x(1))) &&
				(x(0) < (vertices(0,j)-vertices(0,i)) * (x(1)-vertices(1,i)) / (vertices(1,j)-vertices(1,i)) + vertices(0,i)) ){
			c = !c;
		}
	}
	return c;
}


double Utilities::nChooseK(double n, double k){
	double res = 1;
	double multiplier = n;
	double divisor = k<n-k?k:n-k;

    while (divisor > sqrt(std::numeric_limits<double>::epsilon())) {
        res *= multiplier / divisor;
        multiplier--;
        divisor--;
    }
    return res;
}

int Utilities::calculateSign(const std::vector<bool> &state, int i, int j){
	return countParticlesBetween(state, i, j)%2==0?1:-1;
}

int Utilities::countParticlesFromTo(const std::vector<bool> &state, int i, int j)
{
	if (i==j) return 0;
    int num=0;
	for (int p=i; p<j; p++){
	   if (state[p]){
			num++;
	   }
	}
    return num;
}

int Utilities::countParticlesBetween(const std::vector<bool> &state, int i, int j){
    if (j<i){
        uint temp=i;
        i=j;
        j=temp;
    }
    return countParticlesFromTo(state,i+1,j);
}
