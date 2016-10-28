/*
 * File.cpp
 *
 *  Created on: Dec 5, 2008
 *      Author: michael
 */


#include "File.h"

#include <fstream>
#include <strstream>
//#include "/home/maga/C++/boost_1_49_0/boost/filesystem.hpp"
//#include "/home/maga/C++/boost_1_49_0/boost/filesystem/v3/operations.hpp"
//#include "/home/maga/C++/boost_1_49_0/boost/filesystem/v3/path.hpp"
#include <boost/filesystem.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include "../misc/Utilities.h"


/**
 * @subfolder folder with this filename will be created
 */
File::File(std::string pathname, std::string subfolder, std::string filename){
	prec=Utilities::ouputPrecission;
	if (pathname[pathname.size() - 1] != '/') pathname += '/';
	if (boost::filesystem::is_directory(boost::filesystem::path(pathname))){
		if (subfolder.length()>0){
			if (subfolder[subfolder.size() - 1] != '/') subfolder += '/';
		}else{
			int folderNumber=0;
			int maxFolderNumber=0;

			for (boost::filesystem::directory_iterator itr(pathname); itr!=boost::filesystem::directory_iterator(); ++itr){
				if (boost::filesystem::is_directory(itr->path())){
					std::stringstream p(itr->path().leaf().native(),std::istringstream::in);
					p >> folderNumber;
					if (folderNumber>maxFolderNumber){
						maxFolderNumber=folderNumber;
					}
				}
			}
			std::stringstream out;
			out << ++maxFolderNumber << '/';
			subfolder = out.str();
		}
		std::string fullpath=pathname+subfolder;
		boost::filesystem::create_directory(fullpath);
		this->pathname=fullpath;
	}else{
		std::cerr << "Not found: " << pathname << std::endl;
	}
	this->filename=filename;
}

/**
 * @subfolder folder with this filename will be created
 */
File::File(std::string pathname, std::string subfolder){
	prec=Utilities::ouputPrecission;
	std::cout << "path="<<pathname << ", sub=" << subfolder << std::endl;
	if (pathname[pathname.size() - 1] != '/') pathname += "/";
	std::cout << "path="<<pathname << ", sub=" << subfolder << std::endl;
	boost::filesystem::path p; p/=pathname;
	std::cout << boost::filesystem::is_directory(boost::filesystem::path(pathname))
			<< boost::filesystem::exists(boost::filesystem::path("/afs/itp.tugraz.at/user/knami/Data/imp/")) << std::endl;
	if (boost::filesystem::is_directory(p)){
		if (subfolder.length()>0){
			if (subfolder[subfolder.size() - 1] != '/') subfolder += '/';
		}else{
			int folderNumber=0;
			int maxFolderNumber=0;

			for (boost::filesystem::directory_iterator itr(pathname); itr!=boost::filesystem::directory_iterator(); ++itr){
				if (boost::filesystem::is_directory(itr->path())){
					std::stringstream p(itr->path().leaf().native(),std::istringstream::in);
					p >> folderNumber;
					if (folderNumber>maxFolderNumber){
						maxFolderNumber=folderNumber;
					}
				}
			}
			std::stringstream out;
			out << ++maxFolderNumber << '/';
			subfolder = out.str();
		}
		std::cout << "path="<<pathname << ", sub=" << subfolder << std::endl;
		std::string fullpath=pathname+subfolder;
		boost::filesystem::create_directory(fullpath);
		this->pathname=fullpath;
	}else{
		std::cerr << "Not found: " << pathname << std::endl;
	}
}

/**
 * @subfolder folder with this filename will be created
 */
File::File(std::string pathname){
	prec=Utilities::ouputPrecission;
	if (boost::filesystem::exists(pathname)){
		std::string fullpath=pathname;
		this->pathname=fullpath;
	}else{
		std::cerr << "Not found: " << pathname << std::endl;
	}
}

File::~File() {

}

void File::setData(boost::numeric::ublas::matrix<double> mat){
	this->mat=mat;
}

void File::setData(boost::numeric::ublas::vector<double> vec){
	setData(vec, true);
}

void File::setData(boost::numeric::ublas::vector<double> vec, bool col){
	if (col){
		mat.resize(vec.size(),1); mat.clear();
		boost::numeric::ublas::matrix_column<boost::numeric::ublas::matrix<double> >(mat,0)=vec;
	}else{
		mat.resize(1,vec.size()); mat.clear();
		boost::numeric::ublas::matrix_row<boost::numeric::ublas::matrix<double> >(mat,0)=vec;
	}
}

void File::setData(boost::numeric::ublas::vector<std::complex<double> > vec){
	mat.resize(vec.size(),2);
	mat.clear();
	boost::numeric::ublas::matrix_column<boost::numeric::ublas::matrix<double> >(mat,0)=real(vec);
	boost::numeric::ublas::matrix_column<boost::numeric::ublas::matrix<double> >(mat,1)=imag(vec);
}

void File::setData(boost::numeric::ublas::vector<std::pair<double,double> > vec){
	mat.resize(vec.size(),2);
	mat.clear();
	for (uint i=0; i< vec.size(); i++){
		mat(i,0)=vec(i).first;
		mat(i,1)=vec(i).second;
	}
}

void File::setData(std::vector<std::pair<double,double> > vec){
	mat.resize(vec.size(),2);
	mat.clear();
	for (uint i=0; i< vec.size(); i++){
		mat(i,0)=vec[i].first;
		mat(i,1)=vec[i].second;
	}
}

void File::setData(boost::numeric::ublas::matrix<std::complex<double> > mat){
	this->mat.resize(mat.size1(),2*mat.size2());
	this->mat.clear();
	for (uint i=0; i<mat.size2(); i++){
		boost::numeric::ublas::matrix_column<boost::numeric::ublas::matrix<double> >(this->mat,2*i)=real(boost::numeric::ublas::matrix_column<boost::numeric::ublas::matrix<std::complex<double> > >(mat, i));
		boost::numeric::ublas::matrix_column<boost::numeric::ublas::matrix<double> >(this->mat,2*i+1)=imag(boost::numeric::ublas::matrix_column<boost::numeric::ublas::matrix<std::complex<double> > >(mat, i));
	}
}

void File::setData(boost::numeric::ublas::matrix<boost::numeric::ublas::vector<std::complex<double> > > mat){
	this->mat.resize(mat(0,0).size(),2*mat.size1()*mat.size2());
	this->mat.clear();
	for (uint i1=0; i1<mat.size1(); i1++){
		for (uint i2=0; i2<mat.size2(); i2++){
			boost::numeric::ublas::matrix_column<boost::numeric::ublas::matrix<double> >(this->mat,2*(mat.size1()*i1+i2))=real(mat(i1,i2));
			boost::numeric::ublas::matrix_column<boost::numeric::ublas::matrix<double> >(this->mat,2*(mat.size1()*i1+i2)+1)=imag(mat(i1,i2));
		}
	}
}

/** Adds the vector vec as a column vector to mat
 *
 *
 */
void File::addVector(boost::numeric::ublas::vector<double> vec){
	uint matWidth=mat.size2();
	mat.resize(vec.size(),matWidth+1,true);
	boost::numeric::ublas::matrix_column<boost::numeric::ublas::matrix<double> >(mat,matWidth)=vec;
}

void File::addRowVector(boost::numeric::ublas::vector<double> vec){
	uint matHeight=mat.size1();
	mat.resize(matHeight+1,vec.size(),true);
	boost::numeric::ublas::matrix_row<boost::numeric::ublas::matrix<double> >(mat,matHeight)=vec;
}

void File::addColVector(boost::numeric::ublas::vector<double> vec){
	addVector(vec);
}

void File::setFilename(std::string filename){
	this->filename=filename;
}

void File::setPathname(std::string pathname){
	this->pathname=pathname;
}

std::string File::getPathname(){
	return this->pathname;
}

void File::writeToFile(){
	std::string f=pathname+filename;
	std::ofstream out(f.c_str());
	out.precision(prec);

	for (uint i=0; i<mat.size1(); i++){
		for (uint j=0; j<mat.size2(); j++){
			out << mat(i,j) << " ";
		}
		out << std::endl;
	}

}

void File::readFromFile(){
	std::string f=pathname+filename;
	std::ifstream in(f.c_str());
	const int sz=4096;
	char buf[sz];
	double value;

	std::vector<std::vector<double> > data;
	std::vector<double> line;
	// read lines from input file
	while(in.getline(buf, sz)) { // Removes \n
		std::string st=buf;

		while(st.length()!=0){
			// equal to trim left
			st.erase(0, st.find_first_not_of(" "));

			std::istrstream s(st.c_str());
			s >> value;
			line.push_back(value);
			if (st.find(" ")==std::string::npos){
				break;
			}
			st.erase(0, st.find(" ")+1);
		}

		if (line.size()>0){
			data.push_back(line);
		}
	    line.clear();
	}

	if (data.size()>0){
		mat.resize(data.size(), data[0].size());
		mat.clear();
		for (uint i=0; i<mat.size1(); i++){
			line=data[i];
			for (uint j=0; j<mat.size2(); j++){
				mat(i,j)=line[j];
			}
		}
	}
}

boost::numeric::ublas::matrix<double> File::getData(){
	return mat;
}
