/*
 * File.h
 *
 *  Created on: Dec 5, 2008
 *      Author: michael
 */


#ifndef FILE_H_
#define FILE_H_

#include <string>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

class File {
public:
	File(std::string pathname);
	File(std::string pathname, std::string subfolder);
	File(std::string pathname, std::string subfolder, std::string filename);
	virtual ~File();

	std::string getPathname();

	void setFilename(std::string filename);
	void setPathname(std::string pathname);
	void setData(boost::numeric::ublas::matrix<double> mat);
	void setData(boost::numeric::ublas::matrix<std::complex<double> > mat);
	void setData(boost::numeric::ublas::matrix<boost::numeric::ublas::vector<std::complex<double> > > mat);
	void setData(boost::numeric::ublas::vector<double> vec);
	void setData(boost::numeric::ublas::vector<double> vec, bool col);
	void setData(boost::numeric::ublas::vector<std::complex<double> > vec);
	void setData(boost::numeric::ublas::vector<std::pair<double,double> > vec);
	void setData(std::vector<std::pair<double,double> > vec);
	void addVector(boost::numeric::ublas::vector<double> vec);
	void addRowVector(boost::numeric::ublas::vector<double> vec);
	void addColVector(boost::numeric::ublas::vector<double> vec);
	void writeToFile();
	void readFromFile();
	boost::numeric::ublas::matrix<double> getData();
	void setPrec(int prec){this->prec=prec;}


protected:
	std::string pathname;
	std::string filename;

	boost::numeric::ublas::matrix<double> mat;
	int prec;
};

#endif /* FILE_H_ */
