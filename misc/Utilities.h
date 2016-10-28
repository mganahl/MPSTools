/*
 * Utilities.h
 *
 *  Created on: Mar 18, 2009
 *      Author: michael
 */

#ifndef UTILITIES_H_
#define UTILITIES_H_

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <vector>

class Utilities {
public:
	Utilities();
	virtual ~Utilities();

	static const int doublePlacesMap=6;
	static const int ouputPrecission=16;
	static const int fNL=1000;				// factor to multiply with the numeric limit

	static double roundDouble(double n, int places);
	static boost::numeric::ublas::vector<double> roundDoubleVector(boost::numeric::ublas::vector<double> v, int places);
	static std::vector<double> roundDoubleVector(std::vector<double> v, int places);

	static boost::numeric::ublas::matrix<double> convertStringToMatrix(std::string mString);
	static boost::numeric::ublas::vector<double> convertStringToVector(std::string mString);
	static std::vector<std::string> divideStringAtBlanc(std::string mString);
	static std::vector<double> convertRangeToVector(std::string range);

	template<class T>
	static bool equalVectors(boost::numeric::ublas::vector<T> v1, boost::numeric::ublas::vector<T> v2){
		return (boost::numeric::ublas::norm_2(v1-v2) < std::numeric_limits<double>::epsilon()*fNL);
	}

	template<class T>
	static boost::numeric::ublas::matrix<T> convertStdStdVecToBoostMatrix(std::vector<std::vector<T> > v){
		boost::numeric::ublas::matrix<T> res(v.size(),v[0].size());
		for (uint i=0; i<res.size1(); i++){
			for (uint j=0; j<res.size2(); j++){
				res(i,j)=v[i][j];
			}
		}
		return res;
	}

	template<class T>
	static boost::numeric::ublas::matrix<T> convertStdBoostVecToBoostMatrix(std::vector<boost::numeric::ublas::vector<T> > v){
		boost::numeric::ublas::matrix<T> res(v.size(),v[0].size());
		for (uint i=0; i<res.size1(); i++){
			for (uint j=0; j<res.size2(); j++){
				res(i,j)=v[i][j];
			}
		}
		return res;
	}

	template<class T>
	static boost::numeric::ublas::matrix<T> convertPairStdVecToBoostMatrix(std::vector<std::pair<std::vector<T>, std::vector<T> > > v){
		boost::numeric::ublas::matrix<T> res(v.size(),v[0].first.size()+v[0].second.size());
		for (uint i=0; i<res.size1(); i++){
			for (uint j=0; j<v[0].first.size(); j++){
				res(i,j)=v[i].first[j];
			}
			for (uint j=0; j<v[0].second.size(); j++){
				res(i,j+v[0].first.size())=v[i].second[j];
			}
		}
		return res;
	}

	template<class T>
	static std::vector<T> convertBoostToStdVector(boost::numeric::ublas::vector<T> v){
		std::vector<T> res;
		for (uint i=0; i<v.size(); i++){
			res.push_back(v(i));
		}
		return res;
	}

	template<class T>
	static boost::numeric::ublas::vector<T> convertStdToBoostVector(std::vector<T> v){
		boost::numeric::ublas::vector<T> res(v.size());
		for (uint i=0; i<v.size(); i++){
			res(i)=v[i];
		}
		return res;
	}

	static std::vector<double> convertBoostToStdVectorAndRoundValues(boost::numeric::ublas::vector<double> v, int places){
		std::vector<double> res;
		for (uint i=0; i<v.size(); i++){
			res.push_back(roundDouble(v(i), places));
		}
		return res;
	}

	static boost::numeric::ublas::vector<double> convertStdToBoostVectorAndRoundValues(std::vector<double> v, int places){
		boost::numeric::ublas::vector<double> res(v.size());
		for (uint i=0; i<v.size(); i++){
			res(i)=roundDouble(v[i], places);
		}
		return res;
	}

	/** http://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html */
	static int pointInPolygon(boost::numeric::ublas::matrix<double> vertices, boost::numeric::ublas::vector<double> x);


	template<class matrix_T>
	static matrix_T kron(boost::numeric::ublas::matrix_expression<matrix_T> const& mat_a, boost::numeric::ublas::matrix_expression<matrix_T> const& mat_b){
		using namespace boost::numeric::ublas;
		matrix_T res(mat_a().size1()*mat_b().size1(), mat_a().size2()*mat_b().size2()); res.clear();

		for (std::size_t i=0; i<mat_a().size1();i++){
			for (std::size_t j=0; j<mat_a().size2();j++){
				matrix_range<matrix_T>(res, range (mat_b().size1()*i, mat_b().size1()*(i+1)), range(mat_b().size2()*j, mat_b().size2()*(j+1)))=mat_a()(i,j)*mat_b();
			}
		}
		return res;
	}

	static boost::numeric::ublas::compressed_matrix<double> kronSparse(boost::numeric::ublas::compressed_matrix<double>& mat_a, boost::numeric::ublas::compressed_matrix<double>& mat_b){
		using namespace boost::numeric::ublas;
		compressed_matrix<double> res(mat_a.size1()*mat_b.size1(), mat_a.size2()*mat_b.size2()); res.clear();

		for (compressed_matrix<double>::iterator1 i1 = mat_a.begin1(); i1 != mat_a.end1(); ++i1) {
			for (compressed_matrix<double>::iterator2 i2 = i1.begin(); i2 != i1.end(); ++i2){
				int i=i2.index1();
				int j=i2.index2();
				matrix_range<compressed_matrix<double> >(res, range (mat_b.size1()*i, mat_b.size1()*(i+1)), range(mat_b.size2()*j, mat_b.size2()*(j+1)))=*i2*mat_b;
			}
		}

		return res;
	}

	static boost::numeric::ublas::compressed_matrix<double> kronAxISparse(boost::numeric::ublas::compressed_matrix<double>& mat_a){
		using namespace boost::numeric::ublas;
		std::size_t s1=mat_a.size1();
		std::size_t s2=mat_a.size2();
		if (s1!=s2){
			throw std::runtime_error("Utitlities::kronAxISparse s1 != s2");
		}
		compressed_matrix<double> id=identity_matrix<double>(s1);

		return kronSparse(mat_a, id);
	}

	template<class matrix_T>
	static matrix_T kronAxI(boost::numeric::ublas::matrix_expression<matrix_T> const& mat_a){
		using namespace boost::numeric::ublas;
		std::size_t s1=mat_a().size1();
		std::size_t s2=mat_a().size2();
		if (s1!=s2){
			throw std::runtime_error("Utitlities::kronAxI s1 != s2");
		}
		matrix_T id=identity_matrix<double>(s1);

		return kron(mat_a, id);
	}

	template<class matrix_T>
	static matrix_T kronIxA(boost::numeric::ublas::matrix_expression<matrix_T> const& mat_a){
		using namespace boost::numeric::ublas;
		std::size_t s1=mat_a().size1();
		std::size_t s2=mat_a().size2();

		matrix_T res(s1*s1, s2*s2); res.clear();
		std::size_t sMin=s1<s2?s1:s2;
		for (std::size_t i=0; i<sMin; i++){
			matrix_range<matrix_T>(res, range (s1*i, s1*(i+1)), range(s2*i, s2*(i+1)))=mat_a();
		}

		return res;
	}

	template<class matrix_T>
	static boost::numeric::ublas::vector<typename matrix_T::value_type> convertMatrixToVector(boost::numeric::ublas::matrix_expression<matrix_T> const& mat, bool rowOrder=true){
		using namespace boost::numeric::ublas;
		std::size_t s1=mat().size1();
		std::size_t s2=mat().size2();
		matrix_T m(mat());
		boost::numeric::ublas::vector<typename matrix_T::value_type> res(s1*s2);
		if (rowOrder){
			for (std::size_t i=0; i<s1; i++){
				vector_range<boost::numeric::ublas::vector<typename matrix_T::value_type> >(res, range(s2*i,s2*(i+1)))=matrix_row<matrix_T>(m, i);
			}
		}else{
			for (std::size_t i=0; i<s2; i++){
				vector_range<boost::numeric::ublas::vector<typename matrix_T::value_type> >(res, range(s1*i,s1*(i+1)))=matrix_column<matrix_T>(m, i);
			}
		}
		return res;
	}

	template<class vector_T>
	static boost::numeric::ublas::matrix<typename vector_T::value_type> convertVectorToMatrix(boost::numeric::ublas::vector_expression<vector_T> const& vec, std::size_t s1, std::size_t s2, bool rowOrder=true){
		using namespace boost::numeric::ublas;

		if (s1*s2!=vec().size()){
			throw std::runtime_error("Utilities::convertVectorToMatrix s1*s2!=vec().size()");
		}
		vector_T v(vec());
		matrix<typename vector_T::value_type> res(s1,s2);
		if (rowOrder){
			for (std::size_t i=0; i<s1; i++){
				matrix_row<matrix<typename vector_T::value_type> >(res,i)=vector_range<vector_T>(v, range(s2*i,s2*(i+1)));
			}
		}else{
			for (std::size_t i=0; i<s2; i++){
				matrix_column<matrix<typename vector_T::value_type> >(res,i)=vector_range<vector_T>(v, range(s1*i,s1*(i+1)));
			}
		}
		return res;
	}

	static int countParticlesFromTo(const std::vector<bool> &state, int i, int j);
	static int countParticlesBetween(const std::vector<bool> &state, int i, int j);
	static int calculateSign(const std::vector<bool> &state, int i, int j);

	static double convertStringToNumber(std::string element);
	static double convertFullStringToNumber(std::string str);
	static double nChooseK(double n, double k);
protected:
	static double splitDivString(std::string element);
};

#endif /* UTILITIES_H_ */
