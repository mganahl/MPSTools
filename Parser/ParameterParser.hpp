#ifndef PARAMETERPARSER_HPP_
#define PARAMETERPARSER_HPP_

#include <iostream>
#include <string>
#include <fstream>
#include <exception>
#include <stdexcept>
#include <sstream>
#include <map>
#include <vector>
#include <stdlib.h>


class ParameterParser
{
public:
  ParameterParser(std::string filename);
  virtual ~ParameterParser();

  void parse();

  double getValueOrDefault(std::string param_name,double default_val) const;
  int getValueOrDefault(std::string param_name,int default_val) const;
  unsigned int getValueOrDefault(std::string param_name,unsigned int default_val) const;
  std::string getValueOrDefault(std::string param_name,std::string default_val) const;
  bool getValueOrDefault(std::string param_name,bool default_val) const;
  std::vector<double> getValueOrDefault(std::string param_name,std::vector<double> default_val) const;
  std::vector<unsigned int> getValueOrDefault(std::string param_name,std::vector<unsigned int> default_val) const;

  template <typename T>
  void setValue(std::string param_name,T val){
    std::ostringstream strs; strs<<val;
    params_[param_name]=strs.str();
  }
  void setValue(std::string param_name,bool val){
    params_[param_name]=val?"true":"false";
  }

  bool isDefined(std::string name);

  void copyFile(std::string filename);

  void write() const;
private:
  std::string filename_;

  std::map<std::string,std::string> params_;


};

#endif /*PARAMETERPARSER_HPP_*/
