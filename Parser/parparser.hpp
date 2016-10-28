#ifndef __PARPARSER__
#define __PARPARSER__
/*parameter  parser
  author: iztok pizorn
*/

#include <iostream>
#include <cstdlib>
#include <vector>
#include <string.h>
#include <stdio.h>
#include "array.hpp"
enum ParameterType {
  Parameter_Integer,
  Parameter_Float,
  Parameter_String,
  Parameter_Option
};


struct Tparameter {
  std::string name;
  ParameterType type;
  std::string value;
  bool given;
  bool isargument;
};

class ParameterError
{
public:
  ParameterError(const char* message)
  {
    std::cerr << message << std::endl;
  }

};

class parparser
{
public:
  parparser();

  void addargument(const char*, const ParameterType);
  void addoption(const char*, const ParameterType, const char* defval=NULL);

  unsigned int argcount() const;

  void printhelp() const;
  //  void printinfo() const;
        
  // load: returns false on failure.
  bool load(const unsigned int, const char**);

  double getfloat(const char* x) const;
  long getint(const char* x) const;
  const char* getstring(const char*) const;
  bool isgiven(const char*) const;
  double getdouble(const char* x) const  __attribute__((deprecated)) { return getfloat(x); }
        

  unsigned int size() const  { return options.size(); }
  const Tparameter& operator[](const unsigned int i) const { return options[i]; }

private:

  const Tparameter* findoption(const char* name) const;
  const Tparameter* findoption(const std::string& myname) const;
  Tparameter* findoption(const char* name) ;
  Tparameter* findoption(const std::string& myname) ;

  unsigned int nargs;

  void processoption(const char*);
  void processargument(const char*);

  arrayclass::array< Tparameter > options;
};

parparser::parparser(): nargs(0)
{
  addoption("help", Parameter_Option);
}

const Tparameter* parparser::findoption(const std::string& myname) const
{
  int j=0;
  while (j < options.size() && options[j].name != myname) ++j;
  if (j>=options.size() ) 
    return NULL;
  else 
    return &options[j];
}


const Tparameter* parparser::findoption(const char* name) const
{
  return findoption(std::string(name));
}

Tparameter* parparser::findoption(const std::string& myname)
{
  int j=0;
  while (j < options.size() && options[j].name != myname) ++j;
  if (j>=options.size() ) 
    return NULL;
  else 
    return &options[j];
}


Tparameter* parparser::findoption(const char* name)
{
  return findoption(std::string(name));
}


bool parparser::isgiven(const char* option) const
{
  const Tparameter* par = findoption(option);
  if (!par) throw ParameterError("No such option defined.");
  return par->given;
}

long parparser::getint(const char* option) const
{
  const Tparameter* par = findoption(option);
  if (!par) throw ParameterError("No such option defined.");
  if (par->type != Parameter_Integer) {
    std::cerr << "Option " << par->name << " is not an integer.\n";
    throw ParameterError("Wrong type");
  }
  return strtol(par->value.c_str(), NULL, 10 );
}



double parparser::getfloat(const char* option) const
{
  const Tparameter* par = findoption(option);
  if (!par) throw ParameterError("No such option defined.");
  if (par->type != Parameter_Float) {
    std::cerr << "Option " << par->name << " is not a floating point number.\n";
    throw ParameterError("Wrong type");
  }
  return strtod(par->value.c_str(), NULL );
}

const char* parparser::getstring(const char* option) const
{
  const Tparameter* par = findoption(option);
  if (!par) throw ParameterError("No such option defined.");
  if (par->type != Parameter_String) {
    std::cerr << "Option " << par->name << " is not a string.\n";
    throw ParameterError("Wrong type");
  }
  return par->value.c_str();
}


void parparser::addoption(const char* opt, const ParameterType type, const char* defval)
{
  Tparameter par;
  par.name = opt;
  par.type = type;
  if (defval) 
    par.value = std::string(defval);
  else if (type == Parameter_Option) 
    par.value = std::string("false");
  par.given = false;
  par.isargument = false;
  options.push_back(par);
}

void parparser::addargument(const char* opt, const ParameterType type )
{
  Tparameter par;
  par.name = opt;
  par.type = type;
  par.given = false;
  par.isargument = true;
  options.push_back(par);
  ++nargs;
}


bool parparser::load(const unsigned int argc, const char** argv)
{
  // always add a help argument...
    
  bool continuation = false;
  std::string cont_str;
  int i=0;
  while (++i < argc) {
    const char* arg = argv[i];
        
    // check if it is a recognized name.
    if (arg[0]=='-' && arg[1]=='-') {
      processoption(arg);
      continuation = false;
      cont_str.clear();
      continue;
    }

    if (continuation) {
      assert( !cont_str.empty() );
      cont_str += std::string(arg);
      processoption(cont_str.c_str() );
      cont_str.clear();
      continuation = false;
      continue;
    }
        
        
    if (arg[0]=='-' && arg[1]!='-' && !strstr(arg, "=") ) {
      const char* xend = arg + 1;
      bool bfound = false;
      while (*xend) {
	// check if it's a valid name.
	if (findoption(std::string(arg+1).substr(0, (size_t)(xend-arg)).c_str()) ) {
	  // found it!
	  bfound = true;
	  break;
	} else 
	  ++xend;
      }
      std::string my_str = std::string("--") + std::string(arg+1).substr(0, (size_t)(xend-arg)) + std::string("=");
      if (*xend) {
	my_str += std::string(xend+1);
	//std::cerr << "Would send : " << my_str << std::endl;
	processoption(my_str.c_str() );
	continue;
      } else {
	continuation = true;
	cont_str = my_str;
	continue;
      }
    }
    // else, it's an argument.
    processargument(arg);
  }
    
  // just check if help is given...
  if (isgiven("help")) {
    printhelp();
    return false;
  }


  // check the number of arguments...
  unsigned int givenargs=0;
  for(int j=0;j<options.size();++j)
    if (options[j].isargument && options[j].given) ++givenargs;

  if (nargs != givenargs) {
    printhelp();
    throw ParameterError("Insufficient number of arguments given.\n");
  }
    
  return true;
}


void parparser::processoption(const char* opt)
{
  char myopt[128];
  strcpy(myopt, opt);
  char* myname = myopt;
  while (*myname=='-') ++myname;
    
    
  char* name_end = myname;
  while (*name_end != '=' && *name_end) ++name_end;
    
  const char* val = (*name_end =='=' ? name_end + 1 : NULL);
  *name_end = 0;
    
  std::string name(myname);
    
  Tparameter* par = findoption(name);
  // find the one in the list;
  if (!par) {
    std::cerr << "There is no option with the name " << name << std::endl;
    printhelp();
    throw ParameterError("No such option");
  }

  if (par->given) {
    std::cerr << "Option " << name << " has been given more than once!\n";
    throw ParameterError("more than once.");
  }

  if (par->type == Parameter_Option) {
    par->value = std::string("true");
    par->given = true;
    return;
  }
    
  // now process to get the value.
  if (!val)
    ParameterError("no value given");
  par->value = std::string(val);
  par->given = true;
    
  //std::cerr << "Setting option " << name << " to " << par->value << std::endl;
}

void parparser::processargument(const char* opt)
{
  int j=0;
  while (j < options.size() && (!options[j].isargument  || options[j].given )) ++j;
  if (j >= options.size() ) {
    std::cerr << "Too many arguemnts.\n";
    std::cerr << "Failed while processing " << opt << std::endl;
    printhelp();
    throw ParameterError("too many arguments.");
  }
    
  if (!options[j].isargument) throw ParameterError("unexpected error.");
    
  //std::cerr << "setting argument " << arguments[nargs].name << " to " << opt << std::endl;
  options[j].value = std::string(opt);
  options[j].given = true;
}

void parparser::printhelp() const
{
  std::cout << "Arguments:\n   ";
  for(int i=0;i<options.size();++i) {
    if (!options[i].isargument) continue;
    std::cout << "{"<<options[i].name<<"} ";
  }
  std::cout << std::endl;
  std::cout << "Options:\n   ";
  for(int i=0;i<options.size();++i) {
    if (options[i].isargument) continue;
    if (options[i].name == std::string("help")) continue;
    std::cout << "--"<<options[i].name<<"="<<options[i].value << " " ;
  }
  std::cout << std::endl;
}


// void parparser::printinfo() const
// {
//   std::cout << "Arguments:\n   ";
//   for(int i=0;i<options.size();++i) {
//     if (!options[i].isargument) continue;
//     std::cout << "{"<<options[i].name<<"} ";
//   }
//   std::cout << std::endl;
//   std::cout << "Options:\n   ";
//   for(int i=0;i<options.size();++i) {
//     if (options[i].isargument) continue;
//     if (options[i].name == std::string("help")) continue;
//     std::cout << "--"<<options[i].name<<"="<<options[i].value << " " ;
//   }
//   std::cout << std::endl;
// }


std::ostream& operator<<(std::ostream& stream, const parparser& ego)
{
  // determine width:
  unsigned int maxw = 0;
  for(int i=0;i<ego.size();++i)
    maxw = std::max(maxw, (unsigned int)ego[i].name.size() );
  ++maxw;
  for(int i=0;i<ego.size();++i) {
    const Tparameter& opt = ego[i];
    if (opt.name==std::string("help")) continue;
    stream << "# " << opt.name;
    for(int j=0;j<maxw-opt.name.size();++j) stream << " ";
    stream <<"("<<opt.given <<"): ";
    stream << opt.value << std::endl;
  }
  return stream;
}




#endif
