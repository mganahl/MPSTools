#include "ParameterParser.hpp"

#include<algorithm>

ParameterParser::ParameterParser(std::string filename)
{
  filename_=filename;
  parse();
}

ParameterParser::~ParameterParser()
{
}


void ParameterParser::parse()
{
  params_.clear();

  std::ifstream is(filename_.c_str());

  if ( !is.is_open() )
    throw std::runtime_error("ParameterParser::parse could not open file: " +filename_);

  std::string line;

  while ( getline ( is, line ) )
    {
      std::string::size_type i = line.find_first_not_of ( " \t\n\v" );

      if(i==std::string::npos)
	continue;
      if (line[i] == '#' || line[i] == '%' || (line[i] == '/'&&line[i+1] == '/'))
	continue;

      if (line.find("/*")!=std::string::npos){
	line.erase(line.find("/*"),line.find("*/"));
      }
      if (line.find("//")!=std::string::npos){
	line.erase(line.find("//"));
      }

      std::stringstream line_ss(line);

      std::stringstream name_ss;
      std::stringstream val_ss;
      char letter;
      while(line_ss.get(letter))
        {
	  if(letter=='=')
	    break;
	  name_ss<<letter;
        }
      while(line_ss.get(letter))
        {
	  if(letter==';')
	    break;
	  val_ss<<letter;
        }
      std::string name=name_ss.str();
      name.erase(std::find_if(name.rbegin(), name.rend(), std::not1(std::ptr_fun<int,int>(std::isspace))).base(),name.end());
      std::string val_str=val_ss.str();

      if(params_.find(name) != params_.end())
	throw std::runtime_error("ParameterParser::parse parameter "+name+" defined twice");

      params_.insert(std::pair<std::string,std::string>(name,val_str));

    }

  is.close();
}

bool ParameterParser::isDefined(std::string name)
{
  if(params_.find(name) != params_.end())
    return true;
  return false;
}


void ParameterParser::copyFile(std::string filename)
{
  system(("cp "+filename_+" "+filename).c_str());
}

double ParameterParser::getValueOrDefault(std::string param_name,double default_val) const
{
  std::map<std::string,std::string>::const_iterator iter= params_.find(param_name);
  if(iter != params_.end())
    return atof(iter->second.c_str());
  std::cerr<<"ParameterParser WARNING: reverting to default "<<param_name<<"="<<default_val<<std::endl;
  return default_val;
}

int ParameterParser::getValueOrDefault(std::string param_name,int default_val) const
{
  std::map<std::string,std::string>::const_iterator iter= params_.find(param_name);
  if(iter != params_.end())
    return atoi(iter->second.c_str());
  std::cerr<<"ParameterParser WARNING: reverting to default "<<param_name<<"="<<default_val<<std::endl;
  return default_val;
}

unsigned int ParameterParser::getValueOrDefault(std::string param_name,unsigned int default_val) const
{
  std::map<std::string,std::string>::const_iterator iter= params_.find(param_name);
  if(iter != params_.end())
    return (unsigned int)atoi(iter->second.c_str());
  std::cerr<<"ParameterParser WARNING: reverting to default "<<param_name<<"="<<default_val<<std::endl;
  return default_val;
}

std::string ParameterParser::getValueOrDefault(std::string param_name,std::string default_val) const
{
  std::map<std::string,std::string>::const_iterator iter= params_.find(param_name);
  std::string res;
  if(iter != params_.end()){
    res=iter->second;
    if (res.find('\"')!=std::string::npos){
      res=res.substr(res.find('\"')+1,res.size());
    }
    if (res.find('\"')!=std::string::npos){
      res=res.substr(0,res.find('\"'));
      if (res.size()==1){
	res.clear();
      }
    }
    if (res.find('\"')!=std::string::npos){
      std::cerr<<"ParameterParser WARNING: Too many quotes "<<param_name<<"="<<res<<std::endl;
    }
    return res;
  }

  std::cerr<<"ParameterParser WARNING: reverting to default "<<param_name<<"="<<default_val<<std::endl;
  return default_val;
}
std::vector<double> ParameterParser::getValueOrDefault(std::string param_name,std::vector<double> default_val) const
{
  std::map<std::string,std::string>::const_iterator iter= params_.find(param_name);
  if(iter != params_.end())
    {
      std::vector<double> values;
      std::stringstream valss( iter->second);
      char letter;
      valss.get(letter);
      if(letter != '[')
	throw std::runtime_error("ParameterParser::getValueOrDefault could not parse vector<double>");

      while(true)
        {
	  std::stringstream singlevalss;
	  bool finished=false;
	  while(!(finished=!valss.get(letter)))
            {
	      if(letter==';' || letter==','|| letter==' ')
		break;
	      if(letter==']')
                {
		  finished=true;
		  break;
                }
	      singlevalss<<letter;
            }
	  values.push_back(atof(singlevalss.str().c_str()));

	  if(finished)
	    break;
        }
      return values;

    }
  return default_val;
}

std::vector<unsigned int> ParameterParser::getValueOrDefault(std::string param_name,std::vector<unsigned int> default_val) const
{
  std::map<std::string,std::string>::const_iterator iter= params_.find(param_name);
  if(iter != params_.end())
    {
      std::vector<unsigned int> values;
      std::stringstream valss( iter->second);
      char letter;
      valss.get(letter);
      if(letter != '[')
	throw std::runtime_error("ParameterParser::getValueOrDefault could not parse vector<uint>");

      while(true)
        {
	  std::stringstream singlevalss;
	  bool finished=false;
	  while(!(finished=!valss.get(letter)))
            {
	      if(letter==';' || letter==','|| letter==' ')
		break;
	      if(letter==']')
                {
		  finished=true;
		  break;
                }
	      singlevalss<<letter;
            }
	  values.push_back(atoi(singlevalss.str().c_str()));

	  if(finished)
	    break;
        }
      return values;

    }
  return default_val;
}

bool ParameterParser::getValueOrDefault(std::string param_name,bool default_val) const
{
  std::map<std::string,std::string>::const_iterator iter= params_.find(param_name);
  if(iter != params_.end())
    {
      if(iter->second.find("true")!=std::string::npos)
	return true;
      else if(iter->second.find("false")!=std::string::npos)
	return false;
      else
	throw std::runtime_error("ParameterParser::getValueOrDefault could not parse bool parameter "+param_name);
    }
  std::cerr<<"ParameterParser WARNING: reverting to default "<<param_name<<"="<<default_val<<std::endl;
  return default_val;
}


void ParameterParser::write() const
{
  std::cout<<"Parameters:"<<std::endl;
  for(std::map<std::string,std::string>::const_iterator iter=params_.begin();
      iter!=params_.end();++iter)
    {
      std::cout<<"  "<<iter->first<<" = "<<iter->second<<std::endl;
    }
  std::cout<<"Parameters End"<<std::endl;

}


