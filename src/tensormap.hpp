#ifndef __TENSOR_MAP__
#define __TENSOR_MAP__

#include <iostream>
#include <fstream>
#include <map>

#include "../boost_typedefs.hpp"
//#include "keyTypes.hpp"

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

#include <boost/serialization/utility.hpp>
#include <boost/serialization/collections_save_imp.hpp>
#include <boost/serialization/collections_load_imp.hpp>
#include <boost/serialization/split_free.hpp>

using BoostTypedefs::Range;


template<typename TK, typename TV>
class map_container: public std::map<TK,TV>
{
public:
    map_container();
    ~map_container();
    
    typedef typename std::map<TK,TV>::iterator iterator;
    typedef typename std::map<TK,TV>::const_iterator const_iterator;

    size_t memoryusage() const;

    void update(const map_container<TK,TV>& X);
    
    void tobuffer(std::ofstream& buf) const;
    void frombuffer(std::ifstream& buf);
    
    void write(std::ofstream& buf) const { tobuffer(buf); }
    void read(std::ifstream& buf) { frombuffer(buf); }
    
    void tofile(const char* fname) const;
    void fromfile(const char* fname);

};

template<typename TK, typename TV>
map_container<TK,TV>::map_container()
    : std::map<TK,TV>()
{
}


template<typename TK, typename TV>
map_container<TK,TV>::~map_container()
{
}



template<typename TK, typename TV>
size_t map_container<TK,TV>::memoryusage() const
{
    size_t memusg = 0;
    for (const_iterator i=this->begin();i!=this->end();++i) {
        memusg += sizeof(i->first);
        memusg += i->second.memoryusage();
    }
    return memusg;
}

template<typename TK,typename TV>
void map_container<TK,TV>::tobuffer(std::ofstream& buf) const
{
  int n = this->size();
    
  buf.write((char*)&n, 4);
    
  for(const_iterator i=this->begin();i!=this->end();++i) {
    //buf.write( (char*)&i->first, sizeof(TK)); //this is different from iztok's original code; my ukey is not the same as his, so this line of code doesn't work
    i->first.tobuffer(buf);
    i->second.tobuffer(buf);
  }
}

template<typename TK,typename TV>
void map_container<TK,TV>::frombuffer(std::ifstream& buf)
{
    int n;
    buf.read((char*)&n, 4);

    this->clear();

    for(int i=0;i<n;++i) {
        TK ky;
        TV t;
	ky.frombuffer(buf);
	//buf.read((char*)&ky, sizeof(TK));//see one comment above ...
        t.frombuffer(buf);
        (*this)[ky] = t;
    }
    assert( this->size() == n );
}

template<typename TK,typename TV>
void map_container<TK,TV>::tofile(const char* fname) const
{
    std::ofstream fp(fname, std::ios::out | std::ios::binary);
    tobuffer(fp);
    fp.close();
}

template<typename TK,typename TV>
void map_container<TK,TV>::fromfile(const char* fname) 
{
    std::ifstream fp(fname, std::ios::in | std::ios::binary);
    frombuffer(fp);
    fp.close();
}

template<typename TK,typename TV>
void map_container<TK,TV>::update(const map_container<TK,TV>& src)
{
    for(const_iterator i=src.begin();i!=src.end();++i)
        (*this)[i->first] = i->second;    
}



template<typename Kt,typename Mt>
void merge(map_container<Kt,Mt>&m1,const map_container<Kt,Mt>&m2)
{
  for(typename map_container<Kt,Mt>::const_iterator it = m2.begin();it!=m2.end();++it)
    m1[it->first]=it->second;
}




template<typename T1,typename T2>
class valmap: public std::map<T1,T2>
{
public:
    void write(std::ofstream& fp) const;
    void read(std::ifstream& fp);

    void tobuffer(std::ofstream& fp) const { write(fp); }
    void frombuffer(std::ifstream& fp) { read(fp); }

    
    typedef typename std::map<T1,T2>::iterator iterator;
    typedef typename std::map<T1,T2>::const_iterator const_iterator;
};




template<typename T1,typename T2>
void valmap<T1,T2>::write(std::ofstream& fp) const
{
    
    unsigned int n = this->size();
    unsigned int size_key = sizeof(T1);
    unsigned int size_val = sizeof(T2);

    fp.write((char*)&n, sizeof(unsigned int));
    fp.write((char*)&size_key, sizeof(unsigned int));
    fp.write((char*)&size_val, sizeof(unsigned int));
    
    for(const_iterator i=this->begin();i!=this->end();++i) {
        fp.write((char*)&(i->first), size_key );
        fp.write((char*)&(i->second), size_val);
    }
}

template<typename T1,typename T2>
void valmap<T1,T2>::read(std::ifstream& fp)
{
    unsigned int n;
    unsigned int size_key, size_val;
    
    fp.read((char*)&n, sizeof(unsigned int));
    fp.read((char*)&size_key, sizeof(unsigned int));
    fp.read((char*)&size_val, sizeof(unsigned int));
    assert( sizeof(T1) == size_key );
    assert( sizeof(T2) == size_val );
    this->clear();
    for(int i=0;i<n;++i) {
        T1 ky;
        fp.read((char*)&ky, size_key);
        T2 val;
        fp.read((char*)&val, size_val);
        (*this)[ky] = val;
    }
    
    assert( n == this->size() );
}


// // ---------------- TENSOR_MAP ----------------


// template<typename TK,size_t N2, typename T>
// class tensor_map: public map_container<TK, tensor<N2,T> >
// {
// public:
//     tensor_map();
//     ~tensor_map();
    
//     typedef typename map_container<TK, tensor<N2,T> >::iterator iterator;
//     typedef typename map_container<TK, tensor<N2,T> >::const_iterator const_iterator;
    
    
//     tensor_map<TK,N2,T>& operator*=(const double fac);
//     tensor_map<TK,N2,T>& operator+=(const tensor_map<TK,N2,T>& src);
    
//     double normsq() const;
//     double norm() const { return sqrt(normsq()); }
// };

// template<typename TK,size_t N2,typename T>
// tensor_map<TK,N2,T>::tensor_map()
//     : map_container<TK, tensor<N2,T> >()
// {
    
    
// }

// template<typename TK,size_t N2,typename T>
// tensor_map<TK,N2,T>::~tensor_map()
// {
    
// }

// template<typename TK,size_t N2, typename T>
// tensor_map<TK,N2,T>& tensor_map<TK,N2,T>::operator*=(const double fac)
// {
//     for(iterator i=this->begin();i!=this->end();++i)
//         i->second *= fac;
//     return *this;
// }

// template<typename TK,size_t N2, typename T>
// tensor_map<TK,N2,T>& tensor_map<TK,N2,T>::operator+=(const tensor_map<TK,N2,T>& src)
// {
//     for(const_iterator i=src.begin();i!=src.end();++i) {
//         iterator j = this->find(i->first);
//         if (j==this->end())
//             (*this)[i->first] =i->second;
//         else
//             j->second += i->second;
//     }
    
//     return *this;
// }



// template<typename TK,size_t N2, typename T>
// double tensor_map<TK,N2,T>::normsq() const
// {
//     double s=0;
//     for(const_iterator i=this->begin();i!=this->end();++i)
//         s += i->second.normsq();
//     return s;
// }

//-----------------Number table map --------------------

template<class KeyType>
class range_map:public map_container<KeyType,Range>{

};

template<class KeyType>
class numberMaps:public map_container<int,range_map<KeyType> >
{
};

template<class KeyType>
class i_to_key:public map_container<uint,KeyType>
{
public:
  i_to_key(){_basis=1;};
  i_to_key(const i_to_key<KeyType> & i);
  int getSign()const{return _basis;}
  void setSign(const int basistype){assert(basistype);_basis=basistype;}

  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)    {
      ar & boost::serialization::base_object<map_container<uint,KeyType> >(*this) & _basis;
  }
private:
  int _basis;
};

template<class KeyType>
i_to_key<KeyType>::i_to_key(const i_to_key<KeyType> & i):map_container<uint,KeyType>(i){};

#endif
