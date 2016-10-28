#ifndef __ARRAY_H_
#define __ARRAY_H_
/* this is a substitute for std::vector<T>.
 * 
 * In addition, value-arrays (such as double and complex) are in 
 * ivalarray which is just an overloaded array with a couple of extra capabilities
 * applying to statically allocated objects.
 */
 
// this is the very basic class. it doesn't depend on any other class.
 

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <new>
#include <cstring>
#include <fstream>
namespace arrayclass{
  template<class T>
  class array
  {
  public:
    array();
    array(const size_t size);
    array(const size_t size, const T& elem);
    array(const array<T>& arr);
    ~array();
        
    // assign operator
    array<T>& operator=(const array<T>& arr);
        
    // memory management
    void reserve(const size_t capacity, bool preserve=true);
    void resize(const size_t capacity, bool preserve=true);
        
    void clear() { _size = 0; }
        
    // information
    size_t size() const { return _size; }
    size_t capacity() const { return _capacity; }
        
    // element access
    T& operator[](const size_t i) {assert(i<_size); return _data[i]; }
    const T& operator[](const size_t i) const { assert(i<_size); return _data[i]; }
    const T& front() const { return _data[0]; }
    T& front() { return _data[0]; }
    const T& back() const { return _data[_size-1]; }
    T& back() { return _data[_size-1]; }        
    void push_back(const T& a);
    T pop_back();
        
        
    // raw access        
    const T* raw() const  { return _data; }
    T* raw() { return _data; }
        
        
        
    int index(const T& val) const;
        
    // various actions.
    array<T>& extend(const array<T>& addarr);
        
    // I/O
    void write(std::ofstream& fp) const;
    void read(std::ifstream& fp);
        
    size_t memoryusage() const;
        
    void swap(array<T>& x);
        
  protected:
    size_t _capacity, _size;
    T* _data;
  };


  template<typename T>
  inline array<T>::array(): _capacity(0), _size(0), _data(NULL)
  {
  }

  template<typename T>
  inline array<T>::array(const size_t size): _capacity(0), _size(0), _data(NULL)
  {
    resize(size, false);
  }

  template<typename T>
  array<T>::array(const size_t size, const T& elem): _capacity(0), _size(0), _data(NULL)
  {
    resize(size, false);
    for(size_t i=0;i<size;++i)
      _data[i] = elem;
  }

  template<typename T>
  array<T>::array(const array<T>& arr): _capacity(0), _size(0), _data(NULL)
  {
    resize(arr._size, false);
    for(size_t i=0;i<_size;++i) _data[i]=arr._data[i];
  }

  template<typename T>
  inline array<T>::~array()
  {
#ifdef DEBUG
    //std::cerr << "Freeing array " << this << std::endl;
#endif
    if (_data)
      delete [] _data;
  }



  template<typename T>
  void array<T>::reserve(const size_t suggested_capacity, bool preserve)
  {
    if (suggested_capacity <= _capacity ) return;
    
    size_t capacity = ( ((suggested_capacity>>4) + 1) << 4 );
    assert( capacity );
    T* newdata = NULL;

    try {
      newdata = new T[capacity];
    } catch (std::bad_alloc& ba) {
      std::cerr << "bad_alloc caught in array: " << ba.what() << std::endl;
    }
    assert( newdata );
    if (preserve && _size)
      for(size_t i=0;i<_size;++i) newdata[i] = _data[i];
    
    if (_data) {
      delete [] _data;
    }
    _data = newdata;
    _capacity = capacity;
  }

  template<typename T>
  void array<T>::resize(const size_t newsize, bool preserve)
  {
    if (newsize > _capacity)
      reserve(newsize, preserve);
    _size = newsize;
  }


  template<typename T>
  array<T>& array<T>::operator=(const array<T>& arr)
  {
    resize(arr._size, false);
    assert( _size==arr._size);
    for(size_t i=0;i<_size;++i) _data[i]=arr._data[i];
    return (*this);
  }

  template<typename T>
  void array<T>::push_back(const T& x)
  {
    if (_size + 1 > _capacity)
      reserve(_capacity+64, true);
    resize(_size+1, true);
    back() = x;
  }

  template<typename T>
  T array<T>::pop_back()
  {
    const unsigned int n = this->size();
    if (!n) {
      std::cerr << "cant do that, its empty already.\n";
      throw "emptgfdf";
    }
    T x = (*this)[n-1];
    resize(n-1, true);
    return x;
  }



  template<typename T>
  int array<T>::index(const T& val) const
  {
    int pos = -1;
    int i=0;
    const int n = this->size();
    if (n==0) return pos;
    
    while (i < n) {
      if ((*this)[i] == val) {
	pos = i;
	break;
      }
      ++i;
    }
    return pos;
  }



  // for double array, this can be specialized.

  template<>
  array<double>& array<double>::operator=(const array<double>& arr)
  {
    resize(arr._size, false);
    memcpy(_data, arr._data, sizeof(double)*_size);
    return (*this);
  }

  template<>
  void array<double>::reserve(const size_t sug_capacity, bool preserve)
  {
    if (sug_capacity <= _capacity ) return;
    const size_t capacity = (((sug_capacity>>5)+1)<<5);
    
    assert( capacity );
    double* newdata = new double[capacity];
    assert( newdata );
    if (preserve && _size)
      memcpy(newdata, _data, sizeof(double)*_size);
    
    if (_data) delete [] _data;
    _data = newdata;
    _capacity = capacity;
  }


  template<typename T>
  T sum(const array<T>& x)
  {
    T s=0;
    for(int i=0;i<x.size();++i)
      s += x[i];
    return s;
  }


  template<>
  array<double>::array(const array<double>& arr): _capacity(0), _size(0), _data(NULL)
  {
    resize(arr._size, false);
    for(size_t i=0;i<_size;++i) _data[i]=arr._data[i];
  }



  template<typename T>
  std::ostream& operator<<(std::ostream& str, const array<T>& arr)
  {
    for(int i=0;i<arr.size();++i) str<< arr[i] << " " ;
    return str;

  }




  template<typename T>
  array<T>& array<T>::extend(const array<T>& addarr)
  {
    const unsigned int n1 = this->size();
    const unsigned int n2 = addarr.size();
    this->resize( this->size() + n2, true );
    memcpy( this->raw() + n1, addarr.raw(), sizeof(T) * n2);
    return (*this);
  }



  template<typename T>
  void array<T>::write(std::ofstream& fp) const
  {
    unsigned int n = this->size();
    size_t szt = sizeof(T);
    fp.write((char*)&n, sizeof(unsigned int));
    fp.write((char*)&szt, sizeof(size_t));
    for(int i=0;i<n;++i)
      (*this)[i].write(fp);
    //fp.write((char*)this->raw(), n * sizeof(T) );
  }

  template<typename T>
  void array<T>::read(std::ifstream& fp)
  {
    unsigned int n;
    size_t szt;
    fp.read((char*)&n, sizeof(unsigned int));
    fp.read((char*)&szt, sizeof(size_t));
    if (szt != sizeof(T) ) {
      std::cerr << "szt="<<szt << ", sizeof(T)="<<sizeof(T) << std::endl;
      std::cerr << "type mismatch cannot read array<T>.\n";
      throw "seiorsu fdasf";
    }
    
    this->resize(n, false);
    for(int i=0;i<n;++i)
      (*this)[i].read(fp);
    //fp.read((char*)this->raw(), sizeof(T)*n);
  }


  template<typename T>
  size_t array<T>::memoryusage() const
  {
    size_t memusage = _capacity * sizeof(T);
    memusage += 2 * sizeof(size_t) + sizeof(T*);
    return memusage;
  }



  template<typename T>
  void array<T>::swap(array<T>& x)
  {
    const size_t c = _capacity, s = _size;
    T* d = _data;
    
    _data = x._data;
    _capacity = x._capacity;
    _size = x._size;
    
    x._data = d;
    x._capacity = c;
    x._size = s;
  }

}





#endif
