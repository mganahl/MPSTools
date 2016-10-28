#ifndef __VALARRAY__
#define __VALARRAY__

#include "array.hpp"

using arrayclass::array
template<typename T>
class ivalarray: public array<T>
{
    
public:
    ivalarray(): array<T>() {}
    ivalarray(const size_t size): array<T>(size) {}
    ivalarray(const size_t size, const T& elem): array<T>(size, elem) {}
    ivalarray(const ivalarray<T>& arr): array<T>(arr) {}
    
    void zero();
    void sort(bool reverse=false);
    double norm() const { return nrm2(this->size(), this->raw(), 1); }
    double normsq() const { return dotc(this->size(), this->raw(), 1, this->raw(), 1); }
    
    ivalarray<T>& operator+=(const ivalarray<T>& src);
    ivalarray<T>& operator-=(const ivalarray<T>& src);

    T sum() const;
    
    void write(std::ofstream& fp) const;
    void read(std::ifstream& fp);

    void tobuffer(std::ofstream& fp) const { write(fp); }
    void frombuffer(std::ifstream& fp) { read(fp); }
};


template<typename T>
ivalarray<T> quicksort(const ivalarray<T>& a)
{
    // choose a pivot.
    const unsigned int n = a.size();
    
    if (n <= 1)
        return a;
    
    T pivot = a[0];
    ivalarray<T> a_le, a_gt;
    
    for(int i=1;i<n;++i) {
        if (a[i] <= pivot)
            a_le.push_back(a[i]);
        else
            a_gt.push_back(a[i]);
    }
    
    a_le = quicksort(a_le);
    a_gt = quicksort(a_gt);
    
    a_le.push_back(pivot);
    a_le.extend(a_gt);
    return a_le;
}

template<typename T>
void ivalarray<T>::sort(bool reverse)
{
    ivalarray<T> ar2 = quicksort(*this);
    if (!reverse)
        operator=(ar2);
    else {
        const int n = this->size();
        for(int i=0;i<n;++i)
            (*this)[i] = ar2[n-1-i];
        
    }
}



template<typename T>
void ivalarray<T>::zero()
{
    unsigned int n = this->size();
    if (!n) return;
    memset(this->raw(), 0, sizeof(T)*n);
}


template<typename T>
T ivalarray<T>::sum() const
{
    T x=0;
    const T* v = this->raw();
    const int n = this->size();
    for(int i=0;i<n;++i) x += v[i];
    return x;
}


template<typename T>
ivalarray<T> operator-(const ivalarray<T>& x, const ivalarray<T>& y)
{
    assert( x.size() == y.size() );
    ivalarray<T> z(x);
    z -= y;
    return z;
}

template<typename T>
ivalarray<T> operator+(const ivalarray<T>& x, const ivalarray<T>& y)
{
    assert( x.size() == y.size() );
    ivalarray<T> z(x);
    z += y;
    return z;
}


template<typename T>
bool operator==(const ivalarray<T>& x, const ivalarray<T>& y)
{
    if (x.size() != y.size() ) return false;
    return (!memcmp(x.raw(), y.raw(), x.size()*sizeof(T)) );
}


template<typename T>
void ivalarray<T>::write(std::ofstream& fp) const
{
    unsigned int n = this->size();
    size_t szt = sizeof(T);
    fp.write((char*)&n, sizeof(unsigned int));
    fp.write((char*)&szt, sizeof(size_t));
    fp.write((char*)this->raw(), n * sizeof(T) );
}

template<typename T>
void ivalarray<T>::read(std::ifstream& fp)
{
    unsigned int n;
    size_t szt;
    fp.read((char*)&n, sizeof(unsigned int));
    fp.read((char*)&szt, sizeof(size_t));
    if (szt != sizeof(T) ) {
        std::cerr << "type mismatch cannot read array<T>.\n";
        throw "seiorsu fdasf";
    }
    
    this->resize(n, false);
    fp.read((char*)this->raw(), sizeof(T)*n);
}

#endif
