/*
 * saveroutines.hpp
 *
 *  Created on: 20.05.2011
 *      Author: martinjg
 */

#ifndef SAVEROUTINES_HPP_
#define SAVEROUTINES_HPP_
#include"helpers/binary_io.h"
#include "BlockTypes.h"
#include "boost_typedefs.hpp"
#include <iostream>
using namespace BoostTypedefs;



template<typename TKey
void operator<<(std::ostream &file, TKey key){

}



//TODO: write the saveKey routine for all relevant KeyTypes, like UintVector, ...
//TODO: write the saveTensor routine for all relevant types
template<size_t Nb, class KeyType,class MappedType>
void saveBlockTensor(std::ostream &file,const BlockTensor<Nb,KeyType,MappedType> > &BT)
{
  BlockTensor<Nb,KeyType,MappedType> >::const_iterator it;
  uint size = BT.size();
  //write the size of the map:
  write_number_binary(file,size);
  for(it = BT.begin();it!=BT.end();it++){
    //write the phys ind
    write_number_binary(file,it->first._p);
    //save the size of the key ..
    write_number_binary(file,it->first.size());
    //and the key
    for(size_t i = 0;i<it->first.size();i++)
    saveKey(file,it->first[i]);
    //now save the tensor
    saveTensor(file,it->second);
  }
}

template<size_t Nb, class KeyType,class MappedType>
void loadBlockTensor(std::istream &file, BlockTensor<Nb,KeyType,MappedType> > &BT)
{
  BT.clear();
  uint size= read_number_binary<uint>(file);
  uint pind;

  MappedType tens;
  if (!file.good())
    throw BinaryIOException();
  for(uint i = 0;i<size;i++){
    pind = read_number_binary<uint>(file);
    //load the size 0f the key:
    size_t keysize = read_number_binary<size_t>(file);
    TensorKeyType<keysize,KeyType> TKey;
    TKey.setP(pind);

    //load the key ...
    for(size_t i = 0;i<keysize;i++)
      TKey[i] = read_key_binary<KeyType>(file);

    //now load the matrix
    loadTensor(file,tens);
    BT[TKey]=tens;
  }
}




void saveNumberTables(ostream &file,const number_map &nmap)
{
  number_map::const_iterator it;
  write_number_binary(file,(uint)(nmap.size()));
  for(it = nmap.begin();it!=nmap.end();it++)
    {
      write_vector_binary(file,it->first);
      write_number_binary(file,(uint)it->second.start());
      write_number_binary(file,(uint)it->second.size());
    }
}
void loadNumberTables(istream &file,number_map &nmap)
{
  nmap.clear();
  number_map::iterator it;
  uint sizemap,start,size;
  sizemap = read_number_binary<uint>(file);
  UintVectorType q;
  for(uint i = 0;i<sizemap;i++)
    {
      read_vector_binary(file,q);
      start = read_number_binary<uint>(file);
      size = read_number_binary<uint>(file);
      nmap.insert(pair<UintVectorType,Range>(q,Range(start,start+size)));
    }
}

void saveItoQmap(ostream &file,const ItoQmap &lmap)
{
  ItoQmap::const_iterator it;
  write_number_binary(file,(uint)(lmap.size()));
  for(it = lmap.begin();it!=lmap.end();it++)
    {
      write_number_binary(file,it->first);
      write_vector_binary(file,it->second);
    }
}

void loadItoQmap(istream &file, ItoQmap &lmap)
{
  lmap.clear();
  uint size,key;
  size = read_number_binary<uint>(file);
  UintVectorType q;
  for(uint i = 0;i<size;i++)
    {
      key = read_number_binary<uint>(file);
      read_vector_binary(file,q);
      lmap.insert(pair<uint,UintVectorType>(key,q));
    }
}

#endif /* SAVEROUTINES_HPP_ */
