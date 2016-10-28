/*
 * types.hpp
 *
 *  Created on: 01.05.2011
 *      Author: maga
 */

#ifndef TYPES_HPP_
#define TYPES_HPP_
#include <map>
#include "../boost_typedefs.hpp"



typedef std::multimap<UintVectorType, uint> particle_map;
typedef std::map<UintVectorType, Range> number_map;
typedef std::map<uint,UintVectorType> ItoQmap;
typedef std::multimap<UintVectorType, uint> QtoImulti_map;



#endif /* TYPES_HPP_ */
