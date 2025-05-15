/* Copyright (C) 2016-2017 Ludwig Schneider

 This file is part of SOMA.

 SOMA is free software: you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 SOMA is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public License
 along with SOMA.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef MONOMER_H
#define MONOMER_H
//! \file monomer.h
//! \brief Collection of code for Monomer struct

#include "soma_config.h"



/*! \brief Monomer struct contains spatial position and type.
  \warning The bit pattern of the type \a w is an int in a soma_scalar_t variable.*/
typedef struct{
    soma_scalar_t x; /*!<\brief X-coordinate */
    soma_scalar_t y; /*!<\brief Y-coordinate */
    soma_scalar_t z; /*!<\brief Z-coordinate */
    }Monomer;


#endif//MONOMER_H
