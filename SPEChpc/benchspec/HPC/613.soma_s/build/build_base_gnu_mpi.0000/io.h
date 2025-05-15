/* Copyright (C) 2016 Ludwig Schneider
   Copyright (C) 2016 Ulrich Welling
   Copyright (C) 2016 Marcel Langenberg

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
#ifndef SOMA_IO_H
#define SOMA_IO_H
/*! \file io.h
  \brief Header file for all functions, that handle with input and ouput
   operations of SOMA.
*/

/*!
  \brief Read a configuration from disk.

  And initializes all memory fields.
  \param p Phase struct which will define the state of the system
  after the call.
  \param filename Relative or absolute path the configuration file to
  read.
  \pre Uninitialied Phase* \a p struct (except the int init_MPI(Phase*) call).
  \post Initialized \a p points to an fully initialized configuration.
  \return Error code. Returns not equal to zero if an error occured.
*/

#include "specDefs.h"

#include "soma_config.h"

struct Phase;

/*!\brief Writes the current configuration to the disk.
  \param p Pointer to a fully initialized configuration.
  \param filename Relative or absolute path the configuration file to
  read.
  \warning The call will overwrite all previous information contained
  in the file "filename" is existing.
  \post Current state is written to the given file.
  \return Error code. Return not equal to zero if an error occured.
 */
int write_config(const struct Phase*const p,const char*const filename);

int init_cripple(struct Phase*p,const unsigned int n);

int write_state(struct Phase*p);

int compare_state(struct Phase*p,soma_scalar_t epsilon);

int screen_output(struct Phase*const p,const unsigned int Nsteps);
#endif//SOMA_IO_H
