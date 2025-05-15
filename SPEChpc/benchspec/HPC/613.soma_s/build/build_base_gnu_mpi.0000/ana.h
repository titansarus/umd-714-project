/* Copyright (C) 2016 Ludwig Schneider
   Copyright (C) 2016 Ulrich Welling
   Copyright (C) 2016 Marcel Langenberg
   Copyright (C) 2016 Fabien Leonforte
   Copyright (C) 2016 Juan Orozco
   Copyright (C) 2016 Yongzhi Ren
   Copyright (C) 2016 N. Harshavardhan Reddy

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

/*! \file ana.h
  \brief Functions to calculate observables of the configurations.
*/

#ifndef SOMA_ANA_H
#define SOMA_ANA_H

#include "specDefs.h"

#include "soma_config.h"
#include "stdint.h"

struct Phase;

/*!>\brief calculate the end-to-end distance for the polymers of the phase
  \param p Phase configuration to analyze
  \returns global \f$ Re \f$
  \param result pointer to return the result. The length of this array
  should be p->n_poly_types * 4. Order is Re2 Rex Rey Rez for each
  poly type.
*/
void calc_Re(const struct Phase * p, soma_scalar_t *const result);

/*!>\brief calculate the variance between the current density field and the density field at a past time
  \param p Phase configuration to analyze
  \returns global \f$ dvar \f$
  \param dvar pointer to return density variance
*/
void calc_dvar(const struct Phase * p, soma_scalar_t * dvar);

/*!>\brief calculate the gyration radius for the polymers of the phase
  \param p Phase configuration to analyze
  \param result pointer to return the result. The length of this array
  should be p->n_poly_types * 4. Order is Rg2 Rgx Rgy Rgz for each
  poly type.
  \returns global \f$ Rg^2 \f$ and its squared components independently of the polymer type
*/
void calc_Rg(const struct Phase * p, soma_scalar_t *const result);

/*!>\brief calculate the monomer and chain mean square displacement
  \param p Phase configuration to analyze
  \param result array to store the result. Size should be 8*n_poly_types.
  The order is: x-component of MSD, y-component of MSD,z-component of MSD, MSD in 3D and
  x-componenet of mass center MSD, y-componenet of mass center MSD, z-componenet of mass center MSD,
  mass center MSD.
*/
void calc_MSD(const struct Phase * p, soma_scalar_t *const result);

/*! \brief calculate bond anisotropy

  \param p Phase configuration to analyze.
  \param result Array to store the tensor. xx yy zz xy xz yz. Length of array must be 6*n_poly_types.
  \post the memory of where the pointers point is set to the result of calculation.
*/
void calc_anisotropy(const struct Phase * p, soma_scalar_t *const result);

//! \brief calculate the current acceptance ratio.
//! \param p System to analyze
//! \param acc_ratio Pointer to initialized soma_scalar_t where the result is going to be stored.
//! \note After every call the result is reset.
void calc_acc_ratio(struct Phase*const p,soma_scalar_t*const acc_ratio);

//! \brief calculate the non-bonded energy for each particle type
//! \param p System to analyze
//! \param non_bonded_energy Pointer to array of size p->n_types to store the result
void calc_non_bonded_energy(const struct Phase*const p, soma_scalar_t*const non_bonded_energy);

//! \brief calculate the bonded energy for each bond type
//! \param p System to analyze
//! \param bonded_energy Pointer to array of size NUMBER_SOMA_BOND_TYPES to store the result
void calc_bonded_energy(const struct Phase*const p, soma_scalar_t*const bonded_energy);

#endif//SOMA_ANA_H
