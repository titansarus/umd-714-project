/* Copyright (C) 2016-2017 Ludwig Schneider
   Copyright (C) 2016 Ulrich Welling
   Copyright (C) 2016 Marcel Langenberg
   Copyright (C) 2016 Fabien Leonforte
   Copyright (C) 2016 Juan Orozco
   Copyright (C) 2016 Yongzhi Ren

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
#ifndef PHASE_H
#define PHASE_H

struct IndependetSets;

#include <time.h>
#include <assert.h>
#include <stdlib.h>
#include "mpiroutines.h"
#include "cmdline.h"
#include "soma_util.h"
#include "autotuner.h"
#include "polymer.h"
#include "allocator.h"

//! \file phase.h
//! \brief All relevant aspects for the struct Phase

/*! \brief All relevant information for a system configuration.

 * \note Arrays with higher dimensions are usually unrolled as linear
 * arrays for storage. Such unrolled arrays have pointers, which are
 * name with a suffix "_unified". For accessing such pointers, you can
 * use helper function, which are provided by mesh.h.

 * \warning Never access "_unified" pointer arrays, without the helper
 * functions in mesh.h. The memory layout of the array might change
 * with future releases.

 */
typedef struct Phase{
    unsigned int reference_Nbeads; /*!< \brief number of reference beads for the model polymer */
    soma_scalar_t harmonic_normb; //!< Harmonic energy scale (function of spring constant) const. at runtime.
    soma_scalar_t harmonic_normb_variable_scale; //!< different harmonic energy scale (function of spring constant) const. at runtime.
    unsigned int n_types;  /*!<\brief number of monomer types */

    /*! \brief \f$\chi N\f$
      2D matrix with monomer type Flory-Huggins interactions, trace
      contains compressibility \f$\kappa_i\f$.
    */

    soma_scalar_t * xn;
    uint64_t n_polymers; /*!< \brief \#polymers in the configuration. (Local on the MPI node.)*/
    uint64_t n_polymers_storage; /*!< \brief Storage space for polymers.*/
    uint64_t n_polymers_global; /*!< \brief \#polymers in the global configuration.*/
    Polymer * polymers; /*!< \brief pointer to array of polymers */

    //uint16_t **fields; /*!< \brief n_types fields in 3D, mimics the DENSITY NOT normalized to 1, this has to be done in the omega_field calculation*/
    uint16_t * fields_unified; /*!< \brief one pointer that points to the construct of p->n_types * p->n_cells of fields */
    uint16_t * old_fields_unified; /*!< \brief one pointer that points to the construct of p->n_types * p->n_cells of old fields for density variance calculations*/
    uint32_t *fields_32; //!< \brief linear 32 bit version of the fields. This is required for GPU-simulation, because no 16-bit atomic operations are available.

/*! \brief array of shape of field, containing the information if
 * this cell is free space or not. == 0  is a free cell, !=0 is a forbidden cell.
 * If the pointer is set to NULL, the entire simulation box, contains no forbidden cells.
 *
 * \warning area51 is a unified pointer, but is not type specific, so
 * specify for access always type=0.
 * \warning Before any access to area51, check for NULL.
 */
    uint8_t * area51;

    //soma_scalar_t ** omega_field; /*!< \brief calculates the omega fields according to the Hamiltonian*/
    soma_scalar_t * omega_field_unified;  /*!< \brief calculates the omega fields according to the Hamiltonian, unified access*/
    //soma_scalar_t ** external_field; /*!< \brief external fields that act on the polymers, one field per type */
    soma_scalar_t * external_field_unified; /*!< \brief one pointer that points to the construct of p->n_types * p->n_cells of external_fields */
    soma_scalar_t * tempfield; /*!< \brief a temporal storage for intermediate field calculations, used to save the complete density */
    uint64_t *num_bead_type; /*!< \brief stores the number of beads of a specific type*/
    uint64_t *num_bead_type_local; /*!< \brief stores the number of beads of a specific type locally (for this mpi-core)*/

    soma_scalar_t * A; /*!< \brief stores the diffusion constants for each type */
    soma_scalar_t * R; /*!< \brief stores the derived dR for the diffusion constant */
    //!  Mobility of the center of mass for all polymer types.
    //!  Length is p->n_poly_types. NULL if no mobility wanted.
    soma_scalar_t * cm_a;
    //! Tries of the cm moves
    unsigned int n_tries_cm;
    //! Accepted moves of the cm
    unsigned int n_acc_cm;

    soma_scalar_t *field_scaling_type; /*!< \brief stores the scaling factor according to the density */

    unsigned int time; /*!< \brief MC steps into the simulation */
    uint64_t num_all_beads; //!< Number of all monomer/beads in the global system
    uint64_t num_all_beads_local; //!< Number of all monomer/beads on the local MPI core

    unsigned int nx; /*!< \brief x-spatial discretization */
    unsigned int ny; /*!< \brief y-spatial discretization */
    unsigned int nz; /*!< \brief z-spatial discretization */
    uint64_t n_cells; /*!< \brief number of cells in the field */
    soma_scalar_t Lx; /*!< \brief x-spatial dimensions in units of \f$ Re_0 \f$ */
    soma_scalar_t Ly; /*!< \brief y-spatial dimensions in units of \f$ Re_0 \f$ */
    soma_scalar_t Lz; /*!< \brief z-spatial dimensions in units of \f$ Re_0 \f$ */
    soma_scalar_t iLx; /*!< \brief inverse x-spatial dimensions in units of \f$ Re_0 \f$ */
    soma_scalar_t iLy; /*!< \brief inverse y-spatial dimensions in units of \f$ Re_0 \f$ */
    soma_scalar_t iLz; /*!< \brief inverse z-spatial dimensions in units of \f$ Re_0 \f$ */

    // Variables for statistics/ analytics
    unsigned long int n_moves; /*!< \brief total number of moves */
    unsigned long int n_accepts; /*!< \brief accepted moves */

    soma_scalar_t msd_old; /*!< \brief store the MSD from n steps ago */
    //! two dimensional array, (npoly_types x 2) describing the start
    //! and end monomer of each polymer type.
    unsigned int* end_mono;

    Info_MPI info_MPI; /*!< \brief  stores all information regarding MPI */

    //! Number of different polymer architectures.
    unsigned int n_poly_type;
    //!Offset where to start/end reading the architecture data.
    //!
    //! Array size n_poly_type.
    int* poly_type_offset;
    //! Length of the poly_arch array
    unsigned int poly_arch_length;
    //! Array to all store information about the polymer architecture.
    //!
    //! \warning The useage of this array is \a not straight forward. Please read the doc carefully.
    //!
    //! The array contains three different information about the polymer architecture.
    //!  - The number of monomers per polymer. This information is always located at poly_arch[poly_type_offset[poly_type]].
    //!    You can read this information directly:
    //!    \code const unsigned int N = poly_arch[poly_type_offset[poly_type]]; \endcode
    //!    Writing is not recomeneded.
    //!  - The information about the bonded neighbours is stored in a so called bond_offset_lists.<br>
    //!    The second info tells you, where you find the start of such a
    //!    list for your monomer of interest. We might call this the
    //!    Monomer region. And it contains the offset to the bond list
    //!    as well as the Monomer type.
    //!    \note The Monomer region starts for a polymer_type at:
    //!    \code poly_type_offset[poly_type] +1 \endcode
    //!
    //!    So to find the index where to start reading the bond list you do:
    //!    \code
    //!    const unsigned int start = get_bondlist_offset(
    //!                                  p->poly_arch[poly_type_offset[poly_type] + mono_index + 1] );
    //!    \endcode
    //!    If this start is < 0 the meaning is that the corresponding
    //!    Monomer has no bonded neigbours.<br>
    //!    The other information is the particle type of your monomer of interest.
    //!    You can access this information with the following snippet:
    //!    \code
    //!    const unsigned int type = get_particle_type(
    //!                                 p->poly_arch[poly_type_offset[poly_type] + mono_index + 1] );
    //!    \endcode
    //!    If you want to create an element of the monomer region in
    //!    poly_arch, you can call get_info_bl().
    //!  - The final information is the actual bondlist. Each element of
    //!    the bond list region carries 3 bits of information.
    //!    - The first is the \a end flag, which tells you wheter this
    //!    is last element in the bond list. You can access the ith element of your bondlist via:
    //!    \code const unsigned int end = get_end( p->poly_arch[ start + i] );\endcode
    //!    - The second is the \a offset tells you the offset you apply
    //!    to your monomer index to obtain the neighbour index. (intentional int!)
    //!    \code const int offset = get_offset( p->poly_arch[ start + i]); \endcode
    //!    - The last info is the bond type of the corresponding bond.
    //!    \code const unsigned int bond_type = get_bond_type( p->poly_arch[start + i] );\endcode
    //! To see, how you can iterate all bonds of a specific monomer, you
    //! can see in the following example:
    //! \code
    //! const int start = get_bondlist_offset(
    //!             p->poly_arch[p->poly_type_offset[p->polymers[ipoly].type] + ibead + 1]);
    //! const unsigned int mono_type = get_particle_type(p->poly_arch[
    //!                             p->poly_type_offset[p->polymers[ipoly].type] + ibead +1]);
    //! if(start > 0){
    //!   int i = start;
    //!   unsigned int end;
    //!   do{
    //!     const uint32_t info = p->poly_arch[i++];
    //!     end = get_end( info );
    //!     const unsigned int bond_type = get_bond_type(info);
    //!     const int offset = get_offset( info);
    //!     const unsigned int neighbour_id = ibead + offset;
    //!     const unsigned int neighbour_type = get_particle_type(p->poly_arch[
    //!                             p->poly_type_offset[p->polymers[ipoly].type] + neighbour_id +1]);
    //!     const unsigned int jbead = neighbour_id;
    //!     //Do stuff
    //!   }while(end==0);
    //!  \endcode
    uint32_t* poly_arch;

    //! \brief Struct containing the command line arguments.
    struct som_args args;

    //! \brief clock value close to the start of SOMA.
    //!
    //! Not exactly at the start, but at the beginnig of init_values().
    time_t start_clock;
    //! \brief Start time of this simulation.
    unsigned int start_time;

    //! Indiciates, whether bead data have been read from the input file.
    bool bead_data_read;

    //! Maximum distance a particle can move, without accidentically
    //!passing trough an area51 wall.
    //!
    //! This is equivalent to min( p->Lx/p->n_x , p->Ly/p->n_y ,
    //! p->Lz/p->n_z). Some versions may introduce some safety
    //! parameter 0 < s < 1 and multiply it with the min.
    soma_scalar_t max_safe_jump;

    //! Selected hamiltonian for non-bonded interaction.
    enum Hamiltonian hamiltonian;

    //! Array of independet sets. If NULL no sets are stored, the move type is not available.
    //!
    //! Length of the array is the number of polymer types.
    struct IndependetSets*sets;
    //! Number of max members in an independet set. Used for the length of polymer states.
    unsigned int max_set_members;
    //! Max number of sets for all polymer types.
    unsigned int max_n_sets;
    //! Autotuner for the Monte-Carlo kernels.
    Autotuner mc_autotuner;
    //! Autotuner for the Center of Mass Monte-Carlo kernels.
    Autotuner cm_mc_autotuner;
    //! Last measured TPS
    double tps_elapsed_time;
    //! Measured TPS counter
    unsigned int tps_elapsed_steps;

    struct Allocator* allocator;

    }Phase;

/*! \brief Initializes the values additional after the input init by the read*() functions.

  \note This function should be called before any type of simulaton,
  because several parameters from the input are used to fully initialize the Phase.
  \param p Phase to initialize.
  \return error code.
*/
int init_phase(struct Phase*const p);

//! Initialize and copy in all data to DEVICE
//! \param p Pointer to phase to copy in
//! \return Errorcode
int copyin_phase(struct Phase*const p, bool copyMSD);

/*! \brief Copys the data from the device back to the host
  \param p System in which the host data is updated from the device
  \return Errorcode
 */
int update_self_phase(const Phase * const p);

//! Copy all data out from the DEVICE
//! \param p Pointer to phase to copy out
//! \return Errorcode
int copyout_phase(struct Phase*const p);

/*! \brief Deallocate the memory used by p.
  \param p Initialized configuration.
  \return Errorcode.
  \post p points to an invalid configuration. Do not use it for any other call afterwards.
*/
int free_phase(struct Phase*const p);


#endif//PHASE_H
