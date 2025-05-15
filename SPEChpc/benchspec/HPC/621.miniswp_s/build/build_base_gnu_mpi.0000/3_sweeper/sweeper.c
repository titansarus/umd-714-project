/*---------------------------------------------------------------------------*/
/*!
 * \file   sweeper.c
 * \author Wayne Joubert, Veronica G. Vergara Larrea
 * \date   Sep. 3, 2019
 * \brief  Definitions for performing a sweep.
 * \note   Copyright (C) 2019 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#include "sweeper.h"

#ifdef SWEEPER_SIMPLE
#include "sweeper_simple_c.h"
#endif

#ifdef SWEEPER_TILEOCTANTS
#include "sweeper_tileoctants_c.h"
#endif

#ifdef SWEEPER_KBA
#include "sweeper_kba_c.h"
#endif

#ifdef SWEEPER_KBA_ACC
  #include "sweeper_gpu_c.h"
  #include "sweeper_kba_c.h"
#endif

#ifdef SWEEPER_KBA_OPENMP_TARGET
  #include "sweeper_gpu_c.h"
  #include "sweeper_kba_c.h"
#endif

#ifdef SWEEPER_ACC
#include "sweeper_gpu_c.h"
#endif

#ifdef SWEEPER_OPENMP_TARGET
#include "sweeper_gpu_c.h"
#endif

/*---------------------------------------------------------------------------*/
