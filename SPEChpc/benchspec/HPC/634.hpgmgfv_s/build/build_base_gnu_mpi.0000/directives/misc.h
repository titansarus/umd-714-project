/*
# Copyright (c) 2015, NVIDIA CORPORATION. All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
#  * Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
#  * Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in the
#    documentation and/or other materials provided with the distribution.
#  * Neither the name of NVIDIA CORPORATION nor the names of its
#    contributors may be used to endorse or promote products derived
#    from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ``AS IS'' AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
# OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#define max(x,y) ((x) > (y) ? (x) : (y))

// this kernel zeros out the grid
void directives_zero_vector(level_type *level, int component_id)
{
  int blockStart = 0;
  int blockEnd = level->num_my_blocks;
  if (blockEnd <= 0) return;

#if defined (SPEC_OPENMP_TARGET)
# if defined (SPMD_MODE)
#  pragma omp target teams map(to: level[:LEVEL_ITEMS])
   {
#   pragma omp parallel
   {
   TEAM_DISTRIBUTE_LOOP(level->num_my_blocks, blockStart, blockEnd)
# else
#  pragma omp target teams distribute map(to: level[:LEVEL_ITEMS])
# endif
#elif defined (SPEC_OPENACC)
# if defined (SPEC_MANAGED_MEMORY)
#  pragma acc parallel loop gang copyin(level[:1])
# else
#  pragma acc parallel loop gang present(level) default(present)
# endif
#endif
  for(int block = blockStart; block < blockEnd; block++) {

    const int box = level->my_blocks[block].read.box;
          int ilo = level->my_blocks[block].read.i;
          int jlo = level->my_blocks[block].read.j;
          int klo = level->my_blocks[block].read.k;
          int ihi = level->my_blocks[block].dim.i + ilo;
          int jhi = level->my_blocks[block].dim.j + jlo;
          int khi = level->my_blocks[block].dim.k + klo;
    const int jStride = level->my_boxes[box].jStride;
    const int kStride = level->my_boxes[box].kStride;
    const int  ghosts = level->my_boxes[box].ghosts;
    const int     dim = level->my_boxes[box].dim;

    // expand the size of the block to include the ghost zones...
    if(ilo<=  0)ilo-=ghosts;
    if(jlo<=  0)jlo-=ghosts;
    if(klo<=  0)klo-=ghosts;
    if(ihi>=dim)ihi+=ghosts;
    if(jhi>=dim)jhi+=ghosts;
    if(khi>=dim)khi+=ghosts;

    double * __restrict__ grid = level->my_boxes[box].vectors[component_id] + ghosts*(1+jStride+kStride);

#if defined (SPEC_OPENMP_TARGET)
# if defined (SPMD_MODE)
#  pragma omp for collapse(3)
# else
#  pragma omp parallel for collapse(3)
# endif
#elif defined (SPEC_OPENACC)
# pragma acc loop vector collapse(3)
#endif
    for (int k = klo; k < khi; k++) {
      for (int j = jlo; j < jhi; j++) {
	for (int i = ilo; i < ihi; i++) {
	  const int ijk = i + j*jStride + k*kStride;
	  grid[ijk] = 0.0;
	}
      }
    }
  }

#if defined (SPEC_OPENMP_TARGET) && defined(SPMD_MODE)
  } // end of omp parallel
  } // end of omp target teams
#endif
}

// this kernel provides a generic axpy/mul implementation:
// if mul_vectors = 1: c = scale_a * a * b
// if mul_vectors = 0: c = scale_a * a + scale_b * b + shift_a
void directives_axpy_vector(level_type *level, int id_c, double scale_a, double shift_a, double scale_b, int id_a, int id_b, int mul_vectors)
{
  int blockStart = 0;
  int blockEnd = level->num_my_blocks;
#if defined (SPEC_OPENMP_TARGET)
# if defined (SPMD_MODE)
#  pragma omp target teams map(to: level[:LEVEL_ITEMS])
   {
#   pragma omp parallel
   {
   TEAM_DISTRIBUTE_LOOP(level->num_my_blocks, blockStart, blockEnd)
# else
#  pragma omp target teams distribute map(to: level[:LEVEL_ITEMS])
# endif
#elif defined (SPEC_OPENACC)
# if defined (SPEC_MANAGED_MEMORY)
#  pragma acc parallel loop gang copyin(level[:1])
# else
#  pragma acc parallel loop gang present(level) default(present)
# endif
#endif
  for(int block = blockStart; block < blockEnd; block++) {

    const int box = level->my_blocks[block].read.box;
          int ilo = level->my_blocks[block].read.i;
          int jlo = level->my_blocks[block].read.j;
          int klo = level->my_blocks[block].read.k;
          int ihi = level->my_blocks[block].dim.i + ilo;
          int jhi = level->my_blocks[block].dim.j + jlo;
          int khi = level->my_blocks[block].dim.k + klo;
    const int jStride = level->my_boxes[box].jStride;
    const int kStride = level->my_boxes[box].kStride;
    const int  ghosts = level->my_boxes[box].ghosts;
    double * __restrict__ grid_c = level->my_boxes[box].vectors[id_c] + ghosts*(1+jStride+kStride);
    double * __restrict__ grid_a = level->my_boxes[box].vectors[id_a] + ghosts*(1+jStride+kStride);
    double * __restrict__ grid_b = level->my_boxes[box].vectors[id_b] + ghosts*(1+jStride+kStride);

#if defined (SPEC_OPENMP_TARGET)
# if defined (SPMD_MODE)
#  pragma omp for collapse(3)
# else
#  pragma omp parallel for collapse(3)
# endif
#elif defined (SPEC_OPENACC)
# pragma acc loop vector collapse(3)
#endif
    for (int k = klo; k < khi; k++) {
      for (int j = jlo; j < jhi; j++) {
	for (int i = ilo; i < ihi; i++) {
	  const int ijk = i + j*jStride + k*kStride;
	  if (mul_vectors)
	    grid_c[ijk] = scale_a*grid_a[ijk]*grid_b[ijk];
	  else
	    grid_c[ijk] = scale_a*grid_a[ijk] + scale_b*grid_b[ijk] + shift_a;
	}
      }
    }
  }

#if defined (SPEC_OPENMP_TARGET) && defined(SPMD_MODE)
  } // end of omp parallel
  } // end of omp target teams
#endif
}

// simple coloring kernel, see misc.c for details
void directives_color_vector(level_type *level, int id_a, int colors_in_each_dim, int icolor, int jcolor, int kcolor)
{
  int blockStart = 0;
  int blockEnd = level->num_my_blocks;
  if (blockEnd <= 0) return;

#if defined (SPEC_OPENMP_TARGET)
# if defined (SPMD_MODE)
#  pragma omp target teams map(to: level[:LEVEL_ITEMS])
   {
#   pragma omp parallel
   {
   TEAM_DISTRIBUTE_LOOP(level->num_my_blocks, blockStart, blockEnd)
# else
#  pragma omp target teams distribute map(to: level[:LEVEL_ITEMS])
# endif
#elif defined (SPEC_OPENACC)
# if defined (SPEC_MANAGED_MEMORY)
#  pragma acc parallel loop gang copyin(level[:1])
# else
#  pragma acc parallel loop gang present(level) default(present)
# endif
#endif
  for(int block = blockStart; block < blockEnd; block++) {

    const int box = level->my_blocks[block].read.box;
          int ilo = level->my_blocks[block].read.i;
          int jlo = level->my_blocks[block].read.j;
          int klo = level->my_blocks[block].read.k;
          int ihi = level->my_blocks[block].dim.i + ilo;
          int jhi = level->my_blocks[block].dim.j + jlo;
          int khi = level->my_blocks[block].dim.k + klo;
    const int boxlowi = level->my_boxes[box].low.i;
    const int boxlowj = level->my_boxes[box].low.j;
    const int boxlowk = level->my_boxes[box].low.k;
    const int jStride = level->my_boxes[box].jStride;
    const int kStride = level->my_boxes[box].kStride;
    const int  ghosts = level->my_boxes[box].ghosts;
    double * __restrict__ grid = level->my_boxes[box].vectors[id_a] + ghosts*(1+jStride+kStride);

#if defined (SPEC_OPENMP_TARGET)
# if defined (SPMD_MODE)
#  pragma omp for collapse(3)
# else
#  pragma omp parallel for collapse(3)
# endif
#elif defined (SPEC_OPENACC)
# pragma acc loop vector collapse(3)
#endif
    for (int k = klo; k < khi; k++) {
      for (int j = jlo; j < jhi; j++) {
	for (int i = ilo; i < ihi; i++) {
	  double sk=0.0;if( ((k+boxlowk+kcolor)%colors_in_each_dim) == 0 )sk=1.0; // if colors_in_each_dim==1 (don't color), all cells are set to 1.0
	  double sj=0.0;if( ((j+boxlowj+jcolor)%colors_in_each_dim) == 0 )sj=1.0;
	  double si=0.0;if( ((i+boxlowi+icolor)%colors_in_each_dim) == 0 )si=1.0;
	  const int ijk = i + j*jStride + k*kStride;
	  grid[ijk] = si*sj*sk;
	}
      }
    }
  }

#if defined (SPEC_OPENMP_TARGET) && defined(SPMD_MODE)
  } // end of omp parallel
  } // end of omp target teams
#endif
}

// 0: summation, 1: maximum absolute
double directives_reduction(level_type *level, int id, int red_type)
{
#if defined(__PGI) && (__PGIC__ <= 19) && (__PGIC_MINOR__ <= 6)
# warning "The function directives_reduction specially handles old PGI compilers"
#endif
  double res_max = 0.0, res_add = 0.0;

#if defined (SPEC_OPENMP_TARGET)
# pragma omp target teams distribute map(to: level[:LEVEL_ITEMS]) map(tofrom: res_max, res_add) reduction(max:res_max) reduction(+:res_add)
#elif defined (SPEC_OPENACC)
# if defined (SPEC_MANAGED_MEMORY)
#  pragma acc parallel loop gang copyin(level[:1]) copy(res_max, res_add) reduction(max:res_max) reduction(+:res_add)
# else
#  pragma acc parallel loop gang present(level) default(present) copy(res_max, res_add) reduction(max:res_max) reduction(+:res_add)
# endif
#endif
  for(int block = 0; block < level->num_my_blocks; block++) {
    const int box = level->my_blocks[block].read.box;
          int ilo = level->my_blocks[block].read.i;
          int jlo = level->my_blocks[block].read.j;
          int klo = level->my_blocks[block].read.k;
          int ihi = level->my_blocks[block].dim.i + ilo;
          int jhi = level->my_blocks[block].dim.j + jlo;
          int khi = level->my_blocks[block].dim.k + klo;
    const int jStride = level->my_boxes[box].jStride;
    const int kStride = level->my_boxes[box].kStride;
    const int  ghosts = level->my_boxes[box].ghosts;
    double * __restrict__ grid = level->my_boxes[box].vectors[id] + ghosts*(1+jStride+kStride);
    double thread_val_max = 0.0, thread_val_add = 0.0;
    // accumulate per thread first (multiple elements)


#if defined (SPEC_OPENMP_TARGET)
# pragma omp parallel for collapse(3) reduction(max:thread_val_max) reduction(+:thread_val_add)
#elif defined (SPEC_OPENACC)
# if defined(__PGI) && (__PGIC__ <= 19) && (__PGIC_MINOR__ <= 6)
    // worker parallelism leads to inconsistent reduction values
#  pragma acc loop vector collapse(3) reduction(max:thread_val_max) reduction(+:thread_val_add)
# else
#  pragma acc loop worker vector collapse(3) reduction(max:thread_val_max) reduction(+:thread_val_add)
# endif
#endif
   for (int k = klo; k < khi; k++) {
     for (int j = jlo; j < jhi; j++) {
       for (int i = ilo; i < ihi; i++) {
          const int ijk = i + j*jStride + k*kStride;
          double val = grid[ijk];
          switch (red_type) {
          case 0: thread_val_add += val; break;
          case 1: thread_val_max = max(thread_val_max, fabs(val)); break;
          }
        }
      }
   }
  switch (red_type) {
    case 0:
      res_add += thread_val_add;
      break;
    case 1:
      res_max = max(res_max, thread_val_max);
      break;
    }
  } // end of block loop
  if (red_type == 0){
    return res_add;
  } else {
    return res_max;
  }
}

void directives_scale_vector(level_type *d_level, int id_c, double scale_a, int id_a)
{
  int grid = d_level->num_my_blocks;
  if (grid <= 0) return;

  directives_axpy_vector(d_level, id_c, scale_a, 0.0, 0.0, id_a, id_a,0);
}

void directives_shift_vector(level_type *d_level, int id_c, double shift_a, int id_a)
{
  int grid = d_level->num_my_blocks;
  if (grid <= 0) return;

  directives_axpy_vector(d_level, id_c, 1.0, shift_a, 0.0, id_a, id_a,0);
}

void directives_mul_vectors(level_type *d_level, int id_c, double scale, int id_a, int id_b)
{
  int grid = d_level->num_my_blocks;
  if (grid <= 0) return;

  directives_axpy_vector(d_level, id_c, scale, 0.0, 0.0, id_a, id_b,1);
}

void directives_add_vectors(level_type *d_level, int id_c, double scale_a, int id_a, double scale_b, int id_b)
{
  int grid = d_level->num_my_blocks;
  if (grid <= 0) return;

  directives_axpy_vector(d_level, id_c, scale_a, 0.0, scale_b, id_a, id_b, 0);
}

//HA: not used
double directives_sum(level_type *d_level, int id)
{
  int grid = d_level->num_my_blocks;
  if (grid <= 0) return 0.0;

  double h_res;
  h_res = directives_reduction(d_level, id, 0);
  return h_res;
}

double directives_max_abs(level_type *d_level, int id)
{
  int grid = d_level->num_my_blocks;
  if (grid <= 0) return 0.0;

  double h_res;
  h_res = directives_reduction(d_level, id, 1);
  return h_res;
}
