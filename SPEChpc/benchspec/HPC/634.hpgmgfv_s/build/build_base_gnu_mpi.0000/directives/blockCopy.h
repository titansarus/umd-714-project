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

#define READ(i)	read[i]

void directives_copy_block(level_type *level, int id, communicator_type *exchange_ghosts, int block_type)
{
#if defined (SPEC_OPENMP_TARGET) && defined(SPEC_MANAGED_MEMORY)
/* Two part map because "exchange_ghosts" aliases part of "level":
   operators/exchange_boundary.c: device_copy_block(level,id,&level->exchange_ghosts[shape],0);
   operators/exchange_boundary.c: device_copy_block(level, id, &level->exchange_ghosts[shape], 1);
   operators/exchange_boundary.c: device_copy_block(level,id,&level->exchange_ghosts[shape],2);
   operators/restriction.c: device_copy_block(level_c,id_c,&level_c->restriction[restrictionType],2);

   We could alternatively use target map(to:level[:1], exchange_ghosts[:0]), however,
   this would rule out a future use of this function without aliased data. */
# pragma omp target enter data map(to:level[:1])
# pragma omp target enter data map(to:exchange_ghosts[:1])
#endif

  int blockStart = 0;
  int blockEnd = exchange_ghosts->num_blocks[block_type];
#if defined (SPEC_OPENMP_TARGET)
# if defined (SPMD_MODE)
#  pragma omp target teams map(to:level[:0], exchange_ghosts[:0])
   {
#   pragma omp parallel
   {
   TEAM_DISTRIBUTE_LOOP(exchange_ghosts->num_blocks[block_type], blockStart, blockEnd)
# else
#  pragma omp target teams distribute map(to:level[:0], exchange_ghosts[:0])
# endif
#elif defined (SPEC_OPENACC)
# if defined(SPEC_MANAGED_MEMORY)
#  pragma acc parallel loop gang copyin(exchange_ghosts[:1], level[:1])
# else
#  pragma acc parallel loop gang present(exchange_ghosts, level) default(present)
# endif
#endif
  for (int blockIdx=blockStart; blockIdx < blockEnd; ++blockIdx) {

  // one CUDA thread block operates on one HPGMG tile/block
#if defined(FIX_PGI_STRUCT_ASSIGNMENT_BUG)
  // Explicitly assign data for each struct field.
  // This workaround is needed to successfully compile with PGI compiler.
  // The original default code fails with the following error:
  // PGC-S-0155-Compiler failed to translate accelerator region (see -Minfo messages):
  // Unknown variable reference (directives/operators-directives.fv4.c: 45)
  //
  // The bug has been reported:
  // FS#28014: NERSC/HPC2020 ACC fails with Internal accelerator LILI consistency check failure.
  blockCopy_type block;
  block.dim.i = exchange_ghosts->blocks[block_type][blockIdx].dim.i;
  block.dim.j = exchange_ghosts->blocks[block_type][blockIdx].dim.j;
  block.dim.k = exchange_ghosts->blocks[block_type][blockIdx].dim.k;
  block.read.i = exchange_ghosts->blocks[block_type][blockIdx].read.i;
  block.read.j = exchange_ghosts->blocks[block_type][blockIdx].read.j;
  block.read.k = exchange_ghosts->blocks[block_type][blockIdx].read.k;
  block.read.jStride = exchange_ghosts->blocks[block_type][blockIdx].read.jStride;
  block.read.kStride = exchange_ghosts->blocks[block_type][blockIdx].read.kStride;
  block.write.i = exchange_ghosts->blocks[block_type][blockIdx].write.i;
  block.write.j = exchange_ghosts->blocks[block_type][blockIdx].write.j;
  block.write.k = exchange_ghosts->blocks[block_type][blockIdx].write.k;
  block.write.jStride = exchange_ghosts->blocks[block_type][blockIdx].write.jStride;
  block.write.kStride = exchange_ghosts->blocks[block_type][blockIdx].write.kStride;
  block.read.ptr = exchange_ghosts->blocks[block_type][blockIdx].read.ptr;
  block.write.ptr = exchange_ghosts->blocks[block_type][blockIdx].write.ptr;
  block.read.box = exchange_ghosts->blocks[block_type][blockIdx].read.box;
  block.write.box = exchange_ghosts->blocks[block_type][blockIdx].write.box;
#else
  blockCopy_type block = exchange_ghosts->blocks[block_type][blockIdx];
#endif

  // copy 3D array from read_i,j,k of read[] to write_i,j,k in write[]
  int   dim_i       = block.dim.i;
  int   dim_j       = block.dim.j;
  int   dim_k       = block.dim.k;

  int  read_i       = block.read.i;
  int  read_j       = block.read.j;
  int  read_k       = block.read.k;
  int  read_jStride = block.read.jStride;
  int  read_kStride = block.read.kStride;

  int write_i       = block.write.i;
  int write_j       = block.write.j;
  int write_k       = block.write.k;
  int write_jStride = block.write.jStride;
  int write_kStride = block.write.kStride;

  double * __restrict__  read = block.read.ptr;
  double * __restrict__ write = block.write.ptr;

  int  read_box = block.read.box;
  int write_box = block.write.box;
  if(read_box >=0)
    read = level->my_boxes[ read_box].vectors[id] + level->my_boxes[ read_box].ghosts*(1+level->my_boxes[ read_box].jStride+level->my_boxes[ read_box].kStride);
  if(write_box>=0)
    write = level->my_boxes[write_box].vectors[id] + level->my_boxes[write_box].ghosts*(1+level->my_boxes[write_box].jStride+level->my_boxes[write_box].kStride);

#if defined (SPEC_OPENMP_TARGET)
# if defined (SPMD_MODE)
#  pragma omp for
# else
#  pragma omp parallel for
# endif
#elif defined (SPEC_OPENACC)
# pragma acc loop vector
#endif
  for(int gid=0; gid<dim_i*dim_j*dim_k; ++gid){
    // simple linear mapping of 1D threads to 3D indices
    int k=(gid/dim_i)/dim_j;
    int j=(gid/dim_i)%dim_j;
    int i=gid%dim_i;

    int  read_ijk = (i+ read_i) + (j+ read_j)* read_jStride + (k+ read_k)* read_kStride;
    int write_ijk = (i+write_i) + (j+write_j)*write_jStride + (k+write_k)*write_kStride;
    write[write_ijk] = READ(read_ijk);
  } // end of gid loop
  } // end of blockIdx loop

#if defined (SPEC_OPENMP_TARGET) && defined(SPMD_MODE)
  } // end of omp parallel
  } // end of omp target teams
#endif

#if defined (SPEC_OPENMP_TARGET) && defined(SPEC_MANAGED_MEMORY)
# pragma omp target exit data map(release:exchange_ghosts[:1])
# pragma omp target exit data map(delete:level[:1])
#endif
}

void directives_increment_block(level_type *level, int id, double prescale, communicator_type *exchange_ghosts, int block_type)
{
#if defined (SPEC_OPENMP_TARGET) && defined(SPEC_MANAGED_MEMORY)
/* Two part map because "exchange_ghosts" aliases part of "level":
   operators/interpolation_v4.c: device_increment_block(level_f,id_f,prescale_f,&level_f->interpolation,2);
   operators/interpolation_v2.c: device_increment_block(level_f,id_f,prescale_f,&level_f->interpolation,2);

   We could alternatively use target map(to:level[:1], exchange_ghosts[:0]), however,
   this would rule out a future use of this function without aliased data. */
# pragma omp target enter data map(to:level[:1])
# pragma omp target enter data map(to:exchange_ghosts[:1])
#endif

#if defined (SPEC_OPENMP_TARGET)
# pragma omp target teams distribute map(to:level[:0], exchange_ghosts[:0])
#elif defined (SPEC_OPENACC)
# if defined(SPEC_MANAGED_MEMORY)
#  pragma acc parallel loop gang copyin(exchange_ghosts[:1], level[:1])
# else
#  pragma acc parallel loop gang present(exchange_ghosts, level) default(present)
# endif
#endif
  for (int blockIdx=0; blockIdx < exchange_ghosts->num_blocks[block_type]; ++blockIdx) {

  // one CUDA thread block operates on one HPGMG tile/block
#if defined(FIX_PGI_STRUCT_ASSIGNMENT_BUG)
  // Explicitly assign data for each struct field.
  // This workaround is needed to successfully compile with PGI compiler.
  // The original default code fails with the following error:
  // PGC-S-0155-Compiler failed to translate accelerator region (see -Minfo messages):
  // Unknown variable reference (directives/operators-directives.fv4.c: 45)
  //
  // The bug has been reported:
  // FS#28014: NERSC/HPC2020 ACC fails with Internal accelerator LILI consistency check failure.
  blockCopy_type block;
  block.dim.i = exchange_ghosts->blocks[block_type][blockIdx].dim.i;
  block.dim.j = exchange_ghosts->blocks[block_type][blockIdx].dim.j;
  block.dim.k = exchange_ghosts->blocks[block_type][blockIdx].dim.k;
  block.read.i = exchange_ghosts->blocks[block_type][blockIdx].read.i;
  block.read.j = exchange_ghosts->blocks[block_type][blockIdx].read.j;
  block.read.k = exchange_ghosts->blocks[block_type][blockIdx].read.k;
  block.read.jStride = exchange_ghosts->blocks[block_type][blockIdx].read.jStride;
  block.read.kStride = exchange_ghosts->blocks[block_type][blockIdx].read.kStride;
  block.write.i = exchange_ghosts->blocks[block_type][blockIdx].write.i;
  block.write.j = exchange_ghosts->blocks[block_type][blockIdx].write.j;
  block.write.k = exchange_ghosts->blocks[block_type][blockIdx].write.k;
  block.write.jStride = exchange_ghosts->blocks[block_type][blockIdx].write.jStride;
  block.write.kStride = exchange_ghosts->blocks[block_type][blockIdx].write.kStride;
  block.read.ptr = exchange_ghosts->blocks[block_type][blockIdx].read.ptr;
  block.write.ptr = exchange_ghosts->blocks[block_type][blockIdx].write.ptr;
  block.read.box = exchange_ghosts->blocks[block_type][blockIdx].read.box;
  block.write.box = exchange_ghosts->blocks[block_type][blockIdx].write.box;
#else
  blockCopy_type block = exchange_ghosts->blocks[block_type][blockIdx];
#endif

  // copy 3D array from read_i,j,k of read[] to write_i,j,k in write[]
  int   dim_i       = block.dim.i;
  int   dim_j       = block.dim.j;
  int   dim_k       = block.dim.k;

  int  read_i       = block.read.i;
  int  read_j       = block.read.j;
  int  read_k       = block.read.k;
  int  read_jStride = block.read.jStride;
  int  read_kStride = block.read.kStride;

  int write_i       = block.write.i;
  int write_j       = block.write.j;
  int write_k       = block.write.k;
  int write_jStride = block.write.jStride;
  int write_kStride = block.write.kStride;

  double * __restrict__  read = block.read.ptr;
  double * __restrict__ write = block.write.ptr;

  if(block.read.box >=0){
    read = level->my_boxes[ block.read.box].vectors[id] + level->my_boxes[ block.read.box].ghosts*(1+level->my_boxes[ block.read.box].jStride+level->my_boxes[ block.read.box].kStride);
    read_jStride = level->my_boxes[block.read.box ].jStride;
    read_kStride = level->my_boxes[block.read.box ].kStride;
  }
  if(block.write.box>=0){
    write = level->my_boxes[block.write.box].vectors[id] + level->my_boxes[block.write.box].ghosts*(1+level->my_boxes[block.write.box].jStride+level->my_boxes[block.write.box].kStride);
    write_jStride = level->my_boxes[block.write.box].jStride;
    write_kStride = level->my_boxes[block.write.box].kStride;
  }

#if defined (SPEC_OPENMP_TARGET)
# pragma omp parallel for
#elif defined (SPEC_OPENACC)
# pragma acc loop vector
#endif
  for(int gid=0; gid<dim_i*dim_j*dim_k; ++gid){
    // simple linear mapping of 1D threads to 3D indices
    int k=(gid/dim_i)/dim_j;
    int j=(gid/dim_i)%dim_j;
    int i=gid%dim_i;

    int  read_ijk = (i+ read_i) + (j+ read_j)* read_jStride + (k+ read_k)* read_kStride;
    int write_ijk = (i+write_i) + (j+write_j)*write_jStride + (k+write_k)*write_kStride;
    write[write_ijk] = prescale*write[write_ijk] + READ(read_ijk);
  } // end of parallel for
  } // end of teams distribute

#if defined (SPEC_OPENMP_TARGET) && defined(SPEC_MANAGED_MEMORY)
# pragma omp target exit data map(release:exchange_ghosts[:1])
# pragma omp target exit data map(delete:level[:1])
#endif
}
