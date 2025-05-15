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

#undef  READ
#define READ(i) read[i]
void directives_restriction(level_type *level_c, int id_c, level_type *level_f, int id_f, communicator_type *restriction, int restrictionType, int block_type)
{
  int num_blocks = restriction->num_blocks[block_type]; if(num_blocks<=0) return;
  int blockX = min(level_c->box_dim,BLOCKCOPY_TILE_I); // threadDim.x
  int blockY = BLOCKCOPY_TILE_J; // threadDim.y
  // int blockZ = 1;
  int gridX = (BLOCKCOPY_TILE_I+blockX-1)/blockX; // blockDim.x
  int gridY = (BLOCKCOPY_TILE_J+blockY-1)/blockY; // blockDim.y
  // int gridZ = num_blocks;

#if defined (SPEC_OPENMP_TARGET) && defined(SPEC_MANAGED_MEMORY)
  /* Two part map because "restriction" aliases part of "level_f":
     operators/restriction.c: device_restriction(level_c,id_c,level_f,id_f,&level_f->restriction[restrictionType],restrictionType,0);
     operators/restriction.c: device_restriction(level_c, id_c, level_f, id_f, &level_f->restriction[restrictionType], restrictionType, 1);

     We could alternatively use target map(to:level_f[:1], level_c[:1], restriction[:0]), however,
     this would rule out a future use of this function without aliased data. */
# pragma omp target enter data map(to:level_f[:1])
# pragma omp target enter data map(to:level_c[:1])
# pragma omp target enter data map(to:restriction[:1])
#endif

#if defined(SPEC_OPENMP_TARGET)
# pragma omp target teams distribute collapse(3) map(to: level_f[:0], level_c[:0], restriction[:0])
#elif defined(SPEC_OPENACC)
# if defined(SPEC_MANAGED_MEMORY)
#  pragma acc parallel loop gang collapse(3) copyin(level_c[:1], level_f[:1], restriction[:1])
# else
#  pragma acc parallel loop gang collapse(3) present(level_c, level_f, restriction) default(present)
# endif
#endif
  for (int blockIdx=0; blockIdx < num_blocks; ++blockIdx) {
  for (int yTileIdx=0; yTileIdx < gridY; ++yTileIdx) {
  for (int xTileIdx=0; xTileIdx < gridX; ++xTileIdx) {

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
  block.dim.i = restriction->blocks[block_type][blockIdx].dim.i;
  block.dim.j = restriction->blocks[block_type][blockIdx].dim.j;
  block.dim.k = restriction->blocks[block_type][blockIdx].dim.k;
  block.read.i = restriction->blocks[block_type][blockIdx].read.i;
  block.read.j = restriction->blocks[block_type][blockIdx].read.j;
  block.read.k = restriction->blocks[block_type][blockIdx].read.k;
  block.read.jStride = restriction->blocks[block_type][blockIdx].read.jStride;
  block.read.kStride = restriction->blocks[block_type][blockIdx].read.kStride;
  block.write.i = restriction->blocks[block_type][blockIdx].write.i;
  block.write.j = restriction->blocks[block_type][blockIdx].write.j;
  block.write.k = restriction->blocks[block_type][blockIdx].write.k;
  block.write.jStride = restriction->blocks[block_type][blockIdx].write.jStride;
  block.write.kStride = restriction->blocks[block_type][blockIdx].write.kStride;
  block.read.ptr = restriction->blocks[block_type][blockIdx].read.ptr;
  block.write.ptr = restriction->blocks[block_type][blockIdx].write.ptr;
  block.read.box = restriction->blocks[block_type][blockIdx].read.box;
  block.write.box = restriction->blocks[block_type][blockIdx].write.box;
#else
  blockCopy_type block = restriction->blocks[block_type][blockIdx];
#endif

  // restrict 3D array from read_i,j,k of read[] to write_i,j,k in write[]
  int   dim_i       = block.dim.i; // calculate the dimensions of the resultant coarse block
  int   dim_j       = block.dim.j;
  int   dim_k       = block.dim.k;


  // CD: Save the actual number of tiles in total_tile_i. It is possible
  // that the tileIdx loops create more tiles than is needed to cover the
  // iteration space. Therefore, some tiles must be given zero work.
  int total_tile_i = (dim_i + blockX - 1) / blockX; // integer round up
  int total_tile_j = (dim_j + blockY - 1) / blockY; // integer round up
  int tile_start_i, tile_start_j;
  int tile_end_i, tile_end_j;

  if (xTileIdx < total_tile_i) {
    tile_start_i = xTileIdx * blockX;
    tile_end_i = (xTileIdx+1) * blockX;
    // CD: The if condition is needed when write_dim_i is not a
    // multiple of BLOCKCOPY_TILE_I.
    if (tile_end_i > dim_i) tile_end_i = dim_i;
  } else {
    // CD: Do nothing. Equivalent to a return from the loop which is
    // not allowed in OpenACC
    tile_start_i = 0;
    tile_end_i = 0;
  }

  if (yTileIdx < total_tile_j) {
    tile_start_j = yTileIdx * blockY;
    tile_end_j = (yTileIdx+1) * blockY;
    // CD: The if condition is needed when write_dim_j is not a
    // multiple of BLOCKCOPY_TILE_J.
    if (tile_end_j > dim_j) tile_end_j = dim_j;
  } else {
    // CD: Do nothing. Equivalent to a return from the loop which is
    // not allowed in OpenACC
    tile_start_j = 0;
    tile_end_j = 0;
  }

#if defined(SPEC_OPENMP_TARGET)
# pragma omp parallel for collapse(2)
#elif defined(SPEC_OPENACC)
# pragma acc loop worker vector collapse(2)
#endif
  for (int i=tile_start_i; i<tile_end_i; ++i) {
  for (int j=tile_start_j; j<tile_end_j; ++j) {


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
     read = level_f->my_boxes[ block.read.box].vectors[id_f] + level_f->my_boxes[ block.read.box].ghosts*(1+level_f->my_boxes[ block.read.box].jStride+level_f->my_boxes[ block.read.box].kStride);
     read_jStride = level_f->my_boxes[block.read.box ].jStride;
     read_kStride = level_f->my_boxes[block.read.box ].kStride;
  }
  if(block.write.box>=0){
    write = level_c->my_boxes[block.write.box].vectors[id_c] + level_c->my_boxes[block.write.box].ghosts*(1+level_c->my_boxes[block.write.box].jStride+level_c->my_boxes[block.write.box].kStride);
    write_jStride = level_c->my_boxes[block.write.box].jStride;
    write_kStride = level_c->my_boxes[block.write.box].kStride;
  }

  switch(restrictionType){
    case RESTRICT_CELL:
#if defined(SPEC_OPENACC)
# pragma acc loop seq
#endif
         for(int k=0;k<dim_k;k++){
           int write_ijk = ((i   )+write_i) + ((j   )+write_j)*write_jStride + ((k   )+write_k)*write_kStride;
           int  read_ijk = ((i<<1)+ read_i) + ((j<<1)+ read_j)* read_jStride + ((k<<1)+ read_k)* read_kStride;
	   double r11 = READ(read_ijk                            );
	   double r12 = READ(read_ijk+1                          );
	   double r21 = READ(read_ijk  +read_jStride             );
	   double r22 = READ(read_ijk+1+read_jStride             );
	   double r31 = READ(read_ijk               +read_kStride);
	   double r32 = READ(read_ijk+1             +read_kStride);
	   double r41 = READ(read_ijk  +read_jStride+read_kStride);
	   double r42 = READ(read_ijk+1+read_jStride+read_kStride);
	   write[write_ijk] = ( r11+r12 + r21+r22 + r31+r32 + r41+r42 ) * 0.125;
           //write[write_ijk] = ( READ(read_ijk                            )+READ(read_ijk+1                          ) +
           //                     READ(read_ijk  +read_jStride             )+READ(read_ijk+1+read_jStride             ) +
           //                     READ(read_ijk               +read_kStride)+READ(read_ijk+1             +read_kStride) +
           //                     READ(read_ijk  +read_jStride+read_kStride)+READ(read_ijk+1+read_jStride+read_kStride) ) * 0.125;
         }break;
    case RESTRICT_FACE_I:
#if defined(SPEC_OPENACC)
# pragma acc loop seq
#endif
	 for(int k=0;k<dim_k;k++){
           int write_ijk = ((i   )+write_i) + ((j   )+write_j)*write_jStride + ((k   )+write_k)*write_kStride;
           int  read_ijk = ((i<<1)+ read_i) + ((j<<1)+ read_j)* read_jStride + ((k<<1)+ read_k)* read_kStride;
	   double r1 = READ(read_ijk                          );
	   double r2 = READ(read_ijk+read_jStride             );
	   double r3 = READ(read_ijk             +read_kStride);
	   double r4 = READ(read_ijk+read_jStride+read_kStride);
           write[write_ijk] = ( r1 + r2 + r3 + r4 ) * 0.25;
           //write[write_ijk] = ( READ(read_ijk                          ) +
           //                     READ(read_ijk+read_jStride             ) +
           //                     READ(read_ijk             +read_kStride) +
           //                     READ(read_ijk+read_jStride+read_kStride) ) * 0.25;
         }break;
    case RESTRICT_FACE_J:
#if defined(SPEC_OPENACC)
# pragma acc loop seq
#endif
	 for(int k=0;k<dim_k;k++){
           int write_ijk = ((i   )+write_i) + ((j   )+write_j)*write_jStride + ((k   )+write_k)*write_kStride;
           int  read_ijk = ((i<<1)+ read_i) + ((j<<1)+ read_j)* read_jStride + ((k<<1)+ read_k)* read_kStride;
  	   double r1 = READ(read_ijk               );
	   double r2 = READ(read_ijk+1             );
	   double r3 = READ(read_ijk  +read_kStride);
	   double r4 = READ(read_ijk+1+read_kStride);
	   write[write_ijk] = ( r1 + r2 + r3 + r4 ) * 0.25;
           //write[write_ijk] = ( READ(read_ijk               ) +
           //                     READ(read_ijk+1             ) +
           //                     READ(read_ijk  +read_kStride) +
           //                     READ(read_ijk+1+read_kStride) ) * 0.25;
         }break;
    case RESTRICT_FACE_K:
#if defined(SPEC_OPENACC)
# pragma acc loop seq
#endif
	 for(int k=0;k<dim_k;k++){
           int write_ijk = ((i   )+write_i) + ((j   )+write_j)*write_jStride + ((k   )+write_k)*write_kStride;
           int  read_ijk = ((i<<1)+ read_i) + ((j<<1)+ read_j)* read_jStride + ((k<<1)+ read_k)* read_kStride;
	   double r1 = READ(read_ijk               );
	   double r2 = READ(read_ijk+1             );
	   double r3 = READ(read_ijk  +read_jStride);
	   double r4 = READ(read_ijk+1+read_jStride);
	   write[write_ijk] = ( r1 + r2 + r3 + r4 ) * 0.25;
           //write[write_ijk] = ( READ(read_ijk               ) +
           //                     READ(read_ijk+1             ) +
           //                     READ(read_ijk  +read_jStride) +
           //                     READ(read_ijk+1+read_jStride) ) * 0.25;
         }break;
  }
  } // End of j loop
  } // End of i loop
  } // End of x tile loop
  } // End of y tile loop
  } // End of block loop

#if defined (SPEC_OPENMP_TARGET) && defined(SPEC_MANAGED_MEMORY)
# pragma omp target exit data map(release:restriction[:1])
# pragma omp target exit data map(delete:level_c[:1])
# pragma omp target exit data map(delete:level_f[:1])
#endif
}
