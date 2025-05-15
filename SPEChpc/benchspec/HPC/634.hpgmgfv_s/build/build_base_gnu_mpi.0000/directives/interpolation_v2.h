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
#define READ(i)	read[i]


void directives_interpolation_v2(level_type *level_f, int id_f, double prescale_f, level_type *level_c, int id_c, communicator_type *interpolation, int block_type) {
  // one CUDA thread block operates on one HPGMG tile/block
  int num_blocks = interpolation->num_blocks[block_type]; if(num_blocks<=0) return;
  int blockX = min(level_f->box_dim,BLOCKCOPY_TILE_I); // threadDim.x
  // int blockY = BLOCKCOPY_TILE_J; // threadDim.y
  // int blockZ = 1;
  int gridX = (level_f->box_dim+blockX-1)/blockX; // blockDim.x
  // int gridY = 1; // blockDim.y
  // int gridZ = num_blocks;

#if defined (SPEC_OPENMP_TARGET) && defined(SPEC_MANAGED_MEMORY)
  /* Two part map because "interpolation" aliases part of "level_c":
     operators/interpolation_v2.c: device_interpolation_v2(level_f,id_f,0.0,level_c,id_c,&level_c->interpolation,0);
     operators/interpolation_v2.c: device_interpolation_v2(level_f,id_f,prescale_f,level_c,id_c,&level_c->interpolation,1);

     We could alternatively use target map(to:level_f[:1], level_c[:1], interpolation[:0]), however,
     this would rule out a future use of this function without aliased data. */
# pragma omp target enter data map(to:level_f[:1])
# pragma omp target enter data map(to:level_c[:1])
# pragma omp target enter data map(to:interpolation[:1])
#endif

#if defined (SPEC_OPENMP_TARGET)
# pragma omp target teams distribute collapse(2) map(to: level_f[:0], level_c[:0], interpolation[:0])
#elif defined (SPEC_OPENACC)
# if defined (SPEC_MANAGED_MEMORY)
#  pragma acc parallel loop gang collapse(2) copyin(level_f[:1], level_c[:1], interpolation[:1])
# else
#  pragma acc parallel loop gang collapse(2) present(level_f, level_c, interpolation) default(present)
# endif
#endif
  for (int blockIdx=0; blockIdx < num_blocks; ++blockIdx) {
  for (int tileIdx=0; tileIdx < gridX; ++tileIdx) {

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
  block.dim.i = interpolation->blocks[block_type][blockIdx].dim.i;
  block.dim.j = interpolation->blocks[block_type][blockIdx].dim.j;
  block.dim.k = interpolation->blocks[block_type][blockIdx].dim.k;
  block.read.i = interpolation->blocks[block_type][blockIdx].read.i;
  block.read.j = interpolation->blocks[block_type][blockIdx].read.j;
  block.read.k = interpolation->blocks[block_type][blockIdx].read.k;
  block.read.jStride = interpolation->blocks[block_type][blockIdx].read.jStride;
  block.read.kStride = interpolation->blocks[block_type][blockIdx].read.kStride;
  block.write.i = interpolation->blocks[block_type][blockIdx].write.i;
  block.write.j = interpolation->blocks[block_type][blockIdx].write.j;
  block.write.k = interpolation->blocks[block_type][blockIdx].write.k;
  block.write.jStride = interpolation->blocks[block_type][blockIdx].write.jStride;
  block.write.kStride = interpolation->blocks[block_type][blockIdx].write.kStride;
  block.read.ptr = interpolation->blocks[block_type][blockIdx].read.ptr;
  block.write.ptr = interpolation->blocks[block_type][blockIdx].write.ptr;
  block.read.box = interpolation->blocks[block_type][blockIdx].read.box;
  block.write.box = interpolation->blocks[block_type][blockIdx].write.box;
#else
  blockCopy_type block = interpolation->blocks[block_type][blockIdx];
#endif

  // interpolate 3D array from read_i,j,k of read[] to write_i,j,k in write[]
  int write_dim_i   = block.dim.i<<1; // calculate the dimensions of the resultant fine block
  int write_dim_j   = block.dim.j<<1;
  int write_dim_k   = block.dim.k<<1;

  // CD: Save the actual number of tiles in total_tile_i. It is possible
  // that the tileIdx loop creates more tiles than is needed to cover the
  // iteration space. Therefore, some tiles must be given zero work.
  int total_tile_i = (write_dim_i + blockX - 1) / blockX; // integer round up
  int tile_start_i;
  int tile_end_i;

  if (tileIdx < total_tile_i) {
    tile_start_i = tileIdx * blockX;
    tile_end_i = (tileIdx+1) * blockX;
    // CD: The if condition is needed when write_dim_i is not a
    // multiple of BLOCKCOPY_TILE_I.
    if (tile_end_i > write_dim_i) tile_end_i = write_dim_i;
  } else {
    // CD: Do nothing. Equivalent to a return from the loop which is
    // not allowed in OpenACC
    tile_start_i = 0;
    tile_end_i = 0;
  }

#if defined (SPEC_OPENMP_TARGET)
# pragma omp parallel for collapse(2)
#elif defined (SPEC_OPENACC)
# pragma acc loop vector collapse(2)
#endif
  for (int i=tile_start_i; i<tile_end_i; ++i) {
  for (int j=0; j<write_dim_j; j+=2) {

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
     read_jStride = level_c->my_boxes[block.read.box ].jStride;
     read_kStride = level_c->my_boxes[block.read.box ].kStride;
     read = level_c->my_boxes[ block.read.box].vectors[id_c] + level_c->my_boxes[ block.read.box].ghosts*(1+ read_jStride+ read_kStride);
  }
  if(block.write.box>=0){
    write_jStride = level_f->my_boxes[block.write.box].jStride;
    write_kStride = level_f->my_boxes[block.write.box].kStride;
    write = level_f->my_boxes[block.write.box].vectors[id_f] + level_f->my_boxes[block.write.box].ghosts*(1+write_jStride+write_kStride);
  }

  double c1 = 1.0/8.0;
#if defined (SPEC_OPENACC)
# pragma acc loop seq
#endif
  for(int k=0;k<write_dim_k;k+=2){
    double c1i=c1;if(i&0x1){c1i=-c1;}
    double c1j=c1;//if(j&0x1){c1j=-c1;}
    double c1k=c1;//if(k&0x1){c1k=-c1;}
    int write_ijk = ((i   )+write_i) + (((j   )+write_j)*write_jStride) + (((k   )+write_k)*write_kStride);
    int  read_ijk = ((i>>1)+ read_i) + (((j>>1)+ read_j)* read_jStride) + (((k>>1)+ read_k)* read_kStride);
    //
    // |  1/8  |  1.0  | -1/8  | coarse grid
    // |---+---|---+---|---+---|
    // |   |   |???|   |   |   | fine grid
    //
    double r11 = READ(read_ijk-1-read_jStride-read_kStride);
    double r12 = READ(read_ijk-read_jStride-read_kStride  );
    double r13 = READ(read_ijk+1-read_jStride-read_kStride);
    double r21 = READ(read_ijk-1             -read_kStride);
    double r22 = READ(read_ijk             -read_kStride  );
    double r23 = READ(read_ijk+1             -read_kStride);
    double r31 = READ(read_ijk-1+read_jStride-read_kStride);
    double r32 = READ(read_ijk+read_jStride-read_kStride  );
    double r33 = READ(read_ijk+1+read_jStride-read_kStride);
    double r41 = READ(read_ijk-1-read_jStride             );
    double r42 = READ(read_ijk-read_jStride               );
    double r43 = READ(read_ijk+1-read_jStride             );
    double r51 = READ(read_ijk-1                          );
    double r52 = READ(read_ijk                            );
    double r53 = READ(read_ijk+1                          );
    double r61 = READ(read_ijk-1+read_jStride             );
    double r62 = READ(read_ijk+read_jStride               );
    double r63 = READ(read_ijk+1+read_jStride             );
    double r71 = READ(read_ijk-1-read_jStride+read_kStride);
    double r72 = READ(read_ijk-read_jStride+read_kStride  );
    double r73 = READ(read_ijk+1-read_jStride+read_kStride);
    double r81 = READ(read_ijk-1             +read_kStride);
    double r82 = READ(read_ijk             +read_kStride  );
    double r83 = READ(read_ijk+1             +read_kStride);
    double r91 = READ(read_ijk-1+read_jStride+read_kStride);
    double r92 = READ(read_ijk+read_jStride+read_kStride  );
    double r93 = READ(read_ijk+1+read_jStride+read_kStride);
 
    // i  j  k
    write[write_ijk] = prescale_f*write[write_ijk] +
                       + c1k*( + c1j*( c1i*r11 + r12 - c1i*r13 )
                               +     ( c1i*r21 + r22 - c1i*r23 )
                               - c1j*( c1i*r31 + r32 - c1i*r33 ) )
                       +     ( + c1j*( c1i*r41 + r42 - c1i*r43 )
                               +     ( c1i*r51 + r52 - c1i*r53 )
                               - c1j*( c1i*r61 + r62 - c1i*r63 ) )
                       - c1k*( + c1j*( c1i*r71 + r72 - c1i*r73 )
                               +     ( c1i*r81 + r82 - c1i*r83 )
                               - c1j*( c1i*r91 + r92 - c1i*r93 ) );

   // i  j+1  k
   write_ijk = ((i  )+write_i) + (((j+1)+write_j)*write_jStride) + (((k  )+write_k)*write_kStride);  c1j=-c1;c1k=c1;
   write[write_ijk] = prescale_f*write[write_ijk] +
                       + c1k*( + c1j*( c1i*r11 + r12 - c1i*r13 )
                               +     ( c1i*r21 + r22 - c1i*r23 )
                               - c1j*( c1i*r31 + r32 - c1i*r33 ) )
                       +     ( + c1j*( c1i*r41 + r42 - c1i*r43 )
                               +     ( c1i*r51 + r52 - c1i*r53 )
                               - c1j*( c1i*r61 + r62 - c1i*r63 ) )
                       - c1k*( + c1j*( c1i*r71 + r72 - c1i*r73 )
                               +     ( c1i*r81 + r82 - c1i*r83 )
                               - c1j*( c1i*r91 + r92 - c1i*r93 ) );

   // i  j  k+1
   write_ijk = ((i  )+write_i) + (((j  )+write_j)*write_jStride) + (((k+1)+write_k)*write_kStride);  c1j=c1;c1k=-c1;
   write[write_ijk] = prescale_f*write[write_ijk] +
                       + c1k*( + c1j*( c1i*r11 + r12 - c1i*r13 )
                               +     ( c1i*r21 + r22 - c1i*r23 )
                               - c1j*( c1i*r31 + r32 - c1i*r33 ) )
                       +     ( + c1j*( c1i*r41 + r42 - c1i*r43 )
                               +     ( c1i*r51 + r52 - c1i*r53 )
                               - c1j*( c1i*r61 + r62 - c1i*r63 ) )
                       - c1k*( + c1j*( c1i*r71 + r72 - c1i*r73 )
                               +     ( c1i*r81 + r82 - c1i*r83 )
                               - c1j*( c1i*r91 + r92 - c1i*r93 ) );

    // i  j+1  k+1
    write_ijk = ((i  )+write_i) + (((j+1)+write_j)*write_jStride) + (((k+1)+write_k)*write_kStride);  c1j=-c1;c1k=-c1;
    write[write_ijk] = prescale_f*write[write_ijk] +
                       + c1k*( + c1j*( c1i*r11 + r12 - c1i*r13 )
                               +     ( c1i*r21 + r22 - c1i*r23 )
                               - c1j*( c1i*r31 + r32 - c1i*r33 ) )
                       +     ( + c1j*( c1i*r41 + r42 - c1i*r43 )
                               +     ( c1i*r51 + r52 - c1i*r53 )
                               - c1j*( c1i*r61 + r62 - c1i*r63 ) )
                       - c1k*( + c1j*( c1i*r71 + r72 - c1i*r73 )
                               +     ( c1i*r81 + r82 - c1i*r83 )
                               - c1j*( c1i*r91 + r92 - c1i*r93 ) );
  } // End of sequential k loop
  } // End of j loop
  } // End of i loop
  } // End of tile loop
  } // End of block loop

#if defined (SPEC_OPENMP_TARGET) && defined(SPEC_MANAGED_MEMORY)
# pragma omp target exit data map(release:interpolation[:1])
# pragma omp target exit data map(delete:level_c[:1])
# pragma omp target exit data map(delete:level_f[:1])
#endif
}
