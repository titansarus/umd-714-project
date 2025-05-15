/*
  This source file is included by operators-directives.fv4.c
*/

#if defined(SPEC_OPENMP_TARGET) && !defined(SPEC_MANAGED_MEMORY)
// This file is not a header file - despite the .h extension.
// It is therefore OK to define a constant.
#ifdef PASS_NULL_TO_ATTACH
const int attach_null = 1;
#else
const int attach_null = 0;
#endif

#ifdef MINIMIZE_PRESENT_TABLE
// This avoids adding an entry to the present table
# define ATTACH_DEVICE_ADDRESS(_ptr) omp_attach((void**)&_ptr);
#else
// This is the default way to attach a device address
# define MyPragma(a) _Pragma(#a)
# define ATTACH_DEVICE_ADDRESS(_ptr)				\
  {								\
    if (_ptr || attach_null) {					\
      MyPragma(omp target enter data map (alloc:_ptr[:0]))	\
    }								\
  }
#endif

// omp_attach() takes the address of a pointer variable within
// the level data structure. At entry, the device version of
// this pointer points to a host address. We use the use_device_ptr
// clause to identify the device address corresponding
// to this host address. We then update the device version of this
// pointer to point to this newly-identifed device address.
void omp_attach(void **ptr)
{
  void *dptr = *ptr;
  if (dptr || attach_null) {
#pragma omp target data use_device_ptr(dptr)
    {
#pragma omp target is_device_ptr(dptr)
      {
	*ptr = dptr;
      }
    }
  }
}
#endif

void create_level_device_openmp(level_type *level) {
#if defined(SPEC_OPENMP_TARGET) && !defined(SPEC_MANAGED_MEMORY)
int numVectors = level->numVectors;
int num_my_boxes = level->num_my_boxes;
int box_volume = level->box_volume;
int shape;
int box, block;
int i,j,k;

#pragma omp target enter data map(to:level[:1])
#pragma omp target enter data map(to:level->rank_of_box[:level->boxes_in.i*level->boxes_in.j*level->boxes_in.k])


// Create the vectors
#ifdef VECTOR_MALLOC_BULK
#pragma omp target enter data map(to:level->vectors_base_1[:numVectors*num_my_boxes*box_volume])
// alloc is appropriate for level->vectors when using VECTOR_MALLOC_BULK because
// only level->vectors_base_1 is used to update data between CPU and GPU
#pragma omp target enter data map(alloc:level->vectors[:numVectors])
#pragma omp target
   for (i=0; i < numVectors; ++i) {
     level->vectors[i] = level->vectors_base_1 + (i*level->num_my_boxes*level->box_volume);
   }
#else
#pragma omp target enter data map(to:level->vectors[:numVectors])
   for (i=0; i<numVectors; i++) {
#pragma omp target enter data map(to:level->vectors[i][:num_my_boxes*box_volume])
   }
#endif


// Create the Boxes
#pragma omp target enter data map(to:level->my_boxes[:level->num_my_boxes])
  box=0;
  for(k=0;k<level->boxes_in.k;k++){
  for(j=0;j<level->boxes_in.j;j++){
  for(i=0;i<level->boxes_in.i;i++){
    int jStride = level->boxes_in.i;
    int kStride = level->boxes_in.i*level->boxes_in.j;
    int b=i + j*jStride + k*kStride;
    if(level->rank_of_box[b]==level->my_rank){

#ifdef MINIMIZE_PRESENT_TABLE
      // Attach device addresses without adding present table entries
      double **dptr = (double **)omp_target_alloc(numVectors * sizeof(double*), omp_get_default_device());
      size_t offset = (size_t)box*level->box_volume;
#pragma omp target is_device_ptr(dptr)
      {
	level->my_boxes[box].vectors = dptr;
	for (int c=0; c < numVectors; ++c) {
	  level->my_boxes[box].vectors[c] = level->vectors[c] + offset;
	}
      }
#else
      // Attach device addresses - creates unnecessary present table entries
#pragma omp target enter data map(alloc:level->my_boxes[box].vectors[:numVectors])
      for (int c=0; c < numVectors; ++c) {
#pragma omp target enter data map(alloc:level->my_boxes[box].vectors[c][:0])
      }
#endif

      ++box;
    } // End: if(level->rank_of_box[b]==level->my_rank)
  }}}


// Create the Blocks;
#pragma omp target enter data map(to:level->my_blocks[:level->num_my_blocks])
  for (block=0; block < level->num_my_blocks; ++block) {
    ATTACH_DEVICE_ADDRESS(level->my_blocks[block].read.ptr)
    ATTACH_DEVICE_ADDRESS(level->my_blocks[block].write.ptr)
  }

// Create the boundary condition
   if (level->boundary_condition.type != BC_PERIODIC) {
#pragma omp target enter data map(to:level->boundary_condition.blocks[0:STENCIL_MAX_SHAPES])
     for (shape=0;shape<STENCIL_MAX_SHAPES;shape++) {
       #pragma omp target enter data map(to:level->boundary_condition.blocks[shape][:level->boundary_condition.num_blocks[shape]])
     }
   }

// Create RedBlack_FP
  if(level->num_my_boxes){
    #pragma omp target enter data map(to:level->RedBlack_FP[:2*level->my_boxes[0].kStride])
  }
#endif

}


void build_restriction_device_openmp(level_type *level) {
#if defined(SPEC_OPENMP_TARGET) && !defined(SPEC_MANAGED_MEMORY)
int shape,numSendRanks, numRecvRanks, neighbor;
int block;

// Create the restriction
#pragma omp target update to(level->restriction[:4])
   for (shape =0;shape<4;shape++) {
        numSendRanks = level->restriction[shape].num_sends;
        numRecvRanks = level->restriction[shape].num_recvs;

        if (numSendRanks > 0) {
          #pragma omp target enter data map(to:level->restriction[shape].send_ranks[:numSendRanks])
          #pragma omp target enter data map(to:level->restriction[shape].send_sizes[:numSendRanks])
          #pragma omp target enter data map(alloc:level->restriction[shape].send_buffers[:numSendRanks])
	  for(neighbor=0;neighbor<numSendRanks;neighbor++){
            #pragma omp target enter data map(alloc:level->restriction[shape].send_buffers[neighbor][:level->restriction[shape].send_sizes[neighbor]])
	  }
        }
	if (numRecvRanks>0) {
           #pragma omp target enter data map(to:level->restriction[shape].recv_ranks[:numRecvRanks])
           #pragma omp target enter data map(to:level->restriction[shape].recv_sizes[:numRecvRanks])
           #pragma omp target enter data map(alloc:level->restriction[shape].recv_buffers[:numRecvRanks])
	   for(neighbor=0;neighbor<numRecvRanks;neighbor++){
           #pragma omp target enter data map(alloc:level->restriction[shape].recv_buffers[neighbor][:level->restriction[shape].recv_sizes[neighbor]])
	   }
	}
        for (block=0;block<3;++block) {
	  if (level->restriction[shape].blocks[block] != NULL) {
            #pragma omp target enter data map(to:level->restriction[shape].blocks[block][:level->restriction[shape].num_blocks[block]])
	    for (int b=0; b < level->restriction[shape].num_blocks[block]; ++b) {
	      ATTACH_DEVICE_ADDRESS(level->restriction[shape].blocks[block][b].read.ptr)
	      ATTACH_DEVICE_ADDRESS(level->restriction[shape].blocks[block][b].write.ptr)
	    }
	  }
        } // End of loop over blocks
   }
#endif
}

void build_interpolation_device_openmp(level_type *level) {
#if defined(SPEC_OPENMP_TARGET) && !defined(SPEC_MANAGED_MEMORY)
int numSendRanks, numRecvRanks, neighbor;
int block;

// Create interpolation
   numSendRanks = level->interpolation.num_sends;
   numRecvRanks = level->interpolation.num_recvs;

   if (numSendRanks > 0) {
      #pragma omp target enter data map(to:level->interpolation.send_ranks[:numSendRanks])
      #pragma omp target enter data map(to:level->interpolation.send_sizes[:numSendRanks])
      #pragma omp target enter data map(alloc:level->interpolation.send_buffers[:numSendRanks])
      for(neighbor=0;neighbor<numSendRanks;neighbor++){
        #pragma omp target enter data map(alloc:level->interpolation.send_buffers[neighbor][:level->interpolation.send_sizes[neighbor]])
      }
   }

   if (numRecvRanks>0) {
      #pragma omp target enter data map(to:level->interpolation.recv_ranks[:numRecvRanks])
      #pragma omp target enter data map(to:level->interpolation.recv_sizes[:numRecvRanks])
      #pragma omp target enter data map(alloc:level->interpolation.send_buffers[:numRecvRanks])
      for(neighbor=0;neighbor<numRecvRanks;neighbor++){
        #pragma omp target enter data map(alloc:level->interpolation.recv_buffers[neighbor][:level->interpolation.recv_sizes[neighbor]])
      }
   }

   for (block=0;block<3;++block) {
     if (level->interpolation.blocks[block] != NULL) {
       #pragma omp target enter data map(to:level->interpolation.blocks[block][:level->interpolation.num_blocks[block]])
     }
   }

#endif
}

void build_exchange_ghosts_device_openmp(level_type *level) {
#if defined(SPEC_OPENMP_TARGET) && !defined(SPEC_MANAGED_MEMORY)
int shape,numSendRanks, numRecvRanks, neighbor;
int block;

// Create the exchange_ghosts
   for (shape =0;shape<STENCIL_MAX_SHAPES;shape++) {
        numSendRanks = level->exchange_ghosts[shape].num_sends;
        numRecvRanks = level->exchange_ghosts[shape].num_recvs;
        if ( numSendRanks > 0) {
           #pragma omp target enter data map(to:level->exchange_ghosts[shape].send_ranks[:numSendRanks])
           #pragma omp target enter data map(to:level->exchange_ghosts[shape].send_sizes[:numSendRanks])
           #pragma omp target enter data map(to:level->exchange_ghosts[shape].send_buffers[:numSendRanks])
	   for(neighbor=0;neighbor<numSendRanks;neighbor++){
             #pragma omp target enter data map(to:level->exchange_ghosts[shape].send_buffers[neighbor][:level->exchange_ghosts[shape].send_sizes[neighbor]])
	   }
	}

        if ( numRecvRanks > 0) {
           #pragma omp target enter data map(to:level->exchange_ghosts[shape].recv_ranks[:numRecvRanks])
           #pragma omp target enter data map(to:level->exchange_ghosts[shape].recv_sizes[:numRecvRanks])
           #pragma omp target enter data map(to:level->exchange_ghosts[shape].recv_buffers[:numRecvRanks])
           for(neighbor=0;neighbor<numRecvRanks;neighbor++){
             #pragma omp target enter data map(to:level->exchange_ghosts[shape].recv_buffers[neighbor][:level->exchange_ghosts[shape].recv_sizes[neighbor]])
           }
        }

        for (block=0;block<3;++block) {
           if (level->exchange_ghosts[shape].num_blocks[block] > 0) {
             #pragma omp target enter data map(to:level->exchange_ghosts[shape].blocks[block][:level->exchange_ghosts[shape].num_blocks[block]])
           }
        }
   }
#endif
}


void build_box_device_openmp(level_type *level,int flag) {
#if defined(SPEC_OPENMP_TARGET) && !defined(SPEC_MANAGED_MEMORY)
  int shape, block, b;
  for (shape=0;shape<STENCIL_MAX_SHAPES;shape++) {
    for (block=0;block<3;++block) {
      for (b=0;b<level->exchange_ghosts[shape].num_blocks[block];++b) {
	ATTACH_DEVICE_ADDRESS(level->exchange_ghosts[shape].blocks[block][b].read.ptr)
	ATTACH_DEVICE_ADDRESS(level->exchange_ghosts[shape].blocks[block][b].write.ptr)
      }
    }
  }

  if (flag == 1) {
// Retriction Boxes
    for (shape=0;shape<STENCIL_MAX_SHAPES;shape++) {
      for (block=0;block<3;++block) {
	for (b=0; b<level->restriction[shape].num_blocks[block];++b) {
	  ATTACH_DEVICE_ADDRESS(level->restriction[shape].blocks[block][b].read.ptr)
	  ATTACH_DEVICE_ADDRESS(level->restriction[shape].blocks[block][b].write.ptr)
	}
      }
    }
// Interplolation Boxes
    for (block=0;block<3;++block) {
      for (b=0; b<level->interpolation.num_blocks[block];++b) {
	ATTACH_DEVICE_ADDRESS(level->interpolation.blocks[block][b].read.ptr)
	ATTACH_DEVICE_ADDRESS(level->interpolation.blocks[block][b].write.ptr)
      }
    }
  }
#endif
}

void copy_vector_to_device_openmp(level_type *level) {
#if defined(SPEC_OPENMP_TARGET) && !defined(SPEC_MANAGED_MEMORY)
int numVectors = level->numVectors;
int num_my_boxes = level->num_my_boxes;
int box_volume = level->box_volume;
#ifdef VECTOR_MALLOC_BULK
#pragma omp target update to(level->vectors_base_1[:numVectors*num_my_boxes*box_volume])
#else
for (int i=0; i<numVectors; i++) {
  #pragma omp target update to(level->vectors[i][:num_my_boxes*box_volume])
}
#endif
#endif
}

void copy_one_vector_to_device_openmp(level_type *level, int var) {
#if defined(SPEC_OPENMP_TARGET) && !defined(SPEC_MANAGED_MEMORY)
int numVectors = level->numVectors;
int num_my_boxes = level->num_my_boxes;
int box_volume = level->box_volume;
if (var >= 0 && var < numVectors) {
#ifdef VECTOR_MALLOC_BULK
  int start = var*num_my_boxes*box_volume;
  int length = num_my_boxes*box_volume;
#pragma omp target update to(level->vectors_base_1[start:length])
#else
#pragma omp target update to(level->vectors[var][:num_my_boxes*box_volume])
#endif
}
#endif
}

void copy_vector_to_host_openmp(level_type *level) {
#if defined(SPEC_OPENMP_TARGET) && !defined(SPEC_MANAGED_MEMORY)
int numVectors = level->numVectors;
int num_my_boxes = level->num_my_boxes;
int box_volume = level->box_volume;
#ifdef VECTOR_MALLOC_BULK
#pragma omp target update from(level->vectors_base_1[:numVectors*num_my_boxes*box_volume])
#else
for (int i=0; i<numVectors; i++) {
  #pragma omp target update from(level->vectors[i][:num_my_boxes*box_volume])
}
#endif
#endif
}

void copy_one_vector_to_host_openmp(level_type *level, int var) {
#if defined(SPEC_OPENMP_TARGET) && !defined(SPEC_MANAGED_MEMORY)
int numVectors = level->numVectors;
int num_my_boxes = level->num_my_boxes;
int box_volume = level->box_volume;
if (var >= 0 && var < numVectors) {
#ifdef VECTOR_MALLOC_BULK
  int start = var*num_my_boxes*box_volume;
  int length = num_my_boxes*box_volume;
#pragma omp target update from(level->vectors_base_1[start:length])
#else
#pragma omp target update from(level->vectors[var][:num_my_boxes*box_volume])
#endif
}
#endif
}
