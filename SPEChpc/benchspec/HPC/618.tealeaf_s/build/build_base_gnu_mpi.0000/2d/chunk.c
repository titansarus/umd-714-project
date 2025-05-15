#include "chunk.h"

// Initialise the chunk
void initialise_chunk(Chunk* chunk, Settings* settings, int x, int y)
{
  // Initialise the key variables
  chunk->x = x + settings->halo_depth*2;
  chunk->y = y + settings->halo_depth*2;
  chunk->dt_init = settings->dt_init;

  // Allocate the neighbour list
  chunk->neighbours = (int*)malloc(sizeof(int)*NUM_FACES);

  // Allocate the MPI comm buffers
  int lr_len = chunk->y*settings->halo_depth*NUM_FIELDS;
  chunk->left_send = (double*)malloc(sizeof(double)*lr_len);
  chunk->left_recv = (double*)malloc(sizeof(double)*lr_len);
  chunk->right_send = (double*)malloc(sizeof(double)*lr_len);
  chunk->right_recv = (double*)malloc(sizeof(double)*lr_len);

  int tb_len = chunk->x*settings->halo_depth*NUM_FIELDS;
  chunk->top_send = (double*)malloc(sizeof(double)*tb_len);
  chunk->top_recv = (double*)malloc(sizeof(double)*tb_len);
  chunk->bottom_send = (double*)malloc(sizeof(double)*tb_len);
  chunk->bottom_recv = (double*)malloc(sizeof(double)*tb_len);

#if defined(SPEC_OPENACC) && defined(SPEC_ACCEL_AWARE_MPI)
#pragma acc enter data create(                                    \
	chunk->left_send[:lr_len], chunk->left_recv[:lr_len],     \
	chunk->right_send[:lr_len], chunk->right_recv[:lr_len],   \
	chunk->top_send[:tb_len], chunk->top_recv[:tb_len],       \
        chunk->bottom_send[:tb_len], chunk->bottom_recv[:tb_len]) 
#endif
	

  // Initialise the ChunkExtension, which allows composition of extended
  // fields specific to individual implementations
  chunk->ext = (ChunkExtension*)malloc(sizeof(ChunkExtension));
}

// Finalise the chunk
void finalise_chunk(Chunk* chunk)
{
#if defined(SPEC_OPENACC) && defined(SPEC_ACCEL_AWARE_MPI)
#pragma acc exit data delete(                     \
	chunk->left_send, chunk->left_recv,     \
	chunk->right_send, chunk->right_recv,   \
	chunk->top_send, chunk->top_recv,       \
        chunk->bottom_send, chunk->bottom_recv) 
#endif
  free(chunk->neighbours);
  free(chunk->ext);
  free(chunk->left_send);
  free(chunk->left_recv);
  free(chunk->right_send);
  free(chunk->right_recv);
  free(chunk->top_send);
  free(chunk->top_recv);
  free(chunk->bottom_send);
  free(chunk->bottom_recv);
}
