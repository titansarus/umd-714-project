#ifndef OFFLOAD_FNS
#define OFFLOAD_FNS
#include "defines.h"
#include "level.h"

void set_cuda_execution_array(void);
int use_cuda_function(int fid);
void check_consistent_build(void);

//------------------------------------------------------------------------------------------------------------------------------
// access policies for UM
#define UM_ACCESS_CPU  0
#define UM_ACCESS_GPU  1
#define UM_ACCESS_BOTH 2
// custom memory management routines
void* um_malloc(size_t size, int access_policy);
void* um_malloc_pinned(size_t size, int access_policy);
void* um_realloc(void *ptr, size_t size, int access_policy);
void  um_free(void *ptr, int access_policy);
void  um_free_pinned(void *ptr, int access_policy);
//------------------------------------------------------------------------------------------------------------------------------
void do_sync(void);
void do_sync_if_mixed(void);
int device_runtime_init(void);
void device_increment_block(level_type *level, int id, double prescale, communicator_type *exchange_ghosts, int block_type);
void device_copy_block(level_type *level, int id, communicator_type *exchange_ghosts, int block_type);
void device_extrapolate_betas(level_type *level, int shape);
void device_residual(level_type *level, int res_id, int x_id, int rhs_id, double a, double b);
void device_rebuild(level_type *level, int x_id, int Aii_id, int sumAbsAij_id, double a, double b);
void device_color_vector(level_type *d_level, int id_a, int colors_in_each_dim, int icolor, int jcolor, int kcolor);
void device_smooth(level_type *level, int x_id, int rhs_id, double a, double b, int s, double *c, double *d);
void device_restriction(level_type *level_c, int id_c, level_type *level_f, int id_f, communicator_type *restriction, int restrictionType, int block_type);
void device_apply_BCs_v1(level_type *level, int x_id, int shape);
void device_apply_BCs_v2(level_type *level, int x_id, int shape);
void device_apply_BCs_v4(level_type *level, int x_id, int shape);
void device_interpolation_v4(level_type *level_f, int id_f, double prescale_f, level_type *level_c, int id_c, communicator_type *interpolation, int block_type);
void device_interpolation_v2(level_type *level_f, int id_f, double prescale_f, level_type *level_c, int id_c, communicator_type *interpolation, int block_type);
void device_zero_vector(level_type *d_level, int id);
void device_scale_vector(level_type *d_level, int id_c, double scale_a, int id_a);
void device_shift_vector(level_type *d_level, int id_c, double shift_a, int id_a);
void device_mul_vectors(level_type *d_level, int id_c, double scale, int id_a, int id_b);
void device_add_vectors(level_type *d_level, int id_c, double scale_a, int id_a, double scale_b, int id_b);
double device_sum(level_type *d_level, int id);
double device_max_abs(level_type *d_level, int id);
#endif
