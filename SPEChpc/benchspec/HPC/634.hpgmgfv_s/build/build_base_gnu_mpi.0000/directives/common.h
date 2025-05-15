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

// available device functions
int directives_runtime_init(void);
void directives_apply_BCs_v2(level_type *level, int x_id, int shape);
void directives_apply_BCs_v4(level_type *level, int x_id, int shape);
void directives_extrapolate_betas(level_type *level, int shape);
void directives_zero_ghost_region(level_type *level, int x_id, int shape);
void directives_interpolation_v4(level_type *level_f, int id_f, double prescale_f, level_type *level_c, int id_c, communicator_type *interpolation, int block_type);
void directives_copy_block(level_type *level, int id, communicator_type *exchange_ghosts, int block_type);
void directives_increment_block(level_type *level, int id, double prescale, communicator_type *exchange_ghosts, int block_type);
void directives_smooth_unoptimized(level_type *level, int x_id, int rhs_id, double a, double b, int s, double *c, double *d);
void directives_smooth_optimized(level_type *level, int x_id, int rhs_id, double a, double b, int s, double *c, double *d);
void directives_residual(level_type *level, int res_id, int x_id, int rhs_id, double a, double b, int REBUILD);
void directives_zero_vector(level_type *level, int component_id);
void directives_axpy_vector(level_type *level, int id_c, double scale_a, double shift_a, double scale_b, int id_a, int id_b, int mul_vectors);
void directives_scale_vector(level_type *d_level, int id_c, double scale_a, int id_a);
void directives_shift_vector(level_type *d_level, int id_c, double shift_a, int id_a);
void directives_mul_vectors(level_type *d_level, int id_c, double scale, int id_a, int id_b);
void directives_add_vectors(level_type *d_level, int id_c, double scale_a, int id_a, double scale_b, int id_b);
void directives_color_vector(level_type *level, int id_a, int colors_in_each_dim, int icolor, int jcolor, int kcolor);
double directives_reduction(level_type *level, int id, int red_type);
double directives_sum(level_type *d_level, int id);
double directives_max_abs(level_type *d_level, int id);
void directives_restriction(level_type *level_c, int id_c, level_type *level_f, int id_f, communicator_type *restriction, int restrictionType, int block_type);
void directives_interpolation_v2(level_type *level_f, int id_f, double prescale_f, level_type *level_c, int id_c, communicator_type *interpolation, int block_type);
void directives_apply_BCs_v2_unoptimized(level_type *level, int x_id, int shape);
