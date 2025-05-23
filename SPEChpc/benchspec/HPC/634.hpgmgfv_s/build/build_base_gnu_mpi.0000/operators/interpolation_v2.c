//------------------------------------------------------------------------------------------------------------------------------
// Samuel Williams
// SWWilliams@lbl.gov
// Lawrence Berkeley National Lab
//------------------------------------------------------------------------------------------------------------------------------
#include <math.h>
//------------------------------------------------------------------------------------------------------------------------------
static inline void interpolation_v2_block(level_type *level_f, int id_f, double prescale_f, level_type *level_c, int id_c, blockCopy_type *block){
  // interpolate 3D array from read_i,j,k of read[] to write_i,j,k in write[] using volume averaged quadratic prolongation
  int write_dim_i   = block->dim.i<<1; // calculate the dimensions of the resultant fine block
  int write_dim_j   = block->dim.j<<1;
  int write_dim_k   = block->dim.k<<1;

  int  read_i       = block->read.i;
  int  read_j       = block->read.j;
  int  read_k       = block->read.k;
  int  read_jStride = block->read.jStride;
  int  read_kStride = block->read.kStride;

  int write_i       = block->write.i;
  int write_j       = block->write.j;
  int write_k       = block->write.k;
  int write_jStride = block->write.jStride;
  int write_kStride = block->write.kStride;

  double * __restrict__  read = block->read.ptr;
  double * __restrict__ write = block->write.ptr;
  if(block->read.box >=0){
     read_jStride = level_c->my_boxes[block->read.box ].jStride;
     read_kStride = level_c->my_boxes[block->read.box ].kStride;
     read = level_c->my_boxes[ block->read.box].vectors[id_c] + level_c->my_boxes[ block->read.box].ghosts*(1+ read_jStride+ read_kStride);
  }
  if(block->write.box>=0){
    write_jStride = level_f->my_boxes[block->write.box].jStride;
    write_kStride = level_f->my_boxes[block->write.box].kStride;
    write = level_f->my_boxes[block->write.box].vectors[id_f] + level_f->my_boxes[block->write.box].ghosts*(1+write_jStride+write_kStride);
  }
 

  #ifdef USE_NAIVE_INTERP
  // naive 27pt per fine grid cell
  int i,j,k;
  double c1 = 1.0/8.0;
  for(k=0;k<write_dim_k;k++){double c1k=c1;if(k&0x1){c1k=-c1;}
  for(j=0;j<write_dim_j;j++){double c1j=c1;if(j&0x1){c1j=-c1;}
  for(i=0;i<write_dim_i;i++){double c1i=c1;if(i&0x1){c1i=-c1;}
    int write_ijk = ((i   )+write_i) + (((j   )+write_j)*write_jStride) + (((k   )+write_k)*write_kStride);
    int  read_ijk = ((i>>1)+ read_i) + (((j>>1)+ read_j)* read_jStride) + (((k>>1)+ read_k)* read_kStride);
    //
    // |  1/8  |  1.0  | -1/8  | coarse grid
    // |---+---|---+---|---+---|
    // |   |   |???|   |   |   | fine grid
    //
    write[write_ijk] = prescale_f*write[write_ijk] +
                       + c1k*( + c1j*( c1i*read[read_ijk-1-read_jStride-read_kStride] + read[read_ijk-read_jStride-read_kStride] - c1i*read[read_ijk+1-read_jStride-read_kStride] )
                               +     ( c1i*read[read_ijk-1             -read_kStride] + read[read_ijk             -read_kStride] - c1i*read[read_ijk+1             -read_kStride] )
                               - c1j*( c1i*read[read_ijk-1+read_jStride-read_kStride] + read[read_ijk+read_jStride-read_kStride] - c1i*read[read_ijk+1+read_jStride-read_kStride] ) )
                       +     ( + c1j*( c1i*read[read_ijk-1-read_jStride             ] + read[read_ijk-read_jStride             ] - c1i*read[read_ijk+1-read_jStride             ] )
                               +     ( c1i*read[read_ijk-1                          ] + read[read_ijk                          ] - c1i*read[read_ijk+1                          ] )
                               - c1j*( c1i*read[read_ijk-1+read_jStride             ] + read[read_ijk+read_jStride             ] - c1i*read[read_ijk+1+read_jStride             ] ) )
                       - c1k*( + c1j*( c1i*read[read_ijk-1-read_jStride+read_kStride] + read[read_ijk-read_jStride+read_kStride] - c1i*read[read_ijk+1-read_jStride+read_kStride] )
                               +     ( c1i*read[read_ijk-1             +read_kStride] + read[read_ijk             +read_kStride] - c1i*read[read_ijk+1             +read_kStride] )
                               - c1j*( c1i*read[read_ijk-1+read_jStride+read_kStride] + read[read_ijk+read_jStride+read_kStride] - c1i*read[read_ijk+1+read_jStride+read_kStride] ) );
  }}}
  #else
  int i,j,k;
  double c1 = 1.0/8.0;
  for(k=0;k<write_dim_k;k+=2){
  for(j=0;j<write_dim_j;j+=2){
  for(i=0;i<write_dim_i;i+=2){
    int write_ijk = ((i   )+write_i) + (((j   )+write_j)*write_jStride) + (((k   )+write_k)*write_kStride);
    int  read_ijk = ((i>>1)+ read_i) + (((j>>1)+ read_j)* read_jStride) + (((k>>1)+ read_k)* read_kStride);
    //
    // |  1/8  |  1.0  | -1/8  | coarse grid
    // |---+---|---+---|---+---|
    // |   |   |???|   |   |   | fine grid
    //

    // grab all coarse grid points...
    const double c000=read[read_ijk-1-read_jStride-read_kStride], c100=read[read_ijk  -read_jStride-read_kStride], c200=read[read_ijk+1-read_jStride-read_kStride];
    const double c010=read[read_ijk-1             -read_kStride], c110=read[read_ijk               -read_kStride], c210=read[read_ijk+1             -read_kStride];
    const double c020=read[read_ijk-1+read_jStride-read_kStride], c120=read[read_ijk  +read_jStride-read_kStride], c220=read[read_ijk+1+read_jStride-read_kStride];
    const double c001=read[read_ijk-1-read_jStride             ], c101=read[read_ijk  -read_jStride             ], c201=read[read_ijk+1-read_jStride             ];
    const double c011=read[read_ijk-1                          ], c111=read[read_ijk                            ], c211=read[read_ijk+1                          ];
    const double c021=read[read_ijk-1+read_jStride             ], c121=read[read_ijk  +read_jStride             ], c221=read[read_ijk+1+read_jStride             ];
    const double c002=read[read_ijk-1-read_jStride+read_kStride], c102=read[read_ijk  -read_jStride+read_kStride], c202=read[read_ijk+1-read_jStride+read_kStride];
    const double c012=read[read_ijk-1             +read_kStride], c112=read[read_ijk               +read_kStride], c212=read[read_ijk+1             +read_kStride];
    const double c022=read[read_ijk-1+read_jStride+read_kStride], c122=read[read_ijk  +read_jStride+read_kStride], c222=read[read_ijk+1+read_jStride+read_kStride];

    // interpolate in i to create fine i / coarse jk points...
    //
    // +-------+-------+-------+      :.......+---+---+.......:
    // |       |       |       |      :       |   |   |       :
    // |   c   |   c   |   c   |      :       | f | f |       :
    // |       |       |       |      :       |   |   |       :
    // +-------+-------+-------+      :.......+---+---+.......:
    // |       |       |       |      :       |   |   |       :
    // |   c   |   c   |   c   |  ->  :       | f | f |       :
    // |       |       |       |      :       |   |   |       :
    // +-------+-------+-------+      :.......+---+---+.......:
    // |       |       |       |      :       |   |   |       :
    // |   c   |   c   |   c   |      :       | f | f |       :
    // |       |       |       |      :       |   |   |       :
    // +-------+-------+-------+      :.......+---+---+.......:
    //
    const double f0c00 = ( c100 + c1*(c000-c200) ); // same as original 3pt stencil... f0c00 = ( c1*c000 + c100 - c1*c200 );
    const double f1c00 = ( c100 - c1*(c000-c200) );
    const double f0c10 = ( c110 + c1*(c010-c210) );
    const double f1c10 = ( c110 - c1*(c010-c210) );
    const double f0c20 = ( c120 + c1*(c020-c220) );
    const double f1c20 = ( c120 - c1*(c020-c220) );

    const double f0c01 = ( c101 + c1*(c001-c201) );
    const double f1c01 = ( c101 - c1*(c001-c201) );
    const double f0c11 = ( c111 + c1*(c011-c211) );
    const double f1c11 = ( c111 - c1*(c011-c211) );
    const double f0c21 = ( c121 + c1*(c021-c221) );
    const double f1c21 = ( c121 - c1*(c021-c221) );

    const double f0c02 = ( c102 + c1*(c002-c202) );
    const double f1c02 = ( c102 - c1*(c002-c202) );
    const double f0c12 = ( c112 + c1*(c012-c212) );
    const double f1c12 = ( c112 - c1*(c012-c212) );
    const double f0c22 = ( c122 + c1*(c022-c222) );
    const double f1c22 = ( c122 - c1*(c022-c222) );

    // interpolate in j to create fine ij / coarse k points...
    //
    // :.......+---+---+.......:      :.......:.......:.......:
    // :       |   |   |       :      :       :       :       :
    // :       |   |   |       :      :       :       :       :
    // :       |   |   |       :      :       :       :       :
    // :.......+---+---+.......:      :.......+---+---+.......:
    // :       |   |   |       :      :       |   |   |       :
    // :       |   |   |       :  ->  :       +---+---+       :
    // :       |   |   |       :      :       |   |   |       :
    // :.......+---+---+.......:      :.......+---+---+.......:
    // :       |   |   |       :      :       :       :       :
    // :       |   |   |       :      :       :       :       :
    // :       |   |   |       :      :       :       :       :
    // :.......+---+---+.......:      :.......:.......:.......:
    //
    const double f00c0 = ( f0c10 + c1*(f0c00-f0c20) );
    const double f10c0 = ( f1c10 + c1*(f1c00-f1c20) );
    const double f01c0 = ( f0c10 - c1*(f0c00-f0c20) );
    const double f11c0 = ( f1c10 - c1*(f1c00-f1c20) );

    const double f00c1 = ( f0c11 + c1*(f0c01-f0c21) );
    const double f10c1 = ( f1c11 + c1*(f1c01-f1c21) );
    const double f01c1 = ( f0c11 - c1*(f0c01-f0c21) );
    const double f11c1 = ( f1c11 - c1*(f1c01-f1c21) );

    const double f00c2 = ( f0c12 + c1*(f0c02-f0c22) );
    const double f10c2 = ( f1c12 + c1*(f1c02-f1c22) );
    const double f01c2 = ( f0c12 - c1*(f0c02-f0c22) );
    const double f11c2 = ( f1c12 - c1*(f1c02-f1c22) );

    // interpolate in k to create fine ijk points...
    const double f000 = ( f00c1 + c1*(f00c0-f00c2) );
    const double f100 = ( f10c1 + c1*(f10c0-f10c2) );
    const double f010 = ( f01c1 + c1*(f01c0-f01c2) );
    const double f110 = ( f11c1 + c1*(f11c0-f11c2) );
    const double f001 = ( f00c1 - c1*(f00c0-f00c2) );
    const double f101 = ( f10c1 - c1*(f10c0-f10c2) );
    const double f011 = ( f01c1 - c1*(f01c0-f01c2) );
    const double f111 = ( f11c1 - c1*(f11c0-f11c2) );

    // commit to memory...
    write[write_ijk                              ] = prescale_f*write[write_ijk                              ] + f000;
    write[write_ijk+1                            ] = prescale_f*write[write_ijk+1                            ] + f100;
    write[write_ijk  +write_jStride              ] = prescale_f*write[write_ijk  +write_jStride              ] + f010;
    write[write_ijk+1+write_jStride              ] = prescale_f*write[write_ijk+1+write_jStride              ] + f110;
    write[write_ijk                +write_kStride] = prescale_f*write[write_ijk                +write_kStride] + f001;
    write[write_ijk+1              +write_kStride] = prescale_f*write[write_ijk+1              +write_kStride] + f101;
    write[write_ijk  +write_jStride+write_kStride] = prescale_f*write[write_ijk  +write_jStride+write_kStride] + f011;
    write[write_ijk+1+write_jStride+write_kStride] = prescale_f*write[write_ijk+1+write_jStride+write_kStride] + f111;
  }}}
  #endif

}


//------------------------------------------------------------------------------------------------------------------------------
// perform a (inter-level) volumetric quadratic interpolation on vector id_c of the coarse level and increments prescale_f*vector id_f on the fine level by the result
// i.e. id_f = prescale_f*id_f + P*id_c
// prescale_f is nominally 1.0 or 0.0
// quadratic interpolation requires a full ghost zone exchange and boundary condition
// This is a rather bulk synchronous implementation which packs all MPI buffers before initiating any sends
// Similarly, it waits for all remote data before copying any into local boxes.
// It does however attempt to overlap local interpolation with MPI
void interpolation_v2(level_type * level_f, int id_f, double prescale_f, level_type *level_c, int id_c){
    exchange_boundary(level_c,id_c,STENCIL_SHAPE_BOX);
         apply_BCs_v2(level_c,id_c,STENCIL_SHAPE_BOX);

  double _timeCommunicationStart = getTime();
  double _timeStart,_timeEnd;
  int buffer=0;

  #ifdef USE_MPI
  int n;
  int my_tag = (level_f->tag<<4) | 0x7;
  // by convention, level_f allocates a combined array of requests for both level_f recvs and level_c sends...
  int nMessages = level_c->interpolation.num_sends + level_f->interpolation.num_recvs;
  MPI_Request *recv_requests = level_f->interpolation.requests;
  MPI_Request *send_requests = level_f->interpolation.requests + level_f->interpolation.num_recvs;


  // loop through packed list of MPI receives and prepost Irecv's...
  if(level_f->interpolation.num_recvs>0){
    _timeStart = getTime();
    #ifdef USE_MPI_THREAD_MULTIPLE
    #pragma omp parallel for schedule(dynamic,1)
    #endif
    for(n=0;n<level_f->interpolation.num_recvs;n++){
      double * recv_buf = level_f->interpolation.recv_buffers[n];

#if !defined(SPEC_MANAGED_MEMORY) && defined(SPEC_OPENACC) && defined(SPEC_ACCEL_AWARE_MPI)
      double * recv_buf_d = acc_deviceptr(recv_buf);
      if (level_f->use_offload && recv_buf_d != NULL) {
        recv_buf = recv_buf_d;
      }
#endif

#if !defined(SPEC_MANAGED_MEMORY) && defined(SPEC_OPENMP_TARGET) && defined(SPEC_ACCEL_AWARE_MPI)
      double * recv_buf_d = recv_buf;
#pragma omp target data use_device_ptr(recv_buf_d)
      {
	if (level_f->use_offload && recv_buf_d != NULL) {
	  recv_buf = recv_buf_d;
	}
	/* The OpenMP target data region remains open for the MPI_Irecv */
#endif

      MPI_Irecv(recv_buf,
                level_f->interpolation.recv_sizes[n],
                MPI_DOUBLE,
                level_f->interpolation.recv_ranks[n],
                my_tag,
                MPI_COMM_WORLD,
                &recv_requests[n]
      );

#if !defined(SPEC_MANAGED_MEMORY) && defined(SPEC_OPENMP_TARGET) && defined(SPEC_ACCEL_AWARE_MPI)
      } /* Close the OpenMP target data region */
#endif

    }
    _timeEnd = getTime();
    level_f->timers.interpolation_recv += (_timeEnd-_timeStart);
  }


  // pack MPI send buffers...
  if(level_c->interpolation.num_blocks[0]>0){
    _timeStart = getTime();
    if(level_f->use_offload) {
      if(!level_c->use_offload) {
	/* Coarser level is on host. Therefore copy coarse level from host to device */
#if !defined(SPEC_MANAGED_MEMORY)
# if defined(SPEC_OPENACC)
	copy_one_vector_to_device(level_c, id_c);
# elif defined(SPEC_OPENMP_TARGET)
	copy_one_vector_to_device_openmp(level_c, id_c);
# endif
#endif
      }
      device_interpolation_v2(level_f,id_f,0.0,level_c,id_c,&level_c->interpolation,0);
      do_sync();  // synchronize so that CPU can see updated buffers
    }
    else {
    PRAGMA_THREAD_ACROSS_BLOCKS(level_f,buffer,level_c->interpolation.num_blocks[0])
    for(buffer=0;buffer<level_c->interpolation.num_blocks[0];buffer++){
      // !!! prescale==0 because you don't want to increment the MPI buffer
      interpolation_v2_block(level_f,id_f,0.0,level_c,id_c,&level_c->interpolation.blocks[0][buffer]);
    }
    }
    _timeEnd = getTime();
    level_f->timers.interpolation_pack += (_timeEnd-_timeStart);
  }


  // loop through MPI send buffers and post Isend's...
  if(level_c->interpolation.num_sends>0){
    _timeStart = getTime();
    #ifdef USE_MPI_THREAD_MULTIPLE
    #pragma omp parallel for schedule(dynamic,1)
    #endif
    for(n=0;n<level_c->interpolation.num_sends;n++){
      double * send_buf = level_c->interpolation.send_buffers[n];

#if !defined(SPEC_MANAGED_MEMORY) && defined(SPEC_OPENACC)
      double * send_buf_d = acc_deviceptr(send_buf);
      if (level_f->use_offload && send_buf_d != NULL) {
#ifdef SPEC_ACCEL_AWARE_MPI
        send_buf = send_buf_d;
#else
        int isp = acc_is_present(send_buf,level_c->interpolation.send_sizes[n]);
#ifdef DEBUG
	printf("INTER2 send: tag=%d n=%d size=%d send_buf=%p send_buf_d=%p isp=%d\n",level_c->tag,n,level_c->interpolation.send_sizes[n],send_buf,acc_deviceptr(send_buf),isp);
#endif
        #pragma acc update host(send_buf[:level_c->interpolation.send_sizes[n]]) if_present
#endif
      }
#endif

#if !defined(SPEC_MANAGED_MEMORY) && defined(SPEC_OPENMP_TARGET)
      double * send_buf_d = send_buf;
# pragma omp target data use_device_ptr(send_buf_d)
      {
	if (level_f->use_offload && send_buf_d != NULL) {
# ifdef SPEC_ACCEL_AWARE_MPI
	  send_buf = send_buf_d;
# else
	  if (omp_get_num_devices() > 0) {
	    int isp = omp_target_is_present(send_buf, omp_get_default_device());
	    if (isp != 0) {
#  pragma omp target update from(send_buf[:level_c->interpolation.send_sizes[n]])
	    }
	  }
# endif
	} /* End of send_buf_d != NULL */
	/* The OpenMP target data region remains open for the MPI_Isend */
#endif

      MPI_Isend(send_buf,
                level_c->interpolation.send_sizes[n],
                MPI_DOUBLE,
                level_c->interpolation.send_ranks[n],
                my_tag,
                MPI_COMM_WORLD,
                &send_requests[n]
      );

#if !defined(SPEC_MANAGED_MEMORY) && defined(SPEC_OPENMP_TARGET)
      } /* Close the OpenMP target data region */
#endif

    }
    _timeEnd = getTime();
    level_f->timers.interpolation_send += (_timeEnd-_timeStart);
  }
  #endif


  // perform local interpolation... try and hide within Isend latency... 
  if(level_c->interpolation.num_blocks[1]>0){
    _timeStart = getTime();
    if(level_f->use_offload){
      if(!level_c->use_offload) {
	/* Coarser level is on host. Therefore copy coarse level from host to device */
	// printf("%s: Copying coarse level from host to device\n", __FILE__);
#if !defined(SPEC_MANAGED_MEMORY)
# if defined(SPEC_OPENACC)
	copy_one_vector_to_device(level_c, id_c);
# elif defined(SPEC_OPENMP_TARGET)
	copy_one_vector_to_device_openmp(level_c, id_c);
# endif
#endif
      }
      device_interpolation_v2(level_f,id_f,prescale_f,level_c,id_c,&level_c->interpolation,1);
    }
    else {
    PRAGMA_THREAD_ACROSS_BLOCKS(level_f,buffer,level_c->interpolation.num_blocks[1])
    for(buffer=0;buffer<level_c->interpolation.num_blocks[1];buffer++){
      interpolation_v2_block(level_f,id_f,prescale_f,level_c,id_c,&level_c->interpolation.blocks[1][buffer]);
    }
    }
    _timeEnd = getTime();
    level_f->timers.interpolation_local += (_timeEnd-_timeStart);
  }


  // wait for MPI to finish...
  #ifdef USE_MPI 
  if(nMessages>0){
    _timeStart = getTime();
    MPI_Waitall(nMessages,level_f->interpolation.requests,level_f->interpolation.status);
  #ifdef SYNC_DEVICE_AFTER_WAITALL
    do_sync();
  #endif

#if !defined(SPEC_MANAGED_MEMORY) && defined(SPEC_OPENACC) && !defined(SPEC_ACCEL_AWARE_MPI)
    if (level_f->use_offload) {
    for(n=0;n<level_f->interpolation.num_recvs;n++){
       if (level_f->interpolation.recv_sizes[n] > 0) {
	 double * recv_buf = level_f->interpolation.recv_buffers[n];
#ifdef DEBUG
	 printf("INTER2 recv tag=%d n=%d size=%d hptr=%p dptr=%p\n",level_f->tag,n,level_f->interpolation.recv_sizes[n],level_f->interpolation.recv_buffers[n],acc_deviceptr(recv_buf));
#endif
         #pragma acc update device(recv_buf[:level_f->interpolation.recv_sizes[n]]) if_present
       }
    }
    }
#endif

#if !defined(SPEC_MANAGED_MEMORY) && defined(SPEC_OPENMP_TARGET) && !defined(SPEC_ACCEL_AWARE_MPI)
    if (level_f->use_offload && omp_get_num_devices() > 0) {
      for(n=0;n<level_f->interpolation.num_recvs;n++){
	if (level_f->interpolation.recv_sizes[n] > 0) {
	  double * recv_buf = level_f->interpolation.recv_buffers[n];
	  int isp = omp_target_is_present(recv_buf, omp_get_default_device());
	  if (isp != 0) {
# pragma omp target update to(recv_buf[:level_f->interpolation.recv_sizes[n]])
	  }
	}
      }
    }
#endif

    _timeEnd = getTime();
    level_f->timers.interpolation_wait += (_timeEnd-_timeStart);
  }


  // unpack MPI receive buffers 
  if(level_f->interpolation.num_blocks[2]>0){
    _timeStart = getTime();
    if(level_f->use_offload) {
      device_increment_block(level_f,id_f,prescale_f,&level_f->interpolation,2);
    }
    else {
    PRAGMA_THREAD_ACROSS_BLOCKS(level_f,buffer,level_f->interpolation.num_blocks[2])
    for(buffer=0;buffer<level_f->interpolation.num_blocks[2];buffer++){
      IncrementBlock(level_f,id_f,prescale_f,&level_f->interpolation.blocks[2][buffer]);
    }
    }
    _timeEnd = getTime();
    level_f->timers.interpolation_unpack += (_timeEnd-_timeStart);
  }
  #endif 
 
 
  level_f->timers.interpolation_total += (double)(getTime()-_timeCommunicationStart);
}
