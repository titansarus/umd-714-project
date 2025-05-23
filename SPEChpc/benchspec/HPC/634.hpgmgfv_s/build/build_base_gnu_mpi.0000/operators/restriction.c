//------------------------------------------------------------------------------------------------------------------------------
// Samuel Williams
// SWWilliams@lbl.gov
// Lawrence Berkeley National Lab
//------------------------------------------------------------------------------------------------------------------------------
static inline void restriction_pc_block(level_type *level_c, int id_c, level_type *level_f, int id_f, blockCopy_type *block, int restrictionType){
  // restrict 3D array from read_i,j,k of read[] to write_i,j,k in write[] using piecewise constant restriction (cell averaged)
  int   dim_i       = block->dim.i; // calculate the dimensions of the resultant coarse block
  int   dim_j       = block->dim.j;
  int   dim_k       = block->dim.k;

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
     read_jStride = level_f->my_boxes[block->read.box ].jStride;
     read_kStride = level_f->my_boxes[block->read.box ].kStride;
     read = level_f->my_boxes[ block->read.box].vectors[id_f] + level_f->my_boxes[ block->read.box].ghosts*(1+level_f->my_boxes[ block->read.box].jStride+level_f->my_boxes[ block->read.box].kStride);
  }
  if(block->write.box>=0){
    write_jStride = level_c->my_boxes[block->write.box].jStride;
    write_kStride = level_c->my_boxes[block->write.box].kStride;
    write = level_c->my_boxes[block->write.box].vectors[id_c] + level_c->my_boxes[block->write.box].ghosts*(1+level_c->my_boxes[block->write.box].jStride+level_c->my_boxes[block->write.box].kStride);
  }



  int i,j,k;
  switch(restrictionType){
    case RESTRICT_CELL:
         for(k=0;k<dim_k;k++){
         for(j=0;j<dim_j;j++){
         for(i=0;i<dim_i;i++){
           int write_ijk = ((i   )+write_i) + ((j   )+write_j)*write_jStride + ((k   )+write_k)*write_kStride;
           int  read_ijk = ((i<<1)+ read_i) + ((j<<1)+ read_j)* read_jStride + ((k<<1)+ read_k)* read_kStride;
           write[write_ijk] = ( read[read_ijk                            ]+read[read_ijk+1                          ] +
                                read[read_ijk  +read_jStride             ]+read[read_ijk+1+read_jStride             ] +
                                read[read_ijk               +read_kStride]+read[read_ijk+1             +read_kStride] +
                                read[read_ijk  +read_jStride+read_kStride]+read[read_ijk+1+read_jStride+read_kStride] ) * 0.125;
         }}}break;
    case RESTRICT_FACE_I:
         for(k=0;k<dim_k;k++){
         for(j=0;j<dim_j;j++){
         for(i=0;i<dim_i;i++){
           int write_ijk = ((i   )+write_i) + ((j   )+write_j)*write_jStride + ((k   )+write_k)*write_kStride;
           int  read_ijk = ((i<<1)+ read_i) + ((j<<1)+ read_j)* read_jStride + ((k<<1)+ read_k)* read_kStride;
           write[write_ijk] = ( read[read_ijk                          ] +
                                read[read_ijk+read_jStride             ] +
                                read[read_ijk             +read_kStride] +
                                read[read_ijk+read_jStride+read_kStride] ) * 0.25;
         }}}break;
    case RESTRICT_FACE_J:
         for(k=0;k<dim_k;k++){
         for(j=0;j<dim_j;j++){
         for(i=0;i<dim_i;i++){
           int write_ijk = ((i   )+write_i) + ((j   )+write_j)*write_jStride + ((k   )+write_k)*write_kStride;
           int  read_ijk = ((i<<1)+ read_i) + ((j<<1)+ read_j)* read_jStride + ((k<<1)+ read_k)* read_kStride;
           write[write_ijk] = ( read[read_ijk               ] +
                                read[read_ijk+1             ] +
                                read[read_ijk  +read_kStride] +
                                read[read_ijk+1+read_kStride] ) * 0.25;
         }}}break;
    case RESTRICT_FACE_K:
         for(k=0;k<dim_k;k++){
         for(j=0;j<dim_j;j++){
         for(i=0;i<dim_i;i++){
           int write_ijk = ((i   )+write_i) + ((j   )+write_j)*write_jStride + ((k   )+write_k)*write_kStride;
           int  read_ijk = ((i<<1)+ read_i) + ((j<<1)+ read_j)* read_jStride + ((k<<1)+ read_k)* read_kStride;
           write[write_ijk] = ( read[read_ijk               ] +
                                read[read_ijk+1             ] +
                                read[read_ijk  +read_jStride] +
                                read[read_ijk+1+read_jStride] ) * 0.25;
         }}}break;
  }

}


//------------------------------------------------------------------------------------------------------------------------------
// perform a (inter-level) restriction on vector id_f of the fine level and stores the result in vector id_c on the coarse level
// restrictionType specifies whether this is either cell-averaged restriction, or one of three face-averaged restrictions
// piecewise constant restriction requires neither a ghost zone exchange nor a boundary condition
// This is a rather bulk synchronous implementation which packs all MPI buffers before initiating any sends
// Similarly, it waits for all remote data before copying any into local boxes.
// It does however attempt to overlap local restriction with MPI
void restriction(level_type * level_c, int id_c, level_type *level_f, int id_f, int restrictionType){
  double _timeCommunicationStart = getTime();
  double _timeStart,_timeEnd;
  int buffer=0;

  #ifdef USE_MPI
  int n;
  int my_tag = (level_f->tag<<4) | 0x5;
  // by convention, level_f allocates a combined array of requests for both level_f sends and level_c recvs...
  int nMessages = level_c->restriction[restrictionType].num_recvs + level_f->restriction[restrictionType].num_sends;
  MPI_Request *recv_requests = level_f->restriction[restrictionType].requests;
  MPI_Request *send_requests = level_f->restriction[restrictionType].requests + level_c->restriction[restrictionType].num_recvs;

  // loop through packed list of MPI receives and prepost Irecv's...
  if(level_c->restriction[restrictionType].num_recvs>0){
    _timeStart = getTime();
    #ifdef USE_MPI_THREAD_MULTIPLE
    #pragma omp parallel for schedule(dynamic,1)
    #endif
    for(n=0;n<level_c->restriction[restrictionType].num_recvs;n++){
      double * recv_buf = level_c->restriction[restrictionType].recv_buffers[n];

#if !defined(SPEC_MANAGED_MEMORY) && defined(SPEC_OPENACC) && defined(SPEC_ACCEL_AWARE_MPI)
      double * recv_buf_d = acc_deviceptr(recv_buf);
      if (level_c->use_offload && recv_buf_d != NULL) {
        recv_buf = recv_buf_d;
      }
#endif

#if !defined(SPEC_MANAGED_MEMORY) && defined(SPEC_OPENMP_TARGET) && defined(SPEC_ACCEL_AWARE_MPI)
      double * recv_buf_d = recv_buf;
#pragma omp target data use_device_ptr(recv_buf_d)
      {
	if (level_c->use_offload && recv_buf_d != NULL) {
	  recv_buf = recv_buf_d;
	}
	/* The OpenMP target data region remains open for the MPI_Irecv */
#endif

      MPI_Irecv(recv_buf,
                level_c->restriction[restrictionType].recv_sizes[n],
                MPI_DOUBLE,
                level_c->restriction[restrictionType].recv_ranks[n],
                my_tag,
                MPI_COMM_WORLD,
                &recv_requests[n]
      );

#if !defined(SPEC_MANAGED_MEMORY) && defined(SPEC_OPENMP_TARGET) && defined(SPEC_ACCEL_AWARE_MPI)
      } /* Close the OpenMP target data region */
#endif

    }
    _timeEnd = getTime();
    level_f->timers.restriction_recv += (_timeEnd-_timeStart);
  }


  // pack MPI send buffers...
  if(level_f->restriction[restrictionType].num_blocks[0]>0){
    _timeStart = getTime();
    if(level_f->use_offload) {
      /* We restrict from fine data already resident on the device. The
	 resulting coarse data will thus also be on the device. This
	 will be sent directly if SPEC_ACCEL_AWARE_MPI is defined.
	 Otherwise it will be copied to the host and sent from there */
      device_restriction(level_c,id_c,level_f,id_f,&level_f->restriction[restrictionType],restrictionType,0);
      do_sync(); // synchronize so the CPU sees the updated buffers which will be used for MPI transfers
      // CD: I think the do_sync() is only needed if SPEC_ACCEL_AWARE_MPI is not defined.
    }
    else {
    PRAGMA_THREAD_ACROSS_BLOCKS(level_f,buffer,level_f->restriction[restrictionType].num_blocks[0])
    for(buffer=0;buffer<level_f->restriction[restrictionType].num_blocks[0];buffer++){
      restriction_pc_block(level_c,id_c,level_f,id_f,&level_f->restriction[restrictionType].blocks[0][buffer],restrictionType);
    }
    }
    _timeEnd = getTime();
    level_f->timers.restriction_pack += (_timeEnd-_timeStart);
  }

 
  // loop through MPI send buffers and post Isend's...
  if(level_f->restriction[restrictionType].num_sends>0){
    _timeStart = getTime();
    #ifdef USE_MPI_THREAD_MULTIPLE
    #pragma omp parallel for schedule(dynamic,1)
    #endif
    for(n=0;n<level_f->restriction[restrictionType].num_sends;n++){
	double * send_buf = level_f->restriction[restrictionType].send_buffers[n];

#if !defined(SPEC_MANAGED_MEMORY) && defined(SPEC_OPENACC)
     double * send_buf_d = acc_deviceptr(send_buf);
     // CD: I would have thought it should be level_c->use_offload but this
     // causes very obvious failures
     if (level_f->use_offload && send_buf_d != NULL) {
#ifdef SPEC_ACCEL_AWARE_MPI
        send_buf = send_buf_d;
#else
        int isp = acc_is_present(send_buf,level_f->restriction[restrictionType].send_sizes[n]);
#ifdef DEBUG
printf("REST send: tag=%d n=%d size=%d send_buf=%p send_buf_d=%p isp=%d\n",level_f->tag,n,level_f->restriction[restrictionType].send_sizes[n],send_buf,acc_deviceptr(send_buf),isp);
#endif
        #pragma acc update host(send_buf[:level_f->restriction[restrictionType].send_sizes[n]]) if_present
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
#  pragma omp target update from(send_buf[:level_f->restriction[restrictionType].send_sizes[n]])
	    }
	  }
# endif
	} /* End of send_buf_d != NULL */
	/* The OpenMP target data region remains open for the MPI_Isend */
#endif

      MPI_Isend(send_buf,
                level_f->restriction[restrictionType].send_sizes[n],
                MPI_DOUBLE,
                level_f->restriction[restrictionType].send_ranks[n],
                my_tag,
                MPI_COMM_WORLD,
                &send_requests[n]
      );

#if !defined(SPEC_MANAGED_MEMORY) && defined(SPEC_OPENMP_TARGET)
      } /* Close the OpenMP target data region */
#endif

    }
    _timeEnd = getTime();
    level_f->timers.restriction_send += (_timeEnd-_timeStart);
  }
  #endif


  // perform local restriction[restrictionType]... try and hide within Isend latency... 
  if(level_f->restriction[restrictionType].num_blocks[1]>0){
    _timeStart = getTime();
    if (level_f->use_offload) {
      device_restriction(level_c, id_c, level_f, id_f, &level_f->restriction[restrictionType], restrictionType, 1);
      if (!level_c->use_offload) {
	/* If the coarser level is owned by the host we must copy the
	   just restricted coarse data from device to the host */
	do_sync();  // switchover point: must synchronize GPU
#if !defined(SPEC_MANAGED_MEMORY)
# if defined(SPEC_OPENACC)
	copy_one_vector_to_host(level_c, id_c);
# elif defined(SPEC_OPENMP_TARGET)
	copy_one_vector_to_host_openmp(level_c, id_c);
# endif
#endif
      }
    }
    else {
    PRAGMA_THREAD_ACROSS_BLOCKS(level_f,buffer,level_f->restriction[restrictionType].num_blocks[1])
    for(buffer=0;buffer<level_f->restriction[restrictionType].num_blocks[1];buffer++){
      restriction_pc_block(level_c,id_c,level_f,id_f,&level_f->restriction[restrictionType].blocks[1][buffer],restrictionType);
    }
    }
    _timeEnd = getTime();
    level_f->timers.restriction_local += (_timeEnd-_timeStart);
  }

  // wait for MPI to finish...
  #ifdef USE_MPI 
  if(nMessages){
    _timeStart = getTime();
    MPI_Waitall(nMessages,level_f->restriction[restrictionType].requests,level_f->restriction[restrictionType].status);
  #ifdef SYNC_DEVICE_AFTER_WAITALL
    do_sync();
  #endif
    _timeEnd = getTime();
    level_f->timers.restriction_wait += (_timeEnd-_timeStart);
   }

#if !defined(SPEC_MANAGED_MEMORY) && defined(SPEC_OPENACC) && !defined(SPEC_ACCEL_AWARE_MPI)
    if (level_c->use_offload) {
    for(n=0;n<level_c->restriction[restrictionType].num_recvs;n++){
       if (level_c->restriction[restrictionType].recv_sizes[n] > 0) {
         double * recv_buf = level_c->restriction[restrictionType].recv_buffers[n];
         int isp = acc_is_present(recv_buf,level_c->restriction[restrictionType].recv_sizes[n]);
#ifdef DEBUG
printf("REST recv tag=%d n=%d size=%d hptr=%p dptr=%p isp=%d\n",level_c->tag,n,level_c->restriction[restrictionType].recv_sizes[n],recv_buf,acc_deviceptr(recv_buf));
#endif
      #pragma acc update device(recv_buf[:level_c->restriction[restrictionType].recv_sizes[n]]) if_present
       }
    }
    }
#endif

#if !defined(SPEC_MANAGED_MEMORY) && defined(SPEC_OPENMP_TARGET) && !defined(SPEC_ACCEL_AWARE_MPI)
    if (level_c->use_offload && omp_get_num_devices() > 0) {
      for(n=0;n<level_c->restriction[restrictionType].num_recvs;n++){
	if (level_c->restriction[restrictionType].recv_sizes[n] > 0) {
	  double * recv_buf = level_c->restriction[restrictionType].recv_buffers[n];
	  int isp = omp_target_is_present(recv_buf, omp_get_default_device());
	  if (isp != 0) {
# pragma omp target update to(recv_buf[:level_c->restriction[restrictionType].recv_sizes[n]])
	  }
	}
      }
    }
#endif

  // unpack MPI receive buffers 
  if(level_c->restriction[restrictionType].num_blocks[2]>0){
    _timeStart = getTime();
    if(level_c->use_offload) {
      device_copy_block(level_c,id_c,&level_c->restriction[restrictionType],2);
    }
    else {
    PRAGMA_THREAD_ACROSS_BLOCKS(level_f,buffer,level_c->restriction[restrictionType].num_blocks[2])
    for(buffer=0;buffer<level_c->restriction[restrictionType].num_blocks[2];buffer++){
      CopyBlock(level_c,id_c,&level_c->restriction[restrictionType].blocks[2][buffer]);
    }
    }
    _timeEnd = getTime();
    level_f->timers.restriction_unpack += (_timeEnd-_timeStart);
  }
  #endif
 
 
  level_f->timers.restriction_total += (double)(getTime()-_timeCommunicationStart);
}
