//------------------------------------------------------------------------------------------------------------------------------
// Samuel Williams
// SWWilliams@lbl.gov
// Lawrence Berkeley National Lab
//------------------------------------------------------------------------------------------------------------------------------
// perform a (intra-level) ghost zone exchange on vector id
//  NOTE exchange_boundary() only exchanges the boundary.  
//  It will not enforce any boundary conditions
//  BC's are either the responsibility of a separate function or should be fused into the stencil
// The argument shape indicates which of faces, edges, and corners on each box must be exchanged
//  If the specified shape exceeds the range of defined shapes, the code will default to STENCIL_SHAPE_BOX (i.e. exchange faces, edges, and corners)
void exchange_boundary(level_type * level, int id, int shape){
  double _timeCommunicationStart = getTime();
  double _timeStart,_timeEnd;

  if(shape>=STENCIL_MAX_SHAPES)shape=STENCIL_SHAPE_BOX;  // shape must be < STENCIL_MAX_SHAPES in order to safely index into exchange_ghosts[]
  int buffer=0;

  #ifdef USE_MPI
  int my_tag = (level->tag<<4) | shape;
  int n;
  int nMessages = level->exchange_ghosts[shape].num_recvs + level->exchange_ghosts[shape].num_sends;
  MPI_Request *recv_requests = level->exchange_ghosts[shape].requests;
  MPI_Request *send_requests = level->exchange_ghosts[shape].requests + level->exchange_ghosts[shape].num_recvs;

  // TODO: investigate why this is necessary for multi-GPU runs
  if(level->use_offload && (level->num_ranks > 1))
    do_sync();

  // loop through packed list of MPI receives and prepost Irecv's...
  if(level->exchange_ghosts[shape].num_recvs>0){
    _timeStart = getTime();

    #ifdef USE_MPI_THREAD_MULTIPLE
    #pragma omp parallel for schedule(dynamic,1)
    #endif
    for(n=0;n<level->exchange_ghosts[shape].num_recvs;n++){
      double * recv_buf = level->exchange_ghosts[shape].recv_buffers[n];

#if !defined(SPEC_MANAGED_MEMORY) && defined(SPEC_OPENACC) && defined(SPEC_ACCEL_AWARE_MPI)
      double * recv_buf_d = acc_deviceptr(recv_buf);
      if (level->use_offload && recv_buf_d != NULL) {
	recv_buf = recv_buf_d;
      }
#endif

#if !defined(SPEC_MANAGED_MEMORY) && defined(SPEC_OPENMP_TARGET) && defined(SPEC_ACCEL_AWARE_MPI)
      double * recv_buf_d = recv_buf;
#pragma omp target data use_device_ptr(recv_buf_d)
      {
	if (level->use_offload && recv_buf_d != NULL) {
	  recv_buf = recv_buf_d;
	}
	/* The OpenMP target data region remains open for the MPI_Irecv */
#endif

      MPI_Irecv(recv_buf,
                level->exchange_ghosts[shape].recv_sizes[n],
                MPI_DOUBLE,
                level->exchange_ghosts[shape].recv_ranks[n],
                my_tag,
                MPI_COMM_WORLD,
                &recv_requests[n]
      );

#if !defined(SPEC_MANAGED_MEMORY) && defined(SPEC_OPENMP_TARGET) && defined(SPEC_ACCEL_AWARE_MPI)
      } /* Close the OpenMP target data region */
#endif

    }
    _timeEnd = getTime();
    level->timers.ghostZone_recv += (_timeEnd-_timeStart);
  }


  // pack MPI send buffers...
  if(level->exchange_ghosts[shape].num_blocks[0]){
    _timeStart = getTime();
    if(level->use_offload) {
      device_copy_block(level,id,&level->exchange_ghosts[shape],0);
      do_sync();	// synchronize so the CPU sees the updated buffers which will be used for MPI transfers
    }
    else {
    PRAGMA_THREAD_ACROSS_BLOCKS(level,buffer,level->exchange_ghosts[shape].num_blocks[0])
    for(buffer=0;buffer<level->exchange_ghosts[shape].num_blocks[0];buffer++){
      CopyBlock(level,id,&level->exchange_ghosts[shape].blocks[0][buffer]);
    }
    }
    _timeEnd = getTime();
    level->timers.ghostZone_pack += (_timeEnd-_timeStart);
  }

  // loop through MPI send buffers and post Isend's...
  if(level->exchange_ghosts[shape].num_sends>0){
    _timeStart = getTime();
#if defined(USE_MPI_THREAD_MULTIPLE)
    #pragma omp parallel for schedule(dynamic,1)
#endif
    for(n=0;n<level->exchange_ghosts[shape].num_sends;n++){
      double * send_buf = level->exchange_ghosts[shape].send_buffers[n];

#if !defined(SPEC_MANAGED_MEMORY) && defined(SPEC_OPENACC)
      double * send_buf_d = acc_deviceptr(send_buf);
      if (level->use_offload && send_buf_d != NULL) {
# ifdef SPEC_ACCEL_AWARE_MPI
        send_buf = send_buf_d;
# else
        int isp = acc_is_present(send_buf,level->exchange_ghosts[shape].send_sizes[n]);
        #pragma acc update host(send_buf[:level->exchange_ghosts[shape].send_sizes[n]]) if_present
# endif
      }
#endif

#if !defined(SPEC_MANAGED_MEMORY) && defined(SPEC_OPENMP_TARGET)
      double * send_buf_d = send_buf;
# pragma omp target data use_device_ptr(send_buf_d)
      {
	if (level->use_offload && send_buf_d != NULL) {
# ifdef SPEC_ACCEL_AWARE_MPI
	  send_buf = send_buf_d;
# else
	  if (omp_get_num_devices() > 0) {
	    int isp = omp_target_is_present(send_buf, omp_get_default_device());
	    if (isp != 0) {
#  pragma omp target update from(send_buf[:level->exchange_ghosts[shape].send_sizes[n]])
	    }
	  }
# endif
	} /* End of send_buf_d != NULL */
	/* The OpenMP target data region remains open for the MPI_Isend */
#endif

      MPI_Isend(send_buf,
                level->exchange_ghosts[shape].send_sizes[n],
                MPI_DOUBLE,
                level->exchange_ghosts[shape].send_ranks[n],
                my_tag,
                MPI_COMM_WORLD,
                &send_requests[n]
      );

#if !defined(SPEC_MANAGED_MEMORY) && defined(SPEC_OPENMP_TARGET)
      } /* Close the OpenMP target data region */
#endif

    }
    _timeEnd = getTime();
    level->timers.ghostZone_send += (_timeEnd-_timeStart);
  }
  #endif

  // exchange locally... try and hide within Isend latency... 
  if(level->exchange_ghosts[shape].num_blocks[1]){
    _timeStart = getTime();
    if (level->use_offload) {
      device_copy_block(level, id, &level->exchange_ghosts[shape], 1);
    }
    else {
    PRAGMA_THREAD_ACROSS_BLOCKS(level,buffer,level->exchange_ghosts[shape].num_blocks[1])
    for(buffer=0;buffer<level->exchange_ghosts[shape].num_blocks[1];buffer++){
      CopyBlock(level,id,&level->exchange_ghosts[shape].blocks[1][buffer]);
    }
    }
    _timeEnd = getTime();
    level->timers.ghostZone_local += (_timeEnd-_timeStart);
  }


  // wait for MPI to finish...
  #ifdef USE_MPI 
  if(nMessages){
    _timeStart = getTime();
    MPI_Waitall(nMessages,level->exchange_ghosts[shape].requests,level->exchange_ghosts[shape].status);
  #ifdef SYNC_DEVICE_AFTER_WAITALL
    do_sync();
  #endif

#if !defined(SPEC_MANAGED_MEMORY) && defined(SPEC_OPENACC) && !defined(SPEC_ACCEL_AWARE_MPI)
    if (level->use_offload) {
    for(n=0;n<level->exchange_ghosts[shape].num_recvs;n++){
      if (level->exchange_ghosts[shape].recv_sizes[n] > 0) {
	double * recv_buf = level->exchange_ghosts[shape].recv_buffers[n];
	int isp = acc_is_present(recv_buf,level->exchange_ghosts[shape].recv_sizes[n]);
#ifdef DEBUG
	printf("EXCH recv_buf: tag=%d n=%d recv_buf=%p recv_buf_d=%p isp=%d\n",level->tag,n,recv_buf,acc_deviceptr(recv_buf),isp);
#endif
        #pragma acc update device(recv_buf[:level->exchange_ghosts[shape].recv_sizes[n]]) if_present
      }
    }
    }
#endif

#if !defined(SPEC_MANAGED_MEMORY) && defined(SPEC_OPENMP_TARGET) && !defined(SPEC_ACCEL_AWARE_MPI)
    if (level->use_offload && omp_get_num_devices() > 0) {
      for(n=0;n<level->exchange_ghosts[shape].num_recvs;n++){
	if (level->exchange_ghosts[shape].recv_sizes[n] > 0) {
	  double * recv_buf = level->exchange_ghosts[shape].recv_buffers[n];
	  int isp = omp_target_is_present(recv_buf, omp_get_default_device());
	  if (isp != 0) {
# pragma omp target update to(recv_buf[:level->exchange_ghosts[shape].recv_sizes[n]])
	  }
	}
      }
    }
#endif

    _timeEnd = getTime();
    level->timers.ghostZone_wait += (_timeEnd-_timeStart);
  } // end of nMessages


  // unpack MPI receive buffers 
  if(level->exchange_ghosts[shape].num_blocks[2]){
    _timeStart = getTime();
    if(level->use_offload) {
      device_copy_block(level,id,&level->exchange_ghosts[shape],2);
    }
    else {
    PRAGMA_THREAD_ACROSS_BLOCKS(level,buffer,level->exchange_ghosts[shape].num_blocks[2])
    for(buffer=0;buffer<level->exchange_ghosts[shape].num_blocks[2];buffer++){
      CopyBlock(level,id,&level->exchange_ghosts[shape].blocks[2][buffer]);
    }
    }
    _timeEnd = getTime();
    level->timers.ghostZone_unpack += (_timeEnd-_timeStart);
  }
  #endif

 
  level->timers.ghostZone_total += (double)(getTime()-_timeCommunicationStart);
}
