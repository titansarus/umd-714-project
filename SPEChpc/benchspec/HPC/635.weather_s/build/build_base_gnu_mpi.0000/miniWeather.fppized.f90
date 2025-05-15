
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! miniWeather
!! Author: Matt Norman <normanmr@ornl.gov>  , Oak Ridge National Laboratory
!! This code simulates dry, stratified, compressible, non-hydrostatic fluid flows
!! For documentation, please see the attached documentation in the "documentation" folder
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program miniweather
  use mpi
  use specmpitime_mod
  implicit none
  !Declare the precision for the real constants (at least 15 digits of accuracy = double precision)
  integer , parameter :: rp = selected_real_kind(15)
  !Define some physical constants to use throughout the simulation
  real(rp), parameter :: pi        = 3.14159265358979323846264338327_rp   !Pi
  real(rp), parameter :: grav      = 9.8_rp                               !Gravitational acceleration (m / s^2)
  real(rp), parameter :: cp        = 1004._rp                             !Specific heat of dry air at constant pressure
  real(rp), parameter :: cv        = 717._rp                              !Specific heat of dry air at constant volume
  real(rp), parameter :: rd        = 287._rp                              !Dry air constant for equation of state (P=rho*rd*T)
  real(rp), parameter :: p0        = 1.e5_rp                              !Standard pressure at the surface in Pascals
  real(rp), parameter :: C0        = 27.5629410929725921310572974482_rp   !Constant to translate potential temperature into 
									  !  pressure (P=C0*(rho*theta)**gamma)
  real(rp), parameter :: gamma     = 1.40027894002789400278940027894_rp   !gamma=cp/Rd
  !Define domain and stability-related constants
  real(rp) :: xlen !     = 2.e4_rp    !Length of the domain in the x-direction (meters)
  real(rp) :: zlen !     = 1.e4_rp    !Length of the domain in the z-direction (meters)
  real(rp), parameter :: hv_beta   = 0.25_rp     !How strong to diffuse the solution: hv_beta \in [0:1]
  real(rp), parameter :: cfl       = 1.50_rp    !"Courant, Friedrichs, Lewy" number (for numerical stability)
  real(rp), parameter :: max_speed = 450        !Assumed maximum wave speed during the simulation 
	                                        !  (speed of sound + speed of wind) (meter / sec)
  integer , parameter :: hs        = 2          !"Halo" size: number of cells beyond the MPI tasks's domain needed for a full 
						!  "stencil" of information for reconstruction
  integer , parameter :: sten_size = 4          !Size of the stencil used for interpolation

  !Parameters for indexing and flags
  integer , parameter :: NUM_VARS  = 4           !Number of fluid state variables
  integer , parameter :: ID_DENS  = 1           !index for density ("rho")
  integer , parameter :: ID_UMOM  = 2           !index for momentum in the x-direction ("rho * u")
  integer , parameter :: ID_WMOM  = 3           !index for momentum in the z-direction ("rho * w")
  integer , parameter :: ID_RHOT  = 4           !index for density * potential temperature ("rho * theta")
  integer , parameter :: DIR_X = 1              !Integer constant to express that this operation is in the x-direction
  integer , parameter :: DIR_Z = 2              !Integer constant to express that this operation is in the z-direction
  integer , parameter :: DATA__COLLISION       = 1
  integer , parameter :: DATA__THERMAL         = 2
  integer , parameter :: DATA__MOUNTAIN        = 3
  integer , parameter :: DATA__TURBULENCE      = 4
  integer , parameter :: DATA__DENSITY_CURRENT = 5
  integer , parameter :: DATA__INJECTION       = 6

  !Gauss-Legendre quadrature points and weights on the domain [0:1]
  integer , parameter :: nqpoints = 3
  real(rp) :: qpoints (nqpoints) = &
   (/ 0.112701665379258311482073460022E0_rp , 0.500000000000000000000000000000E0_rp , 0.887298334620741688517926539980E0_rp /)
  real(rp) :: qweights(nqpoints) = &
   (/ 0.277777777777777777777777777779E0_rp , 0.444444444444444444444444444444E0_rp , 0.277777777777777777777777777779E0_rp /)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Variables that are initialized but remain static over the coure of the simulation
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real(rp) :: sim_time                          !total simulation time in seconds
  real(rp) :: output_freq                       !frequency to perform output in seconds
  real(rp) :: dt                                !Model time step (seconds)
  integer  :: nx, nz                            !Number of local grid cells in the x- and z- dimensions for this MPI task
  real(rp) :: dx, dz                            !Grid space length in x- and z-dimension (meters)
  integer(8)  :: nx_glob, nz_glob                  !Number of total grid cells in the x- and z- dimensions
  integer(8)  :: i_beg, k_beg                      !beginning index in the x- and z-directions for this MPI task
  integer  :: nranks, myrank                    !Number of MPI ranks and my rank id
  integer  :: left_rank, right_rank             !MPI Rank IDs that exist to my left and right in the global domain
  logical  :: masterproc                        !Am I the master process (rank == 0)?
  real(rp) :: data_spec_int                     !Which data initialization to use
  real(rp), allocatable :: hy_dens_cell      (:)      !hydrostatic density (vert cell avgs).   Dimensions: (1-hs:nz+hs)
  real(rp), allocatable :: hy_dens_theta_cell(:)      !hydrostatic rho*t (vert cell avgs).     Dimensions: (1-hs:nz+hs)
  real(rp), allocatable :: hy_dens_int       (:)      !hydrostatic density (vert cell interf). Dimensions: (1:nz+1)
  real(rp), allocatable :: hy_dens_theta_int (:)      !hydrostatic rho*t (vert cell interf).   Dimensions: (1:nz+1)
  real(rp), allocatable :: hy_pressure_int   (:)      !hydrostatic press (vert cell interf).   Dimensions: (1:nz+1)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Variables that are dynamics over the course of the simulation
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real(rp) :: etime                             !Elapsed model time
  real(rp) :: output_counter                    !Helps determine when it's time to do output
  !Runtime variable arrays
  real(rp), allocatable :: state             (:,:,:)  !Fluid state.              Dimensions: (1-hs:nx+hs,1-hs:nz+hs,NUM_VARS)
  real(rp), allocatable :: state_tmp         (:,:,:)  !Fluid state.              Dimensions: (1-hs:nx+hs,1-hs:nz+hs,NUM_VARS)
  real(rp), allocatable :: flux              (:,:,:)  !Cell interface fluxes.    Dimensions: (nx+1,nz+1,NUM_VARS)
  real(rp), allocatable :: tend              (:,:,:)  !Fluid state tendencies.   Dimensions: (nx,nz,NUM_VARS)
  real(rp), allocatable :: sendbuf_l(:,:,:)     !buffers for MPI data exchanges. Dimensions: (hs,nz,NUM_VARS)
  real(rp), allocatable :: sendbuf_r(:,:,:)     !buffers for MPI data exchanges. Dimensions: (hs,nz,NUM_VARS)
  real(rp), allocatable :: recvbuf_l(:,:,:)     !buffers for MPI data exchanges. Dimensions: (hs,nz,NUM_VARS)
  real(rp), allocatable :: recvbuf_r(:,:,:)     !buffers for MPI data exchanges. Dimensions: (hs,nz,NUM_VARS)
  integer(8) :: t1, t2, rate                    !For CPU Timings
  integer :: argcnt, istat
  character(len=32) :: arg, outfile

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! THE MAIN PROGRAM STARTS HERE
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! BEGIN USER-CONFIGURABLE PARAMETERS
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !The x-direction length is twice as long as the z-direction length
  !So, you'll want to have nx_glob be twice as large as nz_glob
  nx_glob = 400      !Number of total cells in the x-dirction
  nz_glob = 200      !Number of total cells in the z-dirction
  sim_time = 1500      !How many seconds to run the simulation
  output_freq = 10   !How frequently to output data to file (in seconds)
  !Model setup: DATA__THERMAL or DATA__COLLISION
  data_spec_int = DATA__INJECTION
  xlen = 2.e4_rp    !Length of the domain in the x-direction (meters)
  zlen = 1.e4_rp    !Length of the domain in the z-direction (meters)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! END USER-CONFIGURABLE PARAMETERS
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  argcnt = command_argument_count()
  if (argcnt .ne. 8) then 
     print *, "Useage: "
     print *, "  a.out OufFile NX NZ Time Freq Model"
     print *, "  OutFile: Name of the output file"
     print *, "       NX: Number of total cells in the x-dirction "
     print *, "       NZ: Number of total cells in the z-dirction "
     print *, "     Time: How many seconds to run the simulation "
     print *, "     Freq: How frequently to output data to file (in seconds)"
     print *, "    Model:  1 = COLLISION "
     print *, "            2 = THERMAL "
     print *, "            3 = MOUNTAIN "
     print *, "            4 = TURBULENCE "
     print *, "            5 = DENSITY_CURRENT "
     print *, "            6 = INJECTION "
     stop 1
   endif
  
   call get_command_argument(1, outfile)
   call get_command_argument(2, arg)
   read(arg,*,iostat=istat) xlen 
   if (istat.ne.0) then
        print *, "Error: Unable to read xlen";
        stop 1
   endif 
   call get_command_argument(3, arg)
   read(arg,*,iostat=istat) zlen 
   if (istat.ne.0) then
        print *, "Error: Unable to read zlen";
        stop 1
   endif 
   call get_command_argument(4, arg)
   read(arg,*,iostat=istat) nx_glob 
   if (istat.ne.0) then
        print *, "Error: Unable to read NX";
        stop 1
   endif 
   call get_command_argument(5, arg)
   read(arg,*,iostat=istat) nz_glob 
   if (istat.ne.0) then
        print *, "Error: Unable to read NZ";
        stop 1
   endif 
   call get_command_argument(6, arg)
   read(arg,*,iostat=istat) sim_time 
   if (istat.ne.0) then
        print *, "Error: Unable to read Time";
        stop 1
   endif 
   call get_command_argument(7, arg)
   read(arg,*,iostat=istat) output_freq  
   if (istat.ne.0) then
        print *, "Error: Unable to read Freq";
        stop 1
   endif 
   call get_command_argument(8, arg)
   read(arg,*,iostat=istat) data_spec_int
   if (istat.ne.0) then
        print *, "Error: Unable to read Model";
        stop 1
   endif 
   if (data_spec_int.lt.1.or.data_spec_int.gt.6) then
        print *, "Error: Invalid Model ", data_spec_int
        stop 1
   endif 

  !Initialize MPI, allocate arrays, initialize the grid and the data
  call init()
  call spectime_start(0)
  call spectime_start(1)


  !Output the initial state
  call output(state,etime,outfile)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! MAIN TIME STEP LOOP
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call spectime_stop(1)
  call spectime_start(3)
  do while (etime < sim_time)
    !If the time step leads to exceeding the simulation time, shorten it for the last step
    if (etime + dt > sim_time) dt = sim_time - etime
    !Perform a single time step
    call perform_timestep(state,state_tmp,flux,tend,dt)
    !Inform the user
    !Update the elapsed time and output counter
    etime = etime + dt
    output_counter = output_counter + dt
    !If it's time for output, reset the counter, and do output
    if (output_counter >= output_freq) then
      output_counter = output_counter - output_freq
      call output(state,etime,outfile)
    endif
  enddo
  call spectime_stop(3)
  if (masterproc) then
  endif
  call spectime_stop(0)
  call spectime_final(.true.,0.0_8)
  !Deallocate and finialize MPI
  call finalize()


contains


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! HELPER SUBROUTINES AND FUNCTIONS GO HERE
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  !Performs a single dimensionally split time step using a simple low-storate three-stage Runge-Kutta time integrator
  !The dimensional splitting is a second-order-accurate alternating Strang splitting in which the
  !order of directions is alternated each time step.
  !The Runge-Kutta method used here is defined as follows:
  ! q*     = q[n] + dt/3 * rhs(q[n])
  ! q**    = q[n] + dt/2 * rhs(q*  )
  ! q[n+1] = q[n] + dt/1 * rhs(q** )
  subroutine perform_timestep(state,state_tmp,flux,tend,dt)
    implicit none
    real(rp), intent(inout) :: state    (1-hs:nx+hs,1-hs:nz+hs,NUM_VARS)
    real(rp), intent(  out) :: state_tmp(1-hs:nx+hs,1-hs:nz+hs,NUM_VARS)
    real(rp), intent(  out) :: flux     (nx+1,nz+1,NUM_VARS)
    real(rp), intent(  out) :: tend     (nx,nz,NUM_VARS)
    real(rp), intent(in   ) :: dt
    logical, save :: direction_switch = .true.
    if (direction_switch) then
      !x-direction first
      call semi_discrete_step( state , state     , state_tmp , dt / 3 , DIR_X , flux , tend )
      call semi_discrete_step( state , state_tmp , state_tmp , dt / 2 , DIR_X , flux , tend )
      call semi_discrete_step( state , state_tmp , state     , dt / 1 , DIR_X , flux , tend )
      !z-direction second
      call semi_discrete_step( state , state     , state_tmp , dt / 3 , DIR_Z , flux , tend )
      call semi_discrete_step( state , state_tmp , state_tmp , dt / 2 , DIR_Z , flux , tend )
      call semi_discrete_step( state , state_tmp , state     , dt / 1 , DIR_Z , flux , tend )
    else
      !z-direction second
      call semi_discrete_step( state , state     , state_tmp , dt / 3 , DIR_Z , flux , tend )
      call semi_discrete_step( state , state_tmp , state_tmp , dt / 2 , DIR_Z , flux , tend )
      call semi_discrete_step( state , state_tmp , state     , dt / 1 , DIR_Z , flux , tend )
      !x-direction first
      call semi_discrete_step( state , state     , state_tmp , dt / 3 , DIR_X , flux , tend )
      call semi_discrete_step( state , state_tmp , state_tmp , dt / 2 , DIR_X , flux , tend )
      call semi_discrete_step( state , state_tmp , state     , dt / 1 , DIR_X , flux , tend )
    endif
  end subroutine perform_timestep


  !Perform a single semi-discretized step in time with the form:
  !state_out = state_init + dt * rhs(state_forcing)
  !Meaning the step starts from state_init, computes the rhs using state_forcing, and stores the result in state_out
  subroutine semi_discrete_step( state_init , state_forcing , state_out , dt , dir , flux , tend )
    implicit none
    real(rp), intent(in   ) :: state_init   (1-hs:nx+hs,1-hs:nz+hs,NUM_VARS)
    real(rp), intent(inout) :: state_forcing(1-hs:nx+hs,1-hs:nz+hs,NUM_VARS)
    real(rp), intent(  out) :: state_out    (1-hs:nx+hs,1-hs:nz+hs,NUM_VARS)
    real(rp), intent(  out) :: flux         (nx+1,nz+1,NUM_VARS)
    real(rp), intent(  out) :: tend         (nx,nz,NUM_VARS)
    real(rp), intent(in   ) :: dt
    integer , intent(in   ) :: dir
    integer :: i,k,ll

    if     (dir == DIR_X) then
      !Set the halo values for this MPI task's fluid state in the x-direction
      call set_halo_values_x(state_forcing)
      !Compute the time tendencies for the fluid state in the x-direction
      call compute_tendencies_x(state_forcing,flux,tend)
    elseif (dir == DIR_Z) then
      !Set the halo values for this MPI task's fluid state in the z-direction
      call set_halo_values_z(state_forcing)
      !Compute the time tendencies for the fluid state in the z-direction
      call compute_tendencies_z(state_forcing,flux,tend)
    endif

    !Apply the tendencies to the fluid state
    do ll = 1 , NUM_VARS
      do k = 1 , nz
        do i = 1 , nx
          state_out(i,k,ll) = state_init(i,k,ll) + dt * tend(i,k,ll)
        enddo
      enddo
    enddo
  end subroutine semi_discrete_step


  !Compute the time tendencies of the fluid state using forcing in the x-direction
  !Since the halos are set in a separate routine, this will not require MPI
  !First, compute the flux vector at each cell interface in the x-direction (including hyperviscosity)
  !Then, compute the tendencies using those fluxes
  subroutine compute_tendencies_x(state,flux,tend)
    implicit none
    real(rp), intent(in   ) :: state(1-hs:nx+hs,1-hs:nz+hs,NUM_VARS)
    real(rp), intent(  out) :: flux (nx+1,nz+1,NUM_VARS)
    real(rp), intent(  out) :: tend (nx,nz,NUM_VARS)
    integer :: i,k,ll,s
    real(rp) :: r,u,w,t,p, stencil(4), d3_vals(NUM_VARS), vals(NUM_VARS), hv_coef
    !Compute the hyperviscosity coeficient
    hv_coef = -hv_beta * dx / (16*dt)
    !Compute fluxes in the x-direction for each cell
    do k = 1 , nz
      do i = 1 , nx+1
        !Use fourth-order interpolation from four cell averages to compute the value at the interface in question
        do ll = 1 , NUM_VARS
          do s = 1 , sten_size
            stencil(s) = state(i-hs-1+s,k,ll)
          enddo
          !Fourth-order-accurate interpolation of the state
          vals(ll) = -stencil(1)/12 + 7*stencil(2)/12 + 7*stencil(3)/12 - stencil(4)/12
          !First-order-accurate interpolation of the third spatial derivative of the state (for artificial viscosity)
          d3_vals(ll) = -stencil(1) + 3*stencil(2) - 3*stencil(3) + stencil(4)
        enddo

        !Compute density, u-wind, w-wind, potential temperature, and pressure (r,u,w,t,p respectively)
        r = vals(ID_DENS) + hy_dens_cell(k)
        u = vals(ID_UMOM) / r
        w = vals(ID_WMOM) / r
        t = ( vals(ID_RHOT) + hy_dens_theta_cell(k) ) / r
        p = C0*(r*t)**gamma

        !Compute the flux vector
        flux(i,k,ID_DENS) = r*u     - hv_coef*d3_vals(ID_DENS)
        flux(i,k,ID_UMOM) = r*u*u+p - hv_coef*d3_vals(ID_UMOM)
        flux(i,k,ID_WMOM) = r*u*w   - hv_coef*d3_vals(ID_WMOM)
        flux(i,k,ID_RHOT) = r*u*t   - hv_coef*d3_vals(ID_RHOT)
      enddo
    enddo
    !Use the fluxes to compute tendencies for each cell
    do ll = 1 , NUM_VARS
      do k = 1 , nz
        do i = 1 , nx
          tend(i,k,ll) = -( flux(i+1,k,ll) - flux(i,k,ll) ) / dx
        enddo
      enddo
    enddo
  end subroutine compute_tendencies_x


  !Compute the time tendencies of the fluid state using forcing in the z-direction
  !Since the halos are set in a separate routine, this will not require MPI
  !First, compute the flux vector at each cell interface in the z-direction (including hyperviscosity)
  !Then, compute the tendencies using those fluxes
  subroutine compute_tendencies_z(state,flux,tend)
    implicit none
    real(rp), intent(in   ) :: state(1-hs:nx+hs,1-hs:nz+hs,NUM_VARS)
    real(rp), intent(  out) :: flux (nx+1,nz+1,NUM_VARS)
    real(rp), intent(  out) :: tend (nx,nz,NUM_VARS)
    integer :: i,k,ll,s
    real(rp) :: r,u,w,t,p, stencil(4), d3_vals(NUM_VARS), vals(NUM_VARS), hv_coef
    !Compute the hyperviscosity coeficient
    hv_coef = -hv_beta * dx / (16*dt)
    !Compute fluxes in the x-direction for each cell
    do k = 1 , nz+1
      do i = 1 , nx
        !Use fourth-order interpolation from four cell averages to compute the value at the interface in question
        do ll = 1 , NUM_VARS
          do s = 1 , sten_size
            stencil(s) = state(i,k-hs-1+s,ll)
          enddo
          !Fourth-order-accurate interpolation of the state
          vals(ll) = -stencil(1)/12 + 7*stencil(2)/12 + 7*stencil(3)/12 - stencil(4)/12
          !First-order-accurate interpolation of the third spatial derivative of the state
          d3_vals(ll) = -stencil(1) + 3*stencil(2) - 3*stencil(3) + stencil(4)
        enddo

        !Compute density, u-wind, w-wind, potential temperature, and pressure (r,u,w,t,p respectively)
        r = vals(ID_DENS) + hy_dens_int(k)
        u = vals(ID_UMOM) / r
        w = vals(ID_WMOM) / r
        t = ( vals(ID_RHOT) + hy_dens_theta_int(k) ) / r
        p = C0*(r*t)**gamma - hy_pressure_int(k)

        !Compute the flux vector with hyperviscosity
        flux(i,k,ID_DENS) = r*w     - hv_coef*d3_vals(ID_DENS)
        flux(i,k,ID_UMOM) = r*w*u   - hv_coef*d3_vals(ID_UMOM)
        flux(i,k,ID_WMOM) = r*w*w+p - hv_coef*d3_vals(ID_WMOM)
        flux(i,k,ID_RHOT) = r*w*t   - hv_coef*d3_vals(ID_RHOT)
      enddo
    enddo
    !Use the fluxes to compute tendencies for each cell
    do ll = 1 , NUM_VARS
      do k = 1 , nz
        do i = 1 , nx
          tend(i,k,ll) = -( flux(i,k+1,ll) - flux(i,k,ll) ) / dz
          if (ll == ID_WMOM) then
            tend(i,k,ID_WMOM) = tend(i,k,ID_WMOM) - state(i,k,ID_DENS)*grav
          endif
        enddo
      enddo
    enddo
  end subroutine compute_tendencies_z


  !Set this MPI task's halo values in the x-direction. This routine will require MPI
  subroutine set_halo_values_x(state)
    implicit none
    real(rp), intent(inout) :: state(1-hs:nx+hs,1-hs:nz+hs,NUM_VARS)
    integer :: k, ll, s, ierr, req_r(2), req_s(2), status(MPI_STATUS_SIZE,2)
    real(rp) :: z

    !Prepost receives
    call mpi_irecv(recvbuf_l,hs*nz*NUM_VARS,MPI_REAL8, left_rank,0,MPI_COMM_WORLD,req_r(1),ierr)
    call mpi_irecv(recvbuf_r,hs*nz*NUM_VARS,MPI_REAL8,right_rank,1,MPI_COMM_WORLD,req_r(2),ierr)

    !Pack the send buffers
    do ll = 1 , NUM_VARS
      do k = 1 , nz
        do s = 1 , hs
          sendbuf_l(s,k,ll) = state(s      ,k,ll)
          sendbuf_r(s,k,ll) = state(nx-hs+s,k,ll)
        enddo
      enddo
    enddo

    !Fire off the sends
    call mpi_isend(sendbuf_l,hs*nz*NUM_VARS,MPI_REAL8, left_rank,1,MPI_COMM_WORLD,req_s(1),ierr)
    call mpi_isend(sendbuf_r,hs*nz*NUM_VARS,MPI_REAL8,right_rank,0,MPI_COMM_WORLD,req_s(2),ierr)

    !Wait for receives to finish
    call mpi_waitall(2,req_r,status,ierr)



    !Unpack the receive buffers
    do ll = 1 , NUM_VARS
      do k = 1 , nz
        do s = 1 , hs
          state(-hs+s,k,ll) = recvbuf_l(s,k,ll)
          state(nx+s ,k,ll) = recvbuf_r(s,k,ll)
        enddo
      enddo
    enddo

    !Wait for sends to finish
    call mpi_waitall(2,req_s,status,ierr)

    if (data_spec_int == DATA__INJECTION) then
      if (myrank == 0) then
        do k = 1 , nz
          z = (k_beg-1 + k-0.5_rp)*dz
          if (abs(z-3*zlen/4) <= zlen/16) then
            state(-1:0,k,ID_UMOM) = (state(-1:0,k,ID_DENS)+hy_dens_cell(k)) * 50._rp
            state(-1:0,k,ID_RHOT) = (state(-1:0,k,ID_DENS)+hy_dens_cell(k)) * 298._rp - hy_dens_theta_cell(k)
          endif
        enddo
      endif
    endif
  end subroutine set_halo_values_x


  !Set this MPI task's halo values in the z-direction. This does not require MPI because there is no MPI
  !decomposition in the vertical direction
  subroutine set_halo_values_z(state)
    implicit none
    real(rp), intent(inout) :: state(1-hs:nx+hs,1-hs:nz+hs,NUM_VARS)
    integer :: i, ll
    real(rp) :: mnt_width 
    real(rp) :: x, xloc, mnt_deriv
    mnt_width = xlen/8
   do ll = 1 , NUM_VARS
      do i = 1-hs,nx+hs
        if (ll == ID_WMOM) then
          state(i,-1  ,ll) = 0
          state(i,0   ,ll) = 0
          state(i,nz+1,ll) = 0
          state(i,nz+2,ll) = 0
          !Impose the vertical momentum effects of an artificial cos^2 mountain at the lower boundary
          if (data_spec_int == DATA__MOUNTAIN) then
            x = (i_beg-1+i-0.5_rp)*dx
            if ( abs(x-xlen/4) < mnt_width ) then
              xloc = (x-(xlen/4)) / mnt_width
              !Compute the derivative of the fake mountain
              mnt_deriv = -pi*cos(pi*xloc/2)*sin(pi*xloc/2)*10/dx
              !w = (dz/dx)*u
              state(i,-1,ID_WMOM) = mnt_deriv*state(i,1,ID_UMOM)
              state(i,0 ,ID_WMOM) = mnt_deriv*state(i,1,ID_UMOM)
            endif
          endif
        else
          state(i,-1  ,ll) = state(i,1 ,ll)
          state(i,0   ,ll) = state(i,1 ,ll)
          state(i,nz+1,ll) = state(i,nz,ll)
          state(i,nz+2,ll) = state(i,nz,ll)
        endif
      enddo
    enddo
  end subroutine set_halo_values_z


  !Initialize some grid settings, initialize MPI, allocate data, and initialize the full grid and data
  subroutine init()
    implicit none
    integer :: i, k, ii, kk, ll, ierr, i_end
    real(rp) :: x, z, r, u, w, t, hr, ht, nper

    !Set the cell grid size
    dx = xlen / nx_glob
    dz = zlen / nz_glob

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! BEGIN MPI DUMMY SECTION
    !! TODO: (1) GET NUMBER OF MPI RANKS
    !!       (2) GET MY MPI RANK ID (RANKS ARE ZERO-BASED INDEX)
    !!       (3) COMPUTE MY BEGINNING "I" INDEX (1-based index)
    !!       (4) COMPUTE HOW MANY X-DIRECTION CELLS MY RANK HAS
    !!       (5) FIND MY LEFT AND RIGHT NEIGHBORING RANK IDs
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call mpi_init(ierr)
    call mpi_comm_size(mpi_comm_world,nranks,ierr)
    call mpi_comm_rank(mpi_comm_world,myrank,ierr)
    nper = real(nx_glob)/nranks
    i_beg = nint( nper* (myrank)    )+1
    i_end = nint( nper*((myrank)+1) )
    nx = i_end - i_beg + 1
    left_rank  = myrank - 1
    if (left_rank == -1) left_rank = nranks-1
    right_rank = myrank + 1
    if (right_rank == nranks) right_rank = 0
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! END MPI DUMMY SECTION
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



    !Vertical direction isn't MPI-ized, so the rank's local values = the global values
    k_beg = 1
    nz = nz_glob
    masterproc = (myrank == 0)

    !Allocate the model data
    allocate(state             (1-hs:nx+hs,1-hs:nz+hs,NUM_VARS))
    allocate(state_tmp         (1-hs:nx+hs,1-hs:nz+hs,NUM_VARS))
    allocate(flux              (nx+1,nz+1,NUM_VARS))
    allocate(tend              (nx,nz,NUM_VARS))
    allocate(hy_dens_cell      (1-hs:nz+hs))
    allocate(hy_dens_theta_cell(1-hs:nz+hs))
    allocate(hy_dens_int       (nz+1))
    allocate(hy_dens_theta_int (nz+1))
    allocate(hy_pressure_int   (nz+1))
    allocate(sendbuf_l(hs,nz,NUM_VARS))
    allocate(sendbuf_r(hs,nz,NUM_VARS))
    allocate(recvbuf_l(hs,nz,NUM_VARS))
    allocate(recvbuf_r(hs,nz,NUM_VARS))

    !Define the maximum stable time step based on an assumed maximum wind speed
    dt = min(dx,dz) / max_speed * cfl
    !Set initial elapsed model time and output_counter to zero
    etime = 0
    output_counter = 0

    !If I'm the master process in MPI, display some grid information
    if (masterproc) then
      write(*,*) 'nx_glob, nz_glob:', nx_glob, nz_glob
      write(*,*) 'dx,dz: ',dx,dz
      write(*,*) 'dt: ',dt
    endif
    !Want to make sure this info is displayed before further output
    call MPI_Barrier(MPI_COMM_WORLD,ierr)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Initialize the cell-averaged fluid state via Gauss-Legendre quadrature
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    state = 0
    do k = 1-hs , nz+hs
      do i = 1-hs , nx+hs
        !Initialize the state to zero
        !Use Gauss-Legendre quadrature to initialize a hydrostatic balance + temperature perturbation
        do kk = 1 , nqpoints
          do ii = 1 , nqpoints
            !Compute the x,z location within the global domain based on cell and quadrature index
            x = (i_beg-1 + i-0.5_rp) * dx + (qpoints(ii)-0.5_rp)*dx
            z = (k_beg-1 + k-0.5_rp) * dz + (qpoints(kk)-0.5_rp)*dz

            !Set the fluid state based on the user's specification
            if (data_spec_int == DATA__COLLISION      ) call collision      (x,z,r,u,w,t,hr,ht)
            if (data_spec_int == DATA__THERMAL        ) call thermal        (x,z,r,u,w,t,hr,ht)
            if (data_spec_int == DATA__MOUNTAIN       ) call mountain_waves (x,z,r,u,w,t,hr,ht)
            if (data_spec_int == DATA__TURBULENCE     ) call turbulence     (x,z,r,u,w,t,hr,ht)
            if (data_spec_int == DATA__DENSITY_CURRENT) call density_current(x,z,r,u,w,t,hr,ht)
            if (data_spec_int == DATA__INJECTION      ) call injection      (x,z,r,u,w,t,hr,ht)

            !Store into the fluid state array
            state(i,k,ID_DENS) = state(i,k,ID_DENS) + r                         * qweights(ii)*qweights(kk)
            state(i,k,ID_UMOM) = state(i,k,ID_UMOM) + (r+hr)*u                  * qweights(ii)*qweights(kk)
            state(i,k,ID_WMOM) = state(i,k,ID_WMOM) + (r+hr)*w                  * qweights(ii)*qweights(kk)
            state(i,k,ID_RHOT) = state(i,k,ID_RHOT) + ( (r+hr)*(t+ht) - hr*ht ) * qweights(ii)*qweights(kk)
          enddo
        enddo
        do ll = 1 , NUM_VARS
          state_tmp(i,k,ll) = state(i,k,ll)
        enddo
      enddo
    enddo
    !Compute the hydrostatic background state over vertical cell averages
    hy_dens_cell = 0
    hy_dens_theta_cell = 0
    do k = 1-hs , nz+hs
      do kk = 1 , nqpoints
        z = (k_beg-1 + k-0.5_rp) * dz + (qpoints(kk)-0.5_rp)*dz
        !Set the fluid state based on the user's specification
        if (data_spec_int == DATA__COLLISION      ) call collision      (0._rp,z,r,u,w,t,hr,ht)
        if (data_spec_int == DATA__THERMAL        ) call thermal        (0._rp,z,r,u,w,t,hr,ht)
        if (data_spec_int == DATA__MOUNTAIN       ) call mountain_waves (0._rp,z,r,u,w,t,hr,ht)
        if (data_spec_int == DATA__TURBULENCE     ) call turbulence     (0._rp,z,r,u,w,t,hr,ht)
        if (data_spec_int == DATA__DENSITY_CURRENT) call density_current(0._rp,z,r,u,w,t,hr,ht)
        if (data_spec_int == DATA__INJECTION      ) call injection      (0._rp,z,r,u,w,t,hr,ht)
        hy_dens_cell(k)       = hy_dens_cell(k)       + hr    * qweights(kk)
        hy_dens_theta_cell(k) = hy_dens_theta_cell(k) + hr*ht * qweights(kk)
      enddo
    enddo
    !Compute the hydrostatic background state at vertical cell interfaces
    do k = 1 , nz+1
      z = (k_beg-1 + k-1) * dz
      if (data_spec_int == DATA__COLLISION      ) call collision      (0._rp,z,r,u,w,t,hr,ht)
      if (data_spec_int == DATA__THERMAL        ) call thermal        (0._rp,z,r,u,w,t,hr,ht)
      if (data_spec_int == DATA__MOUNTAIN       ) call mountain_waves (0._rp,z,r,u,w,t,hr,ht)
      if (data_spec_int == DATA__TURBULENCE     ) call turbulence     (0._rp,z,r,u,w,t,hr,ht)
      if (data_spec_int == DATA__DENSITY_CURRENT) call density_current(0._rp,z,r,u,w,t,hr,ht)
      if (data_spec_int == DATA__INJECTION      ) call injection      (0._rp,z,r,u,w,t,hr,ht)
      hy_dens_int      (k) = hr
      hy_dens_theta_int(k) = hr*ht
      hy_pressure_int  (k) = C0*(hr*ht)**gamma
    enddo
  end subroutine init


  subroutine injection(x,z,r,u,w,t,hr,ht)
    implicit none
    real(rp), intent(in   ) :: x, z        !x- and z- location of the point being sampled
    real(rp), intent(  out) :: r, u, w, t  !Density, uwind, wwind, and potential temperature
    real(rp), intent(  out) :: hr, ht      !Hydrostatic density and potential temperature
    call hydro_const_theta(z,hr,ht)
    r = 0
    t = 0
    u = 0
    w = 0
  end subroutine injection


  subroutine density_current(x,z,r,u,w,t,hr,ht)
    implicit none
    real(rp), intent(in   ) :: x, z        !x- and z- location of the point being sampled
    real(rp), intent(  out) :: r, u, w, t  !Density, uwind, wwind, and potential temperature
    real(rp), intent(  out) :: hr, ht      !Hydrostatic density and potential temperature
    call hydro_const_theta(z,hr,ht)
    r = 0
    t = 0
    u = 0
    w = 0
    t = t + sample_ellipse_cosine(x,z,-20._rp ,xlen/2,5000._rp,4000._rp,2000._rp)
  end subroutine density_current


  subroutine turbulence(x,z,r,u,w,t,hr,ht)
    implicit none
    real(rp), intent(in   ) :: x, z        !x- and z- location of the point being sampled
    real(rp), intent(  out) :: r, u, w, t  !Density, uwind, wwind, and potential temperature
    real(rp), intent(  out) :: hr, ht      !Hydrostatic density and potential temperature
    call hydro_const_theta(z,hr,ht)
    r = 0
    t = 0
    call random_number(u)
    call random_number(w)
    u = (u-0.5)*20
    w = (w-0.5)*20
  end subroutine turbulence


  subroutine mountain_waves(x,z,r,u,w,t,hr,ht)
    implicit none
    real(rp), intent(in   ) :: x, z        !x- and z- location of the point being sampled
    real(rp), intent(  out) :: r, u, w, t  !Density, uwind, wwind, and potential temperature
    real(rp), intent(  out) :: hr, ht      !Hydrostatic density and potential temperature
    call hydro_const_bvfreq(z,0.02_rp,hr,ht)
    r = 0
    t = 0
    u = 15
    w = 0
  end subroutine mountain_waves


  !Rising thermal
  subroutine thermal(x,z,r,u,w,t,hr,ht)
    implicit none
    real(rp), intent(in   ) :: x, z        !x- and z- location of the point being sampled
    real(rp), intent(  out) :: r, u, w, t  !Density, uwind, wwind, and potential temperature
    real(rp), intent(  out) :: hr, ht      !Hydrostatic density and potential temperature
    call hydro_const_theta(z,hr,ht)
    r = 0
    t = 0
    u = 0
    w = 0
    t = t + sample_ellipse_cosine(x,z, 3._rp ,xlen/2,2000._rp,2000._rp,2000._rp)
  end subroutine thermal


  !Colliding thermals
  subroutine collision(x,z,r,u,w,t,hr,ht)
    implicit none
    real(rp), intent(in   ) :: x, z        !x- and z- location of the point being sampled
    real(rp), intent(  out) :: r, u, w, t  !Density, uwind, wwind, and potential temperature
    real(rp), intent(  out) :: hr, ht      !Hydrostatic density and potential temperature
    call hydro_const_theta(z,hr,ht)
    r = 0
    t = 0
    u = 0
    w = 0
    t = t + sample_ellipse_cosine(x,z, 20._rp,xlen/2,2000._rp,2000._rp,2000._rp)
    t = t + sample_ellipse_cosine(x,z,-20._rp,xlen/2,8000._rp,2000._rp,2000._rp)
  end subroutine collision


  !Establish hydrstatic balance using constant potential temperature (thermally neutral atmosphere)
  subroutine hydro_const_theta(z,r,t)
    implicit none
    real(rp), intent(in   ) :: z  !x- and z- location of the point being sampled
    real(rp), intent(  out) :: r, t  !Density and potential temperature at this point
    real(rp), parameter :: theta0 = 300._rp  !Background potential temperature
    real(rp), parameter :: exner0 = 1._rp    !Surface-level Exner pressure
    real(rp) :: p,exner,rt
    !Establish hydrostatic balance first using Exner pressure
    t = theta0                                  !Potential Temperature at z
    exner = exner0 - grav * z / (cp * theta0)   !Exner pressure at z
    p = p0 * exner**(cp/rd)                     !Pressure at z
    rt = (p / c0)**(1._rp / gamma)              !rho*theta at z
    r = rt / t                                  !Density at z
  end subroutine hydro_const_theta


  !Establish hydrstatic balance using constant Brunt-Vaisala frequency
  subroutine hydro_const_bvfreq( z , bv_freq0 , r , t )
    implicit none
    real(rp), intent(in   ) :: z , bv_freq0
    real(rp), intent(  out) :: r , t
    real(rp) :: theta0 = 300, exner0 = 1
    real(rp) :: p, exner, rt
    t = theta0 * exp( bv_freq0**2 / grav * z )                                  !Pot temp at z
    exner = exner0 - grav**2 / (cp * bv_freq0**2) * (t - theta0) / (t * theta0) !Exner pressure at z
    p = p0 * exner**(cp/rd)                                                     !Pressure at z
    rt = (p / c0)**(1._rp / gamma)                                              !rho*theta at z
    r = rt / t                                                                  !Density at z
  end subroutine hydro_const_bvfreq


  !Sample from an ellipse of a specified center, radius, and amplitude at a specified location
  function sample_ellipse_cosine( x , z , amp , x0 , z0 , xrad , zrad )   result(val)
    implicit none
    real(rp), intent(in) :: x , z         !Location to sample
    real(rp), intent(in) :: amp           !Amplitude of the bubble
    real(rp), intent(in) :: x0 , z0       !Center of the bubble
    real(rp), intent(in) :: xrad , zrad   !Radius of the bubble
    real(rp)             :: val           !The output sampled value
    real(rp) :: dist
    !Compute distance from bubble center
    dist = sqrt( ((x-x0)/xrad)**2 + ((z-z0)/zrad)**2 ) * pi / 2._rp
    !If the distance from bubble center is less than the radius, create a cos**2 profile
    if (dist <= pi / 2._rp) then
      val = amp * cos(dist)**2
    else
      val = 0._rp
    endif
  end function sample_ellipse_cosine


  !Deallocate and call MPI_Finalize
  subroutine finalize()
    implicit none
    integer :: ierr
    deallocate(state             )
    deallocate(state_tmp         )
    deallocate(flux              )
    deallocate(tend              )
    deallocate(hy_dens_cell      )
    deallocate(hy_dens_theta_cell)
    deallocate(hy_dens_int       )
    deallocate(hy_dens_theta_int )
    deallocate(hy_pressure_int   )
    deallocate(sendbuf_l)
    deallocate(sendbuf_r)
    deallocate(recvbuf_l)
    deallocate(recvbuf_r)
    call MPI_Finalize(ierr)
  end subroutine finalize


  !Output the fluid state (state) to a NetCDF file at a given elapsed model time (etime)
  !The file I/O uses parallel-netcdf, the only external library required for this mini-app.
  !If it's too cumbersome, you can comment the I/O out, but you'll miss out on some potentially cool graphics
  subroutine output(state,etime,outfile)
  use mpi
  implicit none
    real(rp), intent(inout) :: state(1-hs:nx+hs,1-hs:nz+hs,NUM_VARS)
    real(rp), intent(in) :: etime
    character(len=32), intent(in) :: outfile
    integer :: ncid, t_dimid, x_dimid, z_dimid, dens_varid, uwnd_varid, wwnd_varid, theta_varid, t_varid
    integer :: i,k,ierr
    integer, save :: num_out = 0
    integer(kind=MPI_OFFSET_KIND) :: len, st1(1),ct1(1),st3(3),ct3(3)
    !Temporary arrays to hold density, u-wind, w-wind, and potential temperature (theta)
    real(rp) :: tmp(4)
    real(rp) :: etimearr(1)


    if (etime .ne. 0 .and. masterproc) then
       open(UNIT=15,FILE=trim(outfile))
       write(15,'(A40)') "-------------- DENS --------------"
       do k = 1 , nz
         do i = 1,4
           tmp(i)= state(i,k,ID_DENS)
         enddo
         write(15,'(I8,A1,F15.10,A1,F15.10,A1,F15.10,A1,F15.10)') k,' ',tmp(1),' ',tmp(2),' ',tmp(3),' ',tmp(4)
       enddo
       write(15,'(A40)') "-------------- UWND --------------"
       do k = 1 , nz
         do i = 1,4
           tmp(i) = state(i,k,ID_UMOM) / ( hy_dens_cell(k) + state(i,k,ID_DENS) )
         enddo
         write(15,'(I8,A1,F15.10,A1,F15.10,A1,F15.10,A1,F15.10)') k,' ',tmp(1),' ',tmp(2),' ',tmp(3),' ',tmp(4)
       enddo
       write(15,'(A40)') "-------------- WWND --------------"
       do k = 1 , nz
         do i = 1,4
           tmp(i) = state(i,k,ID_WMOM) / ( hy_dens_cell(k) + state(i,k,ID_DENS) )
         enddo
         write(15,'(I8,A1,F15.10,A1,F15.10,A1,F15.10,A1,F15.10)') k,' ',tmp(1),' ',tmp(2),' ',tmp(3),' ',tmp(4)
       enddo
       write(15,'(A40)') "-------------- THETA --------------"
       do k = 1 , nz
         do i = 1,4
           tmp(i) = ( state(i,k,ID_RHOT) + hy_dens_theta_cell(k) ) / ( hy_dens_cell(k) + state(i,k,ID_DENS) ) &
                    - hy_dens_theta_cell(k) / hy_dens_cell(k)
         enddo
         write(15,'(I8,A1,F15.10,A1,F15.10,A1,F15.10,A1,F15.10)') k,' ',tmp(1),' ',tmp(2),' ',tmp(3),' ',tmp(4)
       enddo
       close(15)
    endif


  end subroutine output


end program miniweather
