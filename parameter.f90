module para
    use mpi
    implicit none
    
    real, parameter :: kb=0.08617 ! Boltzmann parameter, in meV/K
    integer, parameter :: stdout=100

    integer :: num_atom, num_j
    integer :: num_iterr
    integer, allocatable :: num_nn(:)
    integer :: maxnn
    logical :: afm
    integer, allocatable :: sublat(:,:) ! spin index of each sublattice.
    real :: S ! magnetic momentum, in \muB
    real, allocatable :: spin(:,:)
    integer :: num_ensemble ! ensemble number
    integer :: nem ! ensemble number for short
    real :: latvec(2,3)
    logical :: periodic ! periodic boundary condition
    logical :: randominit ! initialization method, default .true.
    real :: randbias ! Random initialization bias when randominit=.true., 
                     ! Default should be 0.5 but here we set to 0.6 to let converge fast in low temperature.
                     ! When equals to 1, init the system to FM.
    real, allocatable :: atomsite(:,:)
    real, allocatable :: j(:)
    integer, allocatable :: nn(:,:,:)
    real :: temperature_btm, temperature_top ! temperature
    integer :: temperature_num, tnn, nit
    integer :: offset
    integer :: outfile
    integer :: stat
    real, allocatable :: d(:,:) !shape: (num_j, 2), in this d range, set
!corresponding nearst neighbour.
    real :: h ! external filed, default 0. in unit
    logical :: readspin, savespin
    integer :: cpuid  ! CPU id for mpi
    integer :: num_cpu  ! Number of processors for mpi
    integer :: ierr

#if defined (MPI)
    integer, parameter :: mpi_in= mpi_integer
    integer, parameter :: mpi_dp= mpi_double_precision
    integer, parameter :: mpi_rl= mpi_real
    integer, parameter :: mpi_dc= mpi_double_complex
    integer, parameter :: mpi_cmw= mpi_comm_world
#endif


contains
    
end module para
