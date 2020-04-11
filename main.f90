program MonteCarlo
    use para
    implicit none
    integer :: k,l,m,n

    ierr=0
    cpuid=0
    num_cpu=0
#if defined (MPI)
    call mpi_init(ierr)
    call mpi_comm_rank(mpi_cmw, cpuid, ierr)
    call mpi_comm_size(mpi_cmw, num_cpu, ierr)
#endif
    if (cpuid==0) open(unit=stdout, file="mc.out" )

    call readinput()
    call spinflip()
    
    if (cpuid==0) close(stdout)
#if defined(MPI)
    call mpi_finalize(ierr)
#endif
end program MonteCarlo

