  subroutine now(time_now)
     use para
     implicit none
     integer   :: time_new(8)
     real  :: time_now
     call Date_and_time(values=time_new)
     time_now= time_new(3)*24*3600+time_new(5)*3600+&
               time_new(6)*60+time_new(7)+time_new(8)/1000d0

     return
  end subroutine now

subroutine mean(s, shp, avg)
    ! average function
    implicit none
    integer, intent(in) :: shp
    real, intent(in):: s(shp)
    real, intent(inout)::avg
    integer :: m,n,k,l
    
    avg=0e0
    do m=1, shp
        avg=avg+s(m)
    enddo
    avg=avg/real(shp)
end subroutine mean

subroutine init_random
    ! init random function to generate different random number
    implicit none
    integer :: i, n, clock
    integer, allocatable :: seed(:)
    call random_seed(size=n)
    allocate(seed(n))
    call system_clock(count=clock)
    seed=clock+37*(/(i-1, i=1,n)/)
    call random_seed(put=seed)
    deallocate(seed)
end subroutine init_random


subroutine redirect_spin(k)
    ! redirect spin orientation to make average easier
    use para
    implicit none
    integer, intent(in) :: k
    integer::m,n
    integer:: nup(num_ensemble), ndn(num_ensemble)
    
  if (k.eq.num_iterr/2 .or. h>1e-5) then 
    return
  else ! if external field are applied, then turn off redirect spin.
    
    nup=0
    ndn=0

    if (afm) then
      do n=1,num_ensemble
        do m=1, num_atom/2 ! two lattice model, half spin are included.
            if (spin(n, m) > 0) then
                nup(n)=nup(n)+1
            else
                ndn(n)=ndn(n)+1
            endif
        enddo
      enddo
      
    else
      do n=1,num_ensemble
        do m=1, num_atom
            if (spin(n, m) > 0) then
                nup(n)=nup(n)+1
            else
                ndn(n)=ndn(n)+1
            endif
        enddo
      enddo
    endif

    ! for each spin lat, if up<dn, flip all spin
    do n=1,num_ensemble
      if (nup(n)<ndn(n)) then
        do m=1,num_atom
            spin(n,m)=-spin(n,m)
        enddo
      endif
    enddo

  endif   
    
end subroutine
