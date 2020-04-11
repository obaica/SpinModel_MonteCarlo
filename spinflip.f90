subroutine spinflip()
    use para
    implicit none
    integer :: m,n,k,l,iembl
    real :: esite, beta,randx, avg, avg3(3)
    real :: time_start, time_end
    logical::flip
    integer :: sp
    real :: p

    real::t(tnn), e(nit, tnn),mag(nit, tnn),sus(nit, tnn), c(nit, tnn)
    real :: ttp, etp, magtp, sustp, ctp, matp, mbtp, xatp, xbtp
    real::tmpi(tnn),empi(nit, tnn),magmpi(nit, tnn), susmpi(nit, tnn), cmpi(nit, tnn)
    real::ma(nit, tnn),mampi(nit,tnn),mb(nit,tnn),mbmpi(nit,tnn)
    real::xa(nit, tnn),xampi(nit,tnn),xb(nit,tnn),xbmpi(nit,tnn)
! memory explosion
    sp=nit-nit/10
    t=0d0
    e=0d0
    mag=0d0
    sus=0d0
    c=0d0

    ma=0d0
    mampi=0d0
    mb=0d0
    mbmpi=0d0
    xa=0d0
    xampi=0d0
    xb=0d0
    xbmpi=0d0

    tmpi=0d0
    empi=0d0
    magmpi=0d0
    susmpi=0d0
    cmpi=0d0
    time_start=0d0
    call init_random()
    if (cpuid.eq.0) write(stdout, *) "Start spin flipping proccess..."

    do m=1+cpuid, tnn, num_cpu
        tmpi(m)=(temperature_top-temperature_btm)*(m-1)/tnn+temperature_btm
        if (cpuid.eq.0) then
            call now(time_end)
            write(stdout, '(a, f10.2, "/", f10.2 )') 'temperature',  tmpi(m), temperature_top
            time_start= time_end
        endif   
        beta=1/(kb*tmpi(m))
        do k=1,nit
          do iembl=1,nem
            do l=1,offset
            do n=l,num_atom,offset
                esite=0d0
                call energy_site_flip(n, esite, iembl)
                p=1.0/(1.0+exp(esite*beta))
                call random_number(randx)
                if (p.gt.randx) spin(iembl, n)=-spin(iembl, n)
            enddo
            enddo
          enddo
          call redirect_spin(k)
          do iembl=1,nem
              call energy(etp, iembl)
              empi(k,m)=empi(k,m)+etp/real(nem)
              call magnetism(magtp, iembl)
              magmpi(k,m)=magmpi(k,m)+magtp/real(nem)
              call susceptibility(sustp, tmpi(m), iembl)
              susmpi(k,m)=susmpi(k,m)+sustp/real(nem)
              call specific_heat(ctp, tmpi(m), iembl)
              cmpi(k,m)=cmpi(k,m)+ctp/real(nem)
              if (afm) then
                  call submag(matp, mbtp, iembl)
                  mampi(k,m)=mampi(k,m)+matp/real(nem)
                  mbmpi(k,m)=mbmpi(k,m)+mbtp/real(nem)  
                  call subsusceptibility(xatp, xbtp, tmpi(m), iembl)
                  xampi(k,m)=xampi(k,m)+xatp/real(nem)
                  xbmpi(k,m)=xbmpi(k,m)+xbtp/real(nem)
              endif
          enddo
        enddo

    enddo
#if defined(MPI)
    call mpi_allreduce(tmpi, t,size(t),mpi_rl,mpi_sum,mpi_cmw,ierr)
    call mpi_allreduce(empi, e,size(e),mpi_rl,mpi_sum,mpi_cmw,ierr)
    call mpi_allreduce(magmpi, mag,size(mag),mpi_rl,mpi_sum,mpi_cmw,ierr)
    call mpi_allreduce(susmpi, sus,size(sus),mpi_rl,mpi_sum,mpi_cmw,ierr)
    call mpi_allreduce(cmpi, c,size(c),mpi_rl,mpi_sum,mpi_cmw,ierr)
    if (afm) then
        call mpi_allreduce(mampi, ma,size(ma),mpi_rl,mpi_sum,mpi_cmw,ierr)
        call mpi_allreduce(mbmpi, mb,size(ma),mpi_rl,mpi_sum,mpi_cmw,ierr)
        call mpi_allreduce(xampi, xa,size(xa),mpi_rl,mpi_sum,mpi_cmw,ierr)
        call mpi_allreduce(xbmpi, xb,size(xa),mpi_rl,mpi_sum,mpi_cmw,ierr)
    endif
#endif
    if (cpuid.eq.0) write(stdout, *) "End spin flipping proccess!"

if (cpuid==0) then
    outfile=103
    open(unit=outfile,file="energy.dat")
    do m=1, tnn
        avg=0e0
        call mean(e(sp:nit, m), nit-sp+1, avg)
        write(outfile, *)t(m), avg
    enddo
    close(outfile)
  if (afm) then
    outfile=outfile+1
    open(unit=outfile,file="magnetization.dat")
    write(outfile,*)"#temperature, total Mag, Mag A, Mag B, abs.Mag A, -abs.Mag B"
    do m=1, tnn
        avg3(:)=0e0
        call mean(mag(sp:nit, m), nit-sp+1, avg3(1))
        call mean(ma(sp:nit, m), nit-sp+1, avg3(2))
        call mean(mb(sp:nit, m), nit-sp+1, avg3(3))
        write(outfile, '(6F10.3)')t(m), avg3(:), abs(avg3(2)), -abs(avg3(3))
    enddo
    close(outfile)
    open(unit=outfile,file="susceptibility.dat")
    write(outfile,*)"#temperature,      Xa,      Xb"
    do m=1, tnn
        avg3(:)=0e0
        call mean(xa(sp:nit, m), nit-sp+1, avg3(1))
        call mean(xb(sp:nit, m), nit-sp+1, avg3(2))
        write(outfile, '(F10.3,2E15.5)')t(m), avg3(1), avg3(2)
    enddo
    close(outfile)
  else
    outfile=outfile+1
    open(unit=outfile,file="magnetization.dat")
    do m=1, tnn
        avg=0e0
        call mean(mag(sp:nit, m), nit-sp+1, avg)
        write(outfile, *)t(m), avg
    enddo
    close(outfile)
    outfile=outfile+1
    open(unit=outfile,file="susceptibility.dat")
    do m=1, tnn
        avg=0e0
        call mean(sus(sp:nit, m), nit-sp+1, avg)
        write(outfile, '(F10.3,E15.5)')t(m), avg
    enddo
    close(outfile)
    outfile=outfile+1
    open(unit=outfile,file="specific_heat.dat")
    do m=1, tnn
        avg=0e0
        call mean(c(sp:nit, m), nit-sp+1, avg)
        write(outfile, '(F10.3,E15.5)')t(m), avg
    enddo
    close(outfile)
  endif
endif


end subroutine spinflip
