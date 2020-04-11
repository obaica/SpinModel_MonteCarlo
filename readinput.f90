subroutine distance(a, b, td)
    use  para
    implicit none
    real,intent(in) :: a(3), b(3)
    real,intent(inout) :: td(9)
    integer::m,n
    real::v1(3),v2(3),v3(3),v4(3),v5(3),v6(3),v7(3),v8(3)
    v1=a-latvec(1,:)
    v2=a+latvec(1,:)
    v3=a-latvec(2,:)
    v4=a+latvec(2,:)
    v5=a-latvec(1,:)-latvec(2,:)
    v6=a-latvec(1,:)+latvec(2,:)
    v7=a+latvec(1,:)-latvec(2,:)
    v8=a+latvec(1,:)+latvec(2,:)
    td(1)=sqrt((a(1)-b(1))**2+(a(2)-b(2))**2+(a(3)-b(3))**2)
    td(2)=sqrt((v1(1)-b(1))**2+(v1(2)-b(2))**2+(v1(3)-b(3))**2)
    td(3)=sqrt((v2(1)-b(1))**2+(v2(2)-b(2))**2+(v2(3)-b(3))**2)
    td(4)=sqrt((v3(1)-b(1))**2+(v3(2)-b(2))**2+(v3(3)-b(3))**2)
    td(5)=sqrt((v4(1)-b(1))**2+(v4(2)-b(2))**2+(v4(3)-b(3))**2)
    td(6)=sqrt((v5(1)-b(1))**2+(v5(2)-b(2))**2+(v5(3)-b(3))**2)
    td(7)=sqrt((v6(1)-b(1))**2+(v6(2)-b(2))**2+(v6(3)-b(3))**2)
    td(8)=sqrt((v7(1)-b(1))**2+(v7(2)-b(2))**2+(v7(3)-b(3))**2)
    td(9)=sqrt((v8(1)-b(1))**2+(v8(2)-b(2))**2+(v8(3)-b(3))**2)

end subroutine distance

subroutine readinput()
    use para
    implicit none
    
    character*10:: filename="mc.in"
    logical :: fileexist
    integer :: m, n, k, l, nd
    integer, allocatable :: nn_count(:)
    real::randx
    real::dmat(9)
    logical :: lfound
    character*256 :: inline

    dmat=0d0
    fileexist=.false.
    inquire(file=filename, exist=fileexist)
    if (fileexist) then
        if (cpuid==0)write(stdout, *)" "
        if (cpuid==0)write(stdout, *)" Read file from mc.in"
        if (cpuid==0)write(stdout, *)'=================================================='
        if (cpuid==0)write(stdout, *)" Start read input!"
        open(unit=101, file=filename, status="old")
    else 
        if (cpuid==0)write(stdout, *)" Input file mc.in dosn't exist!"
        stop
    endif

    num_atom=0
    num_j=1
    offset=4 ! offset to avoid trend
    h=0
    temperature_btm=0d0
    temperature_top=10d0
    temperature_num=1000
    periodic=.true.
    S=1d0
    randominit=.true.
    randbias=0.6
    num_ensemble=100
    afm=.false.
    namelist /parameters/ num_atom, num_j, h, temperature_btm, temperature_top,&
                        S, num_iterr,temperature_num,offset, afm, &
                        num_ensemble,periodic,randominit,randbias
    if (h>1e-3) then ! when external field applied, strickly random
                     ! initialization method adopted.
        randominit=.true.
        randbias=0.5
    endif

    read(101, parameters, iostat=stat)
    if (afm.and.cpuid==0) write(stdout, *) "only TWO-sublattice AFM is supported!!!"
    if (afm.and.cpuid==0) write(stdout, *) "In atomsite card:"
    if (afm.and.cpuid==0) write(stdout, *) "              1 ~ Num_atom/2 are sublattice A"
    if (afm.and.cpuid==0) write(stdout, *) "     Num_atom/2 ~ Num_atom   are sublattice B"
    tnn=temperature_num
    nit=num_iterr
    nem=num_ensemble
    temperature_btm=temperature_btm+1e-10 ! a small value to avoid devided by 0

    allocate(spin(nem, num_atom+1))
    allocate(atomsite(num_atom, 3))
    allocate(d(num_j,2))
    allocate(j(num_j))
    allocate(num_nn(num_j))
    allocate(nn_count(num_j))
    
    if (afm) then ! set sublattice for afm calculation
        allocate(sublat(2,num_atom/2))
        sublat(:,:)=0
        k=0
        do n=1,num_atom/2
            k=k+1
            sublat(1,k)=n
        enddo
        k=0
        do n=1,num_atom/2
            k=k+1
            sublat(2,k)=n+num_atom/2
        enddo
    endif
  if (randominit)  then 
    call init_random()
    do m=1,num_atom
      do n=1, nem
        call random_number(randx)
        if (randx.ge.randbias) then
            spin(n,m)=S
        else
            spin(n,m)=-S
        endif
        spin(n,num_atom+1)=0e0
      enddo
    enddo
  else
    if (afm) then
      do m=1,num_atom/2
        do n=1, nem
          spin(n,m)=S
          spin(n,m+num_atom/2)=-S
          spin(n,num_atom+1)=0e0
        enddo
      enddo
    else 
      do n=1, nem
        do m=1,num_atom
          spin(n,m)=S
        enddo
        spin(n,num_atom+1)=0e0
      enddo
    endif
  endif

    rewind(101)
    lfound = .false.
    do while (.true.)
       read(101, *, end= 1001)inline
       if (trim(adjustl(inline))=='lattice_vector') then
          lfound= .true.
          exit
       endif
    enddo
 1001 continue
    read(101,*)latvec(1,1),latvec(1,2),latvec(1,3)
    read(101,*)latvec(2,1),latvec(2,2),latvec(3,3)

    rewind(101)
    lfound = .false.
    do while (.true.)
       read(101, *, end= 1005)inline
       if (trim(adjustl(inline))=='atomsite') then
          lfound= .true.
          exit
       endif
    enddo
 1005 continue
    do m=1,num_atom
        read(101,*)atomsite(m,1),atomsite(m,2),atomsite(m,3)
    enddo


    rewind(101)
    lfound = .false.
    do while (.true.)
       read(101, *, end= 1002)inline
       if (trim(adjustl(inline))=='num_nn') then
          lfound= .true.
          exit
       endif
    enddo
 1002 continue
    read(101, *) num_nn
    maxnn=maxval(num_nn)
    allocate(nn(num_atom, num_j, maxnn))
    nn=num_atom+1

    rewind(101)
    lfound = .false.
    do while (.true.)
       read(101, *, end= 1003)inline
       if (trim(adjustl(inline))=='j') then
          lfound= .true.
          exit
       endif
    enddo
 1003 continue
    read(101,*) j
    if(cpuid==0) write(stdout,*) j

    
    rewind(101)
    lfound = .false.
    do while (.true.)
       read(101, *, end= 1004)inline
       if (trim(adjustl(inline))=='d') then
          lfound= .true.
          exit
       endif
    enddo
 1004 continue
    do m=1,num_j
        read(101,*) d(m,1),d(m,2)
    enddo

    !Calculate nearest neighbour site
    do m=1,num_atom
        nn_count(:)=1
        do n=1,num_atom
            call distance(atomsite(m,:),atomsite(n,:),dmat)
            do l=1,num_j
                if (periodic) then ! periodic boundary condition
                    nd=9
                else
                    nd=1
                endif
                do k=1,nd
                    if (dmat(k)>d(l,1) .and. dmat(k)<d(l,2)) then
                        nn(m,l,nn_count(l))=n
                        nn_count(l)=nn_count(l)+1
                        exit
                    endif
                enddo
            enddo
        enddo
    enddo

    if (cpuid==0) write(stdout, *) " Your >>>>>>> mc.in <<<<<< file is:"
    if (cpuid==0) call savemodel(stdout)
    
if (cpuid==0) then
    if (cpuid==0) write(stdout, *) ""
    if (cpuid==0) write(stdout, *) " >>>>>>> Nearst neighbor <<<<<"
    write(stdout,*) "#nearst neighbor atom index"
    write(stdout,*) "#atomN,  j,     nn1,     nn2, ..."
    do m=1,num_atom
        do n=1,num_j
            write(stdout,'(10I7)')m,n,nn(m,n,:)
        enddo
    enddo
endif

    if (cpuid==0) write(stdout, *)' '
    if (cpuid==0) write(stdout, *)' End read input!'
    if (cpuid==0) write(stdout, *)'=================================================='

end subroutine readinput

subroutine savemodel
    ! save model
    use para
    implicit none
    integer :: m,n,k,l
    write(stdout, *)"&parameters"
    write(stdout, *)"num_atom=",num_atom
    write(stdout, *)"num_j=",num_j
    write(stdout, *)"num_iterr=",num_iterr
    write(stdout, *)"num_ensemble=",num_ensemble
    write(stdout, *)"S=", S
    write(stdout, *)"h=", h
    write(stdout, *)"afm=", afm
    write(stdout, *)"randominit=", randominit
    write(stdout, *)"randbias=", randbias
    write(stdout, *)"periodic=", periodic
    write(stdout, *)"temperature_btm=",temperature_btm, "! avoid 0"
    write(stdout, *)"temperature_top=",temperature_top
    write(stdout, *)"temperature_num=",temperature_num
    write(stdout, *)"/"
    
    write(stdout, *)"num_nn"
    write(stdout, '(4I3)')num_nn
    write(stdout, *)""
    write(stdout, *)"d"
    do m=1, num_j
      write(stdout, '(2F10.3)')d(m,1), d(m,2)
    enddo
    write(stdout, *)""
    write(stdout, *)"j"
    do m=1, num_j
      write(stdout, *) j(m)
    enddo
    write(stdout, *)" "
    
    write(stdout, *)"lattice_vector"
    write(stdout, '(3F15.6)')latvec(1,:)
    write(stdout, '(3F15.6)')latvec(2,:)
    
    write(stdout, *)" "
    write(stdout, *) "atomsite"
    do m=1, num_atom
        write(stdout,'(3F15.6)')atomsite(m,:)
    enddo
    
    
end subroutine savemodel
