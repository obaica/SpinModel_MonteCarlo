subroutine energy_site(site, esite, ie)
    use para
    implicit none
    integer, intent(in) :: site, ie
    real, intent(inout) :: esite
    integer :: m, n,k,l
    esite=0d0
    
    do m=1,num_j
        do n=1,num_nn(m)
            l=nn(site,m,n)! nearest site
            esite=esite+1.0/2.0*spin(ie,site)*spin(ie,l)*j(m)
        enddo
        esite=esite-spin(ie,site)*h
    enddo
end subroutine energy_site

subroutine energy_site_flip(site, esite,ie)
    use para
    implicit none
    integer, intent(in) :: site,ie
    real, intent(inout) :: esite
    integer :: m, n,k,l
    esite=0d0
    
    do m=1,num_j
        do n=1,num_nn(m)
            l=nn(site,m,n)! nearest site
            esite=esite-spin(ie,site)*spin(ie,l)*j(m)
        enddo
        esite=esite-spin(ie,site)*h
    enddo
end subroutine energy_site_flip


subroutine energy(E, ie)
    use para
    implicit none
    integer, intent(in):: ie
    real, intent(inout) :: E
    real::temp_E
    integer :: m
    do m=1,num_atom
        call energy_site(m, temp_E, ie)
        E=E+temp_E
    enddo
end subroutine energy

subroutine magnetism(mag,ie)
    use para
    implicit none
    integer, intent(in)::ie
    real, intent(inout) :: mag
    integer ::m,n,k,l
    mag=0e0
    do m=1,num_atom
        mag=mag+spin(ie,m)
    enddo
    mag=mag/num_atom
end subroutine magnetism

subroutine susceptibility(sus,t,ie)
    use para
    implicit none
    real, intent(in) :: t
    integer, intent(in) :: ie
    real,intent(inout) :: sus
    integer ::m,n,k,l
    real::beta
    real:: mag, m2
    sus=0e0

    beta=1/(kb*t)
    call magnetism(mag, ie)
    m2=0
    do m=1,num_atom
        m2=m2+spin(ie, m)**2
    enddo

    sus=beta*num_atom*(m2/num_atom-mag**2)
end subroutine susceptibility

subroutine specific_heat(spe, t,ie)
    use para
    implicit none
    real,intent(in)::t
    integer, intent(in)::ie
    real,intent(inout) :: spe
    integer ::m,n,k,l
    real:: e2, etemp, e
    spe=0e0
    call energy(e,ie)
    do m=1,num_atom
        call energy_site(m,etemp,ie)
        e2=e2+etemp**2
    enddo

    spe=num_atom*(e2/num_atom-(e/num_atom)**2)/(kb*t**2)
end subroutine specific_heat

subroutine submag(ma, mb,ie)
    use para
    implicit none
    integer, intent(in)::ie
    real, intent(inout) :: ma, mb
    integer ::m,n,k,l
    ma=0e0
    do m=1,num_atom/2
        ma=ma+spin(ie,sublat(1,m))
    enddo
    ma=ma/(num_atom/2)

    mb=0e0
    do m=1,num_atom/2
        mb=mb+spin(ie,sublat(2,m))
    enddo
    mb=mb/(num_atom/2)
end subroutine submag

subroutine subsusceptibility(xa, xb,t,ie)
    use para
    implicit none
    real, intent(in) :: t
    integer, intent(in) :: ie
    real,intent(inout) :: xa, xb
    integer ::m,n,k,l
    real::beta
    real:: ma,mb, ma2, mb2
    xa=0e0
    xb=0e0

    beta=1/(kb*t)
    call submag(ma,mb, ie)
    ma2=0
    mb2=0
    do m=1,num_atom/2
        ma2=ma2+spin(ie, sublat(1,m))**2
    enddo

    do m=1,num_atom/2
        mb2=mb2+spin(ie, sublat(2,m))**2
    enddo
    xa=beta*(num_atom/2)*(ma2/(num_atom/2)-ma**2)
    xb=beta*(num_atom/2)*(mb2/(num_atom/2)-mb**2)
end subroutine subsusceptibility
