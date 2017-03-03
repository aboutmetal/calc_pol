program test_ortho

    implicit none

    complex(8), dimension(:,:,:), allocatable :: wfc1, wfc2
    real(8), dimension(3) :: dis
    integer, dimension(3) :: grid

    integer :: ir, itheta, iphi
    real(8) :: r, theta, phi, tmp1, tmp2, tmp3, factor

    integer, parameter :: iread1 = 18, iread2 = 20

    character(len=50) file1, file2

    print *, "File1?"
    read(5,*) file1
    print *, "File2?"
    read(5,*) file2

    open(iread1, file=file1, status='old', action='read', form='unformatted')
    rewind(iread1)
    read(iread1) grid
    print *, grid
    read(iread1) dis
    print *, dis

    open(iread2, file=file2, status='old', action='read', form='unformatted')
    rewind(iread2)
    read(iread2) grid
    print *, grid
    read(iread2) dis
    print *, dis

    allocate(wfc1(0:grid(1),0:grid(2),0:(grid(3)-1)))
    allocate(wfc2(0:grid(1),0:grid(2),0:(grid(3)-1)))

    read(iread1) wfc1
    read(iread2) wfc2

    close(iread1)
    close(iread2)

    tmp1 = 0.d0
    tmp2 = 0.d0
    tmp3 = 0.d0
    factor = dis(1) * dis(2) * dis(3)
    do iphi = 0, grid(3)-1
        phi = dble(iphi) * dis(3)
        do itheta = 0, grid(2)
            theta = dble(itheta) * dis(2)
            do ir = 0, grid(1)
                r = dble(ir) * dis(1)
                tmp1 = tmp1 + abs(wfc1(ir,itheta,iphi))**2 * r**2 * sin(theta)
                tmp2 = tmp2 + abs(wfc2(ir,itheta,iphi))**2 * r**2 * sin(theta)
                tmp3 = tmp3 + real(conjg(wfc1(ir,itheta,iphi)) * wfc2(ir,itheta,iphi)) * r**2 * sin(theta)
            enddo
        enddo
    enddo
    print *, "Spherical normal: ", tmp1, tmp2, tmp3

    deallocate(wfc1)
    deallocate(wfc2)

end
