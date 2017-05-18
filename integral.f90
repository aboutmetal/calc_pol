subroutine integral(V, epsl_in, epsl_out, radius, lamda, grid, dis)

    use module_data
    implicit none

    integer, dimension(3)   :: grid
    real(8), dimension(3)   :: dis
    real(8)                 :: epsl_in, epsl_out, radius, lamda
    complex(8)              :: V

    integer :: i, j, k, l, iepsl
    real(8) :: r, theta, phi
    real(8), dimension(:), allocatable :: Vs

    character(len=20) :: filename, fileindex

    allocate(Vs(0:grid(1)))

    write(fileindex, '(E10.3)') lamda
    filename = 'log.'//trim(adjustl(fileindex))
    open(19, file=filename)
    rewind(19)

    ! Lengendre
    print *, "Calculating polarization potential..."
    print *, "******************************************"

    do iepsl = 0, 190

    epsl_out = 1.d0 + 0.1d0 * dble(iepsl)

    Vs = 0.d0
    do i = 0, grid(1)
        r = i * dis(1)
        if(r.le.radius) then
            do l = 0, 100
                Vs(i) = Vs(i) + (1.d0 - exp(0.d0 - (r-radius)**2 / lamda**2)) * &
                        (dble(l) + 1.d0) * (r/radius)**(2*l) / &
                        (epsl_out + dble(l) * (epsl_in + epsl_out))
            enddo
            Vs(i) = Vs(i) * (epsl_in - epsl_out) / epsl_in / radius / 2.d0
        else
            do l = 1, 100
                Vs(i) = Vs(i) - (1.d0 - exp(0.d0 - (r-radius)**2 / lamda**2)) * &
                        dble(l) * (radius/r)**(2*l+2) / &
                        (epsl_out + dble(l) * (epsl_in + epsl_out))
            enddo
            Vs(i) = Vs(i) * (epsl_in - epsl_out) / epsl_out / radius / 2.d0
        endif
    enddo

!    open(22, file='Vs_plot', status='new', action='write')
!    rewind(22)
!    do i = 0, grid(1)
!        r = i * dis(1)
!        write(22,*) r, Vs(i)
!    enddo
!    close(22)
!
!    stop

    ! Integral
    V = (0.d0, 0.d0)
    do k = 0, grid(3)-1
        phi = dble(k) * dis(3)
        do j = 0, grid(2)
            theta = dble(j) * dis(2)
            do i = 0, grid(1)
                r = dble(i) * dis(1)
                V = V + abs(wfc_sph(i,j,k))**2 * Vs(i) * r**2 * sin(theta)
            enddo
        enddo
    enddo

    write(19, '(2(1X,F17.9))') epsl_out, real(V)*27.211396d0

    enddo

    close(19)
    deallocate(Vs)

end subroutine
