subroutine integral(V, epsl_in, epsl_out, radius, grid, dis)

    use module_data
    implicit none

    integer, dimension(3)   :: grid
    real(8), dimension(3)   :: dis
    real(8)                 :: epsl_in, epsl_out, radius
    complex(8)              :: V

    integer :: i, j, k, l
    real(8) :: lamda, r, theta, phi
    real(8), dimension(:), allocatable :: Vs

    allocate(Vs(0:grid(1)))

    lamda = 1.d-2

    ! Lengendre
    print *, "Calculating polarization potential..."
    print *, "******************************************"
    Vs = 0.d0
    do i = 0, grid(1)
        r = i * dis(1)
        if(r.le.radius) then
            do l = 1, 50
                Vs(i) = Vs(i) + (1.d0 - exp(0.d0 - (r-radius)**2 / lamda**2)) * &
                        (dble(l) + 1.d0) * (r/radius)**(2*l) / &
                        (epsl_out + dble(l) * (epsl_in + epsl_out))
            enddo
            Vs(i) = 0.d0 - Vs(i) * (epsl_in - epsl_out) / epsl_in / 2.d0 / radius
        else
            do l = 1, 50
                Vs(i) = Vs(i) - (1.d0 - exp(0.d0 - (r-radius)**2 / lamda**2)) * &
                        dble(l) * (radius/r)**(2*l+2) / &
                        (epsl_out + dble(l) * (epsl_in + epsl_out))
            enddo
            Vs(i) = 0.d0 - Vs(i) * (epsl_in - epsl_out) / epsl_out / 2.d0 / radius
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
                V = V + conjg(wfc_sph(i,j,k)) * Vs(i) * wfc_sph(i,j,k) * r**2 * sin(theta)
            enddo
        enddo
    enddo

    deallocate(Vs)

end subroutine
