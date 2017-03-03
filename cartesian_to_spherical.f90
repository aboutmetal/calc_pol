subroutine car_to_sph(grid_in, dis_in, grid_out, dis_out)

    ! transfrom from cartesian to spherical
    ! using trilinear interpolation

    use module_data
    implicit none

    integer,    dimension(3)    :: grid_in      ! xyz gridding, input
    real(8),    dimension(3)    :: dis_in       ! xyz grid length, input

    ! spherical coordinate : r, theta, phi
    integer,    dimension(3)    :: grid_out     ! spherical gridding, input
    real(8),    dimension(3)    :: dis_out      ! spherical grid length, output

    integer                     :: ii, jj, kk, ir, itheta, iphi
    integer                     :: x0, x1, y0, y1, z0, z1   ! neighbors
    real(8)                     :: r, theta, phi, factor, tmp, tmp2, tmp3
    complex(8)                  :: c00, c01, c10, c11, c0, c1, c
    real(8),    dimension(3)    :: center, pos

    ! set center point to cell center
    center(1:3) = dis_in(1:3) * dble(grid_in(1:3)) / 2.d0

    ! select the shortest direction, then determine the r gridding length
    dis_out(1) = min(center(1), center(2), center(3)) / dble(grid_out(1))

    dis_out(2) = Pi / dble(grid_out(2))
    dis_out(3) = 2.d0*Pi / dble(grid_out(3))


    print *, "Re-gridding from Cartesian to spherical..."
    print *, "******************************************"
    write(6, '(A,3(1X,F13.5))') "   Grid size (sphe):  ", (dis_out(ii), ii=1,3)
    print *, "******************************************"
    pos = 0.d0
    do iphi = 0, grid_out(3)-1
        phi = dble(iphi) * dis_out(3)
        do itheta = 0, grid_out(2)
            theta = dble(itheta) * dis_out(2)
            do ir = 0, grid_out(1)
                r = dble(ir) * dis_out(1)

                ! calc position
                pos(1) = r * sin(theta) * cos(phi) + center(1)
                pos(2) = r * sin(theta) * sin(phi) + center(2)
                pos(3) = r * cos(theta)            + center(3)

                ! find neighbors
                x0 = floor(pos(1) / dis_in(1))
                x1 = x0 + 1
                if(x1.gt.grid_in(1)) x1 = 1
                pos(1) = (pos(1) - dble(x0) * dis_in(1)) / dis_in(1)
                y0 = floor(pos(2) / dis_in(2))
                y1 = y0 + 1
                if(y1.gt.grid_in(2)) y1 = 1
                pos(2) = (pos(2) - dble(y0) * dis_in(2)) / dis_in(2)
                z0 = floor(pos(3) / dis_in(3))
                z1 = z0 + 1
                if(z1.gt.grid_in(3)) z1 = 1
                pos(3) = (pos(3) - dble(z0) * dis_in(3)) / dis_in(3)

                ! interpolation
                c00 = wfc_car(x0,y0,z0) * (1.d0 - pos(1)) + &
                      wfc_car(x1,y0,z0) * pos(1)
                c01 = wfc_car(x0,y0,z1) * (1.d0 - pos(1)) + &
                      wfc_car(x1,y0,z1) * pos(1)
                c10 = wfc_car(x0,y1,z0) * (1.d0 - pos(1)) + &
                      wfc_car(x1,y1,z0) * pos(1)
                c11 = wfc_car(x0,y1,z1) * (1.d0 - pos(1)) + &
                      wfc_car(x1,y1,z1) * pos(1)
                c0 = c00 * (1.d0 - pos(2)) + c10 * pos(2)
                c1 = c01 * (1.d0 - pos(2)) + c11 * pos(2)
                c = c0 * (1.d0 - pos(3)) + c1 * pos(3)
                wfc_sph(ir, itheta, iphi) = c
            enddo
        enddo
    enddo

    ! check the normalization
    print *, "Checking normalization..."
    print *, "******************************************"
    tmp = 0.d0
    tmp2 = 0.d0
    tmp3 = 0.d0
    factor = dis_in(1) * dis_in(2) * dis_in(3)
    do kk = 0, grid_in(3)-1
        do jj = 0, grid_in(2)-1
            do ii = 0, grid_in(1)-1
                tmp = tmp + abs(wfc_car(ii,jj,kk))**2
                tmp2 = sqrt((dble(ii) * dis_in(1) - center(1))**2 + &
                            (dble(jj) * dis_in(2) - center(2))**2 + &
                            (dble(kk) * dis_in(3) - center(3))**2)
                if(tmp2 .gt. (dble(grid_out(1))*dis_out(1))) then
                    tmp3 = tmp3 + abs(wfc_car(ii,jj,kk))**2
                endif
            enddo
        enddo
    enddo
    print *, "  Cartesian normal: ", tmp, tmp * factor
    print *, "  Cartesian corner: ", tmp3, tmp3 * factor

    tmp = 0.d0
    factor = dis_out(1) * dis_out(2) * dis_out(3)
    do iphi = 0, grid_out(3)-1
        phi = dble(iphi) * dis_out(3)
        do itheta = 0, grid_out(2)
            theta = dble(itheta) * dis_out(2)
            do ir = 0, grid_out(1)
                r = dble(ir) * dis_out(1)
                tmp = tmp + abs(wfc_sph(ir,itheta,iphi))**2 * r**2 * sin(theta)
            enddo
        enddo
    enddo
    print *, "  Spherical normal: ", tmp, tmp * factor
    print *, "******************************************"

    factor = sqrt(1.d0 / tmp)
    wfc_sph = wfc_sph * factor

!    ! test output
!    open(22, file='wfc_sph', status='new', action='write', form='unformatted')
!    rewind(22)
!    write(22) grid_out
!    write(22) dis_out
!    write(22) wfc_sph
!    close(22)
!
!    open(22, file='car_plot', status='new', action='write')
!    rewind(22)
!    do jj = 0, grid_in(2)
!        phi = dble(jj) * dis_in(2)
!        do ii = 0, grid_in(1)
!            r = dble(ii) * dis_in(1)
!            write(22,'(4(1X,E15.8))') r, phi, real(wfc_car(ii, jj, grid_in(3)/2)), imag(wfc_car(ii, jj, grid_in(3)/2))
!        enddo
!    enddo
!    close(22)
!
!    stop
!
!    open(22, file='sph_plot', status='new', action='write')
!    rewind(22)
!    do iphi = 0, grid_out(3)-1
!        phi = dble(iphi) * dis_out(3) * 180.d0 / Pi
!        do ir = 0, grid_out(1)
!            r = dble(ir) * dis_out(1)
!            write(22,'(4(1X,E15.8))') phi, r, real(wfc_sph(ir, grid_out(2)/2, iphi)), imag(wfc_sph(ir, grid_out(2)/2, iphi))
!        enddo
!    enddo
!    close(22)
!
!    stop

end subroutine
