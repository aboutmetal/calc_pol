subroutine read_wfc_CPMD(wfc_file, grid_car)

    ! Read wavefunction (wfc) from external files
    ! The wfc should be read into 3-D Cartesian coordinates
    ! Note that the index should start from 0, not 1

    use module_data
    implicit none

    integer, dimension(3)               :: grid_car
    character(len=50)                   :: wfc_file(2), chtmp

    integer                 :: natom
    integer                 :: i, j, k, ii, jj, inode, ierr
    integer                 :: n1, n2, n3, nnodes, nr, nr_n
    real(8)                 :: rtmp, tmp(6)
    real(8), dimension(3,3) :: AL
    integer, parameter      :: iread1 = 18, iread2 = 20

    ! open wfc files
    ! PEtot outputs wfc in real and imag parts, so read them separately
    open(iread1, file=trim(adjustl(wfc_file(1))), status='old', action='read', iostat=ierr)
    if(ierr.ne.0) then
        print *, "--ERROR: cannot open file ", trim(adjustl(wfc_file(1)))
        stop
    endif

    print *, "Reading wavefunctions..."
    print *, "******************************************"

    rewind(iread1)

    ! Begin to read the file
    read(iread1,*)
    read(iread1,*)
    read(iread1,*) natom, rtmp, rtmp, rtmp
    read(iread1,*) n1, AL(1,1),AL(2,1),AL(3,1)
    read(iread1,*) n2, AL(1,2),AL(2,2),AL(3,2)
    read(iread1,*) n3, AL(1,3),AL(2,3),AL(3,3)
    AL(:,1) = AL(:,1) * n1
    AL(:,2) = AL(:,2) * n2
    AL(:,3) = AL(:,3) * n3

    ! Check if the grid setting here is the same with input parameters
    write(6,*) n1,n2,n3
    write(6,*) AL(1,1),AL(2,1),AL(3,1)
    write(6,*) AL(1,2),AL(2,2),AL(3,2)
    write(6,*) AL(1,3),AL(2,3),AL(3,3)
    if(n1.ne.grid_car(1).or.&
       n2.ne.grid_car(2).or.&
       n3.ne.grid_car(3)) then
        print *, "--ERROR: grid dismatch, stop"
        stop
    endif

    do ii = 1, natom
        read(iread1,*)
    enddo

    nr_n = floor(dble(n3) / 6.d0)

    ! Read wfc
    ! Note again that the index starts from 0, not 1
    do i = 1, n1
        do j = 1, n2
            do k = 1, nr_n
                ii = (k - 1) * 6
                read(iread1,*) (tmp(jj), jj=1,6)
                wfc_car(i,j,(ii+1):(ii+6)) = cmplx(tmp(1:6), 0.d0)
            enddo
            read(iread1,*) (tmp(jj), jj=1,(n3-6*nr_n))
            wfc_car(i,j,(6*nr_n+1):n3) = cmplx(tmp(1:(n3-6*nr_n)), 0.d0)
        enddo
    enddo

    close(iread1)

    ! Set values on the boundaries
    wfc_car(0,:,:) = wfc_car(n1,:,:)
    wfc_car(:,0,:) = wfc_car(:,n2,:)
    wfc_car(:,:,0) = wfc_car(:,:,n3)

!    do k = 0, 3
!    do j = 0, 3
!    do i = 0, n3
!        print *, i, wfc_car(1,1,i)
!    enddo
!    stop
!    enddo
!    enddo

    print *, "******************************************"

end subroutine
