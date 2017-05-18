subroutine read_wfc_PEtot(wfc_file, grid_car)

    ! Read wavefunction (wfc) from external files
    ! The wfc should be read into 3-D Cartesian coordinates
    ! Note that the index should start from 0, not 1

    use module_data
    implicit none

    integer, dimension(3)               :: grid_car
    real(8), dimension(:), allocatable  :: wfc_real, wfc_imag
    character(len=50)                   :: wfc_file(2)

    integer                 :: i, j, k, ii, jj, inode, ierr
    integer                 :: n1, n2, n3, nnodes, nr, nr_n
    real(8), dimension(3,3) :: AL
    integer, parameter      :: iread1 = 18, iread2 = 20

    ! open wfc files
    ! PEtot outputs wfc in real and imag parts, so read them separately
    open(iread1, file=trim(adjustl(wfc_file(1))), status='old', action='read', form='unformatted', iostat=ierr)
    if(ierr.ne.0) then
        print *, "--ERROR: cannot open file ", trim(adjustl(wfc_file(1)))
        stop
    endif
    open(iread2, file=trim(adjustl(wfc_file(2))), status='old', action='read', form='unformatted', iostat=ierr)
    if(ierr.ne.0) then
        print *, "--ERROR: cannot open file ", trim(adjustl(wfc_file(2)))
        stop
    endif

    print *, "Reading wavefunctions..."
    print *, "******************************************"

    rewind(iread1)

    ! PEtot wfc file starts with grid setting and cell size
    read(iread1) n1,n2,n3,nnodes
    write(6,*) n1,n2,n3,nnodes
    read(iread1) AL
    write(6,*) AL(1,1),AL(2,1),AL(3,1)
    write(6,*) AL(1,2),AL(2,2),AL(3,2)
    write(6,*) AL(1,3),AL(2,3),AL(3,3)

    ! Check if the grid setting here is the same with input parameters
    if(n1.ne.grid_car(1).or.&
       n2.ne.grid_car(2).or.&
       n3.ne.grid_car(3)) then
        print *, "--ERROR: grid dismatch, stop"
        stop
    endif

    rewind(iread2)
    read(iread2) n1,n2,n3,nnodes
    write(6,*) n1,n2,n3,nnodes
    read(iread2) AL
    write(6,*) AL(1,1),AL(2,1),AL(3,1)
    write(6,*) AL(1,2),AL(2,2),AL(3,2)
    write(6,*) AL(1,3),AL(2,3),AL(3,3)
    if(n1.ne.grid_car(1).or.&
       n2.ne.grid_car(2).or.&
       n3.ne.grid_car(3)) then
        print *, "--ERROR: grid dismatch, stop"
        stop
    endif

    nr = n1*n2*n3
    nr_n = nr / nnodes

    allocate(wfc_real(nr_n))
    allocate(wfc_imag(nr_n))

    ! Read wfc
    ! Note again that the index starts from 0, not 1
    do inode = 1, nnodes
        read(iread1) wfc_real
        read(iread2) wfc_imag
        do ii = 1, nr_n
            jj = ii + (inode-1) * nr_n
            i = (jj-1) / (n2*n3) + 1
            j = (jj-1 - (i-1) * n2*n3) / n3 + 1
            k = jj - (i-1) * n2*n3 - (j-1) * n3
            wfc_car(i-1,j-1,k-1) = cmplx(wfc_real(ii), wfc_imag(ii))
        enddo
    enddo

    close(iread1)
    close(iread2)

    ! Set values on the boundaries
    wfc_car(n1,:,:) = wfc_car(0,:,:)
    wfc_car(:,n2,:) = wfc_car(:,0,:)
    wfc_car(:,:,n3) = wfc_car(:,:,0)

    ! Output wfc file for mxmat
    open(19,file='graph.mxmat.output',status='new',action='write')
    rewind(19)
    write(19,*) "Wave function NR:           1"
    write(19,*) "(((wr_out_cmplx(it,jt,kt,1),it=1,n1r),jt=1,n2r),kt=1,n3r)"
    do k = 1, n3
    do j = 1, n2
    do i = 1, n1
        write(19,*) wfc_car(i,j,k)
    enddo
    enddo
    enddo
    write(19,*) "Wave function NR:           1"
    write(19,*) "(((wr_out_cmplx(it,jt,kt,1),it=1,n1r),jt=1,n2r),kt=1,n3r)"
    do k = 1, n3
    do j = 1, n2
    do i = 1, n1
        write(19,*) cmplx(0.d0, 0.d0)
    enddo
    enddo
    enddo
    close(19)
    stop

    print *, "******************************************"

    deallocate(wfc_real)
    deallocate(wfc_imag)

end subroutine
