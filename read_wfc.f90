subroutine read_wfc(wfc_file, grid_car)

    use module_data
    implicit none

    integer, dimension(3)               :: grid_car
    real(8), dimension(:), allocatable  :: wfc_real, wfc_imag
    character(len=50)                   :: wfc_file(2)

    integer                 :: i, j, k, ii, jj, inode, ierr
    integer                 :: n1, n2, n3, nnodes, nr, nr_n
    real(8), dimension(3,3) :: AL
    integer, parameter      :: iread1 = 18, iread2 = 20

    open(iread1, file=trim(adjustl(wfc_file(1))), status='old', action='read', iostat=ierr)
    if(ierr.ne.0) then
        print *, "--ERROR: cannot open file ", trim(adjustl(wfc_file(1)))
        stop
    endif
    open(iread2, file=trim(adjustl(wfc_file(2))), status='old', action='read', iostat=ierr)
    if(ierr.ne.0) then
        print *, "--ERROR: cannot open file ", trim(adjustl(wfc_file(2)))
        stop
    endif

    allocate(wfc_real(grid_car(1)*grid_car(2)*grid_car(3)))
    allocate(wfc_imag(grid_car(1)*grid_car(2)*grid_car(3)))

    rewind(iread1)
    read(iread1) n1,n2,n3,nnodes
    write(6,*) n1,n2,n3,nnodes
    read(iread1) AL
    write(6,*) AL(1,1),AL(2,1),AL(3,1)
    write(6,*) AL(1,2),AL(2,2),AL(3,2)
    write(6,*) AL(1,3),AL(2,3),AL(3,3)
    if(n1.ne.grid_car(1).or.&
       n2.ne.grid_car(2).or.&
       n3.ne.grid_car(3)) then
        print *, "--ERROR: grid dismatch, stop"
        stop
    endif
    read(iread1) wfc_real
    close(iread1)

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
    read(iread2) wfc_imag
    close(iread2)

    nr = n1*n2*n3
    nr_n = nr / nnodes

    do inode = 1, nnodes
        do ii = 1, nr_n
            jj = ii + (inode-1) * nr_n
            i = (jj-1) / (n2*n3) + 1
            j = (jj-1 - (i-1) * n2*n3) / n3 + 1
            k = jj - (i-1) * n2*n3 - (j-1) * n3
            wfc_car(i-1,j-1,k-1) = cmplx(wfc_real(jj), wfc_imag(jj))
        enddo
    enddo

    wfc_car(n1,:,:) = wfc_car(0,:,:)
    wfc_car(:,n2,:) = wfc_car(:,0,:)
    wfc_car(:,:,n3) = wfc_car(:,:,0)

    deallocate(wfc_real)
    deallocate(wfc_imag)

end subroutine
