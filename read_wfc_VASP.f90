subroutine read_wfc_VASP(wfc_file, grid_car)

    ! Read wavefunction (wfc) from external files
    ! The wfc should be read into 3-D Cartesian coordinates
    ! Note that the index should start from 0, not 1

    use module_data
    implicit none

    integer, dimension(3)               :: grid_car
    character(len=50)                   :: wfc_file(2), chtmp

    integer                 :: natom
    integer                 :: i, j, k, ii, jj, kk, inode, iline, ierr
    integer                 :: n1, n2, n3, nnodes, nr, nr_n
    real(8)                 :: rtmp, tmp(10)
    real(8), dimension(3,3) :: AL
    integer, parameter      :: iread1 = 18, iread2 = 20
    character(255)          :: line

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

    ! Search for the beginning of wfc array
    iline = 0
    ierr = 1
    do while(ierr.ne.0)
        iline = iline + 1
        read(iread1, '(a)') line
        read(line, *, iostat=ierr) (tmp(i), i=1,10)
    enddo
    print *, "--DEBUG: wfc array starts from line", iline

    n1 = grid_car(1)
    n2 = grid_car(2)
    n3 = grid_car(3)

    ! Read wfc
    nr_n = floor(dble(n1*n2*n3) / 10.d0)

    ! Read wfc
    ! Note again that the index starts from 0, not 1
    do ii = 1, nr_n
        if(ii.ne.1) then
            iline = iline + 1
            read(iread1, '(a)') line
        endif
        read(line, *) (tmp(i), i=1,10)
        do jj = 1, 10
            kk = (ii - 1) * 10 + jj
            k = (kk - 1) / (n1 * n2) + 1
            j = (kk - 1 - (k - 1) * n1 * n2) / n1 + 1
            i = kk - (k - 1) * n1 * n2 - (j - 1) * n1
            wfc_car(i-1,j-1,k-1) = cmplx(tmp(jj), 0.d0)
        enddo
    enddo
    ii = n1*n2*n3 - nr_n*10
    if(ii.ne.0) then
        iline = iline + 1
        read(iread1, '(a)') line
        read(line, *) (tmp(i), i=1, ii)
        do jj = 1, ii
            kk = nr_n * 10 + jj
            k = (kk - 1) / (n1 * n2) + 1
            j = (kk - 1 - (k - 1) * n1 * n2) / n1 + 1
            i = kk - (k - 1) * n1 * n2 - (j - 1) * n1
            wfc_car(i-1,j-1,k-1) = cmplx(tmp(jj), 0.d0)
        enddo
    endif
    print *, "--DEBUG: wfc array ends with line", iline
    close(iread1)

    ! Set values on the boundaries
    wfc_car(n1,:,:) = wfc_car(0,:,:)
    wfc_car(:,n2,:) = wfc_car(:,0,:)
    wfc_car(:,:,n3) = wfc_car(:,:,0)

    write(6,'(10(1X,E10.3))'), (abs(wfc_car(i,0,0)), i=0,9)
    write(6,'(10(1X,E10.3))'), (abs(wfc_car(i,0,0)), i=10,19)
    write(6,'(10(1X,E10.3))'), (abs(wfc_car(i,0,0)), i=20,29)
    write(6,'(10(1X,E10.3))'), (abs(wfc_car(i,0,0)), i=30,39)
    write(6,'(10(1X,E10.3))'), (abs(wfc_car(i,0,0)), i=40,49)
    stop

    print *, "******************************************"

end subroutine
