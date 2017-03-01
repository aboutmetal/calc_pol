subroutine input(input_file, wfc_file, epsl_qd, epsl_sol, qd_radius, cell_size, grid_car, grid_sph)

    use module_data
    implicit none

    integer, dimension(3)   :: grid_car, grid_sph
    real(8)                 :: epsl_qd, epsl_sol, qd_radius
    real(8), dimension(3)   :: cell_size
    character(len=50)       :: input_file, wfc_file(2)

    integer                 :: i, j, ierr, iindex
    integer, parameter      :: iread = 18

    open(iread, file=trim(adjustl(input_file)), status='old', action='read', iostat=ierr)
    if(ierr.ne.0) then
        print *, "--ERROR: cannot open file ", trim(adjustl(input_file))
        stop
    endif

    rewind(iread)

    read(iread, *) iindex, (cell_size(i), i=1,3)
    write(6, '(A,3(1X,F13.5))') "Cell size:         ", (cell_size(i), i=1,3)

    read(iread, *) iindex, qd_radius, epsl_qd, epsl_sol
    write(6, '(A,1X,F13.5)')    "QD radius:         ", qd_radius
    write(6, '(A,2(1X,F13.5))') "epsilon (QD, sol): ", epsl_qd, epsl_sol

    read(iread, *) iindex, (grid_car(i), i=1,3)
    write(6, '(A,3(1X,I5))') "Grid (Cartesian):  ", (grid_car(i), i=1,3)

    read(iread, *) iindex, (grid_sph(i), i=1,3)
    write(6, '(A,3(1X,I5))') "Grid (spherical):  ", (grid_sph(i), i=1,3)

    read(iread, *) iindex, wfc_file(1)
    read(iread, *) iindex, wfc_file(2)

    close(iread)

end subroutine
