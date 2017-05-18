subroutine input(input_file, wfc_type, wfc_file, epsl_qd, epsl_sol, &
                 qd_radius, lamda, cell_size, grid_car, grid_sph)

    use module_data
    implicit none

    integer                 :: wfc_type
    integer, dimension(3)   :: grid_car, grid_sph
    real(8)                 :: epsl_qd, epsl_sol, qd_radius, lamda
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

    print *, "Reading input parameters..."
    print *, "******************************************"
    read(iread, *) iindex, (cell_size(i), i=1,3)
    write(6, '(A,3(1X,F13.5))') "   Cell size:         ", (cell_size(i), i=1,3)

    read(iread, *) iindex, qd_radius, epsl_qd, epsl_sol, lamda
    write(6, '(A,1X,F13.5)')    "   QD radius:         ", qd_radius
    write(6, '(A,2(1X,F13.5))') "   epsilon (QD, sol): ", epsl_qd, epsl_sol
    write(6, '(A,1X,F13.5)')    "   lamda:             ", lamda

    read(iread, *) iindex, (grid_car(i), i=1,3)
    write(6, '(A,3(1X,I5))') "   Grid (Cartesian):  ", (grid_car(i), i=1,3)

    read(iread, *) iindex, (grid_sph(i), i=1,3)
    write(6, '(A,3(1X,I5))') "   Grid (spherical):  ", (grid_sph(i), i=1,3)

    read(iread, *) iindex, wfc_type
    read(iread, *) iindex, wfc_file(1)
    if(wfc_type.eq.1) then
        read(iread, *) iindex, wfc_file(2)
        print *, "  Wfc type:  PEtot"
        print *, "  Wfc files: "
        print *, trim(adjustl(wfc_file(1))), "  ", trim(adjustl(wfc_file(2)))
    elseif(wfc_type.eq.2) then
        print *, "  Wfc type:  CPMD"
        print *, "  Wfc files: "
        print *, trim(adjustl(wfc_file(1)))
    elseif(wfc_type.eq.3) then
        print *, "  Wfc type:  VASP"
        print *, "  Wfc files: "
        print *, trim(adjustl(wfc_file(1)))
    else
        print *, "-ERROR: wfc type wrong"
        stop
    endif

    close(iread)

end subroutine
