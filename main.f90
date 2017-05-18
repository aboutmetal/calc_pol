program polarization

    use module_data
    implicit none

    integer                 :: wfc_type
    integer, dimension(3)   :: grid_car, grid_sph
    real(8)                 :: epsl_qd, epsl_sol, qd_radius, lamda
    real(8), dimension(3)   :: cell_size, dis_car, dis_sph
    complex(8)              :: V

    integer                 :: i
    character(len=50)       :: input_file, wfc_file(2)

    input_file = "sigma_pol.input"

    call input(input_file, wfc_type, wfc_file, epsl_qd, epsl_sol, &
               qd_radius, lamda, cell_size, grid_car, grid_sph)

    call module_allocate_wfc(1, grid_car)
    dis_car(1:3) = cell_size(1:3) / dble(grid_car(1:3))
    write(6, '(A,3(1X,F13.5))') "   Grid size (Cart):  ", (dis_car(i), i=1,3)
    print *, "******************************************"

    if(wfc_type.eq.1) then
        call read_wfc_PEtot(wfc_file, grid_car)
    elseif(wfc_type.eq.2) then
        call read_wfc_CPMD(wfc_file, grid_car)
    elseif(wfc_type.eq.3) then
        call read_wfc_VASP(wfc_file, grid_car)
    endif

    call module_allocate_wfc(2, grid_sph)

    call car_to_sph(grid_car, dis_car, grid_sph, dis_sph)

    call integral(V, epsl_qd, epsl_sol, qd_radius, lamda, grid_sph, dis_sph)

    call module_deallocate_wfc(1)
    call module_deallocate_wfc(2)

end
