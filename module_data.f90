      module module_data

      complex*16, allocatable, dimension(:,:,:) :: wfc_car, wfc_sph
      real(8), parameter :: Pi = 3.1415926534d0

      contains

      subroutine module_allocate_wfc(flag, grid)

      implicit none

      integer :: flag
      integer, dimension(3) :: grid
      if(flag.eq.1) allocate(wfc_car(0:grid(1),0:grid(2),0:grid(3)))
      if(flag.eq.2) allocate(wfc_sph(0:grid(1),0:grid(2),0:(grid(3)-1)))

      end subroutine module_allocate_wfc

      subroutine module_deallocate_wfc(flag)

      implicit none

      integer :: flag
      if(flag.eq.1) deallocate(wfc_car)
      if(flag.eq.2) deallocate(wfc_sph)

      end subroutine module_deallocate_wfc

      end module module_data
