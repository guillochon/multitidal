subroutine com_accel (ptx, pty, ptz, mode, accelx, accely, accelz)
    use Gravity_data, ONLY: grv_factor, orb_t, orb_dt
    use Grid_data, ONLY: gr_myPE
    use Logfile_interface, ONLY: Logfile_stamp
    use gr_mpoleData, ONLY: Xcm, Ycm, Zcm, oXcm, oYcm, oZcm
    use PhysicalConstants_interface, ONLY: PhysicalConstants_get
    use Grid_interface, ONLY: Grid_getMinCellSize

    implicit none
#include "Flash.h"
#include "constants.h"
 
    include "Flash_mpi.h"
    real :: dr32, x, y, z
    real, intent(IN)  :: ptx, pty, ptz
    real :: xcm_accel, ycm_accel
    real :: newton, dx, potr, potl
    integer, intent(IN) :: mode
    real, intent(OUT) :: accelx, accely, accelz
    character(len=256) :: str_buffer
    
    !write(str_buffer, *) 'COM: ', Xcm, Ycm, Zcm, ptx, pty, ptz
    !call Logfile_stamp(gr_myPE, str_buffer, 'Note')
    if (mode .eq. 1) then
        x = Xcm - ptx
        y = Ycm - pty
        z = Zcm - ptz
    else
        x = oXcm - ptx
        y = oYcm - pty
        z = oZcm - ptz
    endif
    dr32 = (x**2.d0 + y**2.d0 + z**2.d0)**(3./2.)

    !call PhysicalConstants_get("Newton", newton)
    !call Grid_getMinCellSize(dx)
    !dx = dx / 2.d0
    !if (mode .eq. 1) then
    !    call gr_zonePotential(-dx, 0.d0, 0.d0, potl)
    !    call gr_zonePotential( dx, 0.d0, 0.d0, potr)
    !    xcm_accel = newton*(potr - potl) / dx / 2.d0
    !    call gr_zonePotential(0.d0, -dx, 0.d0, potl)
    !    call gr_zonePotential(0.d0,  dx, 0.d0, potr)
    !    ycm_accel = newton*(potr - potl) / dx / 2.d0
    !else
    !    call gr_zoneOldPotential(-dx, 0.d0, 0.d0, potl)
    !    call gr_zoneOldPotential( dx, 0.d0, 0.d0, potr)
    !    xcm_accel = newton*(potr - potl) / dx / 2.d0
    !    call gr_zoneOldPotential(0.d0, -dx, 0.d0, potl)
    !    call gr_zoneOldPotential(0.d0,  dx, 0.d0, potr)
    !    ycm_accel = newton*(potr - potl) / dx / 2.d0
    !endif

    accelx = grv_factor*x/dr32! + xcm_accel
    accely = grv_factor*y/dr32! + ycm_accel
    accelz = grv_factor*z/dr32

    !dr32 = ((Xcm - ptx)**2 + (Ycm - pty)**2 + (Zcm - ptz)**2)**(3./2.)
    !accelx = grv_factor*(Xcm - ptx)/dr32
    !accely = grv_factor*(Ycm - pty)/dr32
    !accelz = grv_factor*(Zcm - ptz)/dr32
end
