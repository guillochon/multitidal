subroutine com_accel (ptx, pty, ptz, accelx, accely, accelz)
    use Gravity_data, ONLY: grv_factor
    use Grid_data, ONLY: gr_myPE
    use Logfile_interface, ONLY: Logfile_stamp
    !use Grid_interface, ONLY: Grid_getListOfBlocks
    use gr_mpoleData, ONLY: Xcm, Ycm, Zcm
    use Simulation_data, ONLY: sim_xCenter, sim_yCenter, sim_zCenter
    implicit none
#include "Flash.h"
#include "constants.h"
 
    include "Flash_mpi.h"
    real :: dr32
    real, intent(IN)  :: ptx, pty, ptz
    real, intent(OUT) :: accelx, accely, accelz
    character(len=256) :: str_buffer
    
    !write(str_buffer, *) 'COM: ', Xcm, Ycm, Zcm, ptx, pty, ptz
    !call Logfile_stamp(gr_myPE, str_buffer, 'Note')
    dr32 = ((sim_xCenter - ptx)**2 + (sim_yCenter - pty)**2 + (Zcm - ptz)**2)**(3./2.)
    accelx = grv_factor*(sim_xCenter - ptx)/dr32
    accely = grv_factor*(sim_yCenter - pty)/dr32
    accelz = grv_factor*(Zcm - ptz)/dr32

    !dr32 = ((Xcm - ptx)**2 + (Ycm - pty)**2 + (Zcm - ptz)**2)**(3./2.)
    !accelx = grv_factor*(Xcm - ptx)/dr32
    !accely = grv_factor*(Ycm - pty)/dr32
    !accelz = grv_factor*(Zcm - ptz)/dr32
end
