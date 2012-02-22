subroutine Mass_Loss_Correction
    use Gravity_data, ONLY: grv_obvec, grv_exactvec, grv_oexactvec
    use Driver_data, ONLY: dr_dt

    implicit none 
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"

    grv_obvec(1:3) = grv_obvec(1:3) - (grv_exactvec(1:3) - grv_oexactvec(1:3) - &
        0.5d0*dr_dt*(grv_exactvec(4:6) + grv_oexactvec(4:6)))
end subroutine
