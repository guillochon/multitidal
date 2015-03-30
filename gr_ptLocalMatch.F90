!!****if* source/Grid/GridParticles/GridParticlesMove/gr_ptLocalMatch
!!
!! NAME
!!
!!  gr_ptLocalMatch
!!
!! SYNOPSIS
!!
!!  gr_ptLocalMatch( real,(inout)    :: dataBuf(:,:),
!!                 integer(inout)  :: localCount,
!!                 integer(in)     :: propCount,
!!                 integer (in)    :: maxCount,
!!                 real,(inout)    :: sourceBuf,
!!                 integer,(IN)    :: numSource,
!!                 integer,(OUT)   :: destBuf,
!!                 integer,(OUT)   :: numDest)
!!
!! DESCRIPTION
!!     
!!      In the algorithm to move non stationaly data elements to their correct destination, elements are
!!      passed around to different processors, and each processor examines them and keeps
!!      the ones belonging to itself. It passes on the elements that don't belong to it.
!!      In this routine, the elements that haven't yet found their home are kept in 
!!      "sourceBuf", the elements in their right place are kept in "elements", and 
!!      elements to be passed on are kept in "destBuf". This routine examines all
!!      elements in sourceBuf to find if any of them belong in the section of physical
!!      domain residing on myPE. If they do, they are moved to elements, if not they
!!      are moved to destBuf to be passed on to the next processor.
!!      
!!      This routine also handles the situation if the local number of elements exceeds the 
!!      maximum allowed. If the runtime parameter "grpt_rmElements" is turned on, excess
!!      elements are removed from the simulation based upon the the implementation in
!!      gr_removeElements, otherwise the code aborts. 
!!
!! ARGUMENTS
!!
!!     dataBuf -           A 2 dimensional array of elements and their property
!!                           values
!!     localCount -   The number of elements in the particle array that
!!                           are stored on this processor.
!!     propCount : number of particle attributes
!!     maxCount - The maximum number of elements that can be stored
!!                           on a single processor.
!!     sourceBuf -           temporary storage for elements that may have arrived from 
!!                           another processor
!!     numSource -           number of elements in sourceBuf
!!     destBuf   -           temporary storage for elements that need to move off processor
!!     numDest   -           number of elements in destBuf
!! 
!!
!!***

subroutine gr_ptLocalMatch(dataBuf,localCount,propCount,maxCount,&
     sourceBuf,numSource,destBuf,numDest)

#include "constants.h"
#include "Flash.h"
  use gr_ptData, ONLY : gr_ptBlkCount, gr_ptBlkList,&
       gr_ptBlk,gr_ptPosx,gr_ptPosy,gr_ptPosz, gr_ptProc, gr_ptKeepLostParticles
  use Grid_data, ONLY : gr_meshMe, gr_globalDomain
  use Grid_interface, ONLY :  Grid_outsideBoundBox
  use gr_ptInterface, ONLY : gr_ptHandleExcess
  use gr_interface, ONLY:  gr_findBlock,gr_xyzToBlock

  implicit none

  integer, intent(IN) :: propCount
  integer, intent(IN) :: numSource
  integer, intent(OUT) :: numDest
  integer, intent(INOUT) :: localCount
  integer, intent(IN) :: maxCount
  real,dimension(propCount,maxCount),intent(INOUT) :: dataBuf
  real,dimension(propCount,numSource),intent(INOUT) :: sourceBuf,destBuf
  
  integer ::   blockID,procID,oblockID
  integer :: i,m
  real,dimension(MDIM)::pos
  integer, dimension(MDIM) :: negh, pxyz
  logical :: outside
  
  numDest=0
  m=localCount
  pxyz = (/gr_ptPosx,gr_ptPosy,gr_ptPosz/)
  
#ifdef BITTREE
  do i =1,numSource
     pos(1:NDIM) = sourceBuf(pxyz(1:NDIM),i)
     call gr_xyzToBlock(pos,procID,blockID)
     if(procID==gr_meshMe)then
        m=m+1
        databuf(:,m)=sourceBuf(:,i)
        databuf(gr_ptBlk,m)=blockID
        databuf(gr_ptProc,m)=procID
     else
        numDest=numDest+1
        destbuf(:,numDest)=sourceBuf(:,i)
        destbuf(gr_ptBlk,numDest)=blockID
        destbuf(gr_ptProc,numDest)=procID
     end if
  end do
#else
  do i = 1,numSource
     !! find if the elements belongs to any of the local blocks
     blockID=int(sourceBuf(gr_ptBlk,i))
     if(blockID/=UNKNOWN) then
        procID = int(sourceBuf(gr_ptProc,i))
        if(procID/=gr_meshMe)blockID=NONEXISTENT
     else
        pos(1:NDIM) = sourceBuf(pxyz(1:NDIM),i)
        oblockID = blockID
        call gr_findBlock(gr_ptBlkList,gr_ptBlkCount,pos,blockID)
        if(gr_ptKeepLostParticles) then
           call Grid_outsideBoundBox(pos,gr_globalDomain,outside,Negh)
           !! JFG
           if(outside)blockID=oblockID ! Move to block 1, which should always exist?
        end if
     end if
#ifdef DEBUG_ELEMENTS
993  format(I5,' unknown on ',I2,':(',3G15.8,')',I6)
     print 993,i,gr_meshMe,pos,blockID
#endif
     
     if(blockID==NONEXISTENT) then !! if not then move it to destBuf
        numDest=numDest+1
        destBuf(:,numDest)=sourceBuf(:,i)
        
     else              !! If yes, make sure local number
                       !! of elements is not more than maximum.
        if(m>=maxCount)call gr_ptHandleExcess(dataBuf,propCount,m,maxCount)
        m=m+1   
        dataBuf(:,m)=sourceBuf(:,i)         !! if particle belongs locally, and the
        dataBuf(gr_ptBlk,m) = blockID  !! local maximum is not exceeded, move
        dataBuf(gr_ptProc,m)=gr_meshMe
        !! particle to the primary data structure
     end if
  end do
#endif
  localCount=m   !! this is the new local number of elements on my PE.
end subroutine gr_ptLocalMatch
