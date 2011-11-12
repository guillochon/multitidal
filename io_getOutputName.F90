!!****if* source/IO/IOMain/io_getOutputName
!!
!! NAME
!!  io_getOutputName
!!
!! SYNOPSIS
!!
!!  io_getOutputName(integer(in)         :: fileNum, 
!!                  character(len=*)(in) :: fileTypeStr,
!!                  character(len=*)(in) :: strname, 
!!                  character(len=MAX_STRING_LENGTH)(out) :: filename,
!!                  logical (in)         ::  forced)
!!
!! DESCRIPTION
!!
!!  gets the name of the output file for a checkpoint or plotfile
!!  it takes an integer filenum along with 2 strings describing the
!!  file, ie 'hdf5' and 'chk' or 'ncmpi' and 'part' and combines 
!!  it with the number in
!!  a string format 
!!
!! ARGUMENTS
!!
!!  fileNum - integer representing the checkpoint or plotfile number
!!  fileTypeStr - indicates type of file like 'hdf5' or 'ncmpi' 
!!  strname - string descriptor of file like, 'chk','plt' or 'part'
!!          
!!  filename - returned resulting string name given to the file    
!!  forced - if true, the file name will reflect that this file is considered
!!           a forced output file.  
!!
!!
!! EXAMPLE
!!
!!  call io_getOutputName(io_cpNumber, 'hdf5', '_chk', filename, .false.)
!!  
!!
!!***


subroutine io_getOutputName(fileNum, fileTypeStr, strname, filename, forced)

  use IO_data, ONLY : io_baseName, io_outputDir, io_splitFileNum, &
      io_outputSplitNum

  implicit none
  
#include "constants.h"


  integer, intent(in) :: fileNum
  character (len=MAX_STRING_LENGTH), intent(inout) :: filename
  character(len=*), intent(in) :: strname, fileTypeStr
  logical, intent(in) :: forced

  ! create a character variable to hold the string representation of the block
  ! number.  Note this is set to be 5 characters long (i.e. max = 9999).  
  character (len=5) ::  fnumStr, splitFileNumStr
  integer :: pos

  
  ! if io_outputDir isn't empty, and it doesn't have a 
  ! directory seperator at the end, add the directory seperator.
  pos = index(io_outputDir, ' ')
  if (pos.gt.1) then
     if (io_outputDir(pos-1:pos-1).ne."/") then
        io_outputDir(pos:pos) = "/"
     end if
  end if
  
  ! if io_outputDir is current directory specified 
  ! with a './', just get rid of it.
  if (io_outputDir == "./") then
     io_outputDir = ""
  end if

  filename = " "
  write (fnumStr, '(i5.5)') fileNum

  pos = index(io_baseName, ' ')
  !only if we are splitting files do we adjust the filename.
  !these also have to be adjusted if we are forcing output.
  if (io_outputSplitNum > 1) then
     write (splitFileNumStr, '(i5.5)') io_splitFileNum
     if(forced) then
        filename = io_baseName(:pos-1)// 'forced_' // 's'// splitFileNumStr // '_' // fileTypeStr // strname // fnumStr 
     else
        filename = io_baseName(:pos-1) // 's'// splitFileNumStr // '_' // fileTypeStr // strname // fnumStr 
     end if
  else
     
     if(forced) then
        filename = io_baseName(:pos-1)// 'forced_' // fileTypeStr // strname // fnumStr  
     else
        filename = io_baseName(:pos-1)// fileTypeStr // strname // fnumStr 
     end if
  end if

  pos = index(io_outputDir, ' ')
  filename = io_outputDir(:pos-1) // filename 


end subroutine io_getOutputName
