! This file forms part of the NETWORK_SYNTHESIS project
!
! Software to generate Spice sub-circuit models to reproduce
! impedance functions specified as either s-domain rational functions
! or pole-residue representations. 
!
! Copyright (C) 2018 University of Nottingham
!
! NETWORK_SYNTHESIS is free software: you can redistribute it and/or modify it under the 
! terms of the GNU General Public License as published by the Free Software 
! Foundation, either version 3 of the License, or (at your option) any later 
! version.
! 
! NETWORK_SYNTHESIS is distributed in the hope that it will be useful, but 
! WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
! or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License 
! for more details.
! 
! A copy of the GNU General Public License version 3 can be found in the 
! file COPYING.txt in the root or at <http://www.gnu.org/licenses/>.
! 
! NETWORK_SYNTHESIS uses the EISPACK library. EISPACK is subject to 
! the GNU Lesser General Public License. A copy of the GNU Lesser General Public 
! License version can be found in the file COPPYING.LESSER.txt 
! or at <http://www.gnu.org/licenses/>.
! 
! The University of Nottingham can be contacted at: ggiemr@nottingham.ac.uk
!
! Author C Smartt
!
!
! File Contents:
!
! SUBROUTINE write_program_status
!
! DESCRIPTION
!     
! write a program status message to the program_status file
!
! COMMENTS
! 
!
! HISTORY
!
!     started 20/11/2015 CJS 
!
!
  SUBROUTINE write_program_status()

  IMPLICIT NONE

! START

  OPEN(unit=program_status_file_unit,file=program_status_filename)
  
  write(program_status_file_unit,'(A,A,A)')trim(program_name),':',trim(run_status)
  write(*,'(A,A,A)')trim(program_name),':',trim(run_status)
  
  CLOSE(unit=program_status_file_unit)

  END SUBROUTINE write_program_status

 
