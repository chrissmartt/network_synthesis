!
! This file is part of SACAMOS, State of the Art CAble MOdels in Spice. 
! It was developed by the University of Nottingham and the Netherlands Aerospace 
! Centre (NLR) for ESA under contract number 4000112765/14/NL/HK.
! 
! Copyright (C) 2016-2017 University of Nottingham
! 
! SACAMOS is free software: you can redistribute it and/or modify it under the 
! terms of the GNU General Public License as published by the Free Software 
! Foundation, either version 3 of the License, or (at your option) any later 
! version.
! 
! SACAMOS is distributed in the hope that it will be useful, but 
! WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
! or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License 
! for more details.
! 
! A copy of the GNU General Public License version 3 can be found in the 
! file GNU_GPL_v3 in the root or at <http://www.gnu.org/licenses/>.
! 
! SACAMOS uses the EISPACK library (in /SRC/EISPACK). EISPACK is subject to 
! the GNU Lesser General Public License. A copy of the GNU Lesser General Public 
! License version can be found in the file GNU_LGPL in the root of EISPACK 
! (/SRC/EISPACK ) or at <http://www.gnu.org/licenses/>.
! 
! The University of Nottingham can be contacted at: ggiemr@nottingham.ac.uk
!
! File Contents:
!
! SUBROUTINE read_version
!
!
! DESCRIPTION
!
! read the version tag for the code,the version date and the compilation date
! 
!     
! COMMENTS
!
!
! HISTORY
!
!     started 12/1/2016 CJS 
!
  SUBROUTINE read_version()

  IMPLICIT NONE

! START

  include '../../version_information.inc'  ! this file specifies the version number
  
  include '../compilation_date.inc'        ! this file specifies the compilation date
  
  write(*,*)'Version :',trim(NETWORK_SYNTHESIS_version)
  write(*,*)'Date :',trim(NETWORK_SYNTHESIS_date)
  write(*,*)'Compilation date:',trim(NETWORK_SYNTHESIS_compilation_date)
  
  END SUBROUTINE read_version

 
