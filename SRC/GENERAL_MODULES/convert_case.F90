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
! NAME
!    convert_to_lower_case
!
! DESCRIPTION
!     returns the input string with all characters converted to lower case
!
! HISTORY
!
!     started 10/02/09 CJS
!
! COMMENTS
!     
SUBROUTINE convert_to_lower_case(input_string,input_string_length)

IMPLICIT NONE

! variables passed to subroutine

  integer,intent(IN)                               :: input_string_length
  character(LEN=input_string_length),intent(INOUT) :: input_string
  
! local variables

  character(LEN=input_string_length) ::  output_string   ! local output string
  
  integer :: i                   ! loop variable
  integer :: na,nCA,nz,nCZ,diff  ! character numbers of lower case a, upper case A, lower case z, upper case Z
  integer :: c                   ! character number
  
! START

  na=IACHAR('a')
  nCA=IACHAR('A')
  nz=IACHAR('z')
  nCZ=IACHAR('Z')
  
  diff=na-nCA   ! diff is the difference in character number between lower case and upper case i.e. what has to be added to convert to lower case
  
  do i=1,input_string_length
  
    c=IACHAR(input_string(i:i))
    if ((c.ge.nCA).AND.(c.le.nCZ)) then    ! this is an upper case character which needs conversions
      c=c+diff
      output_string(i:i)=ACHAR(c)
    else
      output_string(i:i)=input_string(i:i) ! this is not an upper case so leave this alone
    end if 
    
  end do

! copy the output string to the input string
  input_string=output_string
  
  RETURN
  
END SUBROUTINE convert_to_lower_case
!
! NAME
!    convert_to_upper_case
!
! DESCRIPTION
!     returns the input string with all characters converted to upper case
!
! HISTORY
!
!     started 26/02/09 CJS
!
! COMMENTS
!     
SUBROUTINE convert_to_upper_case(input_string,input_string_length)

IMPLICIT NONE

! variables passed to subroutine

  integer,intent(IN)                               :: input_string_length
  character(LEN=input_string_length),intent(INOUT) :: input_string
  
! local variables

  character(LEN=input_string_length) ::  output_string   ! local output string
  
  integer :: i                   ! loop variable
  integer :: na,nCA,nz,nCZ,diff  ! character numbers of lower case a, upper case A, lower case z, upper case Z
  integer :: c                   ! character number
  
! START

  na=IACHAR('a')
  nCA=IACHAR('A')
  nz=IACHAR('z')
  nCZ=IACHAR('Z')
  
  diff=nCA-na ! diff is the difference in character number between lower case and upper case i.e. what has to be added to convert to upper case
  
  do i=1,input_string_length
  
    c=IACHAR(input_string(i:i))
    if ((c.ge.na).AND.(c.le.nz)) then      ! this is a lower case character which needs conversions
      c=c+diff
      output_string(i:i)=ACHAR(c)
    else
      output_string(i:i)=input_string(i:i) ! this is not a lower case so leave this alone
    end if
    
  end do
  
! copy the output string to the input string
  input_string=output_string
  
  RETURN
  
END SUBROUTINE convert_to_upper_case
