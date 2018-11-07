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
!
! FILE CONTAINS:
! FUNCTION conjugate_pair
! FUNCTION imaginary_pair
! FUNCTION complex_pair
!
!
! NAME
!     conjugate_pair
!
! DESCRIPTION
!      Test whether the two complex numbers a and b are a complex conjugate pair
!
! HISTORY
!
!     started 15/09/2017 CJS
!

FUNCTION conjugate_pair(a,b) RESULT(res)

USE type_specifications
USE constants

IMPLICIT NONE

! argument types
  complex(dp),intent(IN) :: a,b
! Result type
  logical                :: res
  
!local variables

  complex(dp) :: diff
  real(dp) :: mag_diff,mag1,mag2

! function definition

! Test whether (a*-b)/(|a|+|b|) is zero (very small) 

  mag1=abs(a)
  mag2=abs(b)
  diff=conjg(a)-b     
  mag_diff=abs(diff)
  
  if ( (mag_diff/(mag1+mag2)).GT.zero_test_small) then
    res=.FALSE.
  else
    res=.TRUE.
  end if
  
END FUNCTION conjugate_pair
!
! NAME
!     imaginary_pair
!
! DESCRIPTION
!      Test whether the two complex numbers a and b are an imaginary pair
!      It is assumed that a and b are a complex conjugate pair
!
! HISTORY
!
!     started 15/09/2017 CJS
!

FUNCTION imaginary_pair(a,b) RESULT(res)

USE type_specifications
USE constants

IMPLICIT NONE

! argument types
  complex(dp),intent(IN) :: a,b
! Result type
  logical                :: res
  
! local variables

  complex(dp) :: sum
  real(dp) :: mag_sum,mag1,mag2

! function definition

! Test whether the real part of (a*+b)/(|a|+|b|) is zero (very small) 

  mag1=abs(a)
  mag2=abs(b)
  sum=conjg(a)+b     
  mag_sum=abs(real(sum))
  
  if ( (mag_sum/(mag1+mag2)).GT.zero_test_small) then
    res=.FALSE.
  else
    res=.TRUE.
  end if
  
END FUNCTION imaginary_pair
!
! NAME
!     complex_pair
!
! DESCRIPTION
!      Test whether the two complex numbers a and b are a complex pair
!      i.e. the real parts of a and b are not zero
!      It is assumed that a and b are a complex conjugate pair
!
! HISTORY
!
!     started 15/09/2017 CJS
!


FUNCTION complex_pair(a,b) RESULT(res)

USE type_specifications
USE constants

IMPLICIT NONE

! argument types
  complex(dp),intent(IN) :: a,b
! Result type
  logical                :: res
  
! local variables

  complex(dp) :: sum
  real(dp) :: mag_sum,mag1,mag2

! function definition

! Test whether the real part of (a*+b)/(|a|+|b|) is greater than zero (something very small) 

  mag1=abs(a)
  mag2=abs(b)
  sum=conjg(a)+b     
  mag_sum=abs(real(sum))
  
  if ( (mag_sum/(mag1+mag2)).LT.zero_test_small) then
    res=.FALSE.
  else
    res=.TRUE.
  end if
  
END FUNCTION complex_pair
