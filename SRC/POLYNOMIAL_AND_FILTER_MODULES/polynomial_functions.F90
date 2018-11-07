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
! NAME
!     polynomial_functions
!
! DESCRIPTION
!     definitions of functions operating on polynomial types 
!     to be included in filter_module.F90
!
! HISTORY
!
!     started 23/01/09 CJS
!
! FILE CONTAINS:
! FUNCTION allocate_polynomial
! FUNCTION allocate_complex_polynomial
! FUNCTION evaluate_polynomial
! FUNCTION evaluate_complex_polynomial
! FUNCTION polynomial_is_zero
!

FUNCTION allocate_polynomial(order) RESULT(res)

IMPLICIT NONE

! argument types
  integer,intent(IN)           :: order
! Result type
  type(polynomial)             :: res
  
!local variable
  type(polynomial) :: ans

! function definition

  ans%order=order
  allocate (ans%coeff(0:ans%order))
  ans%coeff(0:ans%order)=0d0
  
  res=ans
  
END FUNCTION allocate_polynomial

!
! NAME
!     allocate_complex_polynomial
!
! DESCRIPTION
!     set up a complex polynomial type with a given order
!
! HISTORY
!
!     started 23/01/09 CJS
!

FUNCTION allocate_complex_polynomial(order) RESULT(res)

USE type_specifications

IMPLICIT NONE

! argument types
  integer,intent(IN)                   :: order
! Result type
  type(complex_polynomial)             :: res
  
! local varialbe
  type(complex_polynomial) :: ans

! function definition

  ans%order=order
  allocate (ans%coeff(0:ans%order))
  ans%coeff(0:ans%order)=(0d0,0d0)
  
  res=ans
  
END FUNCTION allocate_complex_polynomial

!
! NAME
!     evaluate_polynomial
!
! DESCRIPTION
!     evaluate a polynomial with real coefficients and complex argument
!
! HISTORY
!
!     started 26/01/09 CJS
!

FUNCTION evaluate_polynomial(a,s) RESULT(res)

USE type_specifications

IMPLICIT NONE
  
  type(polynomial), intent(IN) :: a
  complex(dp)     , intent(IN) :: s
  
  complex(dp)                  :: res

! local variables  
  complex(dp) :: sn
  integer     :: n

!START

  sn=(1d0,0d0)
  res=(0d0,0d0)
  do n=0,a%order
    res=res+a%coeff(n)*sn
    sn=sn*s
  end do
 
  return
  end function evaluate_polynomial
!
! NAME
!     evaluate_complex_polynomial
!
! DESCRIPTION
!     evaluate a complex polynomial with complex coefficients and complex argument
!
! HISTORY
!
!     started 26/01/09 CJS
!

FUNCTION evaluate_complex_polynomial(a,s) RESULT(res)

USE type_specifications
 
IMPLICIT NONE
 
  type(complex_polynomial), intent(IN) :: a
  complex(dp)     , intent(IN)         :: s
  
  complex(dp)                          :: res

! local variables  
  complex(dp) :: sn
  integer     :: n

!START

  sn=(1d0,0d0)
  res=(0d0,0d0)
  do n=0,a%order
    res=res+a%coeff(n)*sn
    sn=sn*s
  end do
 
  return
  
  END FUNCTION evaluate_complex_polynomial
!
! NAME
!     polynomial_is_zero
!
! DESCRIPTION
!     return TRUE if the polynomial is zero i.e. all its coefficients are zero
!
! HISTORY
!
!     started 26/01/09 CJS
!

FUNCTION polynomial_is_zero(a) RESULT(res)

USE type_specifications
USE constants

IMPLICIT NONE
  
  type(polynomial), intent(IN) :: a
  
  logical                 :: res

! local varaibles
  integer :: i

!START

  res=.TRUE.
    
  do i=0,a%order   
    
    if (abs(a%coeff(i)).GT.zero_test_small) then
      res=.FALSE.
      RETURN      
    end if
    
  end do

  return
  end function polynomial_is_zero
