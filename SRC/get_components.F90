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
! FILE CONTENTS:
!
! SUBROUTINE series_s_term(num,den,value,num_out,remainder_OK,remainder_zero)
! SUBROUTINE series_const_term(num,den,value,num_out,remainder_OK,remainder_zero)
! SUBROUTINE series_pole_at_zero_term(num,den,value,num_out,den_out,remainder_OK,remainder_zero)
!
!
! NAME
!     series_s_term
!
! DESCRIPTION
!       
!
! SEE ALSO
!
!
! HISTORY
!
!     started 25/08/09 CJS
!
SUBROUTINE series_s_term(num,den,value,num_out,remainder_OK,remainder_zero)

! do the division num/den to get the coefficient of s at this stage 

USE type_specifications
USE general_module
USE constants
USE frequency_spec
USE filter_module
USE Sfilter_fit_module
 
IMPLICIT NONE
 
  type(polynomial),intent(IN):: num
  type(polynomial),intent(IN):: den
  real(dp),intent(OUT) :: value
  type(polynomial),intent(INOUT):: num_out
  logical,intent(OUT) :: remainder_OK
  logical,intent(OUT) :: remainder_zero

! local variables

  integer i

!START

  value=num%coeff(num%order)/den%coeff(den%order)
  
  if(verbose) write(*,*)'**** Test_s_term, value=',value
  
  if(verbose) CALL write_polyrat_local(num,den)
   
  remainder_zero=.TRUE.
  
  do i=num_out%order,0,-1
    if (i.GT.0) then
      num_out%coeff(i)=num%coeff(i)-value*den%coeff(i-1)
    else
      num_out%coeff(i)=num%coeff(i)
    end if
    
    if (abs(num_out%coeff(i)).GT.zero_test_small) remainder_zero=.FALSE.
    
  end do
    
  if (verbose) then
    write(*,*)'Remainder:'
    CALL write_poly_local(num_out)
  end if
  
  CALL get_min_order_poly(num_out)
    
  if (verbose) then
    write(*,*)'Min order remainder:'
    CALL write_poly_local(num_out)
  end if
  
! Check whether the remainder rational function is a viable impedance/admittance function

  CALL check_transfer_function_num_den(num_out,den,remainder_OK)

  if((.NOT.remainder_OK).AND.(verbose)) write(*,*)'**** THE REMAINDER IS NOT OK *****'
  
END SUBROUTINE series_s_term
!
! NAME
!     series_const_term
!
! DESCRIPTION
!       
!
! SEE ALSO
!
!
! HISTORY
!
!     started 25/08/09 CJS
!
SUBROUTINE series_const_term(num,den,value,num_out,remainder_OK,remainder_zero)

! do the division num/den to get the constant coefficient at this stage 

USE type_specifications
USE general_module
USE constants
USE frequency_spec
USE filter_module
USE Sfilter_fit_module
 
IMPLICIT NONE
 
  type(polynomial),intent(IN):: num
  type(polynomial),intent(IN):: den
  real(dp),intent(OUT) :: value
  type(polynomial),intent(INOUT):: num_out
  logical,intent(OUT) :: remainder_OK
  logical,intent(OUT) :: remainder_zero

! local variables

  integer i

!START

  value=num%coeff(num%order)/den%coeff(den%order)
  
  if(verbose) write(*,*)'**** Test_const_term_1, value=',value
  
  if(verbose) CALL write_polyrat_local(num,den)
   
  remainder_zero=.TRUE.

! remainder is zero
  if (num%order.EQ.0) then
    remainder_OK=.TRUE.
    RETURN
  end if
  
  do i=num_out%order,0,-1

    num_out%coeff(i)=num%coeff(i)-value*den%coeff(i)
    
    if (abs(num_out%coeff(i)).GT.zero_test_small) remainder_zero=.FALSE.
    
  end do
    
  if (verbose) then
    write(*,*)'Remainder:'
    CALL write_poly_local(num_out)
  end if
  
  CALL get_min_order_poly(num_out)
    
  if (verbose) then
    write(*,*)'Min order remainder:'
    CALL write_poly_local(num_out)
  end if
  
! Check whether the remainder rational function is a viable impedance/admittance function

  CALL check_transfer_function_num_den(num_out,den,remainder_OK)

  if((.NOT.remainder_OK).AND.(verbose)) write(*,*)'**** THE REMAINDER IS NOT OK *****'
  
END SUBROUTINE series_const_term
!
! NAME
!     series_const_term_2
!
! DESCRIPTION
!       
!
! SEE ALSO
!
!
! HISTORY
!
!     started 12/09/09 CJS
!
SUBROUTINE series_const_term_2(num,den,value,num_out,remainder_OK,remainder_zero)

! do the division num/den to get the constant coefficient at this stage 

USE type_specifications
USE general_module
USE constants
USE frequency_spec
USE filter_module
USE Sfilter_fit_module
 
IMPLICIT NONE
 
  type(polynomial),intent(IN):: num
  type(polynomial),intent(IN):: den
  real(dp),intent(OUT) :: value
  type(polynomial),intent(INOUT):: num_out
  logical,intent(OUT) :: remainder_OK
  logical,intent(OUT) :: remainder_zero

! local variables

  integer i

!START

  value=num%coeff(0)/den%coeff(0)
  
  if(verbose) write(*,*)'**** Test_const_term_2, value=',value
  
  if(verbose) CALL write_polyrat_local(num,den)
   
  remainder_zero=.TRUE.

! remainder is zero
  if (num%order.EQ.0) then
    remainder_OK=.TRUE.
    RETURN
  end if
  
  do i=num_out%order,0,-1

    num_out%coeff(i)=num%coeff(i)-value*den%coeff(i)
    
    if (abs(num_out%coeff(i)).GT.zero_test_small) remainder_zero=.FALSE.
    
  end do
    
  if (verbose) then
    write(*,*)'Remainder:'
    CALL write_poly_local(num_out)
  end if
  
  CALL get_min_order_poly(num_out)
    
  if (verbose) then
    write(*,*)'Min order remainder:'
    CALL write_poly_local(num_out)
  end if
  
! Check whether the remainder rational function is a viable impedance/admittance function

  CALL check_transfer_function_num_den(num_out,den,remainder_OK)

  if((.NOT.remainder_OK).AND.(verbose)) write(*,*)'**** THE REMAINDER IS NOT OK *****'
  
END SUBROUTINE series_const_term_2
!
! NAME
!     series_pole_at_zero_term
!
! DESCRIPTION
!       
!
! SEE ALSO
!
!
! HISTORY
!
!     started 31/08/09 CJS
!
SUBROUTINE series_pole_at_zero_term(num,den,value,num_out,den_out,remainder_OK,remainder_zero)

! get the coefficient C in 1/sC term at this stage 

USE type_specifications
USE general_module
USE constants
USE frequency_spec
USE filter_module
USE Sfilter_fit_module
 
IMPLICIT NONE
 
  type(polynomial),intent(IN):: num
  type(polynomial),intent(IN):: den
  real(dp),intent(OUT) :: value
  type(polynomial),intent(INOUT):: num_out
  type(polynomial),intent(INOUT):: den_out
  logical,intent(OUT) :: remainder_OK
  logical,intent(OUT) :: remainder_zero

! local variables

  integer i

!START

  value=den%coeff(1)/num%coeff(0)
  
  if(verbose) write(*,*)'**** Test_s_term, value=',value
  
  if(verbose) CALL write_polyrat_local(num,den)
  
  if (den%order.EQ.1) then
    remainder_zero=.TRUE.
    RETURN
  end if
  
! Denominator terms

  do i=0,den_out%order
    den_out%coeff(i)=den%coeff(i+1)/value
  end do
      
  if (verbose) then
    write(*,*)'Remainder denominator:'
    CALL write_poly_local(den_out)
  end if
  
  CALL get_min_order_poly(den_out)
    
  if (verbose) then
    write(*,*)'Min order denominator remainder:'
    CALL write_poly_local(den_out)
  end if
 
! Numerator terms  
  remainder_zero=.TRUE.
  
  do i=0,num_out%order
  
    if (i.LT.num_out%order) then
      num_out%coeff(i)=(num%coeff(i+1)-den_out%coeff(i+1))/value
    else
      num_out%coeff(i)=num%coeff(i+1)/value    
    end if
    if (abs(num_out%coeff(i)).GT.zero_test_small) remainder_zero=.FALSE.
    
  end do
      
  if (verbose) then
    write(*,*)'Remainder numerator:'
    CALL write_poly_local(num_out)
  end if
  
  CALL get_min_order_poly(num_out)
    
  if (verbose) then
    write(*,*)'Min order remainder:'
    CALL write_poly_local(num_out)
  end if
  
! Check whether the remainder rational function is a viable impedance/admittance function

  CALL check_transfer_function_num_den(num_out,den,remainder_OK)

  if((.NOT.remainder_OK).AND.(verbose)) write(*,*)'**** THE REMAINDER IS NOT OK *****'
  
END SUBROUTINE series_pole_at_zero_term
