!
! This file is part of SACAMOS, State of the Art CAble MOdels for Spice. 
! It was developed by the University of Nottingham and the Netherlands Aerospace 
! Centre (NLR) for ESA under contract number 4000112765/14/NL/HK.
! 
! Copyright (C) 2016-2018 University of Nottingham
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
!
!
! SUBROUTINE R2_test(H_PR,type,CFtype,R,L,C,found,HR,remainder_OK,remainder_zero)
!
!
! NAME
!     R2_test
!
! DESCRIPTION
!       look for a viable R2 branch in a given impedance/admittance function
!       See sections 7.2.2 and 7.2.3 of the Theory manual
!
! SEE ALSO
!
!
! HISTORY
!
!     started 14/09/17 CJS
!
  SUBROUTINE R2_test(H_PR,type,CFtype,R,L,C,found,HR,remainder_OK,remainder_zero)

USE type_specifications
USE general_module
USE constants
USE filter_module

IMPLICIT NONE
  
  type(Sfilter_PR),INTENT(IN) :: H_PR
  integer :: type
  integer :: CFtype
  real(dp):: R,L,C
  logical :: found
  type(Sfilter),INTENT(INOUT) :: HR
  logical :: remainder_OK,remainder_zero

! local variables
  
  real(dp) :: Rdc
  integer :: pole,pole1,pole2
  type(Sfilter_PR) :: HR_PR_local
  type(Sfilter) :: Rdc_filter
  type(Sfilter) :: temp_filter
  logical :: stable
  
  integer :: i,ii
  
  logical  :: positive_R

!START

  if (verbose) write(*,*)'CALLED: R2_test'

  found=.FALSE.
    
! test for whether we have a R2 branch here...
  
! convert to a rational function form  
  HR_PR_local=H_PR
  HR=Convert_filter_S_PR_to_S(HR_PR_local)
  
  if ( (abs(HR%b%coeff(0)).LT.zero_test_small).OR.(abs(HR%a%coeff(0)).LT.zero_test_small) ) then
  
    remainder_OK=.FALSE.
    found=.FALSE.
    remainder_zero=.FALSE.
    RETURN

  end if
  
  Rdc=HR%a%coeff(0)/HR%b%coeff(0)

  positive_R=(Rdc.GT.0d0)
   
  if (verbose) then
    write(*,*)'Tests:'
    write(*,*)'positive R test  : ',positive_R,' ',H_PR%R
  end if
    
! check for stable poles which are not on the imaginary s=jw axis AND positive residues
  if (positive_R) then
    
    if (verbose) write(*,*)'Found possible R2 branch'

! this could be a viable R branch - calculate the remainder when this term is removed

      Rdc_filter=allocate_Sfilter(0,0)
      Rdc_filter%wnorm=HR%wnorm
      Rdc_filter%a%coeff(0)=-Rdc
      Rdc_filter%b%coeff(0)=1d0
      
      temp_filter=Hr+Rdc_filter
      CALL deallocate_Sfilter(HR)
      HR=temp_filter
      
! Check the transfer funcion for stability and for whether it is positive real

      if (verbose) then
        CALL write_Sfilter(HR,0)
      end if

      CALL check_transfer_function(HR,stable)
      
      CALL deallocate_Sfilter_PR(HR_PR_local)
      
      if (verbose) then
        if (stable) then
          write(*,*)'Remainder is stable'
        else
          write(*,*)'Remainder is unstable'        
        end if
      end if  
      
      if (stable) then
        remainder_zero=.FALSE.
        GOTO 8000
      end if
    
! Test whether the remainder is positive real    
    
    end if ! positive R2
   
! we only get here if we have not found a viable R2 branch

  remainder_OK=.FALSE.
  found=.FALSE.
  remainder_zero=.FALSE.
  
  RETURN
  
8000 CONTINUE
! jump here if we have found a viable R2 branch

  remainder_OK=.TRUE.
  found=.TRUE.
  CALL deallocate_Sfilter_PR(HR_PR_local)
  CALL deallocate_Sfilter(Rdc_filter)
  CALL deallocate_Sfilter(temp_filter)
  
  pole=i
  
  if (type.EQ.type_impedance) then
    CFtype=series_R
    R=Rdc
    L=0d0
    C=0d0
    if (verbose) then
      write(*,*)'FOUND VIABLE SERIES R BRANCH'
      write(*,*)'R=',R
      write(*,*)'remainder_OK   :',remainder_OK
      write(*,*)'remainder_zero :',remainder_zero
    end if
  else
    CFtype=shunt_R
    R=1d0/Rdc
    C=0d0
    L=0d0
    if (verbose) then
      write(*,*)'FOUND VIABLE SHUNT R BRANCH'
      write(*,*)'R=',R
      write(*,*)'remainder_OK   :',remainder_OK
      write(*,*)'remainder_zero :',remainder_zero
    end if
  end if
  
  RETURN

  END SUBROUTINE R2_test
