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
! SUBROUTINE R_test(H_PR,type,CFtype,R,L,C,found,HR,remainder_OK,remainder_zero) 
!
!
! NAME
!     R_test
!
! DESCRIPTION
!       look for a viable R branch in a given impedance/admittance function
!       See sections 7.2.2 and 7.2.3 of the Theory manual
!
! SEE ALSO
!
!
! HISTORY
!
!     started 14/09/17 CJS
!
  SUBROUTINE R_test(H_PR,type,CFtype,R,L,C,found,HR,remainder_OK,remainder_zero)

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
  
  integer :: pole,pole1,pole2
  type(Sfilter_PR) :: HR_PR_local
  logical :: stable
  
  integer :: i,ii
  
  logical  :: positive_R

!START

  if (verbose) write(*,*)'CALLED: R_test'

  found=.FALSE.
    
! test for whether we have a R branch here...
    
    positive_R=(dble(H_PR%R).GT.0d0)
   
    if (verbose) then
      write(*,*)'Testing pole ',i
      write(*,*)'Tests:'
      write(*,*)'positive R test  : ',positive_R,' ',H_PR%R
    end if
    
! check for stable poles which are not on the imaginary s=jw axis AND positive residues
    if (positive_R) then
    
      if (verbose) write(*,*)'Found possible R branch'

! this could be a viable L branch - calculate the remainder when this pole is removed

      CALL deallocate_Sfilter(HR)
      
! Test whether the remainder is zero
      
! build a local pole-residue filter without the test pole

! allocate the structure for the local pole-residue form function and copy the 
! required information across            
      HR_PR_local%wnorm=H_PR%wnorm
      HR_PR_local%order=H_PR%order
      HR_PR_local%n_real_poles=H_PR%n_real_poles
      HR_PR_local%n_complex_poles=H_PR%n_complex_poles
      HR_PR_local%n_complex_pole_pairs=H_PR%n_complex_pole_pairs
      HR_PR_local%n_real_poles=H_PR%n_real_poles

! constant term and sL term
      
      HR_PR_local%R=0d0
      HR_PR_local%L=H_PR%L
      
! Test whether the remainder is zero
     if ( (HR_PR_local%order.EQ.0).AND.              &
           (abs(HR_PR_local%R).LT.zero_test_R).AND.   &
           (abs(HR_PR_local%L).LT.zero_test_L) ) then
        remainder_zero=.TRUE.
        GOTO 8000    
      end if

! copy any poles/ residues in the remainder
      
      ALLOCATE( HR_PR_local%complex_pole(HR_PR_local%order) )
      ALLOCATE( HR_PR_local%poles(HR_PR_local%order) )
      ALLOCATE( HR_PR_local%residues(HR_PR_local%order) )

! copy real poles      
      pole1=0
      pole2=0
      do ii=1,HR_PR_local%n_real_poles
        pole1=pole1+1
        pole2=pole2+1
        HR_PR_local%complex_pole(pole1)=.FALSE.
        HR_PR_local%poles(pole1)   =H_PR%poles(pole2)
        HR_PR_local%residues(pole1)=H_PR%residues(pole2)
      end do ! next real pole

! copy complex poles in pairs 
      do ii=1,H_PR%n_complex_pole_pairs
      
        pole1=pole1+1
        pole2=pole2+1
        HR_PR_local%complex_pole(pole1)=.TRUE.
        HR_PR_local%poles(pole1)=H_PR%poles(pole2)
        HR_PR_local%residues(pole1)=H_PR%residues(pole2)
        pole1=pole1+1
        pole2=pole2+1
        HR_PR_local%complex_pole(pole1)=.TRUE.
        HR_PR_local%poles(pole1)=H_PR%poles(pole2)
        HR_PR_local%residues(pole1)=H_PR%residues(pole2)
        
      end do ! next real pole

! convert to a rational function form      
      HR=Convert_filter_S_PR_to_S(HR_PR_local)

! Check the transfer funcion for stability and for whether it is positive real

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
    
    end if ! positive R
   
! we only get here if we have not found a viable C branch

  remainder_OK=.FALSE.
  found=.FALSE.
  remainder_zero=.FALSE.
  
  RETURN
  
8000 CONTINUE
! jump here if we have found a viable RC branch

  remainder_OK=.TRUE.
  found=.TRUE.
  CALL deallocate_Sfilter_PR(HR_PR_local)
  
  pole=i
  
  if (type.EQ.type_impedance) then
    CFtype=series_R
    R=H_PR%R
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
    R=1d0/H_PR%R
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

  END SUBROUTINE R_test
