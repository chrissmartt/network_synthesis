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
! SUBROUTINE RL_test(H_PR,type,CFtype,R,L,C,found,HR,remainder_OK,remainder_zero)
!
!
! NAME
!     RL_test
!
! DESCRIPTION
!       look for a viable RL branch in a given impedance/admittance function
!       See sections 7.2.2 and 7.2.3 of the Theory manual
!
! SEE ALSO
!
!
! HISTORY
!
!     started 14/09/17 CJS
!
  SUBROUTINE RL_test(H_PR,type,CFtype,R,L,C,found,HR,remainder_OK,remainder_zero)

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
  
  real(dp) :: const_term
  
  logical  :: negative_residue,non_zero_pole,const_term_big_enough
  
! function types
  logical :: conjugate_pair
  logical :: imaginary_pair
  logical :: complex_pair

!START

  if (verbose) write(*,*)'CALLED: RL_test'

  found=.FALSE.
  
! loop over real poles
  do i=1,H_PR%n_real_poles
  
    pole=i
    
! test for whether we have an RL branch here...
    
    negative_residue=(dble(H_PR%residues(pole)).LT.0d0)
    non_zero_pole=(abs(H_PR%poles(pole)).GT.zero_test_small)
    
! The RL branch gives a contribution to the constant term, we need to check whether
! the constant in the filter function is big enough to include this
    const_term_big_enough=.FALSE.
    const_term=0d0
    if (negative_residue.AND.non_zero_pole) then
      const_term=dble(H_PR%residues(pole))/dble(H_PR%poles(pole))
      const_term_big_enough=H_PR%R.GE.const_term
    end if
    
    if (verbose) then
      write(*,*)'Testing pole ',i
      write(*,*)'Tests:'
      write(*,*)'negative residue test  : ',negative_residue,' ',H_PR%residues(pole)
      write(*,*)'non zero pole test     : ',non_zero_pole,' ',H_PR%poles(pole)
      write(*,*)'const_term big enough test     : ',const_term_big_enough,' ',const_term
    end if
    
! check for stable poles which are not on the imaginary s=jw axis AND positive residues
    if (negative_residue.AND.non_zero_pole.AND.const_term_big_enough) then
    
      if (verbose) write(*,*)'Found possible RL branch'

! this could be a viable RL branch - calculate the remainder when this pole is removed

      CALL deallocate_Sfilter(HR)
      
! Test whether the remainder is zero
      
! build a local pole-residue filter without the test pole

! allocate the structure for the local pole-residue form function and copy the 
! required information across            
      HR_PR_local%wnorm=H_PR%wnorm
      HR_PR_local%order=H_PR%order-1
      HR_PR_local%n_real_poles=H_PR%n_real_poles-1
      HR_PR_local%n_complex_poles=H_PR%n_complex_poles
      HR_PR_local%n_complex_pole_pairs=H_PR%n_complex_pole_pairs
      HR_PR_local%n_real_poles=H_PR%n_real_poles

! constant term and sL term
      
      HR_PR_local%R=H_PR%R-const_term
      HR_PR_local%L=H_PR%L
      
! Test whether the remainder is zero
      if (     (HR_PR_local%order.EQ.0).AND.              &
           (abs(HR_PR_local%R).LT.zero_test_small).AND.   &
           (abs(HR_PR_local%L).LT.zero_test_small) )        then
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
        if (ii.NE.i) then
          pole1=pole1+1
          pole2=pole2+1
          HR_PR_local%complex_pole(pole1)=.FALSE.
          HR_PR_local%poles(pole1)   =H_PR%poles(pole2)
          HR_PR_local%residues(pole1)=H_PR%residues(pole2)
        else
! this is the pole we wish to remove so just increase the pole2 counter by 2. 
          pole2=pole2+1
        end if
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
    
    end if ! positive residue for this pole 
  
  end do ! next real pole
  
! we only get here if we have not found a viable RL branch

  remainder_OK=.FALSE.
  found=.FALSE.
  remainder_zero=.FALSE.
  
  RETURN
  
8000 CONTINUE
! jump here if we have found a viable RRL branch

  remainder_OK=.TRUE.
  found=.TRUE.
  CALL deallocate_Sfilter_PR(HR_PR_local)
  
  pole=i
  
  if (type.EQ.type_impedance) then
    CFtype=series_RL
    C=0d0
    R=dble(H_PR%residues(pole))/dble(H_PR%poles(pole))
    L=-dble(H_PR%residues(pole))/dble(H_PR%poles(pole)*H_PR%poles(pole))
    L=L/H_PR%wnorm
    if (verbose) then
      write(*,*)'FOUND VIABLE SERIES RL BRANCH'
      write(*,*)'R=',R
      write(*,*)'L=',L
      write(*,*)'remainder_OK   :',remainder_OK
      write(*,*)'remainder_zero :',remainder_zero
    end if
  else
    CFtype=shunt_RC
    L=0d0
    R=dble(H_PR%poles(pole))/dble(H_PR%residues(pole))
    C=-dble(H_PR%residues(pole))/dble(H_PR%poles(pole)*H_PR%poles(pole))
    C=C/H_PR%wnorm
    if (verbose) then
      write(*,*)'FOUND VIABLE SHUNT RC BRANCH'
      write(*,*)'R=',R
      write(*,*)'C=',C
      write(*,*)'remainder_OK   :',remainder_OK
      write(*,*)'remainder_zero :',remainder_zero
    end if
  end if
  
  RETURN

  END SUBROUTINE RL_test
