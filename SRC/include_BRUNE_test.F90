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
!
! SUBROUTINE BRUNE_test(H,type,CFtype,R,L1out,C,L2out,K,found,HR,remainder_OK,remainder_zero)
! SUBROUTINE LC_test_BRUNE(H_PR,L,C,found,HR,remainder_OK,remainder_zero)
! SUBROUTINE L_test_BRUNE(H_PR,L,found,HR,remainder_OK,remainder_zero)
!
! 
! NAME
!     BRUNE_test
!
! DESCRIPTION
!       look for a viable BRUNE branch in a given impedance function
!       See sections 7.2.4 of the Theory manual
!
! SEE ALSO
!
!
! HISTORY
!
!     started 14/09/17 CJS
!
  SUBROUTINE BRUNE_test(H,type,CFtype,R,L1out,C,L2out,K,found,HR,remainder_OK,remainder_zero)

USE type_specifications
USE general_module
USE constants
USE filter_module

IMPLICIT NONE
  
  type(Sfilter),INTENT(INOUT) :: H
  integer :: type
  integer :: CFtype
  real(dp):: R,L1out,L2out,C,K
  logical :: found
  type(Sfilter),INTENT(INOUT) :: HR
  logical :: remainder_OK,remainder_zero

! local variables
  
  real(dp) :: w0,f0
  complex(dp) :: Xw0
  
  type(Sfilter) :: Z2,Y2
  type(Sfilter_PR) :: Y2_PR
  type(Sfilter) :: Y3
  type(Sfilter) :: Z3
  type(Sfilter_PR) :: Z3_PR
  type(Sfilter) :: temp_filter
  real(dp):: L1,L2,L3,M
  

!START

  if (verbose) write(*,*)'CALLED: BRUNE_test'
  
  found=.FALSE.
  
! Calculate the resistance which reduces H_PR to a minimum resistive impedance function
! and the frequency at which the resistance is zero

  CALL calculate_min_resistance_function(H,R,w0)

! evaluate the min_resistance filter at w0. At this frequency the response is purely imaginary.
  f0=w0/(2d0*pi)
  

  f0=w0*H%wnorm/(2d0*pi)
  if (abs(f0).LT.zero_test_small) then
    w0=zero_test_small
    f0=w0*H%wnorm/(2d0*pi)
  end if
  Xw0=evaluate_Sfilter_frequency_response(H,f0)
  
  L1=aimag(Xw0)/w0
  if (verbose) then
    write(*,*)'L1=',L1
  end if

! subtract this inductive impedance from the filter function Z2(s)=H(s)-L1*s

  temp_filter=allocate_Sfilter(1,0)
  temp_filter%wnorm=H%wnorm
  temp_filter%a%coeff(0)=0d0
  temp_filter%a%coeff(1)=-L1
  temp_filter%b%coeff(0)=1d0
  
  Z2=AddSfilter(H,temp_filter)
  CALL get_min_order_poly(Z2%a)
  
  if (verbose) then
    write(*,*)'Z2 filter function:'
    CALL write_Sfilter(Z2,0)
  end if
  
  Y2=reciprocal_Sfilter(Z2)
  
  if (verbose) then
    write(*,*)'Y2 filter function:'
    CALL write_Sfilter(Y2,0)
  end if
  
! now extract an L-C series admittance from Y2

  Y2_PR=Convert_filter_S_to_S_PR(Y2)
  
  if (verbose) CALL write_S_PR_filter(Y2_PR)
  
  CALL LC_test_BRUNE(Y2_PR,L2,C,found,Y3,remainder_OK,remainder_zero)
               
  if (verbose) then
    write(*,*)'CALLED LC_test_BRUNE'
    write(*,*)'L2=',L2
    write(*,*)'C =',C
    write(*,*)'found =',found
    write(*,*)'remainder_OK =',remainder_OK
    write(*,*)'remainder_zero =',remainder_zero
    write(*,*)'remainder pole residue function'
    CALL write_Sfilter(Y3,0)
  end if
  
  Z3=reciprocal_Sfilter(Y3)
  
!  write(*,*)'Z3:'
!  CALL write_Sfilter(Z3,0)
  Z3_PR=Convert_filter_S_to_S_PR(Z3)

! extract a series inductance (which may be negative)  
  CALL L_test_BRUNE(Z3_PR,L3,found,HR,remainder_OK,remainder_zero)
  
  if (verbose) then
    write(*,*)'CALLED L_test_BRUNE'
    write(*,*)'L3=',L3
    write(*,*)'found =',found
    write(*,*)'remainder_OK =',remainder_OK
    write(*,*)'remainder_zero =',remainder_zero
    write(*,*)'remainder pole residue function'
  end if
  
! we can only proceed if the remainder is OK at this point

  if (found.AND.remainder_OK) then 
  
! carry on and work out the mutual impedance form which
! should eliminate the negative inductance
  
    L1out=L1+L2
    L2out=L3+L2
    M=L2
    K=1d0
    
    L1out=L1out/H%wnorm
    L2out=L2out/H%wnorm
    C=C/H%wnorm
    
! check for positive component values
    if (R.LT.0d0) then
      if (verbose) write(*,*)'Error, R<0,  R =',R
      found=.FALSE.
    end if
    if (L1out.LT.0d0) then
      if (verbose) write(*,*)'Error, L1<0, L1=',L1
      found=.FALSE.
    end if
    if (L2out.LT.0d0) then
      if (verbose) write(*,*)'Error, L2<0, L2=',L2
      found=.FALSE.
    end if
    if (K.LT.0d0) then
      if (verbose) write(*,*)'Error, K<0,  K =',K
      found=.FALSE.
    end if
    if (C.LT.0d0) then
      if (verbose) write(*,*)'Error, C<0,  C =',C
      found=.FALSE.
    end if
    
    if (found) GOTO 8000

  end if

  remainder_OK=.FALSE.
  found=.FALSE.
  remainder_zero=.FALSE.
  
  RETURN
  
8000 CONTINUE
! jump here if we have found a viable BRUNE branch
  
  CFtype=series_BRUNE
  found=.TRUE.
    
  if (verbose) then
    write(*,*)'FOUND VIABLE BRUNE BRANCH'
    write(*,*)'R =',R
    write(*,*)'L1=',L1out
    write(*,*)'L2=',L2out
    write(*,*)'C =',C
    write(*,*)'K=',K
    write(*,*)'remainder_OK   :',remainder_OK
    write(*,*)'remainder_zero :',remainder_zero
    write(*,*)'Remainder impedance:'
    CALL write_Sfilter(HR,0)
  end if
  
  RETURN

  END SUBROUTINE BRUNE_test
!
! NAME
!     LC_test_BRUNE
!
! DESCRIPTION
!       look for a viable LC branch in a given impedance/admittance function
!       This is specific to the Brune synthesis as the inductance can be negative
!       so the checks are different 
!
! SEE ALSO
!
!
! HISTORY
!
!     started 14/09/17 CJS
!
  SUBROUTINE LC_test_BRUNE(H_PR,L,C,found,HR,remainder_OK,remainder_zero)

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
  
  integer :: first_complex_pole,pole1,pole2,pole
  type(Sfilter_PR) :: HR_PR_local
  logical :: stable
  
  integer :: i,ii
  real(dp) :: diff,mag_diff,mag1,mag2
  
  real(dp) :: const_term
  logical  :: imaginary_poles,positive_residues,zero_const_term
  
!START

  if (verbose) write(*,*)'CALLED: LC_test_BRUNE'

  found=.FALSE.

  first_complex_pole=H_PR%n_real_poles+1
  
! loop over complex pole pairs
  do i=1,H_PR%n_complex_pole_pairs
  
    pole1=first_complex_pole+(i-1)*2
    pole2=pole1+1
    
! check that we have a conjugate pole pair
    if ( .NOT.conjugate_pair(H_PR%poles(pole1),H_PR%poles(pole2)) ) then
      write(*,*)'ERROR in LC_test_BRUNE: poles are not a conjugate pair'
      write(*,*)'pole 1:',H_PR%poles(pole1)
      write(*,*)'pole 2:',H_PR%poles(pole2)
      STOP
    end if
    
! check that we have a conjugate residue pair
    if ( .NOT.conjugate_pair(H_PR%residues(pole1),H_PR%residues(pole2)) ) then
      write(*,*)'ERROR in LC_test_BRUNE: residues are not a conjugate pair'
      write(*,*)'residues 1:',H_PR%residues(pole1)
      write(*,*)'residues 2:',H_PR%residues(pole2)
      STOP
    end if

! test for whether we have an LC branch here...
    
    imaginary_poles    =imaginary_pair(H_PR%poles(pole1),H_PR%poles(pole2))
    positive_residues=(dble(H_PR%residues(pole1)).GT.0d0)
    
    const_term=dble(H_PR%residues(pole1))*dble(H_PR%poles(pole1))+    &
               aimag(H_PR%poles(pole1))*aimag(H_PR%residues(pole1))
    zero_const_term=(abs(const_term).LT.zero_test_small)
        
    if (verbose) then
      write(*,*)'Testing pole pair',i
      write(*,*)'Tests:'
      write(*,*)'imaginary pole test      : ',imaginary_poles,' ',H_PR%poles(pole1)
      write(*,*)'positive residue test    : ',positive_residues,' ',H_PR%residues(pole1)
      write(*,*)'zero constant term test: : ',zero_const_term,' ', const_term 
    end if
    
! check for stable poles which are not on the imaginary s=jw axis AND positive residues
    if (imaginary_poles.AND.positive_residues.AND.zero_const_term) then
    
      if (verbose) write(*,*)'Found possible LC branch'

! this could be a viable LC branch - calculate the remainder when this pole pair is removed

      CALL deallocate_Sfilter(HR)
      
! Test whether the remainder is zero
      
! build a local pole-residue filter without the test pole pair. 

! allocate the structure for the local pole-residue form function and copy the 
! required information across            
      HR_PR_local%wnorm=H_PR%wnorm
      HR_PR_local%order=H_PR%order-2
      HR_PR_local%n_real_poles=H_PR%n_real_poles
      HR_PR_local%n_complex_poles=H_PR%n_complex_poles-2
      HR_PR_local%n_complex_pole_pairs=H_PR%n_complex_pole_pairs-1
      HR_PR_local%n_real_poles=H_PR%n_real_poles

! constant term and sL term
      
      HR_PR_local%R=H_PR%R
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
      
        if (ii.NE.i) then
! this is not the pole pair we are trying to remove so copy them over
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
        else
! this is the pair of poles we wish to remove so just increase the pole2 counter by 2. 
          pole2=pole2+2 
        end if
        
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
      
      remainder_zero=.FALSE.
      GOTO 8000
        
    end if ! positive residues for this complex pole pair
  
  end do ! next complex pole pair
  
! we only get here if we have not found a viable LC branch

  remainder_OK=.FALSE.
  found=.FALSE.
  remainder_zero=.FALSE.
  
  RETURN
  
8000 CONTINUE
! jump here if we have found a viable LC branch

  remainder_OK=.TRUE.
  found=.TRUE.
  CALL deallocate_Sfilter_PR(HR_PR_local)
  pole1=first_complex_pole+(i-1)*2
  pole2=pole1+1
  L=1d0/dble(H_PR%residues(pole1)+H_PR%residues(pole2))
  C= dble( (H_PR%residues(pole1)+H_PR%residues(pole2))/(H_PR%poles(pole1)*H_PR%poles(pole2)) )
  if (verbose) then
    write(*,*)'FOUND VIABLE LC BRANCH'
    write(*,*)'L=',L
    write(*,*)'C=',C
    write(*,*)'remainder_OK   :',remainder_OK
    write(*,*)'remainder_zero :',remainder_zero
  end if
  
  RETURN

  END SUBROUTINE LC_test_BRUNE
!
! NAME
!     L_test_BRUNE
!
! DESCRIPTION
!       look for a viable series inductance in the BRUNE model which may be negative.
!
! SEE ALSO
!
!
! HISTORY
!
!     started 14/09/17 CJS
!
  SUBROUTINE L_test_BRUNE(H_PR,L,found,HR,remainder_OK,remainder_zero)

USE type_specifications
USE general_module
USE constants
USE filter_module

IMPLICIT NONE
  
  type(Sfilter_PR),INTENT(IN) :: H_PR
  real(dp):: L
  logical :: found
  type(Sfilter),INTENT(INOUT) :: HR
  logical :: remainder_OK,remainder_zero

! local variables
  
  integer :: pole,pole1,pole2
  type(Sfilter_PR) :: HR_PR_local
  logical :: stable
  
  integer :: i,ii

!START

  if (verbose) write(*,*)'CALLED: L_test_BRUNE'

  found=.FALSE.
          
  if (verbose) write(*,*)'Found possible L branch'

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
      
      HR_PR_local%R=H_PR%R
      HR_PR_local%L=0d0
      
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
      
! we only get here if we have not found a viable L branch

  remainder_OK=.FALSE.
  found=.FALSE.
  remainder_zero=.FALSE.
  
  RETURN
  
8000 CONTINUE
! jump here if we have found a viable RC branch

  remainder_OK=.TRUE.
  found=.TRUE.
  CALL deallocate_Sfilter_PR(HR_PR_local)
   
  L=H_PR%L
  if (verbose) then
    write(*,*)'FOUND VIABLE SERIES L BRANCH'
    write(*,*)'L=',L
    write(*,*)'remainder_OK   :',remainder_OK
    write(*,*)'remainder_zero :',remainder_zero
  end if
  
  RETURN

  END SUBROUTINE L_test_BRUNE
