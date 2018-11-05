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
! Description
!   Local filter functions used specifically in the network synthesis processes
!
! Comments:
!   These subroutines could be rationalised and merged into the POLYNOMIAL_AND_FILTER_MODULES
!
! History
!
!     16/11/2017 CJS Include network synthesis process to replace s-domain transfer functions
!
! SUBROUTINE write_Sfilter_local(s1)
! SUBROUTINE write_CF_local(CFterm,CFtype,CF_dim,max_order)
! SUBROUTINE check_transfer_function(Z,stable)
! SUBROUTINE check_transfer_function_num_den(a,b,stable)
! SUBROUTINE pole_zero_cancel(H)
!
!
!
! NAME
!     write_Sfilter_local
!
! DESCRIPTION
!       write Laplace domain filter coefficients to screen or file 
!
! SEE ALSO
!
!
! HISTORY
!
!     started 23/01/09 CJS
!
  SUBROUTINE write_Sfilter_local(s1)


USE type_specifications
USE general_module
USE constants
USE frequency_spec
USE filter_module
USE Sfilter_fit_module
 
IMPLICIT NONE
 
  type(Sfilter),intent(IN):: s1
  
  integer i

!START

! write to screen
  write(*,*)'Laplace domain filter'
  write(*,*)'a order=',s1%a%order
  write(*,*)'b order=',s1%b%order
  write(*,*)'wnorm=',s1%wnorm
  write(*,*)''
  
  write(*,'(A)',advance='NO')'H(s) = '
  
  do i=s1%a%order,1,-1
    write(*,8000,advance='NO')s1%a%coeff(i),' s^',i,' + '
  end do
  write(*,8010)s1%a%coeff(0)
  
  write(*,'(A)',advance='NO')'       '
  do i=s1%a%order,1,-1
    write(*,'(A)',advance='NO')'------------'
  end do
  write(*,'(A)')'-----'
  
  write(*,'(A)',advance='NO')'       '
  do i=s1%b%order,1,-1
    write(*,8000,advance='NO')s1%b%coeff(i),' s^',i,' + '
  end do
  write(*,8010)s1%b%coeff(0)
      
8000 format(F5.2,A3,I1,A3)
8010 format(F5.2)
  
  RETURN
  END SUBROUTINE write_Sfilter_local

!
! NAME
!     write_CF_local
!
! DESCRIPTION
!       write continued fraction
!
! SEE ALSO
!
!
! HISTORY
!
!     started 23/01/09 CJS
!
SUBROUTINE write_CF_local(CFterm,CFtype,CF_dim,max_order)

USE type_specifications
USE general_module
USE constants
USE frequency_spec
USE filter_module
USE Sfilter_fit_module

IMPLICIT NONE

integer CF_dim,max_order

real(dp) :: CFterm(1:CF_dim)
integer  :: CFtype(1:CF_dim)

! local variables

real(dp):: term
integer :: type
integer :: last_type

integer :: loop,i

logical :: first_term_on_line

!START

! write the continued fraction to the screen
    
  write(*,*)'Continued fraction form:'
  
  last_type=1 ! assume this is an impedance function
  first_term_on_line=.TRUE.
  
  do loop=1,max_order
  
    if (loop.Eq.1) then
      write(*,'(A)',advance='NO')'H(s) = '
    end if
      
    term=CFterm(loop)
    type=CFtype(loop)
            
!    if(type*last_type.LT.0) then
!      write(*,'(A)')' + 1/ '
!      write(*,'(A)',advance='NO')'       '    
!      do i=1,loop-1
!        write(*,'(A)',advance='NO')'                 '
!      end do
!      first_term_on_line=.TRUE.
!    end if
            
    if(type*last_type.LT.0) then
      write(*,'(A)',advance='NO')' + 1__'
      do i=1,max_order-1
        write(*,'(A)',advance='NO')'_________________'
      end do
      write(*,*)
      
      write(*,'(A)',advance='NO')'       '    
      do i=1,loop-1
        write(*,'(A)',advance='NO')'                 '
      end do
      first_term_on_line=.TRUE.
    end if
    
    if (.NOT.first_term_on_line) then
      write(*,'(A)',advance='NO')' + '       
    end if
    
    if (abs(type).EQ.1) then
    
      write(*,8000,advance='NO')term,' s ' 
      first_term_on_line=.FALSE.
     
    else if (abs(type).EQ.2) then
    
      write(*,8000,advance='NO')term,'   ' 
      first_term_on_line=.FALSE.
   
    else if (abs(type).EQ.3) then
    
      write(*,8020,advance='NO')'1/(',term,'s) '     
      first_term_on_line=.FALSE.
      
    end if

8000 format(F5.2,A)
8010 format(A,F5.2)
8020 format(A,F5.2,A)
    
    last_type=type
   
  end do

  write(*,*)' '
  
  RETURN
  END SUBROUTINE write_CF_local

!
! NAME
!     check_transfer_function
!
! DESCRIPTION
!       check that a transfer function describing an impedance is 'physical' 
!       i.e. that the tansfer function is stable 
!       and that the impedance described by the transfer function is dissipative for all frequencies
!
! SEE ALSO
!
!
! HISTORY
!
!     started 23/01/09 CJS
!
  SUBROUTINE check_transfer_function(Z,stable)


USE type_specifications
USE general_module
USE constants
USE frequency_spec
USE filter_module
USE Sfilter_fit_module
 
IMPLICIT NONE
 
  type(Sfilter),intent(INOUT):: Z
  logical,intent(OUT) :: stable
  
! local variables

  type(Sfilter):: Y
  integer :: i
  logical :: check_passed
  
  integer :: lpn,lpd

!START

  if(verbose) write(*,*)'Check Impedance function:'
  
! Initial process to remove any leading zero coefficients
  
  CALL get_min_order_poly(Z%a)
  CALL get_min_order_poly(Z%b)
  
! assume stability
  stable=.TRUE.

! 1. Check that the numerator order does not exceed the denominator order by more than 1
  
! check that the restrictions on the polynomial orders are satisfied
  if (Z%a%order.GT.Z%b%order+1) then
    if (verbose) write(*,*)'Numerator order greater than denominator order by more than 1: ' &
    ,Z%a%order,Z%b%order
    stable=.FALSE.
  end if
  
! 2. Check that poles are in the left half of the s-plane
  CALL test_filter_pole_stability(Z,check_passed)
  
  if (.NOT.check_passed) then
    if (verbose) write(*,*)'The impedance filter function has unstable poles'
    stable=.FALSE.
  else
    if (verbose) write(*,*)'The impedance filter function has stable poles'    
  end if

! 3. check that we only have simple poles on the jw axis

  CALL test_filter_simple_poles(Z,check_passed)
  
  if (.NOT.check_passed) then
    if (verbose) write(*,*)'The impedance filter function has multiple poles on the jw axis'
    stable=.FALSE.
  else
    if (verbose) write(*,*)'The impedance filter function does not have multiple poles on the jw axis'    
  end if

! 4 check that the lowest powers of numerator and denominator differ by 1 at the most
  lpn=0
  do i=0,Z%a%order
    if ( (Z%a%coeff(i).NE.0d0).AND.(lpn.Eq.0) ) then
      lpn=i
    end if
  end do
  lpd=0
  do i=0,Z%b%order
    if ( (Z%b%coeff(i).NE.0d0).AND.(lpd.Eq.0) ) then
      lpd=i
    end if
  end do
  
  if (abs(lpn-lpd).GT.1) then
    if (verbose) write(*,*)'The impedance filter function lowest powers of numerator and denominator differ by > 1'
    stable=.FALSE.
  else
    if (verbose) write(*,*)'The impedance filter function lowest powers of numerator and denominator differ by <= 1'    
  end if

! check that all coefficients are positive
  do i=0,Z%a%order
    if (Z%a%coeff(i).LT.0d0) then
      if (verbose) write(*,*)'Negative numerator coefficient'
      stable=.FALSE.     
    end if
  end do
  
  do i=0,Z%b%order
    if (Z%b%coeff(i).LT.0d0) then
      if (verbose) write(*,*)'Negative denominator coefficient'
      stable=.FALSE.     
    end if
  end do

! 6 check for positive real impedance function
  CALL test_filter_positive_real(Z,check_passed,.FALSE.) ! LAST VARIABLE IS LOCAL VERBOSE FLAG
!  CALL test_filter_positive_real(Z,check_passed,.TRUE.) ! LAST VARIABLE IS LOCAL VERBOSE FLAG
  if (.NOT.check_passed) then
    if (verbose) write(*,*)'The impedance filter function is not positive real'
    stable=.FALSE.
  else
    if (verbose) write(*,*)'The impedance filter function is positive real'    
  end if

!  write(*,*)'Check admittance function:'

  Y=allocate_Sfilter(Z%b%order,Z%a%order)
  Y%wnorm=Z%wnorm
  do i=0,Z%b%order
    Y%a%coeff(i)=Z%b%coeff(i)
  end do
  do i=0,Z%a%order
    Y%b%coeff(i)=Z%a%coeff(i)
  end do

! 1b. Check that the numerator order does not exceed the denominator order by more than 1
  
! check that the restrictions on the polynomial orders are satisfied
  if (Y%a%order.GT.Y%b%order+1) then
    if (verbose) write(*,*)'Numerator order greater than denominator order by more than 1: ' &
    ,Y%a%order,Y%b%order
    stable=.FALSE.
  end if
  
! 2b. Check that poles are in the left half of the s-plane
  CALL test_filter_pole_stability(Y,check_passed)
  
  if (.NOT.check_passed) then
    if (verbose) write(*,*)'The admittance filter function has unstable poles'
    stable=.FALSE.
  else
    if (verbose) write(*,*)'The admittance filter function has stable poles'    
  end if

! 3b. check that we only have simple poles on the jw axis

  CALL test_filter_simple_poles(Y,check_passed)
  
  if (.NOT.check_passed) then
    if (verbose) write(*,*)'The filter function has multiple poles on the jw axis'
    stable=.FALSE.
  else
    if (verbose) write(*,*)'The admittance filter function has no multiple poles on the jw axis'    
  end if
  
  CALL deallocate_Sfilter(Y)

  RETURN
  END SUBROUTINE check_transfer_function
!
! NAME
!     check_transfer_function_num_den
!
! DESCRIPTION
!       check that a transfer function describing an impedance is 'physical'
!       i.e. that the tansfer function is stable 
!       and that the impedance described by the transfer function is dissipative for all frequencies
!       Here the transfer function is specified as separate numerator and denominator polynomials
!       The subroutine assembles the rational function and calls check_transfer_function.
!
! SEE ALSO
!       SUBROUTINE check_transfer_function above.
!
! HISTORY
!
!     started 23/01/09 CJS
!
  SUBROUTINE check_transfer_function_num_den(a,b,stable)


USE type_specifications
USE general_module
USE constants
USE frequency_spec
USE filter_module
USE Sfilter_fit_module
 
IMPLICIT NONE
 
  type(polynomial),intent(IN):: a,b       ! input numerator and denominator polynomials
  logical,intent(OUT) :: stable           ! output flag to indicate stability
  
! local variables

  type(Sfilter):: Z    ! impedance filter function
  integer :: i

!START
 
! construct the impedance filter

  Z=allocate_Sfilter(a%order,b%order)
  Z%wnorm=1d0
  do i=0,a%order
    Z%a%coeff(i)=a%coeff(i)
  end do
  do i=0,b%order
    Z%b%coeff(i)=b%coeff(i)
  end do

  CALL check_transfer_function(Z,stable)

  CALL deallocate_Sfilter(Z)
  
  RETURN
  END SUBROUTINE check_transfer_function_num_den

!
! NAME
!     pole_zero_cancel
!
! DESCRIPTION
!       Find the poles and zeros of the transfer function and do any cancellation 
!       that is possible when a pole is equal to a zero
!
! SEE ALSO
!
!
! HISTORY
!
!     started 25/9/2017 CJS
!
  SUBROUTINE pole_zero_cancel(H)


USE type_specifications
USE general_module
USE constants
USE frequency_spec
USE filter_module
USE Sfilter_fit_module
 
IMPLICIT NONE
 
  type(Sfilter),intent(INOUT):: H      ! transfer function input and output
  
! local variables

  integer :: i
  integer      :: aorder,border
  complex(dp),allocatable  :: roots(:)
  complex(dp),allocatable  :: poles(:)
  real(dp) :: AK,BK
  
  integer :: zero,pole
  integer :: zero2,pole2
  
  type(complex_polynomial) :: function,term,temp_poly
  


!START

  if (verbose) then
    write(*,*)'CALLED: pole_zero_cancel'
    write(*,*)'H='
    CALL write_Sfilter(H,0)  
  end if

1000 CONTINUE
 
  aorder=H%a%order  
  ALLOCATE( roots(1:aorder) )  
  AK=H%a%coeff(aorder)
  CALL findroots(H%a,roots,aorder)
  if (verbose) then
    do zero=1,aorder
      if (verbose) write(*,*)'Found a zero',roots(zero)
    end do
  end if
  
  border=H%b%order  
  ALLOCATE( poles(1:border) )  
  BK=H%b%coeff(border)
  CALL findroots(H%b,poles,border)
  if (verbose) then
    do pole=1,border
      if (verbose) write(*,*)'Found a pole',poles(pole)
    end do
  end if
  
! loop over zeros
  do zero=1,aorder
  
! loop over poles
    do pole=1,border
    
      if ( (abs(roots(zero)-poles(pole)).LT.zero_test_small) ) then
        
! this is a repeated root set the mutliple roots flag
          
        if (verbose) write(*,*)'Found a pole - zero pair to cancel',roots(zero)
        
        if (verbose) write(*,*)'numerator'
        term=allocate_complex_polynomial(1)
        term%coeff(1)=(1d0,0d0)
        function=allocate_complex_polynomial(0)
        function%coeff(0)=(1d0,0d0)
        
        do zero2=1,aorder
          if (verbose) write(*,*)'zero',zero2
          if (zero2.NE.zero) then
            if (verbose) write(*,*)'Include this zero:',roots(zero2)
            term%coeff(0)=-roots(zero2)
            temp_poly=function*term
            function=temp_poly        
          end if
        end do
        function%coeff(:)=function%coeff(:)*AK
        
        deallocate(H%a%coeff)
        H%a=allocate_polynomial(function%order)
        H%a%coeff(0:H%a%order)=dble(function%coeff(0:H%a%order))
        if (verbose) CALL write_poly_local2('N','x',H%a)  
        
        CALL deallocate_complex_poly(function)
        function=allocate_complex_polynomial(0)
        function%coeff(0)=(1d0,0d0)
        
        if (verbose) write(*,*)'denominator'
        do pole2=1,border
          if (pole2.NE.pole) then
            if (verbose) write(*,*)'Include this pole:',poles(pole2)
            term%coeff(0)=-poles(pole2)
            temp_poly=function*term
            function=temp_poly        
          end if
        end do
        function%coeff(:)=function%coeff(:)*BK
        
        deallocate(H%b%coeff)
        H%b=allocate_polynomial(function%order)
        H%b%coeff(0:H%b%order)=dble(function%coeff(0:H%b%order))
        if (verbose) CALL write_poly_local2('D','x',H%b)  
        
        CALL deallocate_complex_poly(function)
        CALL deallocate_complex_poly(term)
        CALL deallocate_complex_poly(temp_poly)

! Start the process again with the reduced filter
        DEALLOCATE(roots)
        DEALLOCATE(poles)
        GOTO 1000  
          
      end if
    
    end do ! next pole
    
  end do ! next root
  
  DEALLOCATE(roots)
  DEALLOCATE(poles)

  if (verbose) then
    write(*,*)'FINISHED: pole_zero_cancel'
    write(*,*)'H='
    CALL write_Sfilter(H,0)  
  end if
  
  RETURN
  END SUBROUTINE pole_zero_cancel

