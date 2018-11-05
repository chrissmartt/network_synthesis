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
! File Contents:
! SUBROUTINE calculate_min_resistance_function(H,R,w0)
! SUBROUTINE calculate_min_resistance_value(H,R,w0)
!
! NAME
!     calculate_min_resistance_value
!
! DESCRIPTION
!     This subroutine calculates the minimum resistance of a rational transfer function
!     and the angular frequency at which it occurs
!     See F. F. Kuo, "Network Analysis and Synthesis" section 10.3
!
!     The algorithm first obtains an expression for the real part of the rational function as
!     a function of frequency. This is then differentiated and the minima calculated.
!     The lowest minimum, or the resistance as f-> infinty if is smaller is then the minimum resistance
!
! SEE ALSO
!
!
! HISTORY
!
!     started 05/10/2017 CJS
!
  SUBROUTINE calculate_min_resistance_value(H,R,w0)

USE type_specifications
USE general_module
USE constants
USE frequency_spec
USE filter_module
USE Sfilter_fit_module

IMPLICIT NONE

! variables passed to subroutine  
  type(Sfilter),intent(INOUT)       :: H        ! Input Rational transfer function 
  real(dp),intent(OUT)              :: R        ! Output Minimum resistance (real part of the transfer function)
  real(dp),intent(OUT)              :: w0       ! Output Angular frequency at which the minimum occurs
  
! local variables

  type(Sfilter)    :: Hnew
  type(polynomial) :: m1,m2,n1,n2,t1,t2,A2,A,Ap,B,B2,Bp,Q
  integer  :: i,ii,order,m_order,n_order,a_order
  
  integer      :: poly_order
  complex(dp),allocatable :: roots(:)
  complex(dp),allocatable :: rroots(:)
  complex(dp),allocatable :: croots(:)
  integer                 :: nreal
  integer                 :: ncomplex
  
  integer  :: root
  real(dp) :: min_R,min_w0,w,f,R_w
  
  logical :: pole_at_zero
  logical :: local_verbose

    
!START

!  local_verbose=.TRUE.
  local_verbose=.FALSE.

  if (local_verbose) then
    write(*,*)'CALLED calculate_min_resistance_value'
    write(*,*)'Input filter function:'
    CALL write_Sfilter(H,0)
  end if
  
  if (abs(H%b%coeff(0)).LT.zero_test_small) then
    pole_at_zero=.TRUE.
  else
    pole_at_zero=.FALSE.
  end if
  
! assemble the odd and even polynomials m1 and n1 from the numerator of H 

  order=H%a%order
  if (mod(order,2).EQ.0) then
    m_order=order
    n_order=max(order-1,0)
  else
    m_order=max(order-1,0)
    n_order=order 
  end if
  
  m1=allocate_polynomial(m_order)
  n1=allocate_polynomial(n_order)
  
  m1%coeff(:)=0d0
  do i=0,m_order,2
    m1%coeff(i)=H%a%coeff(i)
  end do
  
  n1%coeff(:)=0d0
  do i=1,n_order,2
    n1%coeff(i)=H%a%coeff(i)
  end do
  
  if (local_verbose) then
    write(*,*)'Numerator:'
    CALL write_poly_local3(H%a)  
    write(*,*)'m1: even order terms'
    CALL write_poly_local3(m1)  
    write(*,*)'n1: odd order terms'
    CALL write_poly_local3(n1)  
  end if

! assemble the odd and even polynomials m2 and n2 from the denominator of H 

  order=H%b%order
  if (mod(order,2).EQ.0) then
    m_order=order
    n_order=max(order-1,0)
  else
    m_order=max(order-1,0)
    n_order=order 
  end if
  
  m2=allocate_polynomial(m_order)
  n2=allocate_polynomial(n_order)
  
  m2%coeff(:)=0d0
  do i=0,m_order,2
    m2%coeff(i)=H%b%coeff(i)
  end do
  
  n2%coeff(:)=0d0
  do i=1,n_order,2
    n2%coeff(i)=H%b%coeff(i)
  end do
  
  if (local_verbose) then
    write(*,*)'Denominator:'
    CALL write_poly_local3(H%b)  
    write(*,*)'m2: even order terms'
    CALL write_poly_local3(m2)  
    write(*,*)'n2: odd order terms'
    CALL write_poly_local3(n2)  
  end if

! Calculate the polynomial A2(jw)=m1(jw)m2(jw)-n1(jw)n2(jw)

  t1=m1*m2
  t2=n1*n2
  A2=t1-t2
  
  CALL get_min_order_poly(A2)
  
  if (local_verbose) then
    write(*,*)'A2(jw)=m1m2-n1n2'
    CALL write_poly_local3(A2)  
   end if
   
! Calculate the coefficients of A(w)

  order=A2%order
  if (mod(order,2).EQ.0) then
    a_order=order
  else
    write(*,*)'Error in calculate_min_resistance_function. Order of A2(jw) is not an even number'
  end if
  
  A=allocate_polynomial(a_order)
  
  do i=0,a_order
    A%coeff(i)=A2%coeff(i)
  end do

! At the moment A is a function of (jw) not w with zero odd order coefficients
! We get the function of w by making the jw**2n order coefficients negative
  do i=2,a_order,4
    A%coeff(i)=-A%coeff(i)
  end do
    
  CALL get_min_order_poly(A)
  
  if (local_verbose) then
    write(*,*)'A(w)='
    CALL write_poly_local3(A)  
   end if

! Calculate the polynomial B2(jw)=m2(jw)m2(jw)-n2(jw)n2(jw)

  t1=m2*m2
  t2=n2*n2
  B2=t1-t2
  
  CALL get_min_order_poly(B2)
  
  if (local_verbose) then
    write(*,*)'B2(jw)=m1m2-n1n2'
    CALL write_poly_local3(B2)  
   end if
   
! Calculate the coefficients of B(w)

  order=B2%order
  if (mod(order,2).EQ.0) then
    a_order=order
  else
    write(*,*)'Error in calculate_min_resistance_function. Order of B2(jw) is not an even number'
  end if
  
  B=allocate_polynomial(a_order)
  
  do i=0,a_order
    B%coeff(i)=B2%coeff(i)
  end do

! Bt the moment B is a function of (jw) not w with zero odd order coefficients
! We get the function of w by making the jw**2n order coefficients negative
  do i=2,a_order,4
    B%coeff(i)=-B%coeff(i)
  end do
    
  CALL get_min_order_poly(B)
  
  if (local_verbose) then
    write(*,*)'B(w)='
    CALL write_poly_local3(B)  
   end if

! Calculate the derivative of A(w)/B(w) wrt w

  order=max(0,A%order-1)
  Ap=allocate_polynomial(order)
  
  if (A%order.GT.0) then
    do i=0,Ap%order
      Ap%coeff(i)=A%coeff(i+1)*dble(i+1)
    end do
  end if
  
  if (local_verbose) then
    write(*,*)"Deivative function A'(w)="
    CALL write_poly_local3(Ap)  
   end if

  order=max(0,B%order-1)
  Bp=allocate_polynomial(order)

  if (B%order.GT.0) then
    do i=0,Bp%order
      Bp%coeff(i)=B%coeff(i+1)*dble(i+1)
    end do
  end if
  
  if (local_verbose) then
    write(*,*)"Deivative function B'(w)="
    CALL write_poly_local3(Bp)  
   end if
  
! Calculate the zeros of the derivative function d/dw(A(w)/B(w))=
  t1=Ap*B
  t2=A*Bp
  Q=t1-t2
   
  if (local_verbose) then
    write(*,*)"we require Q(w)=0 where Q="
    CALL write_poly_local3(Q)  
   end if
 
  CALL get_min_order_poly(Q)
  
  if (local_verbose) then
    write(*,*)"we require Q(w)=0 where Q="
    CALL write_poly_local3(Q)  
   end if

  poly_order=Q%order
  ALLOCATE( roots(1:poly_order) )
  CALL findroots(Q,roots,poly_order)
  
  if (local_verbose) then
    write(*,*)'Roots of the derivative function are'
    do root=1,poly_order
      write(*,*)roots(root)
    end do
  end if
  
  ALLOCATE( rroots(1:poly_order) )
  ALLOCATE( croots(1:poly_order) )
  CALL rootsort(poly_order,roots,rroots,             &
                croots,nreal,ncomplex,poly_order)
                     
  if (local_verbose) write(*,*)'Getting the high frequency resistance' 
  
  
  write(*,*)'Getting the high frequency resistance of H:' 
  CALL write_Sfilter(H,0)
  write(*,*)'aorder=',H%a%order
  write(*,*)'border=',H%b%order
  
  min_w0=1D30
  min_R=evaluate_Sfilter_high_frequency_limit(H)
  if (local_verbose) write(*,*)'High frequency resistance=',min_R 
  
  if (local_verbose) write(*,*)'Checking real roots: number of roots=',nreal,poly_order 
  do root=1,nreal
    
    w=-dble(rroots(root))
    f=w*H%wnorm/(2d0*pi)  
    
    if ((pole_at_zero).AND.(abs(f).LT.zero_test_small)) then
! perturb the frequency slightly away from the pole at zero
      write(*,*)'Special case, f=',f
      write(*,*)'(pole_at_zero).AND.(abs(f).LT.zero_test_small),zero_test_small=',zero_test_small
      f=f+zero_test_small
      write(*,*)'new f=',f
    end if

    R_w=evaluate_Sfilter_frequency_response(H,f)
    
    if (local_verbose) then
      write(*,*)'root',rroots(root),'w=',w,' R=',R_w
    end if
   
    if (R_w.EQ.min_R) then
! if we have roots at + and - p, choose the one with positive frequency
      min_R=R_w
      min_w0=max(w,min_w0)
    else if (R_w.LT.min_R) then
      min_R=R_w
      min_w0=w
    end if
    
  end do
  
  R=min_R
  w0=min_w0

  if (local_verbose) then
    write(*,*)'FINISHED calculate_min_resistance_value'
    write(*,*)
    write(*,*)'R0=',R
    write(*,*)'w0=',w0
    write(*,*)
  end if
  
! finish up
  
  DEALLOCATE( roots )
  DEALLOCATE( rroots )
  DEALLOCATE( croots )

  RETURN
  END SUBROUTINE calculate_min_resistance_value
!
! NAME
!     calculate_min_resistance_function
!
! DESCRIPTION
!     This subroutine calculates the minimum resistance value of a PR function
!     and then subtracts it from the input function to give a function
!     whose minimum resistance is zero. This is required in the Brune Synthesis. 
!
! SEE ALSO
!
!
! HISTORY
!
!     started 04/03/09 CJS
!
  SUBROUTINE calculate_min_resistance_function(H,R,w0)

USE type_specifications
USE general_module
USE constants
USE frequency_spec
USE filter_module
USE Sfilter_fit_module

IMPLICIT NONE

! variables passed to subroutine  
  type(Sfilter),intent(INOUT)    :: H            ! Input/ Output Rational transfer function 
  real(dp),intent(OUT)           :: R            ! Output Minimum resistance (real part of the transfer function)
  real(dp),intent(OUT)           :: w0           ! Output Angular frequency at which the minimum occurs
  
! local variables

  type(Sfilter)    :: Hnew
  type(polynomial) :: t1,t2
  real(dp) :: min_R,min_w0,w,f,R_w
  
  logical ::local_verbose

    
!START

  local_verbose=.TRUE.

  if (local_verbose) then
    write(*,*)'CALLED calculate_min_resistance_function'
    write(*,*)'Input filter function:'
    CALL write_Sfilter(H,0)
  end if
  
  CALL calculate_min_resistance_value(H,min_R,min_w0)

! Deal with the numerics around small values of min_R
  if (min_R.LT.-zero_test_small) then
    write(*,*)'calculate_min_resistance_function called with non positive real function'
    STOP
  else if (abs(min_R).LT.zero_test_small) then
    min_R=0d0
  end if
  
  if (local_verbose) then
    write(*,*)'Minimum R=',min_R,' at w=',min_w0 
  end if
  
! now subtract R_min from the input function H(jw)

  Hnew=H
  
  t1=H%a
  t2=H%b
  t2%coeff(:)=t2%coeff(:)*min_R
  
  Hnew%a=t1-t2
  CALL get_min_order_poly(Hnew%a)
  
  H=Hnew
  
  R=min_R
  w0=min_w0

  if (local_verbose) then
    write(*,*)'FINISHED calculate_min_resistance_function'
    write(*,*)'Output filter function:'
    CALL write_Sfilter(H,0)
    write(*,*)
  end if

  RETURN
  END SUBROUTINE calculate_min_resistance_function
