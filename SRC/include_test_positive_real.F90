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
! SUBROUTINE test_filter_positive_real(H,stable,local_verbose)
!
!
! NAME
!     test_filter_positive_real
!
! DESCRIPTION
!     test whether a rational function is positive real. The definition of positive real
!     is described in the Theory Manual.
!     This implements the Sturm test as described in Kim and Meadows section 9.1.3
!
! SEE ALSO
!
!
! HISTORY
!
!     started 9/2017 CJS
!     24/10/2017 CJS Use the minimum resistance value rather than the Sturm test method.
!
  SUBROUTINE test_filter_positive_real(H,stable,local_verbose)

USE type_specifications
USE general_module
USE constants
USE frequency_spec
USE filter_module

IMPLICIT NONE

! variables passed to subroutine  
  type(Sfilter),intent(IN)    :: H
  logical,intent(OUT)         :: stable
  logical,intent(IN)          :: local_verbose
  
! local variables

  type(polynomial) :: m1,m2,n1,n2,P,Q,t1,t2,A2,A,Ap,Am,C,R
  integer  :: i,ii,order,m_order,n_order,a_order
  
  real(dp),allocatable :: value_inf(:)
  real(dp),allocatable :: value_0(:)

  integer  :: n_sturm,sturm_coeff
  
  integer :: S_inf,S_0,n_real_zeros
  real(dp) :: last_value_inf,last_value_0
  
  logical :: odd_multiplicity_flag
  
  real(dp) :: Rmin,Wmin
  type(Sfilter)    :: Hlocal
   
!START

  if (local_verbose) then
    write(*,*)'CALLED test_filter_positive_real'
    write(*,*)'Input filter function:'
    CALL write_Sfilter(H,0)
  end if
  
  Hlocal=H
  CALL calculate_min_resistance_value(Hlocal,Rmin,wmin)
  
  if (Rmin.LT.0d0) then
    stable=.FALSE.
  else
    stable=.TRUE.
  end if
  
  CALL deallocate_Sfilter(Hlocal)

RETURN

! assume the filter is stable initially
  stable=.TRUE.

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
    CALL write_poly_local(H%a)  
    write(*,*)'m1: even order terms'
    CALL write_poly_local(m1)  
    write(*,*)'n1: odd order terms'
    CALL write_poly_local(n1)  
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
    CALL write_poly_local(H%b)  
    write(*,*)'m2: even order terms'
    CALL write_poly_local(m2)  
    write(*,*)'n2: odd order terms'
    CALL write_poly_local(n2)  
  end if

! Calculate the polynomial A2(jw)=m1(jw)m2(jw)-n1(jw)n2(jw)
! This is in fact a function of w**2 only (no odd order terms in jw)

  t1=m1*m2
  t2=n1*n2
  A2=t1-t2
  
  CALL get_min_order_poly(A2)
  
  if (local_verbose) then
    write(*,*)'A2(jw)=m1m2-n1n2'
    CALL write_poly_local(A2)  
   end if
   
! Calculate the coefficients of A(w**2)

  order=A2%order
  if (mod(order,2).EQ.0) then
    a_order=order/2
  else
    write(*,*)'Error in test_filter_positive_real. Order of A2(jw) is not an even number'
  end if
  
  A=allocate_polynomial(a_order)
  
  do i=0,a_order
    A%coeff(i)=A2%coeff(2*i)
  end do

! At the moment A is a function of (jw)**2 not w**2.
! We get the function of w**2 by making the odd order coefficients negative
  do i=1,a_order,2
    A%coeff(i)=-A%coeff(i)
  end do
  
  if (local_verbose) then
    write(*,*)'A(w**2)'
    CALL write_poly_local(A)  
    write(*,*)'Full precision coefficients'
    do i=0,A%order
      write(*,*)i,A%coeff(i)
    end do
  end if
    
  CALL get_min_order_poly(A)
  
  if (local_verbose) then
    write(*,*)'minimum order: A(w**2)'
    CALL write_poly_local(A)  
    write(*,*)'Full precision coefficients'
    do i=0,A%order
      write(*,*)i,A%coeff(i)
    end do
  end if

! Get the high frequency value of the function. 
! If it is less than zero the funciton is not positive real so return
  if (A%coeff(A%order).LT.-zero_test_small) then
    stable=.FALSE.
    if (local_verbose) then
      write(*,*)'Real part of function is less than zero at high frequency'
      write(*,*)'test value:',A%coeff(A%order)
      CALL write_poly_local(A)  
    end if
    RETURN
  end if
  
!! **** QUESTION: If we calculate all the roots to remove any even multiple roots then 
!! maybe we should just test for odd multiplicity real roots directly...
!
  CALL remove_even_multiple_zeros(A,odd_multiplicity_flag,local_verbose)
  
  if (local_verbose) then
    write(*,*)'remove even multiple zeros: A(w**2)'
    CALL write_poly_local(A)  
    write(*,*)'Full precision coefficients'
    do i=0,A%order
      write(*,*)i,A%coeff(i)
    end do
  end if
  
! Checks for a valid at this point
  
  if (a%coeff(a%order).EQ.0d0) then
    write(*,*)'Error in test_filter_positive_real. Highest order A coefficient =0d0'
  end if
  
! Start to assemble the terms of the Sturm sequence

  n_sturm=a_order+1

  allocate( value_inf(1:n_sturm) )
  allocate( value_0(1:n_sturm) )

! We already have the first term which is A

  Ap=A
  
  sturm_coeff=1

! Get a function sign for x->infinity
  value_inf(sturm_coeff)=0d0
  do ii=Ap%order,0,-1
    if (abs(Ap%coeff(ii)).GT.zero_test_small) then
      value_inf(sturm_coeff)=Ap%coeff(ii)
      exit
    end if
  end do
  
! Get a function sign for x->0
  value_0(sturm_coeff)=0d0
  do ii=0,Ap%order
    if (abs(Ap%coeff(ii)).GT.zero_test_small) then
      value_0(sturm_coeff)=Ap%coeff(ii)
      exit
    end if
  end do
  
  if (local_verbose) then
    write(*,*)'A0(x)='
    CALL write_poly_local2('A','x',Ap)  
  end if
  
  if (Ap%order.Eq.0) then
! this is a zero order function so exit
    if (local_verbose) then
      write(*,*)'Zero order function, exiting'
    end if
    RETURN
  end if

! second term is the derivative of Ap wrt x

  Am=allocate_polynomial(ap%order-1)

  do i=0,am%order
    am%coeff(i)=ap%coeff(i+1)*dble(i+1)
  end do
  
  sturm_coeff=sturm_coeff+1

! Get a function sign for x->infinity
  value_inf(sturm_coeff)=0d0
  do ii=Am%order,0,-1
    if (abs(Am%coeff(ii)).GT.zero_test_small) then
      value_inf(sturm_coeff)=Am%coeff(ii)
      exit
    end if
  end do
  
! Get a function sign for x->0
  value_0(sturm_coeff)=0d0
  do ii=0,Am%order
    if (abs(Am%coeff(ii)).GT.zero_test_small) then
      value_0(sturm_coeff)=Am%coeff(ii)
      exit
    end if
  end do
  
  if (local_verbose) then
    write(*,*)'A1(x)='
    CALL write_poly_local2('A','x',Am)  
  end if
  
! now iterate to get the remaining Sturm coefficients

  do i=3,n_sturm
  
! we have the two previous terms Ap and Am. 
! The next term is found as the remainder when Ap is divided by Am

!   R(x)=Ap(x)-Am(x)(C(x)) where C(x)=k1x+k2
!   New Ap=Am
!   New Am=-R

    CALL divide_poly(Ap,Am,C,R,.FALSE.)
    
    if (polynomial_is_zero(R)) then
      n_sturm=i-1
      if (local_verbose) then
        write(*,*)'Zero remainder in Sturm sequence calculation'
        write(*,*)'Number of Sturm coefficients=',n_sturm
      end if
      exit
    end if
    
    deallocate( Ap%coeff )
    Ap=Am
    deallocate( Am%coeff )
    Am=R
    do ii=0,Am%order
      Am%coeff(ii)=-Am%coeff(ii)
    end do
    deallocate( R%coeff )
    deallocate( C%coeff )
  
    sturm_coeff=sturm_coeff+1

! Get a function sign for x->infinity
    value_inf(sturm_coeff)=0d0
    do ii=Am%order,0,-1
      if (abs(Am%coeff(ii)).GT.zero_test_small) then
        value_inf(sturm_coeff)=Am%coeff(ii)
        exit
      end if
    end do
    
! Get a function sign for x->0
    value_0(sturm_coeff)=0d0
    do ii=0,Am%order
      if (abs(Am%coeff(ii)).GT.zero_test_small) then
        value_0(sturm_coeff)=Am%coeff(ii)
        exit
      end if
    end do
    
    if (local_verbose) then
      write(*,'(A1,I1,A4)')'A',i-1,'(x)='
      CALL write_poly_local2('A','x',Am)  
    end if
   
  end do
     
! Work out the number of sign changes in value_inf and value_0 in the sequence

  if (local_verbose) then
    write(*,*)'Number of Sturm coefficients=',n_sturm
    write(*,*)'Sturm coefficients'
    write(*,*)'           A(0)            A(infinity)       '
    do i=1,n_sturm
      write(*,*)value_0(i),value_inf(i)
    end do
  end if

  S_inf=0
  S_0=0
  
  last_value_inf=value_inf(1)
  last_value_0=value_0(1)
  
  do i=2,n_sturm
    
    if (value_0(i)*last_value_0.LT.0d0) S_0=S_0+1
    last_value_0=value_0(i)
  
    if (value_inf(i)*last_value_inf.LT.0d0) S_inf=S_inf+1
    last_value_inf=value_inf(i)
    
  end do
  
  n_real_zeros=abs(S_inf-S_0)
  
  if (local_verbose) then
    write(*,8000)'n_sign changes_0=',S_0,'  n_sign_changes_inf=',S_inf
8000 format(A,I3,A,I3)
    write(*,*)'Number of real zeros =',n_real_zeros
    write(*,*)'Odd_multiplicity_flag=',odd_multiplicity_flag
  end if
  
  if (n_real_zeros.GT.0) stable=.FALSE.
  
! check the consistency between the Sturm test and 
! the previously calculated odd_multiplicity_flag

  if ( (n_real_zeros.GT.0).AND.(odd_multiplicity_flag) ) then
! consistent result
    if (local_verbose) write(*,*)'Sturm test and explicit zero analysis agree: stable=.FALSE.'
  else if ( (n_real_zeros.EQ.0).AND.(.NOT.odd_multiplicity_flag) ) then
! consistent result
    if (local_verbose) write(*,*)'Sturm test and explicit zero analysis agree: stable=.TRUE.'
  else
! inconsistent result
    write(*,*)'Sturm test and explicit zero analysis do NOT agree'
    write(*,8000)'n_sign changes_0=',S_0,'  n_sign_changes_inf=',S_inf
    write(*,*)'Number of real zeros =',n_real_zeros
    write(*,*)'Odd_multiplicity_flag=',odd_multiplicity_flag
!    STOP
  end if
  
! Deallocate arrays and polynomials

  deallocate( value_inf )
  deallocate( value_0 )
  
  RETURN
  END SUBROUTINE test_filter_positive_real
