!
! This file is part of SACAMOS, State of the Art CAble MOdels in Spice. 
! It was developed by the University of Nottingham and the Netherlands Aerospace 
! Centre (NLR) for ESA under contract number 4000112765/14/NL/HK.
! 
! Copyright (C) 2016-2017 University of Nottingham
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
! File Contents:
! PROGRAM create_test_case
!
! NAME
!     create_test_case
!
! AUTHORS
!     Chris Smartt
!
! DESCRIPTION
!     This program creates a transfer impedance model as in SACAMOS.
!     The diffusion impedance comes from a L-R ladder network transfer impedance model
!     then braid/ hole inductance is added to this to give Zt(s).
!
!     The SACAMOS model usually requires the subtraction of the d.c. resistance (optional here)
!     before requiring the implementation of H(s)=(1/s)(Zt(s)-Rdc)*H(s) where H(s) is
!     the propagation correction.
!
!
! COMMENTS
!     
!
! HISTORY
!
!     started 26/09/2017 CJS 
!     29/09/2017 CJS add a hole/ braid inductance term
!     2-4/10/2017 CJS include optional d.c. resistance removal, multiplication by 1/s
!                     then provide two impedances Za and Zb which are both positive-real
!                     and can be implemented in a bridge circuit to implement Zt(s)
!                     Also add a propagation correction Hp(s)
!
!
!
PROGRAM create_ZT_test_case

USE type_specifications
USE general_module
USE constants
USE frequency_spec
USE filter_module
USE Sfilter_fit_module

IMPLICIT NONE

! local variables

! command line argument value and length

! command line argument value and length
character(len=filename_length)    :: argument1
integer                           :: argument1_length

character(len=filename_length)    :: filename1   ! filename for impedance Za
character(len=filename_length)    :: filename2   ! filename for impedance Zb
character(len=filename_length)    :: filename3   ! filename for impedance Z=Za-Zb

integer     :: N   ! number of R-L sections
real(dp)    :: R,L,R0,L0,Rsection,Lbh,b0,R_min,w_R_min

complex(dp):: jw,num,den,res

type(Sfilter) :: a11,a12,a21,a22
type(Sfilter) :: m11,m12,m21,m22

type(Sfilter) :: T1,T2,T3,T4,T5,T6,T7,T8 

type(Sfilter) :: Y
type(Sfilter) :: Zd
type(Sfilter) :: Zbh
type(Sfilter) :: Zt
type(Sfilter) :: Zt2
type(Sfilter) :: Za
type(Sfilter) :: Zb
type(Sfilter) :: integrate
type(Sfilter) :: Rdc_filter
type(Sfilter) :: F
type(Sfilter) :: Hp

logical :: subtract_Rdc
logical :: integrate_function
logical :: apply_propagation_correction

integer :: i

! START

  program_name="create_ZT_test_case"
  run_status='Started'
  CALL write_program_status()
  
  CALL get_command_argument(1 , argument1, argument1_length)

! Set the number of sections
! required for the impedance calculation

  filename1='test_ZA'
  filename2='test_ZB'
  filename3='test_Z'

  N=2

! Set the component values
  R0=0.1d0      ! diffusion resistance
  L0=1d0        ! diffusion inductance
  Lbh=-0.1d0    ! braid+hole transfer inductance (can be -ve)

  subtract_Rdc=.FALSE.
  integrate_function=.FALSE.
  apply_propagation_correction=.FALSE.
  
! comment out the following to set to .FALSE. 

  subtract_Rdc=.TRUE.
!  apply_propagation_correction=.TRUE.
  
  integrate_function=.TRUE.

! automatically set the configurations
  if (argument1_length.NE.0) then
  
    if(argument1(1:1).EQ.'1') then
      subtract_Rdc=.FALSE.
      apply_propagation_correction=.FALSE.
      Lbh=+abs(Lbh)
      
    else if (argument1(1:1).EQ.'2') then
      subtract_Rdc=.FALSE.
      apply_propagation_correction=.FALSE. 
      Lbh=-abs(Lbh)
      
    else if (argument1(1:1).EQ.'3') then
      subtract_Rdc=.FALSE.
      apply_propagation_correction=.TRUE.  
      Lbh=+abs(Lbh)
      
    else if (argument1(1:1).EQ.'4') then
      subtract_Rdc=.FALSE.
      apply_propagation_correction=.TRUE.  
      Lbh=-abs(Lbh)
      
    else if (argument1(1:1).EQ.'5') then
      subtract_Rdc=.TRUE.
      apply_propagation_correction=.FALSE. 
      Lbh=+abs(Lbh)
      
    else if (argument1(1:1).EQ.'6') then
      subtract_Rdc=.TRUE.
      apply_propagation_correction=.FALSE.  
      Lbh=-abs(Lbh)
      
    else if (argument1(1:1).EQ.'7') then
      subtract_Rdc=.TRUE.
      apply_propagation_correction=.TRUE.  
      Lbh=+abs(Lbh)
      
    else if (argument1(1:1).EQ.'8') then
      subtract_Rdc=.TRUE.
      apply_propagation_correction=.TRUE.  
      Lbh=-abs(Lbh)
  
  
    end if
  end if


! first order propagation correction function
  Hp=allocate_Sfilter(1,1)
  Hp%a%coeff(0)=1d0
  Hp%a%coeff(1)=0.2349292019D+01
  Hp%b%coeff(0)=1d0
  Hp%b%coeff(1)=0.2385758553D+01

! get the R and L values in each section of the ladder network  
  R=R0*N
  L=L0/N
    
! Calculate the four elements of the M matrix

  if(N.EQ.1) then
! final section has R*2
    Rsection=R*2d0
  else
    Rsection=R
  end if

! a11=(R+sL)/R
  a11=allocate_Sfilter(1,0)
  a11%a%coeff(0)=Rsection
  a11%a%coeff(1)=L
  a11%b%coeff(0)=Rsection

! a12=1/R
  a12=allocate_Sfilter(0,0)
  a12%a%coeff(0)=1d0
  a12%b%coeff(0)=Rsection

! a21=sL
  a21=allocate_Sfilter(1,0)
  a21%a%coeff(0)=0d0
  a21%a%coeff(1)=L
  a21%b%coeff(0)=1d0

! a22=1
  a22=allocate_Sfilter(0,0)
  a22%a%coeff(0)=1d0
  a22%b%coeff(0)=1d0
  
! Matrix for a single section, m=a

  m11=a11
  m12=a12
  m21=a21
  m22=a22
  
  do i=N-1,1,-1

    if(i.EQ.1) then
! final section has R*2
      Rsection=R*2d0
    else
      Rsection=R
    end if

! a11=(R+sL)/R
    a11=allocate_Sfilter(1,0)
    a11%a%coeff(0)=Rsection
    a11%a%coeff(1)=L
    a11%b%coeff(0)=Rsection

! a12=1/R
    a12=allocate_Sfilter(0,0)
    a12%a%coeff(0)=1d0
    a12%b%coeff(0)=Rsection

! a21=sL
    a21=allocate_Sfilter(1,0)
    a21%a%coeff(0)=0d0
    a21%a%coeff(1)=L
    a21%b%coeff(0)=1d0

! a22=1
    a22=allocate_Sfilter(0,0)
    a22%a%coeff(0)=1d0
    a22%b%coeff(0)=1d0
  
    write(*,*)'m11'
    CALL write_Sfilter(m11,0)
    write(*,*)'m12'
    CALL write_Sfilter(m12,0)
    write(*,*)'m21'
    CALL write_Sfilter(m21,0)
    write(*,*)'m22'
    CALL write_Sfilter(m22,0)
  
    T1=a11*m11
    T2=a12*m21
    T3=a11*m12
    T4=a12*m22
  
    T5=a21*m11
    T6=a22*m21
    T7=a21*m12
    T8=a22*m22
  
! calculate the next M matrix
    m11=T1+T2
    m12=T3+T4
    m21=T5+T6
    m22=T7+T8
 
    write(*,*)'new m11'
    CALL write_Sfilter(m11,0)
    write(*,*)'new m12'
    CALL write_Sfilter(m12,0)
    write(*,*)'new m21'
    CALL write_Sfilter(m21,0)
    write(*,*)'new m22'
    CALL write_Sfilter(m22,0)
  
  end do ! next section

! M is now the matrix for the complete ladder network
! evaluate the impedance function from the termination condition Vn=R*In
! note the final section has R*2
  m11%a%coeff(0:m11%a%order)=m11%a%coeff(0:m11%a%order)/(R*2d0)
  Y=m11+m12
  Zd=reciprocal_Sfilter(Y)
  
! Zd is the diffusion impedance which may have -ve Real and imaginary parts
  
! Calculate the contribution from the braid/hole inductance (which may be negative)
  Zbh=allocate_Sfilter(1,0)
  Zbh%a%coeff(0)=0d0
  Zbh%a%coeff(1)=Lbh
  Zbh%b%coeff(0)=1d0
  
! Zt is now the complete transfer impedance fuction as it would be
! specified in a SACAMOS cable model

  Zt=Zd+Zbh
  
  write(*,*)'Zt=Zd+Zbh'
  CALL write_Sfilter(Zt,0)
   
!  write(*,*)'integrate_filter'
  integrate=allocate_Sfilter(0,1)
  integrate%a%coeff(0)=1d0
  integrate%b%coeff(0)=0d0
  integrate%b%coeff(1)=1d0

!  write(*,*)'Rdc_filter'
  Rdc_filter=allocate_Sfilter(0,0)
  Rdc_filter%a%coeff(0)=-R0
  Rdc_filter%b%coeff(0)=1d0

! The normal transfer impedance model requires the subtraction of the d.c. resistance
! i.e. it is extracted to the termination of the shield conductor
! though it may be left in in specific circumstances

  if(subtract_Rdc) then
  
    CALL deallocate_Sfilter(T1)
    T1=Zt+Rdc_filter
    Zt2=T1
    
  else
  
    Zt2=Zt
    
  end if

! Calculate impedances Za and Zb such that Zt2=Za-Zb and Za and Zb are both positive-real

! once we integrate, this simply becomes a negative resistance which is dealt with below.
! If Lbh is -ve then set Za=Zt2+abs(Zbh) and Zb=+abs(Zbh) so Zt=Za-Zb is unchanged
  if ( (Lbh.LT.0d0).AND.(.NOT.integrate_function ) )then  
    Zbh%a%coeff(1)=abs(Lbh)
    Za=Zt2+Zbh
    Zb=Zbh
  else
    Za=Zt2
    Zb=0d0
  end if
   
  if (integrate_function) then
! multiply by 1/s i.e. time integration term
          
    write(*,*)'Apply Time integration function, integral='
    CALL write_Sfilter(integrate,0)
    
    CALL deallocate_Sfilter(T1)
    T1=Za*integrate
    CALL deallocate_Sfilter(Za)
    Za=T1
    
    CALL deallocate_Sfilter(T1)
    T1=Zb*integrate
    CALL deallocate_Sfilter(Zb)
    Zb=T1
    
    CALL deallocate_Sfilter(T1)
    T1=Zt2*integrate
    CALL deallocate_Sfilter(Zt2)
    Zt2=T1

  end if
    
  if (apply_propagation_correction) then
! multiply by Hp(s) i.e. apply the propagation correction
          
    write(*,*)'Apply Propagation correction, Hp='
    CALL write_Sfilter(Hp,0)
    
    write(*,*)'Za='
    CALL write_Sfilter(Za,0)

    CALL deallocate_Sfilter(T1)
    T1=Za*Hp
    CALL deallocate_Sfilter(Za)
    Za=T1
    
    write(*,*)'ZaHp='
    CALL write_Sfilter(Za,0)
    
    CALL deallocate_Sfilter(T1)
    T1=Zb*Hp
    CALL deallocate_Sfilter(Zb)
    Zb=T1
    
    CALL deallocate_Sfilter(T1)
    T1=Zt2*Hp
    CALL deallocate_Sfilter(Zt2)
    Zt2=T1
      
  end if

! Try to ensure that all the functions are in their minimum form.  
  CALL get_min_order_poly(Za%a)
  CALL get_min_order_poly(Za%b)
  CALL pole_zero_cancel(Za)
   
  CALL get_min_order_poly(Zb%a)
  CALL get_min_order_poly(Zb%b)
  CALL pole_zero_cancel(Zb)
 
  CALL get_min_order_poly(Zt2%a)
  CALL get_min_order_poly(Zt2%b)
  CALL pole_zero_cancel(Zt2)

! If the minimum resistance, R_min, is negative then
! add 2*R_min to Za make the function positive-real, also add this to Zb so Zt=Za-Zb is unchanged

! Calculate R_min directly from the filter function
   
  CALL deallocate_Sfilter(T1)
  T1=Za
  write(*,*)'Calculate minimum resistance of function:'
  CALL write_Sfilter(T1,0)
  CALL calculate_min_resistance_value(T1,R_min,w_R_min)
  
  write(*,*)'Minimum Resistance value is',R_min
  write(*,*)'at w=',w_R_min,' f=',w_R_min/(6.28)
  
  if (R_min.LT.0d0) then
  
    Rdc_filter%a%coeff(0)=2d0*abs(R_min)
  
    CALL deallocate_Sfilter(T1)
    T1=Za+Rdc_filter
    CALL deallocate_Sfilter(Za)
    Za=T1
  
    CALL deallocate_Sfilter(T1)
    T1=Zb+Rdc_filter
    CALL deallocate_Sfilter(Zb)
    Zb=T1

  end if

! Clean up the filter functions
  CALL get_min_order_poly(Za%a)
  CALL get_min_order_poly(Za%b)
  CALL pole_zero_cancel(Za)

  b0=0d0
  do i=0,Za%b%order
    b0=max( abs(Za%b%coeff(i)) , b0 )
  end do
  if (b0.NE.0d0) then
     Za%a%coeff(:)=Za%a%coeff(:)/b0
     Za%b%coeff(:)=Za%b%coeff(:)/b0
  end if
      
  write(*,*)'Za='
  CALL write_Sfilter(Za,0)
  
  CALL get_min_order_poly(Zb%a)
  CALL get_min_order_poly(Zb%b)
  CALL pole_zero_cancel(Zb)

  b0=0d0
  do i=0,Zb%b%order
    b0=max( abs(Zb%b%coeff(i)) , b0 )
  end do
  if (b0.NE.0d0) then
     Zb%a%coeff(:)=Zb%a%coeff(:)/b0
     Zb%b%coeff(:)=Zb%b%coeff(:)/b0
  end if
      
  write(*,*)'Zb='
  CALL write_Sfilter(Zb,0)

  CALL get_min_order_poly(Zt2%a)
  CALL get_min_order_poly(Zt2%b)
  CALL pole_zero_cancel(Zt2)
  
  b0=0d0
  do i=0,Zt2%b%order
    b0=max( abs(Zt2%b%coeff(i)) , b0 )
  end do
  if (b0.NE.0d0) then
     Zt2%a%coeff(:)=Zt2%a%coeff(:)/b0
     Zt2%b%coeff(:)=Zt2%b%coeff(:)/b0
  end if
      
  write(*,*)'Zt2='
  CALL write_Sfilter(Zt2,0)

! Open the input file for the impedance function
  write(*,*)'Open the input file for the impedance function Za'

  open(unit=10,file=filename1)
  
  write(10,*)Za%a%order,'     # aorder'
  write(10,*)Za%b%order,'     # border'
  write(10,*)Za%wnorm  ,'     # wnorm'

  write(10,*)'# a coefficients ' 
  
  do i=0,Za%a%order
    write(10,*)Za%a%coeff(i)
  end do

  write(10,*)'# b coefficients ' 
  
  do i=0,Za%b%order
    write(10,*)Za%b%coeff(i)
  end do
  
  write(10,*)'# wmin,wmax,nw'
  write(10,*)'0.01 4.0 200'
  
  close(unit=10)
 
! Open the input file for the impedance function
  write(*,*)'Open the input file for the impedance function Zb'

  open(unit=10,file=filename2)
  
  write(10,*)Zb%a%order,'     # aorder'
  write(10,*)Zb%b%order,'     # border'
  write(10,*)Zb%wnorm  ,'     # wnorm'

  write(10,*)'# a coefficients ' 
  
  do i=0,Zb%a%order
    write(10,*)Zb%a%coeff(i)
  end do

  write(10,*)'# b coefficients ' 
  
  do i=0,Zb%b%order
    write(10,*)Zb%b%coeff(i)
  end do
  
  write(10,*)'# wmin,wmax,nw'
  write(10,*)'0.01 4.0 200'
  
  close(unit=10)
  
! Open the input file for the impedance function
  write(*,*)'Open the input file for the impedance function Zt2'

  open(unit=10,file=filename3)
  
  write(10,*)Zt2%a%order,'     # aorder'
  write(10,*)Zt2%b%order,'     # border'
  write(10,*)Zt2%wnorm  ,'     # wnorm'

  write(10,*)'# a coefficients ' 
  
  do i=0,Zt2%a%order
    write(10,*)Zt2%a%coeff(i)
  end do

  write(10,*)'# b coefficients ' 
  
  do i=0,Zt2%b%order
    write(10,*)Zt2%b%coeff(i)
  end do
  
  write(10,*)'# wmin,wmax,nw'
  write(10,*)'0.01 4.0 200'
  
  close(unit=10)

! deallocate arrays
  write(*,*)'deallocate arrays'
  
  CALL deallocate_Sfilter(a11)
  CALL deallocate_Sfilter(a12)
  CALL deallocate_Sfilter(a21)
  CALL deallocate_Sfilter(a22)
  CALL deallocate_Sfilter(m11)
  CALL deallocate_Sfilter(m12)
  CALL deallocate_Sfilter(m21)
  CALL deallocate_Sfilter(m22)
  CALL deallocate_Sfilter(T1)
  CALL deallocate_Sfilter(T2)
  CALL deallocate_Sfilter(T3)
  CALL deallocate_Sfilter(T4)
  CALL deallocate_Sfilter(T5)
  CALL deallocate_Sfilter(T6)
  CALL deallocate_Sfilter(T7)
  CALL deallocate_Sfilter(T8)
  CALL deallocate_Sfilter(Y)
  CALL deallocate_Sfilter(Zd)
  CALL deallocate_Sfilter(Zt)
  CALL deallocate_Sfilter(Rdc_filter)
  CALL deallocate_Sfilter(integrate)
  CALL deallocate_Sfilter(F)
  CALL deallocate_Sfilter(Zt2)
  CALL deallocate_Sfilter(Za)
  CALL deallocate_Sfilter(Zb)
  CALL deallocate_Sfilter(Hp)
        
  run_status='Finished_Correctly'
  CALL write_program_status()

END PROGRAM create_ZT_test_case

include 'include_local_filter_subroutines.F90'

INCLUDE 'include_test_positive_real.F90'

INCLUDE 'include_minimum_resistance_function.F90'
