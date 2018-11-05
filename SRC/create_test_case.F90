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
!     Create input files for circuit synthesis test cases as rational functions of jw
!
!
! COMMENTS
!     
!
! HISTORY
!
!     started 18/09/2017 CJS 
!!
!
PROGRAM create_test_case

USE type_specifications
USE general_module
USE constants
USE frequency_spec
USE filter_module
USE Sfilter_fit_module

IMPLICIT NONE

! local variables

! command line argument value and length

character(len=filename_length)    :: filename   ! filename for the transfer function

type(Sfilter)    :: H                     ! rational transfer function model 

integer :: NR,NL,NC
real(dp),allocatable    :: R(:),L(:),C(:)

integer :: NZ,NY
type(Sfilter),allocatable  :: Z(:),Y(:)   ! rational impedance and admittance functions

integer :: i

integer :: test_case

! START

! set the test case number here

  test_case=2
  
  select case (test_case)
  
  case(1)

    include 'TEST_CASES/include_RL_C_R_parallel.F90'
  
  case(2)

    include 'TEST_CASES/include_RLRCR.F90'
    
  end select

! Open the input file for the impedance function
  write(*,*)'Open the input file for the impedance function'

  open(unit=10,file=filename)
  
  write(10,*)H%a%order,'     # aorder'
  write(10,*)H%b%order,'     # border'
  write(10,*)H%wnorm  ,'     # wnorm'

  write(10,*)'# a coefficients ' 
  
  do i=0,H%a%order
    write(10,*)H%a%coeff(i)
  end do

  write(10,*)'# b coefficients ' 
  
  do i=0,H%b%order
    write(10,*)H%b%coeff(i)
  end do
  
  write(10,*)'# wmin,wmax,nw'
  write(10,*)'0.01 4.0 200'
  
  close(unit=10)
  
! deallocate arrays
  write(*,*)'deallocate arrays'

  if (allocated( R )) deallocate( R )  
  if (allocated( L )) deallocate( L )  
  if (allocated( C )) deallocate( C )  
  
  do i=1,NZ
    if (allocated( Z(i)%a%coeff )) deallocate( Z(i)%a%coeff )
    if (allocated( Z(i)%b%coeff )) deallocate( Z(i)%b%coeff )
  end do
  if (allocated( Z )) deallocate( Z )
  
  do i=1,NY
    if (allocated( Y(i)%a%coeff )) deallocate( Y(i)%a%coeff )
    if (allocated( Y(i)%b%coeff )) deallocate( Y(i)%b%coeff )
  end do
  if (allocated( Y )) deallocate( Y )
  
  if (allocated( H%a%coeff )) deallocate( H%a%coeff )
  if (allocated( H%b%coeff )) deallocate( H%b%coeff )
       
END PROGRAM create_test_case


