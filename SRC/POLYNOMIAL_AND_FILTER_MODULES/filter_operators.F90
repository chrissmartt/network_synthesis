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
!+
! Name
!     filter_operators -  overloading operators for filter functions
!
! Description
!     This file defines arithmetic operations with filter functions
!     to be included in filter_module.F90 along with the associated INTERFACE 
!
! Public
!     The following are defined in this include file
!     
!     +
!     *
!     /
!     =
!
! History
!
!     started 22/01/09 CJS
!     include ADDfilter stuff 1/5/2015 CJS
!
!  FILE CONTAINS:
! FUNCTION AddSfilter
! FUNCTION MultiplySfilter
! FUNCTION ConstantMultiplySfilter
! FUNCTION DivideSfilter
!
! NAME
!     AddSfilter
!
! DESCRIPTION
!     Add Sfilter types
! 
!
! HISTORY
!
!     started 1/05/2015 CJS
!

FUNCTION AddSfilter(s1_in,s2_in) RESULT(res)

USE type_specifications

IMPLICIT NONE

! argument types
type(Sfilter), intent(IN)  :: s1_in,s2_in
! Result type

type(Sfilter)              :: res

type(polynomial) :: num1,num2

! local variables

type(Sfilter) :: s1,s2

integer aorder1,aorder2
integer aorder,border
integer order
real(dp)    term

real(dp) Wnorm_scale,new_wnorm

type(polynomial) :: num_term1,num_term2

! function definition

! copy the input argument filters

S1=S1_in
S2=S2_in

! the order of the denominator is the sum of the orders of the argument denominators
border=s1%b%order+s2%b%order

! the order of the numerator is the maximum of a1order+b2order and a2order+b1order
aorder1=s1%a%order+s2%b%order
aorder2=s2%a%order+s1%b%order
aorder=max(aorder1,aorder2)

res%a%order=aorder
res%b%order=border

if (s1%wnorm.ne.s2%wnorm) then

! we need to re-scale the w normalisation of the filters before they are combined
! choose the normalisation of the combined filter to be the maximum wnorm value of the two input filters
  new_wnorm=max(s1%wnorm,s2%wnorm)

! scale the s1 filter to this normalisation value
  wnorm_scale=1d0
  do order=0,s1%a%order
    s1%a%coeff(order)=s1%a%coeff(order)*wnorm_scale
    wnorm_scale=wnorm_scale*new_wnorm/s1%wnorm
  end do
  wnorm_scale=1d0
  do order=0,s1%b%order
    s1%b%coeff(order)=s1%b%coeff(order)*wnorm_scale
    wnorm_scale=wnorm_scale*new_wnorm/s1%wnorm
  end do

! scale the s2 filter to this normalisation value
  wnorm_scale=1d0
  do order=0,s2%a%order
    s2%a%coeff(order)=s2%a%coeff(order)*wnorm_scale
    wnorm_scale=wnorm_scale*new_wnorm/s2%wnorm
  end do
  wnorm_scale=1d0
  do order=0,s2%b%order
    s2%b%coeff(order)=s2%b%coeff(order)*wnorm_scale
    wnorm_scale=wnorm_scale*new_wnorm/s2%wnorm
  end do
  
  s1%wnorm=new_wnorm
  s2%wnorm=new_wnorm
      
end if

res=allocate_Sfilter(res%a%order,res%b%order)
res%wnorm=s1%wnorm

! multiply numerator polynomials to give the numerator result
res%b=s1%b*s2%b

! calculate the two numerator terms
num_term1=s1%a*s2%b
num_term2=s2%a*s1%b

! sum the two numerator polynomials to give the final numerator result

do order=0,aorder
  term=0d0
  if (order.LE.aorder1) then
! add the a1 term
    term=term+num_term1%coeff(order)
  end if
  if (order.LE.aorder2) then
! add the a1 term
    term=term+num_term2%coeff(order)
  end if
  
  res%a%coeff(order)=term
end do ! next term of numerator

CALL deallocate_Sfilter(s1)
CALL deallocate_Sfilter(s2)

END FUNCTION AddSfilter

!*
! NAME
!     MultiplySfilter
!
! DESCRIPTION
!     Multiply Sfilter types
!
! HISTORY
!
!     started 22/01/09 CJS
!     allow different wnorm values by renormalising one of the filters 21/04/16 CJS
!

FUNCTION MultiplySfilter(s1_in,s2) RESULT(res)

USE type_specifications

IMPLICIT NONE

! argument types

type(Sfilter), intent(IN)  :: s1_in,s2
! Result type
type(Sfilter)              :: res

! local variables

type(Sfilter) :: s1
real(dp) :: renorm_factor
integer  :: i

! function definition

s1=s1_in  ! copy filter s1_in to a local filter as it may have to be renormalised

res%a%order=s1%a%order+s2%a%order
res%b%order=s1%b%order+s2%b%order

if (s1%wnorm.ne.s2%wnorm) then
! wnorm discrepancy here we choose to renormalise s1 to wnorm of S2
  
  renorm_factor=s2%wnorm/s1%wnorm
  
  do i=1,s1%a%order
    s1%a%coeff(i)=s1%a%coeff(i)*(renorm_factor)**i
  end do
  
  do i=1,s1%b%order
    s1%b%coeff(i)=s1%b%coeff(i)*(renorm_factor)**i
  end do
  
  s1%wnorm=s2%wnorm

end if

res=allocate_Sfilter(res%a%order,res%b%order)
res%wnorm=s1%wnorm

! multiply polynomials top and bottom to give the result
res%a=s1%a*s2%a
res%b=s1%b*s2%b

CALL deallocate_Sfilter(s1)

END FUNCTION MultiplySfilter


!*
! NAME
!     ConstantMultiplySfilter
!
! DESCRIPTION
!     Multiply Sfilter type by a real constant
!
! HISTORY
!
!     started 21/04/16 CJS
!

FUNCTION ConstantMultiplySfilter(c1,s1) RESULT(res)

USE type_specifications

IMPLICIT NONE

! argument types

real(dp), intent(IN)       :: c1
type(Sfilter), intent(IN)  :: s1

! Result type

type(Sfilter)              :: res

! local variables

integer :: i

! function definition

res%a%order=s1%a%order
res%b%order=s1%b%order

res=allocate_Sfilter(res%a%order,res%b%order)
res%wnorm=s1%wnorm

! multiply numerator polynomial by the constant to give the result
do i=0,res%a%order
  res%a%coeff(i)=c1*s1%a%coeff(i)
end do
res%b=s1%b

END FUNCTION ConstantMultiplySfilter

! ---------------------------------------------------------
!*
! NAME
!     DivideSfilter
!
! DESCRIPTION
!     Divide Sfilter types
!
! HISTORY
!
!     started 26/01/09 CJS
!

FUNCTION DivideSfilter(s1,s2) RESULT(res)

USE type_specifications

IMPLICIT NONE

! argument types
type(Sfilter), intent(IN)  :: s1,s2
! Result type
type(Sfilter)              :: res

! function definition

res%a%order=s1%a%order+s2%b%order
res%b%order=s1%b%order+s2%a%order

if (s1%wnorm.ne.s2%wnorm) then
  write(*,*)'wnorm discrepancy in DivideSfilter'
  stop
end if

res=allocate_Sfilter(res%a%order,res%b%order)
res%wnorm=s1%wnorm

! Divide polynomials top and bottom to give the result
res%a=s1%a*s2%b
res%b=s1%b*s2%a

END FUNCTION DivideSfilter


! ---------------------------------------------------------

SUBROUTINE AssignSfilter(res,s)

USE type_specifications

IMPLICIT NONE

! Result type
type(Sfilter)              :: res
! argument types
real(dp)     , intent(IN)  :: s

! SUBROUTINE definition
  
  if (ALLOCATED(res%a%coeff)) DEALLOCATE( res%a%coeff )
  if (ALLOCATED(res%b%coeff)) DEALLOCATE( res%b%coeff )

  res%wnorm=1d0
  res%a%order=0
  allocate (res%a%coeff(0:res%a%order))
  res%b%order=0
  allocate (res%b%coeff(0:res%b%order))
  res%a%coeff(0:res%a%order)=s
  res%b%coeff(0:res%b%order)=1d0

END SUBROUTINE AssignSfilter

! ---------------------------------------------------------

SUBROUTINE AssignSfilter2(res,s)

USE type_specifications

IMPLICIT NONE

! Result type
type(Sfilter) , intent(OUT) :: res
! argument types
type(Sfilter) , intent(IN)  :: s

! SUBROUTINE definition
  
  if (ALLOCATED (res%a%coeff)) DEALLOCATE (res%a%coeff)
  ALLOCATE (res%a%coeff(0:s%a%order))

  if (ALLOCATED (res%b%coeff)) DEALLOCATE (res%b%coeff)
  ALLOCATE (res%b%coeff(0:s%b%order))
  
  
! set values  
  res%wnorm=s%wnorm
  res%a%order=s%a%order
  res%a%coeff(0:res%a%order)=s%a%coeff(0:s%a%order)
  res%b%order=s%b%order
  res%b%coeff(0:res%b%order)=s%b%coeff(0:s%b%order)

END SUBROUTINE AssignSfilter2


