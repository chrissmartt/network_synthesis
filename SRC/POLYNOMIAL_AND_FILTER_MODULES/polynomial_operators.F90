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
!     Operators -  overloading operators for polynomials
!
! Description
!     This file defines arithmetic operations for polynomials and complex polynomials
!     to be included in filter_module.F90 along with the associated INTERFACE 
!
!     The following are defined in this module
!
!     +
!     -
!     *
!     =
!
! History
!
!     started 23/01/09 CJS
!
!
! FILE CONTAINS:
!FUNCTION Addpoly(p1,p2)
!FUNCTION Addpoly_complex(p1,p2)
!FUNCTION Subtractpoly(p1,p2)
!FUNCTION Subtractpoly_complex(p1,p2)
!FUNCTION Multiplypoly(p1,p2)
!FUNCTION Multiplypoly_complex(p1,p2)
!SUBROUTINE Assignpoly(res,b)
!SUBROUTINE Assignpoly2(res,a)
!SUBROUTINE Assignpoly_complex(res,b)
!SUBROUTINE Assignpoly_complex2(res,a)
!
!+
! NAME
!     Addpoly
!
! DESCRIPTION
!     Add polynomial types
!
! HISTORY
!
!     started 23/01/09 CJS
!

FUNCTION Addpoly(p1,p2) RESULT(res)

USE type_specifications

IMPLICIT NONE

! argument types
type(polynomial), intent(IN)  :: p1,p2
! Result type
type(polynomial)              :: res

! local variable
type(polynomial) :: answer

integer max_order,i,o1,o2

! function definition
o1=p1%order
o2=p2%order
max_order=o1
if (o2.gt.o1) then
  max_order=o2
end if

answer=allocate_polynomial(max_order)
answer%order=max_order

! reset result
do i=0,max_order
  if (i.le.o1) answer%coeff(i)=answer%coeff(i)+p1%coeff(i)
  if (i.le.o2) answer%coeff(i)=answer%coeff(i)+p2%coeff(i)
end do

res=answer

DEALLOCATE( answer%coeff )
RETURN

END FUNCTION Addpoly

!+
! NAME
!     Addpoly_complex
!
! DESCRIPTION
!     Add complex polynomial types
!
! HISTORY
!
!     started 4/01/13 CJS
!

FUNCTION Addpoly_complex(p1,p2) RESULT(res)

USE type_specifications

IMPLICIT NONE

! argument types
type(complex_polynomial), intent(IN)  :: p1,p2
! Result type
type(complex_polynomial)              :: res

! local variable
type(complex_polynomial) :: answer

integer max_order,i,o1,o2

! function definition
o1=p1%order
o2=p2%order
max_order=o1
if (o2.gt.o1) then
  max_order=o2
end if

answer=allocate_complex_polynomial(max_order)
answer%order=max_order

! reset result
do i=0,max_order
  if (i.le.o1) answer%coeff(i)=answer%coeff(i)+p1%coeff(i)
  if (i.le.o2) answer%coeff(i)=answer%coeff(i)+p2%coeff(i)
end do

res=answer

DEALLOCATE( answer%coeff )
RETURN

END FUNCTION Addpoly_complex

! ---------------------------------------------------------


!+
! NAME
!     Subtractpoly
!
! DESCRIPTION
!     Subtract polynomial types
!
! HISTORY
!
!     started 23/01/09 CJS
!

FUNCTION Subtractpoly(p1,p2) RESULT(res)

USE type_specifications

IMPLICIT NONE

! argument types
type(polynomial), intent(IN)  :: p1,p2
! Result type
type(polynomial)              :: res

! local variable
type(polynomial) :: answer

integer max_order,i,o1,o2

! function definition
o1=p1%order
o2=p2%order
max_order=o1
if (o2.gt.o1) then
  max_order=o2
end if

answer=allocate_polynomial(max_order)
answer%order=max_order

! reset result
do i=0,max_order
  if (i.le.o1) answer%coeff(i)=answer%coeff(i)+p1%coeff(i)
  if (i.le.o2) answer%coeff(i)=answer%coeff(i)-p2%coeff(i)
end do

res=answer

DEALLOCATE( answer%coeff )
RETURN

END FUNCTION Subtractpoly


!+
! NAME
!     Subtractpoly_complex
!
! DESCRIPTION
!     Subtract complex polynomial types
!
! HISTORY
!
!     started 4/1/2013 CJS
!

FUNCTION Subtractpoly_complex(p1,p2) RESULT(res)

USE type_specifications

IMPLICIT NONE

! argument types
type(complex_polynomial), intent(IN)  :: p1,p2
! Result type
type(complex_polynomial)              :: res

! local variable
type(complex_polynomial) :: answer

integer max_order,i,o1,o2

! function definition
o1=p1%order
o2=p2%order
max_order=o1
if (o2.gt.o1) then
  max_order=o2
end if

answer=allocate_complex_polynomial(max_order)
answer%order=max_order

! reset result
do i=0,max_order
  if (i.le.o1) answer%coeff(i)=answer%coeff(i)+p1%coeff(i)
  if (i.le.o2) answer%coeff(i)=answer%coeff(i)-p2%coeff(i)
end do

res=answer

DEALLOCATE( answer%coeff )
RETURN

END FUNCTION Subtractpoly_complex


! ---------------------------------------------------------


!+
! NAME
!     Multiplypoly
!
! DESCRIPTION
!     Multiply polynomial types
!
! HISTORY
!
!     started 23/01/09 CJS
!

FUNCTION Multiplypoly(p1,p2) RESULT(res)

USE type_specifications

IMPLICIT NONE

! argument types
type(polynomial), intent(IN)  :: p1,p2
! Result type
type(polynomial)              :: res

! local variables
type(polynomial) :: answer

integer i1,j1,order

! function definition

order=p1%order+p2%order
answer=allocate_polynomial(order)
answer%order=order

  do i1=0,p1%order
    do j1=0,p2%order
! product of p1(i1) p2(j1)
      answer%coeff(i1+j1)=answer%coeff(i1+j1)+p1%coeff(i1)*p2%coeff(j1)
    end do
  end do  

res=answer

DEALLOCATE( answer%coeff )
RETURN

END FUNCTION Multiplypoly
!+
! NAME
!     Multiplypoly_complex
!
! DESCRIPTION
!     Multiply complex polynomial types
!
! HISTORY
!
!     started 4/1/2013 CJS
!

FUNCTION Multiplypoly_complex(p1,p2) RESULT(res)

USE type_specifications

IMPLICIT NONE

! argument types
type(complex_polynomial), intent(IN) :: p1,p2
! Result type
type(complex_polynomial)             :: res

! local variables
type(complex_polynomial) :: answer

integer i1,j1,order

! function definition

order=p1%order+p2%order
answer=allocate_complex_polynomial(order)
answer%order=order

  do i1=0,p1%order
    do j1=0,p2%order
! product of p1(i1) p2(j1)
      answer%coeff(i1+j1)=answer%coeff(i1+j1)+p1%coeff(i1)*p2%coeff(j1)
    end do
  end do  

res=answer

DEALLOCATE( answer%coeff )
RETURN

END FUNCTION Multiplypoly_complex
 
! ---------------------------------------------------------
! NAME
!     assignpoly
!
! DESCRIPTION
!     assign polynomial types: set to real
!
! COMMENTS
!     
! HISTORY
!
!     started 23/01/09 CJS
!

SUBROUTINE Assignpoly(res,b)

USE type_specifications

IMPLICIT NONE

! argument types
real(dp)        , intent(IN) :: b

! Result type
type(polynomial)             ::  res

integer order

! SUBROUTINE definition

  order=0
  if (.NOT. ALLOCATED (res%coeff) ) then
! this filter not yet allocated so allocate here
    res%order=order
    ALLOCATE (res%coeff(0:order))
  else if (res%order.ne.order) then
! this filter has already been allocated but to the wrong order so reallocate  
    res%order=order
    DEALLOCATE (res%coeff)
    ALLOCATE (res%coeff(0:order))
  end if
  
  res%order=0

  res%coeff(0)=b

END SUBROUTINE Assignpoly

! ---------------------------------------------------------
! NAME
!     assignpoly2
!
! DESCRIPTION
!     assign polynomial types: set to existing polynomial
!
! COMMENTS
!     
! HISTORY
!
!     started 23/01/09 CJS
!

SUBROUTINE Assignpoly2(res,a)

USE type_specifications

IMPLICIT NONE

! argument types
type(polynomial), intent(IN) :: a

! Result type
type(polynomial)             :: res

! SUBROUTINE definition

  if (.NOT. ALLOCATED (res%coeff)) then
! this filter not yet allocated so allocate here
    ALLOCATE (res%coeff(0:a%order))
  else if (res%order.ne.a%order) then
! this filter has already been allocated but to the wrong order so reallocate  
    DEALLOCATE (res%coeff)
    ALLOCATE (res%coeff(0:a%order))
  end if
  
  res%order=a%order
  res%coeff(0:res%order)=a%coeff(0:a%order)

END SUBROUTINE Assignpoly2

! ---------------------------------------------------------
! NAME
!     assignpoly_complex
!
! DESCRIPTION
!     assign complex polynomial types: set to real
!
! COMMENTS
!     
! HISTORY
!
!     started 4/1/2013 CJS
!

SUBROUTINE Assignpoly_complex(res,b)

USE type_specifications

IMPLICIT NONE

! argument types
complex(dp)             , intent(IN) :: b

! Result type
type(complex_polynomial)             :: res

integer order

! SUBROUTINE definition

  order=0
  if (.NOT. ALLOCATED (res%coeff) ) then
! this filter not yet allocated so allocate here
    res%order=order
    ALLOCATE (res%coeff(0:order))
  else if (res%order.ne.order) then
! this filter has already been allocated but to the wrong order so reallocate  
    res%order=order
    DEALLOCATE (res%coeff)
    ALLOCATE (res%coeff(0:order))
  end if
  
  res%order=0

  res%coeff(0)=b

END SUBROUTINE Assignpoly_complex

! ---------------------------------------------------------
! NAME
!     assignpoly_complex2
!
! DESCRIPTION
!     assign complex polynomial types: set to existing polynomial
!
! COMMENTS
!     
! HISTORY
!
!     started 4/1/2013 CJS
!

SUBROUTINE Assignpoly_complex2(res,a)

USE type_specifications

IMPLICIT NONE

! argument types
type(complex_polynomial), intent(IN) :: a

! Result type
type(complex_polynomial)             :: res

! SUBROUTINE definition

  if (.NOT. ALLOCATED (res%coeff)) then
! this filter not yet allocated so allocate here
    ALLOCATE (res%coeff(0:a%order))
  else if (res%order.ne.a%order) then
! this filter has already been allocated but to the wrong order so reallocate  
    DEALLOCATE (res%coeff)
    ALLOCATE (res%coeff(0:a%order))
  end if
  
  res%order=a%order
  res%coeff(0:res%order)=a%coeff(0:a%order)

END SUBROUTINE Assignpoly_complex2

