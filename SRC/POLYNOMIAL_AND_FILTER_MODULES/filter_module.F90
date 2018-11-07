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
! Name
!    filter_module
!
! Description
!     module which combines the operators, functions and subroutines relating to polynomials and rational function filters
!
! Comments:
!     The contents of the module are in include files
!
! History
!
!     started 12/01/16 CJS
!

MODULE filter_module

USE type_specifications

IMPLICIT NONE

! OVERLOAD + to operate on polynomial data type and filter types

INTERFACE OPERATOR(+)
  MODULE PROCEDURE Addpoly, Addpoly_complex, AddSfilter
END INTERFACE

! OVERLOAD - to operate on polynomial data type

INTERFACE OPERATOR(-)
  MODULE PROCEDURE Subtractpoly, Subtractpoly_complex
END INTERFACE

! OVERLOAD * to operate on polynomial data type and filter types

INTERFACE OPERATOR(*)
  MODULE PROCEDURE Multiplypoly, Multiplypoly_complex, MultiplySfilter, ConstantMultiplySfilter
END INTERFACE

! OVERLOAD / to operate on filter data type

INTERFACE OPERATOR(/)
  MODULE PROCEDURE DivideSfilter
END INTERFACE

! OVERLOAD = to assign a value to polynomial types and filter types

INTERFACE ASSIGNMENT(=)
  MODULE PROCEDURE Assignpoly, Assignpoly2, Assignpoly_complex, Assignpoly_complex2, &
                   AssignSfilter, AssignSfilter2
                   
END INTERFACE

INCLUDE 'polynomial_types.F90'

INCLUDE 'filter_types.F90'

CONTAINS

INCLUDE 'polynomial_functions.F90'

INCLUDE 'polynomial_operators.F90'

INCLUDE 'polynomial_subroutines.F90'

INCLUDE 'filter_functions.F90'

INCLUDE 'filter_operators.F90'

INCLUDE 'filter_subroutines.F90'

INCLUDE 'findroots.F90'

END MODULE filter_module

