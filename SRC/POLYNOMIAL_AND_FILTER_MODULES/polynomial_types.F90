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
!    polynomial_types
!
! Description
!     type definitions for polynomials 
!     to be included in filter_module.F90
!
! Comments:
!      
!
! History
!
!     started 23/01/09 CJS
!

TYPE::polynomial
  INTEGER         :: order          ! order of the polynomial
  real(dp),allocatable :: coeff(:)  ! real polynomial coefficients: f(x)=coeff(0)+coeff(1)*x+coeff(2)*(x**2)+...+coeff(order)*(x**order)
END TYPE polynomial

TYPE::complex_polynomial
  INTEGER         :: order             ! order of the polynomial
  complex(dp),allocatable :: coeff(:)  ! complex polynomial coefficients: f(x)=coeff(0)+coeff(1)*x+coeff(2)*(x**2)+...+coeff(order)*(x**order)
END TYPE complex_polynomial

