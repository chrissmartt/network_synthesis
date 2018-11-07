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
! Description
!     type definitions for rational function filter function representations
!     to be included in filter_module.F90
!
! Comments:
!      
!
! History
!
!     started 22/01/09 CJS
!

TYPE::Sfilter

  real(dp)  :: wnorm        ! frequency normalisation constants
  TYPE(polynomial) :: a     ! numerator polynomial function 
  TYPE(polynomial) :: b     ! denominator polynomial function
  
END TYPE Sfilter

TYPE::Sfilter_PR

  real(dp)  :: wnorm
  INTEGER   :: order
  INTEGER   :: n_complex_poles
  INTEGER   :: n_complex_pole_pairs
  INTEGER   :: n_real_poles
  real(dp)  :: R
  real(dp)  :: L
  LOGICAL,allocatable     :: complex_pole(:)
  complex(dp),allocatable :: poles(:)
  complex(dp),allocatable :: residues(:)
  
END TYPE Sfilter_PR


TYPE:: Sfilter_matrix

! matrix of rational function filters. These are used for frequency dependent admittance and impedance functions

  integer :: dim
  type(Sfilter),allocatable :: sfilter_mat(:,:)

END TYPE Sfilter_matrix

