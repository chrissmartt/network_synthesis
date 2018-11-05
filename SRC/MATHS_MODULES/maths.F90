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
! Name
!    maths
!
! Description
!     module which combines the modules relating to mathematical processes required
!     for the creation of spice cable bundle models and the solution of multi-conductor
!     transmission line equations.
!     This module does not include eigenvalue/eigenvector solutions- these are in the EISPACK directory
!
!     The subroutines within this module are contained in include files as follows:
!
! dmatrix.F90: SUBROUTINE dwrite_matrix
! dmatrix.F90: SUBROUTINE dread_matrix
! dmatrix.F90: SUBROUTINE dinvert_Gauss_Jordan
!
! cmatrix.F90: SUBROUTINE write_cmatrix
! cmatrix.F90: SUBROUTINE write_cmatrix_re
! cmatrix.F90: SUBROUTINE write_cmatrix_im
! cmatrix.F90: SUBROUTINE cinvert_Gauss_Jordan
! cmatrix.F90: SUBROUTINE c_condition_number
!
! FFT.F90: SUBROUTINE FFT_TIME_TO_FREQ
! FFT.F90: SUBROUTINE FFT_FREQ_TO_TIME
! FFT.F90: SUBROUTINE FFT

!
! Comments:
!     
!
! History
!
!     started 8/01/16 CJS
!

MODULE maths

USE type_specifications

IMPLICIT NONE

CONTAINS

include 'dmatrix.F90'

include 'cmatrix.F90'

include 'FFT.F90'

END MODULE maths
