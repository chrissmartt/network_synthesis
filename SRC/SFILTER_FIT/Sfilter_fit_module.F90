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
! FILE CONTENTS (in include files)
!
! SUBROUTINE Calculate_Sfilter
! SUBROUTINE Weiner_Hopf
! SUBROUTINE filter_MSE
!
! Name
!    Sfilter_fit_module
!
! Description
!     module which provides a best fit rational function to the provided
!     complex frequency domain data 
!
! Comments:
!      
!
! History
!
!     started 28/04/16 CJS based on FILTER_FIT from the GGI_TLM project 
!                                   www.github.com/ggiemr/GGI_TLM
!

MODULE Sfilter_fit_module

IMPLICIT NONE

CONTAINS

include 'Calculate_Sfilter.F90'

include 'Weiner_Hopf.F90'

include 'filter_MSE.F90'

END MODULE Sfilter_fit_module
