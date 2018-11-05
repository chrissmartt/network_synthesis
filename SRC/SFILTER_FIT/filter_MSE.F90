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
! File contents:
! SUBROUTINE filter_MSE
!
! NAME
!     filter_MSE
!
! DESCRIPTION
!
! Given a complex frequency domain dataset and a filter function return the mean square error between the two
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 28/04/2016 CJS 
!
!
SUBROUTINE filter_MSE(Hp,f,n_frequencies,Hfilter,MSE)

USE type_specifications
USE general_module
USE constants
USE filter_module
USE maths

IMPLICIT NONE

! variables passed to subroutine

integer,intent(IN)                 :: n_frequencies           ! number of frequencies at which we have data points
complex(dp),intent(IN)             :: Hp(1:n_frequencies)     ! list of complex function values at each frequency
real(dp),intent(IN)                :: f(1:n_frequencies)      ! list of frequencies at which we have data
type(Sfilter),intent(IN)           :: Hfilter                 ! rational (filter) function
real(dp),intent(OUT)               :: MSE                     ! Mean Square Error to output
  
! local variables

  complex(dp)             :: filter_response
  integer :: i

! START

  MSE=0d0
  
  if (n_frequencies.EQ.0) RETURN

! loop over frequency
  do i=1,n_frequencies

! evaluate the filter function at the current frequency  
    filter_response=evaluate_Sfilter_frequency_response(Hfilter,f(i))

! add the contribution to the square error at the current frequency
    MSE=MSE+(abs(filter_response-Hp(i)))**2
       
  end do ! next row

! divide by the number of samples to return the mean square error  
  MSE=MSE/n_frequencies
                   
  RETURN

END SUBROUTINE filter_MSE
