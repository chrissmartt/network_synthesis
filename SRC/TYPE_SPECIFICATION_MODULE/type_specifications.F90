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
! MODULE type_specifications
!
! NAME
!     type_specifications
!
! AUTHORS
!     Chris Smartt
!
! DESCRIPTION
!     This module defines the types used everywhere in this software for the following quantities:
!     single and long integer types
!     single and double precision real and complex types
!     string lengths for filenames, lines to read from files, spice component names
!     matrix type
!
!     The specifications are designed to be portable between compiilers and conform to the Fortran 2008 standard
!
!
! COMMENTS
!     
!
! HISTORY
!
!     
MODULE type_specifications

  IMPLICIT NONE

! single precision real/complex
  integer,parameter        :: sp=selected_real_kind(6,30) 
  
! double precision real/complex
  integer,parameter        :: dp=selected_real_kind(15,300) 

! normal integer    
  integer,parameter        :: i4=selected_int_kind(5)

! long integer    
  integer,parameter        :: i8=selected_int_kind(15)

! string lengths    
  integer,parameter        :: filename_length=256
  integer,parameter        :: line_length=256
  integer,parameter        :: spice_name_length=30
  integer,parameter        :: error_message_length=200

! matrix type
  TYPE:: matrix

    integer                :: dim
    real(dp),allocatable   :: mat(:,:)

  END TYPE matrix

END MODULE type_specifications
