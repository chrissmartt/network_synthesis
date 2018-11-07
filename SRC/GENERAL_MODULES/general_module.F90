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
! File Contents:
! MODULE general_module
!
! and within include files the following:
!
!
! program_status.F90
! SUBROUTINE write_program_status
!
! read_version.F90
! SUBROUTINE read_version
!
! read_write_subroutines.F90
! SUBROUTINE write_license
! SUBROUTINE path_format
! SUBROUTINE check_and_make_path
! SUBROUTINE write_long_node_list

! NAME
!     general_module
!
! AUTHORS
!     Chris Smartt
!
! DESCRIPTION
!     This module contains varaibles and constants which are used generally within the SACAMOS software
!
!     Also general subroutines are contained within the module in the included files
!
!     add_integer_to_string.F90
!     convert_case.F90
!     read_version.F90
!     program_status.F90
!     read_write_subroutines.F90
!     build_component_name_strings.F90
!     plot_geometry.F90

! COMMENTS
!     
!
! HISTORY
!
MODULE general_module

USE type_specifications

IMPLICIT NONE

SAVE

! version information 
  
character(LEN=line_length)    :: NETWORK_SYNTHESIS_version

character(LEN=line_length)    :: NETWORK_SYNTHESIS_date
  
character(LEN=line_length)    :: NETWORK_SYNTHESIS_compilation_date

! program run status information 

character(LEN=line_length)       :: program_name
character(LEN=line_length)       :: run_status
character(LEN=10),parameter      :: program_status_filename='run_status'

logical :: verbose=.FALSE.

! There are some differences between windows and unix relating to
! file separators and making directories. The different formats for
! the two ooperating systems supported are included below 
! NOTE: Uses mkdir -p in unix which can build the whole path as opposed to mkdir which can only build the bottom level directory
! mkdir -p works in ubuntu 14.04 - not sure that this is a portable solution though...

#ifdef WINDOWS
  character        :: file_separator='\'
  character(LEN=6) :: mkdir_command='mkdir '
#else
  character        :: file_separator='/'
  character(LEN=9) :: mkdir_command='mkdir -p '
#endif


integer,parameter     :: program_status_file_unit=80

integer,parameter     :: local_file_unit=90

integer,parameter     :: temp_file_unit=99

CONTAINS

include 'read_version.F90'

include 'program_status.F90'

include 'read_write_subroutines.F90'

END MODULE general_module
