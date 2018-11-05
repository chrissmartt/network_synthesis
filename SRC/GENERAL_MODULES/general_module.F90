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
! MODULE general_module
!
! and within include files the following:
!
! add_integer_to_string.F90
! SUBROUTINE add_integer_to_string
!
! convert_case.F90
! SUBROUTINE convert_to_lower_case
! SUBROUTINE convert_to_upper_case
!
! build_component_name_strings.F90
! SUBROUTINE build_name_with_conductor_and_end
! SUBROUTINE build_name_with_conductor_domainConductor_and_end     
! SUBROUTINE build_name_with_domain_conductor_and_end     
! SUBROUTINE build_name_with_domain_mode_and_end       
! SUBROUTINE build_name_with_domain_and_mode
! SUBROUTINE build_name_with_domain_conductor_mode_and_end       
!
! generate_shapes.F90
! SUBROUTINE generate_circle_points
! SUBROUTINE generate_rectangle_points
! SUBROUTINE generate_Dshape_points
! SUBROUTINE generate_arc_points
!
! plot_geometry.F90
! SUBROUTINE write_circle
! SUBROUTINE write_rectangle
! SUBROUTINE write_Dshape
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
  
character(LEN=line_length)    :: SPICE_CABLE_MODEL_BUILDER_version

character(LEN=line_length)    :: SPICE_CABLE_MODEL_BUILDER_date
  
character(LEN=line_length)    :: SPICE_CABLE_MODEL_BUILDER_compilation_date

! program run status information 

character(LEN=line_length)       :: program_name
character(LEN=line_length)       :: run_status
character(LEN=10),parameter      :: program_status_filename='run_status'

! Library of cable models, directory information 

character(len=line_length)    :: MOD_cable_lib_dir 
character(len=line_length)    :: MOD_bundle_lib_dir 
character(len=line_length)    :: MOD_spice_bundle_lib_dir 
character(len=line_length)    :: spice_symbol_dir 

! Target spice version information
  
integer,parameter :: ngspice=1  
integer,parameter :: LTspice=2  
integer,parameter :: Pspice =3  

! Shape information 

integer,parameter :: circle =1 
integer,parameter :: rectangle=2
integer,parameter :: Dshape=3

integer :: spice_version 

! Maximum line length for Spice input files

integer :: max_spice_line_length=100

! flags. The default values are set here but they may be changed in the *_spec files.

logical :: use_s_xfer=.TRUE.

logical :: use_Laplace=.FALSE.

logical :: use_Xie=.FALSE.

logical :: plot_potential=.FALSE.

logical :: plot_mesh=.TRUE.

logical :: plot_real=.FALSE.

logical :: high_freq_Zt_model=.FALSE.

logical :: run_validation_test=.TRUE.

logical :: run_validation_test_Vbased=.FALSE.

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

! filename and unit information  
  
integer,parameter     :: cable_spec_file_unit=10
character(LEN=11),parameter:: cable_spec_file_extn='.cable_spec'

integer,parameter     :: cable_file_unit=11
character(LEN=6),parameter :: cable_file_extn='.cable'

integer,parameter     :: bundle_spec_file_unit=12
character(LEN=12),parameter:: bundle_spec_file_extn='.bundle_spec'

integer,parameter     :: bundle_file_unit=13
character(LEN=7),parameter :: bundle_file_extn='.bundle'

integer,parameter     :: spice_model_spec_file_unit=14
character(LEN=17),parameter:: spice_model_spec_file_extn='.spice_model_spec'

integer,parameter     :: spice_model_file_unit=15
integer,parameter     :: test_circuit_file_unit=16

character(LEN=12) :: test_circuit_file_extn
character(LEN=12) :: spice_model_file_extn

character(LEN=12),parameter :: ngspice_test_circuit_file_extn='_NGspice.cir'
character(LEN=12),parameter :: ngspice_spice_model_file_extn='_NGspice.lib'

character(LEN=12),parameter :: LTspice_test_circuit_file_extn='_LTspice.cir'
character(LEN=12),parameter :: LTspice_spice_model_file_extn='_LTspice.lib'

character(LEN=11),parameter :: Pspice_test_circuit_file_extn='_Pspice.cir'
character(LEN=11),parameter :: Pspice_spice_model_file_extn='_Pspice.lib'

integer,parameter           :: analytic_soln_file_unit=17
character(LEN=21),parameter :: analytic_soln_filename='analytic_solution.dat'
 
integer,parameter           :: gmsh_geometry_file_unit=20
character(LEN=4),parameter  :: gmsh_geometry_file_extn='.geo'

integer,parameter           :: mesh_file_unit=21
character(LEN=4),parameter  :: mesh_file_extn='.msh'

integer,parameter           :: boundary_file_unit=22
character(LEN=14),parameter :: boundary_file_name='laplace_in.bnd'

integer,parameter           :: symbol_file_unit=30
character(LEN=4)            :: symbol_file_extn
character(LEN=4),parameter  :: symbol_file_extn_NGspice='.sym'
character(LEN=4),parameter  :: symbol_file_extn_LTspice='.asy'
  
integer,parameter           :: braid_spec_file_unit=40
character(LEN=11),parameter :: braid_spec_file_extn='.braid_spec'
  
integer,parameter           :: shield_model_file_unit=41
character(LEN=13),parameter :: shield_model_file_extn='.shield_model'

integer,parameter           :: dielectric_geometry_file_unit=51
character(LEN=15),parameter :: dielectric_geometry_filename='dielectrics.dat'

integer,parameter           :: conductor_geometry_file_unit=52
character(LEN=14),parameter :: conductor_geometry_filename='conductors.dat'

integer,parameter           :: frame_geometry_file_unit=53
character(LEN=9),parameter  :: frame_geometry_filename='frame.dat'

integer,parameter     :: program_status_file_unit=80

integer,parameter     :: local_file_unit=90

integer,parameter     :: temp_file_unit=99

CONTAINS

include 'add_integer_to_string.F90'

include 'convert_case.F90'

include 'read_version.F90'

include 'program_status.F90'

include 'read_write_subroutines.F90'

include 'build_component_name_strings.F90'

include 'generate_shapes.F90'

include 'plot_geometry.F90'

END MODULE general_module
