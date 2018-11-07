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
! File Contents:
! MODULE constants
!
! NAME
!     constants
!
! AUTHORS
!     Chris Smartt
!
! DESCRIPTION
!     This module defines the constants and parameters used in the SACAMOS software
!     The parameters are specified here and cannot be changed when running the code
!     There are also constants which are set here which may be given different values
!     by the user when runnning the software.
!
! COMMENTS
!     
!
! HISTORY
!
MODULE constants

USE type_specifications

IMPLICIT NONE

! These parameters are fixed and cannot be changed by the user when running the software

  real(dp),parameter ::  small=1d-14

  real(dp),parameter ::  large=1d14

  real(dp),parameter ::  pi=3.14159265358979323846264338d0
  
  real(dp),parameter ::  eps0=8.8541878176D-12
  
  real(dp),parameter ::  mu0=pi*4D-7
  
  real(dp),parameter ::  c0=2.99792458D8
  
  real(dp),parameter ::  Z0=376.73031346177d0
  
  real(dp),parameter ::  Y0=1D0/Z0
  
  complex(dp),parameter :: j=(0d0,1d0)

  real(dp),parameter ::  zero_test_small=1d-6

  real(dp),parameter ::  zero_test_tiny=1d-12
  
  real(dp),parameter ::  zero_test_R=1d-6   ! 1d-6

  real(dp),parameter ::  zero_test_L=1d-12   ! 1d-12

  real(dp),parameter ::  zero_test_C=1d-14   ! 1d-14

  real(dp),parameter ::  filter_fit_err=1d-8

! These could be constants which can be changed from the input files. 
  real(dp),parameter ::  Laplace_ground_plane_ratio=0.8D0     ! the size of the ground plane as a proportion of the radius of the outer boundary used in Laplace
  
  real(dp),parameter ::  Laplace_ground_plane_thickness_ratio=0.1d0 ! the thickness of the ground plane as a proportion of its width for Laplce
  
! The following constants may be changed from the input files. 

! minimum transmission line delay parameters. These may be changed in the spice_model_spec file  
  real(dp) ::  ZT_min_delay=1D-12   ! minimum delay in seconds allowed in the transfer impedance model Tvictim-Tsource delay line

  real(dp) ::  Einc_min_delay=1D-12 ! minimum delay in seconds allowed in the incident field excitation Tvictim-Tz delay line

  real(dp) ::  Tz_min_delay=1D-12 ! minimum delay in seconds allowed in the incident field excitation Tz delay line

! mesh generation parameters. These may be changed in the .cable_spec and bundle_spec files.
  real(dp) ::  Laplace_boundary_constant=3d0 ! the distance to the outer boundary in the Laplace solver is calculated as
                                             ! Laplace_boundary_constant*(radius of circle surrounding all conductors)
                                           
  real(dp) ::  Laplace_surface_mesh_constant=3d0 ! the mesh edge length on boundaries is calculated as radius/Laplace_surface_mesh_constant 
                                                 ! if Laplace_surface_mesh_constant=5 we have just over 30 elements on the circumference
                                               
  real(dp) ::  Twisted_pair_equivalent_radius=1.5d0 ! The twisted pair commmon mode interaction is calculated by treating the
                                                    ! common mode as being carried on an 'equivalent cylindrical conductor'
                                                    ! whose radius is twisted_pair_equivalent_radius*wire radius

! Minimum resistance value. This may be changed in the spice_model_spec file                                  
  real(dp) ::  Rsmall=1D-8   ! A small resistance value to be used in place of 0 as this is not allowed in LTspice

! test value for poles not on the imaginary (s=jw) axis
  real(dp),parameter ::  pole_small=1d-8
                      
  integer  :: first_subcircuit_node_number=1                                     

! note that the series elements are +ve and the shunt elements negative
integer,parameter :: series_RLC=  1            
integer,parameter :: series_LC=   2     
integer,parameter :: series_RC=   3     
integer,parameter :: series_RL=   4    
integer,parameter :: series_C=    5     
integer,parameter :: series_L=    6     
integer,parameter :: series_R=    7     
integer,parameter :: series_BRUNE=8   
integer,parameter :: series_RLRCR=9   

integer,parameter :: shunt_RLC=  -1            
integer,parameter :: shunt_LC=   -2     
integer,parameter :: shunt_RC=   -3     
integer,parameter :: shunt_RL=   -4    
integer,parameter :: shunt_C=    -5     
integer,parameter :: shunt_L=    -6     
integer,parameter :: shunt_R=    -7     
integer,parameter :: shunt_RLRCR=-9  

  
integer,parameter :: type_impedance  = 1
integer,parameter :: type_admittance =-1
  
  
END MODULE constants
