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
! MODULE network_synthesis
! PROGRAM network_synthesis
!

MODULE network_synthesis_global

USE type_specifications
USE general_module
USE constants
USE frequency_spec
USE filter_module
USE Sfilter_fit_module

IMPLICIT NONE

CONTAINS

INCLUDE 'include_local_filter_subroutines.F90'

INCLUDE 'include_test_positive_real.F90'

INCLUDE 'include_write_spice_model.F90'

INCLUDE 'include_pole_residue_test_functions.F90'

INCLUDE 'include_minimum_resistance_function.F90'

INCLUDE 'include_RLC_test.F90'

INCLUDE 'include_LC_test.F90'

INCLUDE 'include_RC_test.F90'

INCLUDE 'include_RL_test.F90'

INCLUDE 'include_C_test.F90'

INCLUDE 'include_L_test.F90'

INCLUDE 'include_R_test.F90'

INCLUDE 'include_R2_test.F90'

INCLUDE 'include_BRUNE_test.F90'

END MODULE network_synthesis_global

!
! NAME
!     network_synthesis
!
! AUTHORS
!     Chris Smartt
!
! DESCRIPTION
!     Implementation of rational transfer functions as
!     passive equivalent circuits
!
!     PROCESS:
!
!  1. Test that the rational function does constitute a physical impedance function
!  2. Generate a ladder network which synthesises this impedance function where the
!     elements of the ladder network consist of RLC elements
!  3. Evaluate the frequency response of the original function
!  4. Evaluate the frequency response of the continued fraction form
!  5. Write a Spice file for the network
!
!
! COMMENTS
!     
!
! HISTORY
!
!     started 22/08/2017 CJS 
!     14/9/2017 CJS Generalise the circuit elements to RLC combinations 
!                   using partial fraction expansions
!
!
PROGRAM network_synthesis

USE type_specifications
USE general_module
USE constants
USE frequency_spec
USE filter_module
USE Sfilter_fit_module
USE network_synthesis_global

IMPLICIT NONE

! local variables

! command line argument value and length
character(len=filename_length)    :: argument1
integer                           :: argument1_length
character(len=filename_length)    :: argument2
integer                           :: argument2_length

character(len=filename_length)    :: filename   ! filename for the transfer function

type(Sfilter)    :: H                     ! rational transfer function model 
type(Sfilter)    :: HR                    ! rational transfer function model 
type(Sfilter_PR) :: H_PR                  ! pole-residue transfer function model 
type(Sfilter)    :: T1                    ! temporary filter function
type(Sfilter)    :: Rdc_filter            ! temporary filter function

integer        :: aorder,border,max_order,CF_dim
integer        :: n_branches

integer        :: test

integer        :: i,loop,type_loop                          ! loop variables

type(Polynomial)  :: num
type(Polynomial)  :: den

real(dp) :: norm

real(dp) :: value

real(dp),allocatable :: CFterm(:,:)
integer,allocatable  :: CFtype(:)

real(dp)    :: wmin,wmax,wstep
integer     :: nw
complex(dp) :: s,sn,num_fs,den_fs,CF_term,last_CF_term
complex(dp) :: H_rational,H_CF
real(dp)    :: R,L,C

real(dp)    :: R_min     ! minimum value of the resistance of the input function
real(dp)    :: w_R_min     ! angular frequency for minimum value of the resistance of the input function
real(dp)    :: R_add     ! additional resistance required to make the function positive real

logical :: stable
logical :: found
logical :: remainder_OK
logical :: remainder_zero
logical :: multiple_poles

integer           :: type,last_type
real(dp)          :: term

integer :: on,od
logical :: pole_at_zero
logical :: zero_at_zero

logical :: renormalise_input_filter

real(dp):: wnorm_save

character(LEN=2) :: ch2
character(LEN=256) :: line
real(dp) :: re_p,im_p,re_r,im_r

! START

!  renormalise_input_filter=.TRUE.
  renormalise_input_filter=.FALSE.
  
! Open the input file containing the complex frequency domain data

  program_name="network_synthesis"
  run_status='Started'
  CALL write_program_status()
  
!  CALL read_version()
    
!  CALL write_license()

INCLUDE 'include_read_and_plot_function.F90'

  if (renormalise_input_filter) then
  
    write(*,*)'Renormalise the filter:'
    CALL write_Sfilter(H,0)
    
    T1=renormalise_Sfilter(H)
    CALL deallocate_Sfilter(H)
    H=T1
  
    write(*,*)'Calculate minimum resistance of function:'
    CALL write_Sfilter(H,0)
       
  end if
  
  CALL calculate_min_resistance_value(H,R_min,w_R_min)
  
  write(*,*)'Minimum Resistance value is',R_min
  write(*,*)'at w=',w_R_min,' f=',w_R_min/(6.283185)
  
  
  if (R_min.LT.0d0) then
  
    write(*,*)'Adding 1.5*abs(Minimum Resistance value) to H'
  
    R_add=1.5d0*abs(R_min)
    Rdc_filter=allocate_Sfilter(0,0)
    Rdc_filter%a%coeff(0)=R_add
    Rdc_filter%b%coeff(0)=1d0
  
    CALL deallocate_Sfilter(T1)
    T1=H+Rdc_filter
    CALL deallocate_Sfilter(H)
    H=T1
    CALL deallocate_Sfilter(T1)
    CALL deallocate_Sfilter(Rdc_filter)

  else
  
    write(*,*)'No resistance added to H'
    R_add=0d0

  end if
  
  write(*,*)'Revised H:'
  CALL write_Sfilter(H,0)

! now we have ensured that Re(H)>=0 at all frequencies, do the checks again.
! Check the transfer funcion for stability and for whether it is positive real
  CALL check_transfer_function(H,stable) 
  
  if (stable) then  
    write(*,*)'INPUT FUNCTION IS A STABLE, PHYSICAL IMPEDANCE'
  else  
    write(*,*)'INPUT FUNCTION IS NOT STABLE'
    run_status='INPUT FUNCTION IS NOT STABLE, even after adding a stabilising d.c. resistance'
    CALL write_program_status()
    STOP   
  end if
    
! Max_order is the maximum number of components, including  
! resistive and reactive components

  max_order=2*max(H%a%order,H%b%order) +1
  
  if (verbose) then
    write(*,*)'Maximum order is estimated to be ',max_order
  end if
  
! allocate the continued fraction data  
  CF_dim=max_order
  allocate( CFterm(1:max_order,1:5) ) ! note 5 terms
  CFterm(:,:)=0d0
  allocate( CFtype(1:max_order) )
  CFtype(:)=0
  
  n_branches=0
  type=type_impedance  ! H(s) is an impedance function to start with
  
  do loop=1,max_order
  
    remainder_OK=.FALSE.
    remainder_zero=.FALSE.

! Loop for trying both impedance and admittance functions   
    do type_loop=1,2
    
      if (verbose) then
        write(*,*)'Stage ',loop,' of ',max_order       
        if (type.EQ.type_impedance) then
          write(*,*)'TRYING TO CALCULATE AN IMPEDANCE'        
        else
          write(*,*)'TRYING TO CALCULATE AN ADMITTANCE'
         end if
      end if

! check the function for multiple poles. If there are no multiple poles then we can
! go on and look for viable branches

      CALL check_for_multiple_roots(H%b,multiple_poles,.FALSE.)
      
      if (.NOT.multiple_poles) then

! Calculate the partial fraction expansion of the function H(s)
        H_PR=Convert_filter_S_to_S_PR(H)
        if (verbose) CALL write_S_PR_filter(H_PR)
      
        do test=1,8
           
! Test number 1: looking for RLC branch
          select case (test)
        
          case(1)
! Test number 1: looking for RLC branch
            CALL RLC_test(H_PR,type,CFtype(loop),CFterm(loop,1),CFterm(loop,2),CFterm(loop,3),found, &
                          HR,remainder_OK,remainder_zero)
                              
          case(2)
! Test number 2: looking for LC branch
            CALL LC_test(H_PR,type,CFtype(loop),CFterm(loop,1),CFterm(loop,2),CFterm(loop,3),found, &
                          HR,remainder_OK,remainder_zero)
       
          case(3)
! Test number 3: looking for RC branch
            CALL RC_test(H_PR,type,CFtype(loop),CFterm(loop,1),CFterm(loop,2),CFterm(loop,3),found, &
                          HR,remainder_OK,remainder_zero)
      
          case(4)
! Test number 4: looking for RL branch
            CALL RL_test(H_PR,type,CFtype(loop),CFterm(loop,1),CFterm(loop,2),CFterm(loop,3),found, &
                         HR,remainder_OK,remainder_zero)
     
          case(5)
! Test number 5: looking for C branch
            CALL C_test(H_PR,type,CFtype(loop),CFterm(loop,1),CFterm(loop,2),CFterm(loop,3),found, &
                         HR,remainder_OK,remainder_zero)
     
          case(6)
! Test number 6: looking for L branch
            CALL L_test(H_PR,type,CFtype(loop),CFterm(loop,1),CFterm(loop,2),CFterm(loop,3),found, &
                         HR,remainder_OK,remainder_zero)
     
          case(7)
! Test number 7: looking for R branch
            CALL R_test(H_PR,type,CFtype(loop),CFterm(loop,1),CFterm(loop,2),CFterm(loop,3),found, &
                         HR,remainder_OK,remainder_zero)
    
          case(8)
! Test number 8: looking for R branch, second form
            CALL R2_test(H_PR,type,CFtype(loop),CFterm(loop,1),CFterm(loop,2),CFterm(loop,3),found, &
                         HR,remainder_OK,remainder_zero)

          end select
                    
          if (found.AND.remainder_OK) then ! Adopt this as a viable branch
            n_branches=n_branches+1
            if (remainder_zero) then
              GOTO 2000
            else
              GOTO 1000
            end if
          end if

        end do ! next test
      
      end if

! If we have been unsuccessful finding a circuit element
! change from impedance->admittance of vice vers    

      if (type.EQ.type_impedance) then
        type=type_admittance
      else
        type=type_impedance
      end if
      
! Calculate the reciprocal filter function

      num=H%a
      den=H%b
      CALL deallocate_poly(H%a)
      CALL deallocate_poly(H%b)
      H%a=den
      H%b=num
      CALL deallocate_poly(num)
      CALL deallocate_poly(den)
          
    end do ! next type_loop
    
! if we get here then we have not found a viable way to proceed with
! building the circuit
    write(*,*)'CANNOT FIND VIABLE COMPONENT'
    write(*,*)'TRYING BRUNE CYCLE'

! convert to an impedance function if required
    if (type.EQ.type_admittance) then
      type=type_impedance     
! Calculate the reciprocal filter function
      num=H%a
      den=H%b
      CALL deallocate_poly(H%a)
      CALL deallocate_poly(H%b)
      H%a=den
      H%b=num
      CALL deallocate_poly(num)
      CALL deallocate_poly(den)
    end if
    
    CALL BRUNE_test(H,type,CFtype(loop),CFterm(loop,1),CFterm(loop,2),CFterm(loop,3), &
                                           CFterm(loop,4),CFterm(loop,5),found,          &
                                           HR,remainder_OK,remainder_zero)
                                           
    if (found.AND.remainder_OK) then ! Adopt this as a viable branch
      n_branches=n_branches+1
      if (remainder_zero) then
        GOTO 2000
      else
        GOTO 1000
      end if
    end if
                       
    run_status='ERROR: CANNOT FIND VIABLE COMPONENT'
    CALL write_program_status()
    STOP
    
! jump here when we have found the next circuit element
1000 CONTINUE
        
    if (verbose) write(*,*)'Prepare for the next stage'
    
    CALL deallocate_Sfilter(H)
    H=HR
               
  end do  ! next circuit element
  
! If we end up here then there is a problem because the remainder is not zero
! and we are supposed to have worked out all the circuit elements by now.


! jump here when the continued fraction truncates with a zero remainder
2000 continue

!  CALL write_CF_local(CFterm,CFtype,CF_dim,n_branches)

  write(*,*)''
  write(*,*)'Re-scale normailsed filter components'
  write(*,*)''

INCLUDE 'include_write_frequency_response.F90'
  
  do loop=1,max_order 
    CFterm(loop,2)=CFterm(loop,2)/wnorm_save      ! L term
    CFterm(loop,3)=CFterm(loop,3)/wnorm_save      ! C term
    CFterm(loop,4)=CFterm(loop,4)/wnorm_save      ! L2 term
  end do

! Write a spice circuit model for the ladder network derived from the
! continued fraction expansion

  CALL write_ladder_network(CFterm,CFtype,CF_dim,n_branches,R_add,nw,wmin,wmax,wnorm_save)

! deallocate memory and finish up
  
  CALL deallocate_poly(num)
  CALL deallocate_poly(den)

  CALL deallocate_Sfilter(H)
  CALL deallocate_Sfilter(HR)
  CALL deallocate_Sfilter(T1)
  CALL deallocate_Sfilter(Rdc_filter)
  deallocate( CFterm )
  
  run_status='Finished_Correctly'
  CALL write_program_status()
  
  STOP
  
9000 write(*,*)'Error reading transfer function file'
  STOP
  
END PROGRAM network_synthesis
