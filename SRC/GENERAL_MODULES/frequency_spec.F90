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
! MODULE frequency_spec
! SUBROUTINE read_and_set_up_frequency_specification
! SUBROUTINE reset_frequency_specification
! SUBROUTINE set_up_frequency_specification
! SUBROUTINE deallocate_frequency_specification
!
! NAME
!     frequency_spec
!
! AUTHORS
!     Chris Smartt
!
! DESCRIPTION
!     Module to deal with all aspects of frequency range specifications
!     The module defines a struncture for a frequency range
!     subroutine to read a frequency specification and set up the structure
!     subroutine to reset the frequency_specification
!     subroutine to set up the frequency_specification
!     subroutine to deallocate the frequency_specification
!
! COMMENTS
!      
!
! HISTORY
!
!     started 22/08/2016 CJS 
!     Add some checks etc for robustness 15/12/2016 CJS 
!
!
MODULE frequency_spec

USE type_specifications

IMPLICIT NONE

TYPE :: frequency_specification

character(LEN=3):: freq_range_type ! frequency range type, 'log' or 'lin'
real(dp)        :: fmin            ! minimum frequency
real(dp)        :: fmax            ! maximum frequency
integer         :: n_frequencies   ! number of frequencies
real(dp),allocatable ::freq_list(:)! list of frequencies

END TYPE frequency_specification

CONTAINS
 
! NAME
!     read_and_set_up_frequency_specification
!
! AUTHORS
!     Chris Smartt
!
! DESCRIPTION
!     Set up frequency range specifications based on information read from file
!
! COMMENTS
!      
!
! HISTORY
!
!     started 22/08/2016 CJS 
!
SUBROUTINE read_and_set_up_frequency_specification(freq_data,unit)

USE general_module

IMPLICIT NONE

! variables passed to subroutine

TYPE(frequency_specification),intent(INOUT) :: freq_data
integer,intent(IN)                          :: unit
  
! local variables  
  
real(dp)        :: fstep               ! frequency step

real(dp)        :: log_fmin            ! log of minimum frequency
real(dp)        :: log_fmax            ! log of maximum frequency
real(dp)        :: log_fstep           ! log of frequency step
real(dp)        :: log_f               ! log of frequency

integer         :: ierr

! START

  read(unit,'(A3)',IOSTAT=ierr)freq_data%freq_range_type  
  if (ierr.NE.0) then
    run_status='ERROR reading .spice_model_spec file:  '// &
               'The frequency range type needs to be specified (log or lin)'
    CALL write_program_status()
    STOP 1
  end if
  
  read(unit,*,IOSTAT=ierr)freq_data%fmin,   &
                          freq_data%fmax,   &
                          freq_data%n_frequencies
  if (ierr.NE.0) then
    run_status='ERROR: reading .spice_model_spec file: '// &
               'The frequency range  needs to be specified (fmin fmax n_frequenceis)'
    CALL write_program_status()
    STOP 1
  end if

! Some checks that the specified frequency range is sensible

  if (freq_data%fmin.LT.0d0) then
    run_status='ERROR: read_and_set_up_frequency_specification.   fmin < 0'
    CALL write_program_status()
    STOP 1
  end if

  if (freq_data%fmax.LT.freq_data%fmin) then
    run_status='ERROR: read_and_set_up_frequency_specification.   fmax < fmin'
    CALL write_program_status()
    STOP 1
  end if

  if (freq_data%n_frequencies.LT.0d01) then
    run_status='ERROR: read_and_set_up_frequency_specification.  number of frequencies < 1'
    CALL write_program_status()
    STOP 1
  end if
  
  CALL set_up_frequency_specification(freq_data)

  RETURN

END SUBROUTINE read_and_set_up_frequency_specification
!
! NAME
!     reset_frequency_specification
!
! AUTHORS
!     Chris Smartt
!
! DESCRIPTION
!     reset frequency range specifications
!
! COMMENTS
!      
!
! HISTORY
!
!     started 15/12/2016 CJS 
!
SUBROUTINE reset_frequency_specification(freq_data)

IMPLICIT NONE

! variables passed to subroutine

TYPE(frequency_specification),intent(INOUT) :: freq_data
  
! local variables  

! START

  freq_data%freq_range_type='lin'
  freq_data%fmin=0d0
  freq_data%fmax=0d0
  freq_data%n_frequencies=1
  
  RETURN

END SUBROUTINE reset_frequency_specification
!
! NAME
!     set_up_frequency_specification
!
! AUTHORS
!     Chris Smartt
!
! DESCRIPTION
!     Set up frequency range specifications i.e. given the frequency range specification, 
!     allocate and set the frequency list array.
!
! COMMENTS
!      
!
! HISTORY
!
!     started 22/08/2016 CJS 
!     16/09/2016 CJS process separated from read_and_set_up_frequency_specification
!
SUBROUTINE set_up_frequency_specification(freq_data)

USE general_module

IMPLICIT NONE

! variables passed to subroutine

TYPE(frequency_specification),intent(INOUT) :: freq_data
  
! local variables  
  
real(dp)        :: fstep               ! frequency step

real(dp)        :: log_fmin            ! log of minimum frequency
real(dp)        :: log_fmax            ! log of maximum frequency
real(dp)        :: log_fstep           ! log of frequency step
real(dp)        :: log_f               ! log of frequency

integer         :: floop               ! frequency loop variable

integer         :: ierr

! START

  if (freq_data%freq_range_type.EQ.'log') then

    log_fmin=log10(freq_data%fmin)
    log_fmax=log10(freq_data%fmax)
    log_fstep=0d0                  ! this is the value used if freq_data%n_frequencies=1
    if (freq_data%n_frequencies.ne.1) then
      log_fstep=(log_fmax-log_fmin)/dble(freq_data%n_frequencies-1)
    end if
    
  else if (freq_data%freq_range_type.EQ.'lin') then

    fstep=0d0                      ! this is the value used if freq_data%n_frequencies=1
    if (freq_data%n_frequencies.ne.1) then
      fstep=(freq_data%fmax-freq_data%fmin)/dble(freq_data%n_frequencies-1)
    end if
        
  else 

    run_status='ERROR in set_up_frequency_specification: frequency range type should be log or lin'
    CALL write_program_status()
    STOP 1
  
  end if
  
  if ( ALLOCATED(freq_data%freq_list) ) CALL deallocate_frequency_specification(freq_data)

  ALLOCATE(freq_data%freq_list(1:freq_data%n_frequencies))
  
  do floop=1,freq_data%n_frequencies
    
    if (freq_data%freq_range_type.EQ.'log') then

      log_f=log_fmin+(floop-1)*log_fstep
      freq_data%freq_list(floop)=10d0**(log_f)
    
    else if (freq_data%freq_range_type.EQ.'lin') then

      freq_data%freq_list(floop)=freq_data%fmin+(floop-1)*fstep

    end if
    
  end do ! next frequency

END SUBROUTINE set_up_frequency_specification
 
! NAME
!     deallocate_frequency_specification
!
! AUTHORS
!     Chris Smartt
!
! DESCRIPTION
!     deallocate the frequency range specification
!
! COMMENTS
!      
!
! HISTORY
!
!     started 22/08/2016 CJS 
!

SUBROUTINE deallocate_frequency_specification(freq_data)

IMPLICIT NONE

! variables passed to subroutine

TYPE(frequency_specification),intent(INOUT) :: freq_data

! START

  if (ALLOCATED(freq_data%freq_list)) DEALLOCATE(freq_data%freq_list)

END SUBROUTINE deallocate_frequency_specification

END MODULE frequency_spec
