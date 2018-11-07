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
! SUBROUTINE FFT_TIME_TO_FREQ
! SUBROUTINE FFT_FREQ_TO_TIME
! SUBROUTINE FFT
!
! NAME
!    FFT_TIME_TO_FREQ
!
! DESCRIPTION
!     wrapper for Forward Fast Fourier Transform routine
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 24/4/2015 CJS
!
!

SUBROUTINE FFT_TIME_TO_FREQ(n,time,f_time,freq,f_freq)

USE type_specifications
USE constants

IMPLICIT NONE

! variables passed to subroutine

  integer,intent(IN)      :: n         ! number of samples
  real(dp),intent(IN)     :: time(n)   ! input time values
  real(dp),intent(OUT)    :: freq(n)   ! output frequency values
  complex(dp),intent(IN)  :: f_time(n) ! input function of time to transform
  complex(dp),intent(OUT) :: f_freq(n) ! output function of frequency
  
! local variables

  complex(dp)    :: x(n)
  integer            :: n2
  complex(dp),allocatable    :: xe(:)
  complex(dp),allocatable    :: xo(:)
  complex(dp)    :: const
  
  real(dp)    :: dt
  real(dp)    :: fmin,fmax,fstep
  integer    :: i

! START

! No action required for n=1  
  if(n .LE. 1) RETURN

! get the timestep
  dt=time(2)-time(1)  
  
! set the frequency values

  fmin=0d0
  fstep=1d0/(n*dt)
  fmax=fstep*(n-1)
  
! Generate the frequency list
  
  do i=1,n
      
    freq(i)=fmin+(i-1)*fstep
  
  end do
  
! copy input function of time to x
  x(:)=f_time(:)
  
! Calculate FFT
  CALL FFT(x,N)
   
! copy x to the output function of frequency
  f_freq(:)=x(:)
 
  RETURN
 
END SUBROUTINE FFT_TIME_TO_FREQ
!
! NAME
!    INVERSE_FFT_FREQ_TO_TIME
!
! DESCRIPTION
!     Inverse Fast Fourier Transform routine 
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 24/4/2015 CJS
!
!
SUBROUTINE FFT_FREQ_TO_TIME(n,time,f_time,freq,f_freq)

USE type_specifications
USE constants

IMPLICIT NONE

! variables passed to subroutine

  integer,intent(IN)       :: n         ! number of samples
  real(dp),intent(OUT)     :: time(n)   ! output time values
  real(dp),intent(IN)      :: freq(n)   ! input frequency values
  complex(dp),intent(OUT)  :: f_time(n) ! output function of time to transform
  complex(dp),intent(IN)   :: f_freq(n) ! input function of frequency
  
! local variables

  complex(dp)    :: x(n)
  integer            :: n2
  complex(dp),allocatable    :: xe(:)
  complex(dp),allocatable    :: xo(:)
  complex(dp)    :: const

  complex(dp)     :: swap
  
  real(dp)    :: df
  real(dp)    :: tmin,tmax,tstep
  integer    :: i

! START

! get the frequency step
  df=freq(2)-freq(1)  
  
! set the time values

  tmin=0d0
  tstep=1d0/(n*df)
  tmax=tstep*(n-1)
  
! Generate the time list
  
  do i=1,n
      
    time(i)=tmin+(i-1)*tstep
  
  end do

! copy the dataset into a local array
  do i=1,n
    x(i)=f_freq(i)
  end do

! Call the forward FFT routine
  CALL FFT(x,n)
  
  do i=1,n
    x(i)=x(i)/n
  end do
  
! loop over the result dividing by n and reversing the imaginary part
  f_time(1)=x(1)
  do i=2,n
    f_time(i)=x(n-i+2)
  end do
  
  RETURN
 
END SUBROUTINE FFT_FREQ_TO_TIME
!
! SUBROUTINE FFT
!
! NAME
!    FFT
!
! DESCRIPTION
!     Fast Fourier Transform routine
!     This is a simple recursive implementation of the fast fourier transform
!     not the most efficient but very compact.
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 20/11/2013 CJS
!
!
RECURSIVE SUBROUTINE FFT(x,N)

USE type_specifications
USE constants

IMPLICIT NONE

! variables passed to subroutine

  integer,intent(IN)           :: n    ! number of samples
  complex(dp),intent(INOUT)    :: x(n) ! sample values

! local variables

  integer            :: n2
  complex(dp),allocatable    :: xe(:)
  complex(dp),allocatable    :: xo(:)
  complex(dp)    :: const
  
  integer    :: i

! START

! No action required for n=1  
  if(n .LE. 1) RETURN
  
  n2=n/2
 
  ALLOCATE(xe(1:n2))
  ALLOCATE(xo(1:n2))
 
! fill odd and even data
 
  do i=1,n2
    xo(i)=x(2*i-1)
    xe(i)=x(2*i)
  end do

! FFT odd and even sequences
  CALL FFT(xo,n2)
  CALL FFT(xe,n2)
 
! combine odd and even FFTs

  do i=1,n2
 
    const=exp(-2d0*pi*j*(i-1)/n)

    x(i)   = xo(i)+xe(i)*const
    x(i+N2)= xo(i)-xe(i)*const
    
  end do
 
  DEALLOCATE(xo)
  DEALLOCATE(xe)
  
  RETURN
 
END SUBROUTINE FFT
