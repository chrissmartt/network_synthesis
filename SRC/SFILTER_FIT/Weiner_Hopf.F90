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
! SUBROUTINE Weiner_Hopf
!
! NAME
!     Weiner_Hopf
!
! DESCRIPTION
!
! Fit rational function model to input data using the Wiener Hopf method
! as described in Dawson paper [**** REFERENCE REQUIRED ****]
!     
! COMMENTS
!     If fit_type=1 then we impose Hfilter(f=0)=1.0. In this case we fit to a finction Hp-1 to
!     calculate a filter function with a(0)=0, then add 1 to this filter to give the result.
!     If fit_type=0 then we fit directly to the given function values.
!
! HISTORY
!
!     started 28/04/2016 CJS based on FILTER_FIT from the GGI_TLM project 
!                                   www.github.com/ggiemr/GGI_TLM
!     8/5/2017         CJS: Include references to Theory_Manual
!!
!
SUBROUTINE Weiner_Hopf(Hp,f,n_frequencies,Hfilter,aorder,border,fit_type)

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
type(Sfilter),intent(OUT)          :: Hfilter                 ! Output 'best fit' rational (filter) function
integer,intent(IN)                 :: aorder                  ! Numerator model order to calculate
integer,intent(IN)                 :: border                  ! Denominator model order to calculate
integer,intent(IN)                 :: fit_type                ! Type of model fit. fit_type=0 has no restrictions on f=0 value
                                                              !                    fit_type=1 imposes Hfilter(f=0)=1.0
! local variables

  integer :: coeff_dim               ! number of unknown filter coefficients to caluclate
  integer :: freq_dim                ! number of frequencies for which we have function values 
  integer :: first_a                 ! First value of a to be calculated (fit_type=1 imposes a%coeff(0)=1)
  integer :: N_a                     ! Number of a values to be calculated 

! matrices used in the calculations
  complex(dp),allocatable :: A(:,:)
  complex(dp),allocatable :: AH(:,:)
  complex(dp),allocatable :: AHA(:,:)
  complex(dp),allocatable :: AHAI(:,:)
  complex(dp),allocatable :: AHAIAH(:,:)

! vectors used in the calculations  
  complex(dp),allocatable :: VC(:)
  complex(dp),allocatable :: VH(:)

! Laplace variable  
  complex(dp) :: s

! normalisation constant  
  real(dp)    :: wnorm

! loop variables   
  integer :: i,row,col

! error code for EISPACK subroutine calls  
  integer :: ierr

! START

  if (verbose) write(*,*)'CALLED: Weiner_Hopf'

! allocate the filter structure for the specified model numerator and denominator order  
  Hfilter=allocate_Sfilter(aorder,border)
  
! Set the angular frequency normalisation to the maximum angular frequency i.e. the last frequency sample
  wnorm=2d0*pi*f(n_frequencies)
  Hfilter%wnorm=wnorm
  
  if ((aorder.eq.0).AND.(border.eq.0)) then
! this is a special case which can be dealt with simply

    Hfilter%b%coeff(0)=1d0
   
    if (fit_type.EQ.1) then   ! this type has the requirement that Hfilter(f=0)=1.0 so return  Hfilter=1 

      Hfilter%a%coeff(0)=1d0
      
      RETURN
      
    else
! Hfilter is effectively a real constant and can be calculated as the average of the real part of the input function values

      Hfilter%a%coeff(0)=0d0
      
      do row=1,n_frequencies
        Hfilter%a%coeff(0)=Hfilter%a%coeff(0)+dble(Hp(row))
      end do
      Hfilter%a%coeff(0)=Hfilter%a%coeff(0)/n_frequencies
      
      RETURN
      
    end if
   
  end if  ! aorder and border are both equal to zero

! Set the system dimensions and set the known filter coefficient values

  freq_dim=2*n_frequencies  ! twice the total number of frequencies for which we have data as we include negative frequencies to make the coefficients real.
  
  if (fit_type.EQ.0) then  ! No restrictions on Hfilter(f=0), include an extra unknown for a%coeff(0)
    coeff_dim=aorder+1+border 
    first_a=0
    N_a=aorder+1
  else                     ! Hfilter(f=0)=1.0 thus a(0)=1.0 and b(0)=1.0
    coeff_dim=aorder+border     
    first_a=1
    N_a=aorder
  end if
  
! Allocate memory for the matrices and vectors used in the calculations
  ALLOCATE( A(freq_dim,coeff_dim) )

  ALLOCATE( AH(coeff_dim,freq_dim) )
  
  ALLOCATE( AHA(coeff_dim,coeff_dim) )
  ALLOCATE( AHAI(coeff_dim,coeff_dim) )
  ALLOCATE( AHAIAH(coeff_dim,freq_dim) )
  
  ALLOCATE( VH(freq_dim) )
  ALLOCATE( VC(coeff_dim) )
  
  if (verbose) then
    write(*,*)'n_frequencies=',n_frequencies
    write(*,*)'fmin=',f(1)
    write(*,*)'fmax=',f(n_frequencies)
  end if
  
  if (freq_dim.LT.coeff_dim) then
    if(verbose) then
      write(*,*)'Number of frequency samples must be greater than the number of unknowns'
      write(*,*)'n_frequencies=',n_frequencies,' number of unknowns=',coeff_dim
    end if
    run_status='ERROR in Weiner_Hopf: Number of frequency samples must be greater than the number of unknowns'
    CALL write_program_status()
    STOP 1
  end if
                	      
! fill matrix A (equation [A](VC)=(VH)                            
! Include fs at negative frequencies to ensure coefficients are real
!  Theory_Manual_Eqn 5.4, 5.5

  do row=1,n_frequencies   ! loop over frequencies

    s=2d0*pi*j*f(row)/wnorm   ! normalised laplace variable s=jw/wnorm
    
! fill columns related to the unknown a coefficients
    col=0
    do i=first_a,aorder
      col=col+1
      A(row,col)=(s**i)
      A(row+n_frequencies,col)=conjg(A(row,col))
    end do

! fill columns related to the unknown b coefficients
    do i=1,border
      col=col+1
      if (fit_type.eq.1) then            ! fit to Hp(f)-1
        A(row,col)=-(Hp(row)-1d0)*(s**i)
      else
        A(row,col)=-(Hp(row))*(s**i)
      end if
      A(row+n_frequencies,col)=conjg( A(row,col) )
    end do  
       
  end do ! next row

! fill vector VH     Theory_Manual_Eqn 5.4, 5.5

  do row=1,n_frequencies
    if (fit_type.eq.1) then            ! fit to Hp(f)-1
      VH(row)=(Hp(row)-1d0)
    else
      VH(row)=Hp(row)
    endif
    VH(row+n_frequencies)=conjg( VH(row) )
  end do

! Invert using the Moore Penrose method    Theory_Manual_Eqn 5.15

! Calculate the conjugate transpose of A
  do row=1,freq_dim
    do col=1,coeff_dim
      AH(col,row)=conjg( A(row,col) )
    end do
  end do
  
  AHA=matmul(AH,A)
  
  ierr=0   ! set ierr=0 on input to matrix inverse to cause the program to stop if we have a singular matrix
  CALL cinvert_Gauss_Jordan(AHA,coeff_dim,AHAI,coeff_dim,ierr) 

  AHAIAH=matmul(AHAI,AH)                        

! calculate the vector of a and b coefficients  
  VC=matmul(AHAIAH,VH)

! add the a and b coefficients to the filter structure
! the denomiator is the same whether we fit to H-1 or H
  
  Hfilter%b%coeff(0)=1d0             ! b(0) is always 1
  
  do i=1,border 
    row=i+N_a
    Hfilter%b%coeff(i)=VC(row)  
  end do
  
  if (fit_type.eq.1) then            ! we have a fit to Hp(f)-1 so calculate the coefficients of H=1+A(jomega)/B(jomega)
    Hfilter%a%coeff(0)=1d0
    do i=1,aorder
      row=i
      Hfilter%a%coeff(i)=VC(row)+Hfilter%b%coeff(i)    
    end do
  else
    do i=0,aorder
      row=i+1
      Hfilter%a%coeff(i)=VC(row)
    end do
  end if
  
  if (verbose) then
    write(*,*)' ********** FILTER FUNCTION **********'
    CALL write_Sfilter(Hfilter,0)
    write(*,*)' ********** FILTER FUNCTION **********'
  end if  
   
! Deallocate memory

  DEALLOCATE( A )

  DEALLOCATE( AH )
  
  DEALLOCATE( AHA )
  DEALLOCATE( AHAI )
  DEALLOCATE( AHAIAH )
  
  DEALLOCATE( VH )
  DEALLOCATE( VC )

  if (verbose) write(*,*)'FINISHED: Weiner_Hopf'
                   
  RETURN

END SUBROUTINE Weiner_Hopf
