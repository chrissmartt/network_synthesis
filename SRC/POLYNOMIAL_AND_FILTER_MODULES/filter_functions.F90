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
!
! Description
!     functions related to rational function filter function representations
!     to be included in filter_module.F90
!
! Comments:
!      
!
! History
!
!     started 22/01/09 CJS
!
! FILE CONTAINS
!FUNCTION allocate_Sfilter(aorder,border) RESULT(res)
!FUNCTION jwA_filter(A) RESULT(res)
!FUNCTION evaluate_Sfilter_frequency_response(s1,f) RESULT(res)
!FUNCTION evaluate_Sfilter_high_frequency_limit(s1) RESULT(res)
!FUNCTION Convert_filter_S_to_S_PR(S)
!
! NAME
!     allocate_Sfilter
!
! DESCRIPTION
!     set up a Sfilter type with a given order
!
! HISTORY
!
!     started 22/01/09 CJS
!

FUNCTION allocate_Sfilter(aorder,border) RESULT(res)

USE type_specifications

IMPLICIT NONE

! argument types
  integer,intent(IN)        :: aorder,border
  
! Result type
  type(Sfilter)             :: res
  
! local variables
  type(Sfilter) :: ans

! function definition

  ans%wnorm=1d0
  ans%a%order=aorder
  allocate (ans%a%coeff(0:ans%a%order))
  ans%b%order=border
  allocate (ans%b%coeff(0:ans%b%order))
  ans%a%coeff(0:ans%a%order)=0d0
  ans%b%coeff(0:ans%b%order)=0d0

  res=ans
  
  CALL deallocate_Sfilter(ans)

END FUNCTION allocate_Sfilter
!
!
! NAME
!     jwA_filter
!
! DESCRIPTION
!     set up a Sfilter type of the form jwA where A is a real argument to the function
!
! HISTORY
!
!     started 21/04/2016 CJS
!

FUNCTION jwA_filter(A) RESULT(res)

USE type_specifications

IMPLICIT NONE

! argument types
  real(dp),intent(IN)       :: A
  
! Result type
  type(Sfilter)             :: res

! function definition

  res%wnorm=1d0
  res%a%order=1
  allocate (res%a%coeff(0:res%a%order))
  res%b%order=0
  allocate (res%b%coeff(0:res%b%order))
  res%a%coeff(0)=0d0
  res%a%coeff(1)=A
  res%b%coeff(0)=1d0

END FUNCTION jwA_filter
!
! NAME
!      evaluate_Sfilter_frequency_response
!     
! DESCRIPTION
!       Evaluate the frequency resaponse of a Laplace domain filter at
!     a specified frequency. Filter format is Rational Function
!
! SEE ALSO
!
!
! HISTORY
!
!     started 23/01/09 CJS
!
  FUNCTION evaluate_Sfilter_frequency_response(s1,f) RESULT(res)

USE type_specifications
USE constants
  
  IMPLICIT NONE
  
  type(Sfilter)  , intent(IN)  :: s1
  real(dp)       , intent(IN)  :: f
  
  complex(dp)                  :: res

! local variables  
  real(dp)   :: w
  complex(dp):: jw,num,den
  integer    :: n

!START

  w=2d0*pi*f/s1%wnorm
  jw=j*w
  
  num=evaluate_polynomial(s1%a,jw)
  den=evaluate_polynomial(s1%b,jw)
  
  res=num/den
  
  RETURN
  END FUNCTION evaluate_Sfilter_frequency_response
!
! NAME
!      evaluate_Sfilter_high_frequency_limit
!     
! DESCRIPTION
!       Evaluate the frequency resaponse of a Laplace domain filter at
!     high frequency i.e. as f-> infinity. Filter format is Rational Function
!
! SEE ALSO
!
!
! HISTORY
!
!     started 21/04/2016 CJS
!
  FUNCTION evaluate_Sfilter_high_frequency_limit(s1) RESULT(res)

USE type_specifications
USE general_module
USE constants
  
  IMPLICIT NONE
  
  type(Sfilter), intent(IN) :: s1
  
  real(dp)                  :: res

! local variables  
  integer :: aorder,border

!START

    aorder=s1%a%order
    border=s1%b%order
    
    if (aorder.LT.border) then
      res=0d0
    else if (aorder.EQ.border) then
      res=s1%a%coeff(aorder)/s1%b%coeff(border)
    else if (aorder.EQ.border+1) then
      res=s1%a%coeff(aorder-1)/s1%b%coeff(border)
    else
      run_status='ERROR in evaluate_Sfilter_high_frequency_limit. Numerator order is more than one greater than denominator order'
      CALL write_program_status()
      STOP 1
    end if
  
  RETURN
  END FUNCTION evaluate_Sfilter_high_frequency_limit
!
! NAME
!      evaluate_Sfilter_high_frequency_resistance
!     
! DESCRIPTION
!       Evaluate the real part of the frequency response of a Laplace domain filter at
!     high frequency i.e. as f-> infinity. Filter format is Rational Function
!
! SEE ALSO
!
!
! HISTORY
!
!     started 21/04/2016 CJS
!
  FUNCTION evaluate_Sfilter_high_frequency_resistance(s1) RESULT(res)

USE type_specifications
USE general_module
USE constants
  
  IMPLICIT NONE
  
  type(Sfilter), intent(IN) :: s1
  
  real(dp)                  :: res

! local variables  
  integer :: aorder,border

!START

    aorder=s1%a%order
    border=s1%b%order
    
    if (aorder.LT.border) then
      res=0d0
    else if (aorder.EQ.border) then
      res=s1%a%coeff(aorder)/s1%b%coeff(border)
    else
      run_status='ERROR in evaluate_Sfilter_high_frequency_limit. Numerator order is greater than denominator order'
      CALL write_program_status()
      STOP 1
    end if
  
  RETURN
  END FUNCTION evaluate_Sfilter_high_frequency_resistance

! NAME
!      Convert_filter_S_to_S_PR
!     
! DESCRIPTION
!       Convert filter from rational function format to format is Pole-Residue format
!       Remove any poles with zero residues
! SEE ALSO
!
!
! HISTORY
!
!     started 4/01/13 CJS
!
  FUNCTION Convert_filter_S_to_S_PR(S) RESULT(res)

USE type_specifications
USE general_module
USE constants
USE maths
  
  IMPLICIT NONE
  
  type(Sfilter),intent(in)      :: S
  	
  type(Sfilter_PR)              :: res

! local variables  
  
  integer                 :: order
  complex(dp),allocatable :: local_poly_roots(:)
  complex(dp),allocatable :: real_roots(:)
  complex(dp),allocatable :: complex_roots(:)
  integer                 :: nreal
  integer                 :: ncomplex
  type(Sfilter_PR)        :: PR_filter
  type(Sfilter_PR)        :: PR_filter2
  
  type(polynomial)        :: Q
  type(polynomial)        :: R
    
  integer :: i
  integer :: row,col
  
  integer                 :: matrix_dim
  integer                 :: ierr
  
  complex(dp),allocatable :: M(:,:)
  complex(dp),allocatable :: MI(:,:)
  complex(dp),allocatable :: LHS(:)
  complex(dp),allocatable :: RHS(:)
  
  type(complex_polynomial):: P1,P_term
  
  integer :: n_zero_real_residues,n_zero_complex_residues,pole,pole_count

!START
  if (verbose) write(*,*)'CALLED Convert_filter_S_to_S_PR'

  if (S%a%order.GT.S%b%order+1) then
    write(*,*)'Error in Convert_filter_S_to_S_PR'
    write(*,*)'numerator order should be less than or equal to denominator order+1'
    STOP
  end if
  
! The rational function order is the number of poles
  order=S%b%order
  PR_filter%order=order
  PR_filter%wnorm=S%wnorm
  
  PR_filter%n_complex_pole_pairs=0
  PR_filter%n_complex_poles=0
  PR_filter%n_real_poles=0
  PR_filter%R=0D0
  PR_filter%L=0D0
    
! Divide the numerator polynomial by the denominator polynomial to extract R+sL terms if they exist
! The new numerator polynomial is now R and has order one less than the denominator
!  CALL divide_poly(S%a,S%b,Q,R,.TRUE.)
  CALL divide_poly(S%a,S%b,Q,R,.FALSE.)
  
  if (Q%order.EQ.1) then
    PR_filter%R=Q%coeff(0)
    PR_filter%L=Q%coeff(1)
    if (verbose) write(*,*)'Setting R=',PR_filter%R,' L=',PR_filter%L
  else if(Q%order.EQ.0) then
    PR_filter%R=Q%coeff(0)
    if (verbose) write(*,*)'Setting R=',PR_filter%R
  end if
  
! special case for zero remainder: no further action required
  if ((R%order.eq.0).AND.(S%b%order.EQ.0))then  
    write(*,*)'Remainder is zero: returning'
    res=PR_filter
    RETURN    
  end if
  
! calculate poles 
  
  ALLOCATE ( local_poly_roots(1:order) )
  ALLOCATE ( real_roots(1:order) )
  ALLOCATE ( complex_roots(1:order) )

  CALL findroots(S%b,local_poly_roots,order)
  CALL rootsort(order,local_poly_roots,real_roots,  &
                complex_roots,nreal,ncomplex,order)
   		
  PR_filter%n_complex_pole_pairs=ncomplex
  PR_filter%n_complex_poles=ncomplex*2
  PR_filter%n_real_poles=nreal
  
  if (verbose) then
    write(*,*)'N complex poles=',PR_filter%n_complex_poles
    write(*,*)'N real poles=',PR_filter%n_real_poles
    do i=1,order
      write(*,*)'pole',i,local_poly_roots(i)
    end do
  end if
  
  if (ALLOCATED( PR_filter%poles )) DEALLOCATE(PR_filter%poles)
  ALLOCATE( PR_filter%poles(1:order) )
  PR_filter%poles(1:nreal)=real_roots(1:nreal)
  PR_filter%poles(nreal+1:nreal+ncomplex*2)=complex_roots(1:ncomplex*2)
  
  if (ALLOCATED( PR_filter%complex_pole )) DEALLOCATE(PR_filter%complex_pole)
  ALLOCATE( PR_filter%complex_pole(1:order) )
  PR_filter%complex_pole(1:nreal)=.FALSE.
  PR_filter%complex_pole(nreal+1:nreal+ncomplex*2)=.TRUE.
  
  DEALLOCATE ( local_poly_roots )
  DEALLOCATE ( real_roots )
  DEALLOCATE ( complex_roots )

! Allocate the matrices M  and MI plus LHS and RHS for matrix equation
  matrix_dim=order
  ALLOCATE( M(matrix_dim,matrix_dim) )
  ALLOCATE( MI(matrix_dim,matrix_dim) )
  ALLOCATE( LHS(matrix_dim) )
  ALLOCATE( RHS(matrix_dim) )
  
  M(1:matrix_dim,1:matrix_dim)=(0d0,0d0)
  MI(1:matrix_dim,1:matrix_dim)=(0d0,0d0)
  LHS(1:matrix_dim)=(0d0,0d0)
  RHS(1:matrix_dim)=(0d0,0d0)
  
! get the LHS as the numerator polynomial of the rational function filter 
  do i=0,order-1
    if (i.LE.R%order) then
      LHS(1+i)=R%coeff(i)/S%b%coeff(order)
    else
      LHS(1+i)=(0d0,0d0)
    end if
  end do
  
! set up P_term as a first order polynomial (used in the form -pole+x )
  P_term=allocate_complex_polynomial(1)
  P_term%coeff(1)=(1d0,0d0)
  
  do col=1,order
  
    P1=(1d0,0d0)   ! initialise polynomial
    
    do i=1,order
      if (i.ne.col) then ! don't include the pole associated with this residue
        P_term%coeff(0)=-PR_filter%poles(i)
        P1=P1*P_term
      end if
    end do
    
! put the coefficients of the polynomial into the matrix M
    do row=1,order
      M(row,col)=P1%coeff(row-1)
    end do
    
  end do ! next row of the matrix M
  
! Invert the matrix M and solve for the coefficients of the pole/ residue expansion
  
  CALL cinvert_Gauss_Jordan(M,matrix_dim,MI,matrix_dim,ierr)
  
  if (ierr.NE.0) then
    write(*,*)'Error in cinvert_Gauss_Jordan. Singular matrix'
    write(*,*)'This may be caused by a multiple pole,'
    write(*,*)'this case has not been allowed for at the moment'
    STOP
  end if
  
  RHS=matmul(MI,LHS)
  
  if (ALLOCATED( PR_filter%residues )) DEALLOCATE(PR_filter%residues)
  ALLOCATE( PR_filter%residues(1:order) )
  PR_filter%residues(1:order)=RHS(1:order)

  DEALLOCATE( M )
  DEALLOCATE( MI )
  DEALLOCATE( LHS )
  DEALLOCATE( RHS )
 
! look for zero residues and remove these poles
  n_zero_real_residues=0
  n_zero_complex_residues=0
  pole=0
  do i=1,PR_filter%n_real_poles
    pole=pole+1
    if (abs(PR_filter%residues(pole)).LT.zero_test_small) n_zero_real_residues= n_zero_real_residues+1
  end do
  do i=1,PR_filter%n_complex_pole_pairs
    pole=pole+1
    if (abs(PR_filter%residues(pole)).LT.zero_test_small) n_zero_complex_residues= n_zero_complex_residues+2
    pole=pole+1  ! skip over the conjugate pole
  end do
  
  if (verbose) then
    write(*,*)'Number of real residues to remove   :',n_zero_real_residues
    write(*,*)'Number of complex residues to remove:',n_zero_complex_residues
  end if
  
! PR_filter2 has the terms with these zero residues removed

  PR_filter2%wnorm=PR_filter%wnorm
  PR_filter2%order=PR_filter%order-n_zero_real_residues-n_zero_complex_residues
  PR_filter2%R    =PR_filter%R
  PR_filter2%L    =PR_filter%L
  
  PR_filter2%n_complex_poles     =PR_filter%n_complex_poles-n_zero_complex_residues
  PR_filter2%n_complex_pole_pairs=PR_filter%n_complex_pole_pairs-n_zero_complex_residues/2
  PR_filter2%n_real_poles        =PR_filter%n_real_poles-n_zero_real_residues

  order=PR_filter2%order
  ALLOCATE( PR_filter2%complex_pole(1:order) ) 
  ALLOCATE( PR_filter2%poles(1:order) ) 
  ALLOCATE( PR_filter2%residues(1:order) ) 
  
  pole=0
  pole_count=0
  do i=1,PR_filter%n_real_poles
    pole_count=pole_count+1
    if (abs(PR_filter%residues(pole_count)).GE.zero_test_small) then
      pole=pole+1
      PR_filter2%poles(pole)   =PR_filter%poles(pole_count)
      PR_filter2%residues(pole)=PR_filter%residues(pole_count)
    else
      if (verbose) then
        write(*,*)'Removing real pole   :',PR_filter%poles(pole_count)
        write(*,*)'Removing real residue:',PR_filter%residues(pole_count)
      end if
    end if
  end do
  do i=1,PR_filter%n_complex_pole_pairs
    pole_count=pole_count+1
    if (abs(PR_filter%residues(pole_count)).GE.zero_test_small) then
      pole=pole+1
      PR_filter2%poles(pole)     =PR_filter%poles(pole_count)
      PR_filter2%residues(pole)  =PR_filter%residues(pole_count) 
      PR_filter2%poles(pole+1)   =PR_filter%poles(pole_count+1)
      PR_filter2%residues(pole+1)=PR_filter%residues(pole_count+1) 
      pole=pole+1  ! skip over the conjugate pole
    else
      if (verbose) then
        write(*,*)'Removing complex pole   :',PR_filter%poles(pole_count)
        write(*,*)'Removing complex pole   :',PR_filter%poles(pole_count+1)
        write(*,*)'Removing real residue   :',PR_filter%residues(pole_count)
        write(*,*)'Removing real residue   :',PR_filter%residues(pole_count+1)
      end if
    end if
    pole_count=pole_count+1  ! skip over the conjugate pole
  end do
  
  res=PR_filter2
  
  CALL deallocate_Sfilter_PR(PR_filter)
  CALL deallocate_Sfilter_PR(PR_filter2)
  
  RETURN
  END FUNCTION Convert_filter_S_to_S_PR
!
! NAME
!      Convert_filter_S_PR_to_S
!     
! DESCRIPTION
!       Convert filter from Pole-Residue format to rational function format 
!
! SEE ALSO
!
!
! HISTORY
!
!     started 4/01/13 CJS
!
  FUNCTION Convert_filter_S_PR_to_S(S_PR) RESULT(res)

USE type_specifications
USE general_module
USE constants
USE maths
  
  IMPLICIT NONE
  
  type(Sfilter_PR), intent(in)  :: S_PR
  	
  type(Sfilter)                 :: res

! local variables  
  
  integer                       :: order,aorder,border

  type(Sfilter)                 :: S_filter
  
  integer                       :: i,term
  
  type(complex_polynomial)      :: P1,P2,P_term
  complex(dp)                   :: Ctemp
  
  real(dp)                      :: norm

!START

! order is the number of poles in the pole-residue function
  order=S_PR%order
  
  border=order

  aorder=order-1
  if (S_PR%R.NE.0d0) aorder=order
  if (S_PR%L.NE.0d0) aorder=order+1
  
  S_filter=allocate_Sfilter(aorder,border)
  S_filter%a%order=aorder
  S_filter%b%order=border
  S_filter%wnorm=S_PR%wnorm
    
! special case for zero order, no further action required
  if ( (aorder.eq.0).AND.(border.eq.0) ) then  
    S_filter%a%coeff(0)=S_PR%R
    S_filter%b%coeff(0)=1d0
    res=S_filter
    RETURN    
  end if
  
! set up P_term as a first order polynomial (used in the form s-pole+ )
  P_term=allocate_complex_polynomial(1)
  P_term%coeff(1)=(1d0,0d0)

! numerator function
  P2=(0d0,0d0)   ! initialise polynomial
  
  do term=1,order
  
    P1=S_PR%residues(term)   ! initialise polynomial
    
    do i=1,order
      if (i.ne.term) then ! don't include the pole associated with this residue
        P_term%coeff(0)=-S_PR%poles(i)
        P1=P1*P_term
      end if
    end do
    
    P2=P2+P1
    
  end do ! next residue term

  if (S_PR%L.NE.0d0) then
! sL term  
    P1=allocate_complex_polynomial(1) ! initialise polynomial p1(s)=sL
    P1%coeff(0)=0d0
    Ctemp=cmplx(S_PR%L,kind=dp)
    P1%coeff(1)=Ctemp
    do i=1,order
      P_term%coeff(0)=-S_PR%poles(i)
      P1=P1*P_term
    end do
  
    P2=P2+P1
  end if     ! constant term

  if (S_PR%R.NE.0d0) then
! constant term  
    Ctemp=cmplx(S_PR%R,kind=dp)
    P1=Ctemp           ! initialise polynomial
    do i=1,order
      P_term%coeff(0)=-S_PR%poles(i)
      P1=P1*P_term
    end do
  
    P2=P2+P1
  end if     ! constant term
  
! put the coefficients of the polynomial into the numerator of S_filter
  S_filter%a%coeff(0:aorder)=P2%coeff(0:aorder)

! denominator function
    
  P1=(1d0,0d0)    ! initialise polynomial
  do i=1,order
    P_term%coeff(0)=-S_PR%poles(i)
    P1=P1*P_term
  end do
  
! put the coefficients of the polynomial into the denominator of S_filter
  S_filter%b%coeff(0:order)=P1%coeff(0:order)

! normalise such that S_filter%b%coeff(0)=1d0 if possible
  norm=S_filter%b%coeff(0)
  if (norm.ne.0d0) then
    S_filter%a%coeff(0:S_filter%a%order)=S_filter%a%coeff(0:S_filter%a%order)/norm
    S_filter%b%coeff(0:S_filter%b%order)=S_filter%b%coeff(0:S_filter%b%order)/norm
  end if

  res=S_filter
  
  RETURN
  END FUNCTION Convert_filter_S_PR_to_S
!
! NAME
!     reciprocal_Sfilter
!
! DESCRIPTION
!     set up a Sfilter type with a given order
!
! HISTORY
!
!     started 22/01/09 CJS
!

FUNCTION reciprocal_Sfilter(H) RESULT(res)

USE type_specifications

IMPLICIT NONE

! argument types
  type(Sfilter),intent(IN)        :: H
  
! Result type
  type(Sfilter)             :: res
  
! local variables
  type(Sfilter)  :: ans

! function definition

  ans%wnorm=1d0
  
  ans%a%order=H%b%order
  allocate (ans%a%coeff(0:ans%a%order))
  ans%a%coeff(0:ans%a%order)=H%b%coeff(0:H%b%order)
  
  ans%b%order=H%a%order
  allocate (ans%b%coeff(0:ans%b%order))
  ans%b%coeff(0:ans%b%order)=H%a%coeff(0:H%a%order)

  res=ans
  
  CALL deallocate_Sfilter(ans)

END FUNCTION reciprocal_Sfilter
!
! NAME
!     renormalise_Sfilter
!
! DESCRIPTION
!     renormalise a filter so it's coefficients are more equal
!
! HISTORY
!
!     started 10/10/09 CJS
!

FUNCTION renormalise_Sfilter(H) RESULT(res)

USE type_specifications

IMPLICIT NONE

! argument types
  type(Sfilter),intent(IN)        :: H
  
! Result type
  type(Sfilter)             :: res
  
! local variables
  type(Sfilter)             :: ans
  integer   :: aorder,border
  real(dp)  :: anorm,bnorm,norm,p
  
  integer :: i

! function definition

  write(*,*)'CALLED renormalise_Sfilter'
  
  write(*,*)'CALLED with function:'
  CALL write_Sfilter(H,0)

  aorder=H%a%order  
  p=1d0/dble(aorder)
  anorm=(H%a%coeff(0)/H%a%coeff(aorder))**p

  border=H%b%order
  p=1d0/dble(border)
  bnorm=(H%b%coeff(0)/H%b%coeff(border))**p

! choose the normalisation to be the average of anorm and bnorm

  norm=(anorm+bnorm)/2d0
  
  write(*,*)'anorm=',anorm
  write(*,*)'bnorm=',bnorm
  write(*,*)'norm =',norm
  
  ans%wnorm=H%wnorm*norm
  
  ans%a%order=H%a%order
  allocate (ans%a%coeff(0:ans%a%order))
  
  anorm=1d0
  do i=0,ans%a%order
    ans%a%coeff(i)=H%a%coeff(i)*anorm
    anorm=anorm*norm
  end do
  
  ans%b%order=H%b%order
  allocate (ans%b%coeff(0:ans%b%order))
  
  bnorm=1d0
  do i=0,ans%b%order
    ans%b%coeff(i)=H%b%coeff(i)*bnorm
    bnorm=bnorm*norm
  end do
  
  res=ans
  
  write(*,*)'RETURNING with function:'
  CALL write_Sfilter(res,0)
  
  CALL deallocate_Sfilter(ans)

END FUNCTION renormalise_Sfilter
