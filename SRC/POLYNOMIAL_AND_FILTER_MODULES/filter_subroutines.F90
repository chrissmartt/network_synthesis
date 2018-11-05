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
!     subroutines related to rational function filter function representations
!     to be included in filter_module.F90
!
! Comments:
!      
!
! History
!
!     started 22/01/09 CJS
!
!
!  SUBROUTINE read_Sfilter(s1,ip_unit)
!  SUBROUTINE write_Sfilter(s1,op_unit)
!  SUBROUTINE test_filter_pole_stability(filter,stable)
!  SUBROUTINE deallocate_Sfilter
!  SUBROUTINE read_Sfilter_matrix(s1,ip_unit)
!  SUBROUTINE write_Sfilter_matrix(s1,op_unit)
!  SUBROUTINE deallocate_Sfilter_matrix
!  SUBROUTINE output_Sfilter_frequency_response(s1,fmin,fmax,fstep,op_unit)
!
! NAME
!     read_Sfilter
!
! DESCRIPTION
!       read Laplace domain filter coefficients from file
!
! SEE ALSO
!
!
! HISTORY
!
!     started 02/03/09 CJS
!
  SUBROUTINE read_Sfilter(s1,ip_unit)

USE type_specifications
  
  type(Sfilter),intent(OUT) :: s1
  integer,intent(IN)        :: ip_unit
  
  integer i

!START

! read filter coefficients from the given file unit

  read(ip_unit,*)s1%wnorm
  read(ip_unit,*)s1%a%order
  if (allocated (s1%a%coeff) ) deallocate(s1%a%coeff)
  allocate( s1%a%coeff(0:s1%a%order) )
  read(ip_unit,*)(s1%a%coeff(i),i=0,s1%a%order)  
  read(ip_unit,*)s1%b%order
  if (allocated (s1%b%coeff) ) deallocate(s1%b%coeff)
  allocate( s1%b%coeff(0:s1%b%order) )
  read(ip_unit,*)(s1%b%coeff(i),i=0,s1%b%order)  
  
  RETURN
  END SUBROUTINE read_Sfilter
!
! NAME
!     write_Sfilter
!
! DESCRIPTION
!       write Laplace domain filter coefficients to screen or file 
!
! SEE ALSO
!
!
! HISTORY
!
!     started 23/01/09 CJS
!
  SUBROUTINE write_Sfilter(s1,op_unit)

USE type_specifications
  
  type(Sfilter),intent(IN):: s1
  integer,intent(IN)      :: op_unit 
  
  integer i

!START

  if (op_unit.eq.0) then
! write to screen
    write(*,*)'Laplace domain filter'
    write(*,*)'wnorm=',s1%wnorm
    write(*,*)'a order=',s1%a%order
    write(*,*)(s1%a%coeff(i),i=0,s1%a%order)  
    write(*,*)'b order=',s1%b%order
    write(*,*)(s1%b%coeff(i),i=0,s1%b%order)  
    write(*,*)' '
  else
    write(op_unit,*)s1%wnorm,'  # w normalisation constant'
    write(op_unit,*)s1%a%order,'  # a order, a coefficients follow below:'
    write(op_unit,*)(s1%a%coeff(i),i=0,s1%a%order)  
    write(op_unit,*)s1%b%order,'  # b order, b coefficients follow below:'
    write(op_unit,*)(s1%b%coeff(i),i=0,s1%b%order)  
  end if
  
  RETURN
  END SUBROUTINE write_Sfilter
!
! NAME
!      output_Sfilter_frequency_response
!     
! DESCRIPTION
!       Evaluate the frequency resaponse of a Laplace domain filter and
!       write the response to file
!
! SEE ALSO
!
!
! HISTORY
!
!     started 23/01/09 CJS
!
  SUBROUTINE output_Sfilter_frequency_response(s1,fmin,fmax,fstep,op_unit)

USE type_specifications
  
  type(Sfilter),intent(IN) :: s1
  real(dp),intent(IN)      :: fmin,fmax,fstep
  integer,intent(IN)       :: op_unit 

! local variables  
  real(dp) f,w
  complex(dp) jw,num,den,response
  integer n
  integer n_frequencies,frequency_loop

!START

!  do f=fmin,fmax,fstep

  n_frequencies=int( (fmax-fmin)/fstep )+1  ! number of frequency samples

!loop over frequency
  do frequency_loop=1,n_frequencies

    f=fmin+(frequency_loop-1)*fstep   ! calculate frequency of this sample
  
    response=evaluate_Sfilter_frequency_response(s1,f)
        
    write(op_unit,1000)f,real(response),aimag(response)
1000 format(3E16.7)

  end do ! next frequency
  
  RETURN
  END SUBROUTINE output_Sfilter_frequency_response
!
! NAME
!     test_filter_pole_stability
!
! DESCRIPTION
!     work out the poles of a filter and check that they are in the LHS of the s-plane
!
! SEE ALSO
!
!
! HISTORY
!
!     started 02/03/09 CJS
!
  SUBROUTINE test_filter_pole_stability(filter,stable)

USE type_specifications
USE general_module
USE constants
USE eispack

! variables passed to subroutine  
  type(Sfilter),intent(IN)    :: filter
  logical,intent(OUT)         :: stable
  
! local variables

  integer   :: poly_order
  integer  :: n
  integer  :: low
  integer  :: igh
  real(dp),allocatable  :: h(:,:)
  real(dp),allocatable  :: wi(:)
  real(dp),allocatable  :: wr(:)
  integer  :: i
  integer  :: ierr
   
!START

! assume the filter is stable initially
  stable=.TRUE.

  poly_order=filter%b%order     ! order of the denominator polynomial function 
  
  if (poly_order.EQ.0) RETURN

! CALL EISPACK subroutine to calculate the eigenvalues of an upper Hessenberg matrix
! which give the eigenvalues of the polynomial

  ALLOCATE( h(1:poly_order,1:poly_order) )
  ALLOCATE( wr(1:poly_order) )
  ALLOCATE( wi(1:poly_order) )

! Fill the upper Hessenberg matrix
  
  h(1:poly_order,1:poly_order)=0d0
  
  do i=1,poly_order
    h(1,i)=-filter%b%coeff(poly_order-i)/filter%b%coeff(poly_order)  ! first row is calculated from the polynomial coefficients
    if (i.ne.poly_order) then
      h(i+1,i)=1d0
    end if
  end do

  n=poly_order
  low=1
  igh=n
  
! EISPACK subroutine to calculate eigenvalues of upper Hessenberg matrix  
! note that the eigenvalues are not ordered (apart from in complex pairs)
  CALL hqr ( n, low, igh, h, wr, wi, ierr )  
  
  if (ierr.ne.0) then
    run_status='ERROR in test_filter_pole_stability. Not all roots of the polynomial could be found correctly'
    CALL write_program_status()
    STOP 1
  end if

! loop over the poles
  do i=1,n
  
! check the real part of each of the poles. The filter is unstable if the poles are on the RHS of the s-plane ie Re{pole}.GT.0 
    if (wr(i).GT.zero_test_small) then   
! unstable pole
      stable=.FALSE.
    end if

  end do
  
  if (.NOT.stable) then
  
    write(*,*)'Poles are:'
    do i=1,n
    
      write(*,'(I4,ES12.5,A4,ES12.5)')i,wr(i),' +j ',wi(i)
    
    end do
  end if
  
  DEALLOCATE( h )
  DEALLOCATE( wr )
  DEALLOCATE( wi )
  
  RETURN
  END SUBROUTINE test_filter_pole_stability
!
! NAME
!     test_filter_simple_poles
!
! DESCRIPTION
!     work out the poles of a filter and check that there are only simple poles on the s=jw axis
!
! SEE ALSO
!
!
! HISTORY
!
!     started 02/03/09 CJS
!
  SUBROUTINE test_filter_simple_poles(filter,stable)

USE type_specifications
USE general_module
USE constants
USE eispack

! variables passed to subroutine  
  type(Sfilter),intent(IN)    :: filter
  logical,intent(OUT)         :: stable
  
! local variables

  integer   :: poly_order
  integer  :: n
  integer  :: low
  integer  :: igh
  real(dp),allocatable  :: h(:,:)
  real(dp),allocatable  :: wi(:)
  real(dp),allocatable  :: wr(:)
  integer  :: i,ii
  integer  :: ierr
   
!START

! assume the filter is stable initially
  stable=.TRUE.

  poly_order=filter%b%order     ! order of the denominator polynomial function 
  
  if (poly_order.EQ.0) RETURN

! CALL EISPACK subroutine to calculate the eigenvalues of an upper Hessenberg matrix
! which give the eigenvalues of the polynomial

  ALLOCATE( h(1:poly_order,1:poly_order) )
  ALLOCATE( wr(1:poly_order) )
  ALLOCATE( wi(1:poly_order) )

! Fill the upper Hessenberg matrix
  
  h(1:poly_order,1:poly_order)=0d0
  
  do i=1,poly_order
    h(1,i)=-filter%b%coeff(poly_order-i)/filter%b%coeff(poly_order)  ! first row is calculated from the polynomial coefficients
    if (i.ne.poly_order) then
      h(i+1,i)=1d0
    end if
  end do

  n=poly_order
  low=1
  igh=n
  
! EISPACK subroutine to calculate eigenvalues of upper Hessenberg matrix  
! note that the eigenvalues are not ordered (apart from in complex pairs)
  CALL hqr ( n, low, igh, h, wr, wi, ierr )  
  
  if (ierr.ne.0) then
    run_status='ERROR in test_filter_pole_stability. Not all roots of the polynomial could be found correctly'
    CALL write_program_status()
    STOP 1
  end if

! loop over the poles
  do i=1,n-1
  
! look for poles on the imaginary axis ie Re{pole}.EQ.0 
    if (wr(i).EQ.0d0) then   
! this pole is on the jw axis
! now check all other poles to see if this is a simple pole or not

      do ii=i+1,n
    
        if ( (wr(ii).EQ.0d0).AND.(abs(wi(i)-wi(ii)).LT.zero_test_small) ) then
!the pole we are checking is on the jw axis and is very close to the pole we are checking against
          stable=.FALSE.
        end if
        
      end do ! next pole to check
      
    end if

  end do
  
  if (.NOT.stable) then
  
    write(*,*)'Poles are:'
    do i=1,n
    
      write(*,'(I4,ES12.5,A4,ES12.5)')i,wr(i),' +j ',wi(i)
    
    end do
  end if
  
  DEALLOCATE( h )
  DEALLOCATE( wr )
  DEALLOCATE( wi )
  
  RETURN
  END SUBROUTINE test_filter_simple_poles
  
!
! NAME
!      deallocate_Sfilter(ft)
!      
! DESCRIPTION
!       deallocate Sfilter type
!     
!
! SEE ALSO
!
!
! HISTORY
!
!     started 26/09/12 CJS
!
  SUBROUTINE deallocate_Sfilter(s1)

USE type_specifications

IMPLICIT NONE
  
  type(Sfilter),intent(INOUT)  :: s1

! local variables  
  integer i

!START

  if (allocated(  s1%a%coeff )) then
    DEALLOCATE( s1%a%coeff )
  end if 
  if (allocated(  s1%b%coeff )) then
    DEALLOCATE( s1%b%coeff )
  end if 
  
  RETURN
  END SUBROUTINE deallocate_Sfilter
!
! NAME
!      deallocate_Sfilter_PR(ft)
!      
! DESCRIPTION
!       deallocate Sfilter_PR type
!     
!
! SEE ALSO
!
!
! HISTORY
!
!     started 14/09/17 CJS
!
  SUBROUTINE deallocate_Sfilter_PR(s1)

USE type_specifications

IMPLICIT NONE
  
  type(Sfilter_PR),intent(INOUT)  :: s1

! local variables  
  integer i

!START

  if (allocated(  s1%complex_pole )) then
    DEALLOCATE( s1%complex_pole )
  end if 
  if (allocated(  s1%poles )) then
    DEALLOCATE( s1%poles )
  end if 
  if (allocated(  s1%residues )) then
    DEALLOCATE( s1%residues )
  end if 
  
  RETURN
  END SUBROUTINE deallocate_Sfilter_PR
!
! NAME
!     read_Sfilter_matrix
!
! DESCRIPTION
!       read a matrix of Laplace domain filter coefficients from file
!
! SEE ALSO
!
!
! HISTORY
!
!     started 21/04/2016 CJS
!
  SUBROUTINE read_Sfilter_matrix(s1,ip_unit)

USE type_specifications

IMPLICIT NONE

! variables passed to subroutine  
  type(Sfilter_matrix),intent(OUT) :: s1
  integer ,intent(IN)              :: ip_unit
  
! local variables
  integer    :: row,col
  integer    :: dim

! START
  
  read(ip_unit,*)     ! comment line
  read(ip_unit,*)dim  ! dimension of matrix
  
  s1%dim=dim
  
  if (dim.GT.0) then
  
    ALLOCATE( s1%sfilter_mat(1:dim,1:dim) )
  
    do row=1,dim
      do col=1,dim

        read(ip_unit,*)  ! comment line
      
        CALL read_Sfilter(s1%sfilter_mat(row,col),ip_unit)
        
      end do
    end do
  
  end if ! matrix dimension.GT.0
  
  RETURN
  END SUBROUTINE read_Sfilter_matrix
!
! NAME
!     write_Sfilter_matrix
!
! DESCRIPTION
!       write Laplace domain filter coefficients to screen or file 
!
! SEE ALSO
!
!
! HISTORY
!
!     started 21/04/2016 CJS
!
  SUBROUTINE write_Sfilter_matrix(s1,op_unit)

USE type_specifications

IMPLICIT NONE
  
  type(Sfilter_matrix),intent(IN) :: s1
  integer,intent(IN)              :: op_unit 
  
! local variables  

  integer :: row,col

!START

  if (op_unit.eq.0) then
    write(*,*)'Matrix of Laplace domain filters'
    write(*,*)s1%dim,'   # dimension of Sfilter matrix'
  else
    write(op_unit,*)'Matrix of Laplace domain filters'
    write(op_unit,*)s1%dim,'   # dimension of Sfilter matrix'
  end if
  
  do row=1,s1%dim
    do col=1,s1%dim

      if (op_unit.eq.0) then
        write(*,8000)'Element:',row,col
      else
        write(op_unit,8000)'Element:',row,col
      end if
8000 format(A8,2I4)
      
      CALL write_Sfilter(s1%sfilter_mat(row,col),op_unit)
        
    end do
  end do
  
  RETURN
  END SUBROUTINE write_Sfilter_matrix
!
! NAME
!      deallocate_Sfilter_matrix(ft)
!      
! DESCRIPTION
!       deallocate Sfilter_matrix type
!     
!
! SEE ALSO
!
!
! HISTORY
!
!     started 26/09/12 CJS
!
  SUBROUTINE deallocate_Sfilter_matrix(s1)

USE type_specifications

IMPLICIT NONE
  
  type(Sfilter_matrix),intent(INOUT)  :: s1

! local variables  

  integer :: row,col

!START

  if (allocated(  s1%sfilter_mat )) then
  
    do row=1,s1%dim
      do col=1,s1%dim
       
        CALL deallocate_Sfilter(s1%sfilter_mat(row,col))
        
      end do
    end do
    
   DEALLOCATE( s1%sfilter_mat )
   
  end if 
  
  RETURN
  END SUBROUTINE deallocate_Sfilter_matrix
!
! NAME
!     write_S_PR_filter
!
! DESCRIPTION
!       write Laplace domain filter in pole/ residue format to screen 
!
! SEE ALSO
!
!
! HISTORY
!
!     started 4/01/13 CJS
!
  SUBROUTINE write_S_PR_filter(s1)

USE type_specifications

IMPLICIT NONE
  
  type(Sfilter_PR) :: s1
  
  integer i

!START

! write to screen
  write(*,*)'Laplace domain filter, Pole Residue format'
  write(*,*)'wnorm=',s1%wnorm
  write(*,*)'order=',s1%order
  write(*,*)'n_complex_poles     =',s1%n_complex_poles
  write(*,*)'n_complex_pole_pairs=',s1%n_complex_pole_pairs
  write(*,*)'n_real_poles        =',s1%n_real_poles
  
! top line of function (numerator terms)

  write(*,'(A6)',advance='NO')'H(s) ='
  
  if (s1%R.NE.0d0) then
    write(*,8000,advance='NO')s1%R
  end if
  
  if (s1%L.NE.0d0) then
    write(*,8010,advance='NO')' + ',s1%L,'s '
  end if
  
  do i=1,s1%order
    write(*,8020,advance='NO')' + ',s1%residues(i),'     '
  end do

  write(*,*)

! write horizontal lines as required

  write(*,'(A6)',advance='NO')'      '
  
  if (s1%R.NE.0d0) then
    write(*,'(A5)',advance='NO')'     '
  end if
  
  if (s1%L.NE.0d0) then
    write(*,'(A10)',advance='NO')'          '
  end if
  
  do i=1,s1%order
    write(*,'(A19)',advance='NO')'---------------    '
  end do

  write(*,*)
  
! write denominator terms

  write(*,'(A6)',advance='NO')'      '
  
  if (s1%R.NE.0d0) then
    write(*,'(A5)',advance='NO')'     '
  end if
  
  if (s1%L.NE.0d0) then
    write(*,'(A10)',advance='NO')'          '
  end if
  
  do i=1,s1%order
    write(*,8020,advance='NO')'s-(',s1%poles(i),')    '
  end do

  write(*,*)
  
8000 format(F5.2)
8010 format(A3,F5.2,A2)
8020 format(A3,F5.2,SP,F5.2,"i",A5)

  
  RETURN
  END SUBROUTINE write_S_PR_filter
