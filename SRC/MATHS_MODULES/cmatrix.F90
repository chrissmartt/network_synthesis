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
!       SUBROUTINE write_cmatrix(Mat,dim,unit)
!       SUBROUTINE write_cmatrix_re(Mat,dim,unit)
!       SUBROUTINE write_cmatrix_im(Mat,dim,unit)
!       SUBROUTINE cinvert_Gauss_Jordan(A,n,AI,dim) 
!       SUBROUTINE c_condition_number (A,n,condition_number,dim)
!
! NAME
!    write_cmatrix
!
! DESCRIPTION
!    write complex(dp) matrix to file or screen
!
!
! HISTORY
!
!     started 2/12/15 CJS
!
! COMMENTS
! 
SUBROUTINE write_cmatrix(Mat,dim,unit)

USE type_specifications

IMPLICIT NONE

! variables passed to subroutine

integer,intent(IN)      :: dim          ! matrix dimension
integer,intent(IN)      :: unit         ! unit to write to. Set to zero for screen output.

complex(dp),intent(IN) :: Mat(dim,dim)  ! matrix to write

! local variables

integer row,col

! START

do row=1,dim

  if (unit.EQ.0) then
    write(*,*)(Mat(row,col),col=1,dim)
  else
    write(unit,*)(Mat(row,col),col=1,dim)
  end if
  
end do
8000 format(1000ES16.6)


END SUBROUTINE write_cmatrix
!
! NAME
!    write_cmatrix_re
!
! DESCRIPTION
!   write the real part of a complex(dp) matrix to file or screen
!
! HISTORY
!
!     started 2/12/15 CJS
!
! COMMENTS
! 
SUBROUTINE write_cmatrix_re(Mat,dim,unit)

USE type_specifications

IMPLICIT NONE

! variables passed to subroutine

integer,intent(IN)      :: dim          ! matrix dimension
integer,intent(IN)      :: unit         ! unit to write to. Set to zero for screen output.

complex(dp),intent(IN) :: Mat(dim,dim)  ! matrix to write

! local variables

integer row,col

! START

do row=1,dim

  if (unit.EQ.0) then
    write(*,*)(dble(Mat(row,col)),col=1,dim)
  else
    write(unit,*)(dble(Mat(row,col)),col=1,dim)
  end if
  
end do
8000 format(20ES16.6)


END SUBROUTINE write_cmatrix_re
!
! NAME
!   write_cmatrix_im
!
! DESCRIPTION
!   write the real part of a complex(dp) matrix to file or screen
!
!
! HISTORY
!
!     started 2/12/15 CJS
!
! COMMENTS
! 
SUBROUTINE write_cmatrix_im(Mat,dim,unit)

USE type_specifications

IMPLICIT NONE

! variables passed to subroutine

integer,intent(IN)      :: dim          ! matrix dimension
integer,intent(IN)      :: unit         ! unit to write to. Set to zero for screen output.

complex(dp),intent(IN) :: Mat(dim,dim)  ! matrix to write

! local variables

integer row,col

! START

do row=1,dim

  if (unit.EQ.0) then
    write(*,*)(AIMAG(Mat(row,col)),col=1,dim)
  else
    write(unit,*)(AIMAG(Mat(row,col)),col=1,dim)
  end if
  
end do
8000 format(20ES16.6)


END SUBROUTINE write_cmatrix_im
!
! NAME
!    cinvert_Gauss_Jordan
!
! DESCRIPTION
!  
! Invert the complex matrix A using Gauss Jordan method with pivoting and return the result in AI
! ierr=0 on the successful calculation of the inverse
! if a singular matrix is found then
! if ierr.EQ.0 on input then the program stops
! if ierr.NE.0 on input then the program returns with ierr=1
!
! HISTORY
!
!     started 2/12/15 CJS
!
! COMMENTS
! 
  SUBROUTINE cinvert_Gauss_Jordan(A,n,AI,dim,ierr) 

USE type_specifications
USE general_module

IMPLICIT NONE

! variables passed to subroutine

  integer,intent(IN)      :: dim       ! matrix dimension
  integer,intent(IN)      :: n         ! size of matrix to invert

  complex(dp),intent(IN)  :: A(dim,dim)  ! matrix to invert
  complex(dp),intent(OUT) :: AI(dim,dim) ! inverse matrix
       
  integer,intent(INOUT)  :: ierr        ! error code

! local variables
                 
  integer    :: row,col,reduce_col,i
  
  real(dp)    :: max_element
  complex(dp)    :: pivot_element
  integer    :: max_row
  
  integer    :: pivot_row
  
  integer    :: pivot_row_save(dim)
  
  complex(dp)    :: row_multiplier
  complex(dp)    :: swap

! START
  
! copy A to AI 
  AI(1:n,1:n)= A(1:n,1:n)

  pivot_row_save(1:dim)=0
  
! loop over columns of the matrix and reduce each column in turn to identity matrix column     
  do reduce_col=1,n
  
! find the largest element in this column and use as the pivot element
    max_element=0d0
    max_row=0
    do row=reduce_col,n
      if (abs(AI(row,reduce_col)).GT.max_element) then
        max_element=abs(AI(row,reduce_col))
        max_row=row
      end if
    end do  
    
    if (max_row.eq.0) then
! all elements are zero so singular matrix
      if(verbose) write(*,*)'Singular matrix found in cinvert_Gauss_Jordan'
      if (ierr.NE.0) then
        run_status='ERROR: Singular matrix in cinvert_Gauss_Jordan'
        CALL write_program_status()
        STOP 1
      else
        ierr=1
        RETURN
      end if
    end if
    
    pivot_row=max_row
    pivot_row_save(reduce_col)=pivot_row
    
! swap pivot row with the row reduce_col

    if (pivot_row.ne.reduce_col) then
      do col=1,n
        swap=AI(reduce_col,col)
    AI(reduce_col,col)=AI(pivot_row,col)
    AI(pivot_row,col)=swap
      end do 
    end if
    
    pivot_row=reduce_col   
    pivot_element=AI(reduce_col,reduce_col)
    
! operate on pivot row    
    do col=1,n
      if (col.ne.reduce_col) then
        AI(pivot_row,col) = AI(pivot_row,col)/pivot_element
      else    
        AI(pivot_row,col) = (1d0,0d0)/pivot_element
      end if
    end do

! operate on rows other than the pivot row   
    do row=1,n
    
      if (row.ne.pivot_row) then
      
        row_multiplier=AI(row,reduce_col)
    
        do col=1,n
          if (col.ne.reduce_col) then
            AI(row,col) = AI(row,col)- AI(pivot_row,col)*row_multiplier
          else    
            AI(row,reduce_col) =-AI(pivot_row,reduce_col)*row_multiplier
          end if
        end do
    
      end if ! not pivot row
    
    end do ! next row
    
  end do ! next column of the matrix to reduce
  
  do reduce_col=n,1,-1
  
    if (reduce_col.ne.pivot_row_save(reduce_col)) then
! rows were swapped so must swap the corresponding columns

      do row=1,n
        swap=AI(row,pivot_row_save(reduce_col))
        AI(row,pivot_row_save(reduce_col))=AI(row,reduce_col)
    AI(row,reduce_col)=swap
      end do
      
    end if
    
  end do
  
  ierr=0

  RETURN
  
  END SUBROUTINE cinvert_Gauss_Jordan
!
! NAME
!    
!
! DESCRIPTION
!
!
! HISTORY
!
!     started 2/12/15 CJS
!
! COMMENTS
! 
  SUBROUTINE c_condition_number(A,n,condition_number,dim) 
  
! Calculate the condition number of the complex matrix A 

USE type_specifications
USE eispack

IMPLICIT NONE

! variables passed to subroutine

  integer,intent(IN)      :: dim       ! matrix dimension
  integer,intent(IN)      :: n         ! size of matrix to process

  complex(dp),intent(IN) :: A(dim,dim)  ! input matrix
              
  real(dp),intent(OUT)   :: condition_number  ! output condition number

! local variables

  complex(dp)    :: AH(dim,dim)
  complex(dp)    :: AHA(dim,dim)
  real(dp)       :: Real_AHA(dim,dim)
  
  real(dp)       :: singular_values(dim)
                 
  integer     :: row,col
  
  real(dp)    :: max_eigenvalue
  real(dp)    :: min_eigenvalue
  
  logical :: matu,matv
  
  integer :: ierr

! START
  
! calculate the Hermitian conjugate of A

  do row=1,n
    do col=1,n
      AH(row,col)=conjg(A(col,row))
    end do
  end do
  
  AHA=matmul(AH,A)
  Real_AHA=dble(AHA)
  
! calculate the Singular Value Decomposition of AHA using Eispack
  matu=.FALSE.
  matv=.FALSE. ! we don't need the matrices U or V
  CALL svd ( n, n, Real_AHA, singular_values, matu, Real_AHA, matv, Real_AHA, ierr )
  
! find the maximum and minimum magnitude of singular values
  max_eigenvalue=sqrt(abs(singular_values(1)))
  min_eigenvalue=sqrt(abs(singular_values(1)))
  
  do row=2,n
   
! Note that the singular values of A are equal to the square root of the singular values of AHA
    max_eigenvalue=max( max_eigenvalue,sqrt(abs(singular_values(row))) )
    min_eigenvalue=min( min_eigenvalue,sqrt(abs(singular_values(row))) )
  
  end do
  
! calculate the condition number
  if (min_eigenvalue.NE.0D0) then
    condition_number=max_eigenvalue/min_eigenvalue
  else
! set the condition number to something large ***** SHOULD PROBABLY USE A PARAMETER FROM MODULE constants HERE *****
    condition_number=1D100
  end if

  RETURN
  
  END SUBROUTINE c_condition_number
