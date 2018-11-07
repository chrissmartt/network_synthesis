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
! SUBROUTINE dwrite_matrix(a,ar,ac,dim,unit)
! SUBROUTINE dread_matrix(a,ar,ac,dim,unit)
! SUBROUTINE dinvert_Gauss_Jordan(A,ar,AI,dim) 
!_____________________________________________________________________
!
! NAME
!    dwrite_matrix
!
! DESCRIPTION
!    write real(dp) matrix to file or screen
!
!
! SEE ALSO
!   
!
! HISTORY
!
!     started 26/04/15 CJS
!
! COMMENTS
!   
!
  SUBROUTINE dwrite_matrix(a,ar,ac,dim,unit)

! Modules used

USE type_specifications

IMPLICIT NONE

! variables passed to subroutine

  integer,intent(IN)  :: dim             ! matrix dimension
  integer,intent(IN)  :: unit               ! output unit. Set to zero to output to screen
  real(dp),intent(IN) :: a(dim,dim)   ! matrix to write
  integer,intent(IN)  :: ar,ac              ! number of rows and columns to write

! Local variables

  integer row,i
 
! START

  do row=1,ar
    if (unit.NE.0) then
      write(unit,8000)(a(row,i),i=1,ac)
    else
      write(*,8000)(a(row,i),i=1,ac) 
    end if
  end do
  
8000 format (1000ES16.8)

! END
 
  RETURN
  
  END SUBROUTINE dwrite_matrix
!
!_____________________________________________________________________
!
! NAME
!    dread_matrix
!
! DESCRIPTION
!    read real(dp) matrix from file
!
! SEE ALSO
!   
!
! HISTORY
!
!     started 2/12/15 CJS
!
! COMMENTS
!   
!
  SUBROUTINE dread_matrix(a,ar,ac,dim,unit)

! Modules used

USE type_specifications

IMPLICIT NONE

! variables passed to subroutine

  integer,intent(IN)  :: dim             ! matrix dimension
  integer,intent(IN)  :: unit               ! input unit
  real(dp),intent(OUT) :: a(dim,dim)   ! matrix to write
  integer,intent(IN)  :: ar,ac              ! number of rows and columns to read

! Local variables

  integer row,i
 
! START

  do row=1,ar
    if (unit.NE.0) then
      read(unit,*)(a(row,i),i=1,ac)
    else
      read(*,*)(a(row,i),i=1,ac) 
    end if
  end do

! END
 
  RETURN
  
  END SUBROUTINE dread_matrix
!
! NAME
!    dread_matrix
!
! DESCRIPTION
!
! Invert the real matrix A using Gauss Jordan method with pivoting and return the result in AI
! ierr=0 on the successful calculation of the inverse
! If a singular matrix is found then:
! if ierr.EQ.0 on input then the program stops
! if ierr.NE.0 on input then the program returns with ierr=1
!
! HISTORY
!
!     started 2/12/15 CJS
!
! COMMENTS
! 
  SUBROUTINE dinvert_Gauss_Jordan(A,n,AI,dim,ierr) 

USE type_specifications
USE general_module

IMPLICIT NONE

! variables passed to subroutine
       
  integer,intent(IN)    :: dim            ! matrix dimension
  integer,intent(IN)    :: n                 ! size iof matrix to invert
  real(dp),intent(IN)   :: A(dim,dim) ! matrix to invert
  real(dp),intent(OUT)  :: AI(dim,dim) ! matrix inverse
  integer,intent(INOUT) :: ierr

! local variables
                 
  integer    :: row,col,reduce_col,i
  
  real(dp)    :: max_element
  real(dp)    :: pivot_element
  integer     :: max_row
  
  integer    :: pivot_row
  
  integer    :: pivot_row_save(dim)
  
  real(dp)    :: row_multiplier
  real(dp)    :: swap

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
      if (verbose) write(*,*)'Singular matrix found in dinvert_Gauss_Jordan'
      if (ierr.NE.0) then
        run_status='ERROR: Singular matrix in dinvert_Gauss_Jordan'
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
        AI(pivot_row,col) = 1d0/pivot_element
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
  
  END SUBROUTINE dinvert_Gauss_Jordan
