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
! FILE CONTENTS:
!SUBROUTINE write_circle
!SUBROUTINE write_rectangle
!SUBROUTINE write_Dshape
!
!SUBROUTINE write_circle
!
! NAME
!     SUBROUTINE write_circle
!
! DESCRIPTION
!     write a circle with specified x,y centre and radius to file for plotting with gnuplot
!
!     
! COMMENTS
!     return the extent of the plotting area...
!
! HISTORY
!
!     started 10/05/2013 CJS
!     using generate_shapes.F90 20/4/2017 CJS
!
!
SUBROUTINE write_circle(x,y,r,unit,xmin,xmax,ymin,ymax)
   
USE type_specifications
USE constants

IMPLICIT NONE

  real(dp),intent(IN)    :: x,y,r                 ! centre x and y coordinates and radius
  real(dp),intent(INOUT) :: xmin,xmax,ymin,ymax   ! extent of the plotting area
  integer, intent(IN)    :: unit                  ! unit to write to  
  
! local variables  
  
  integer :: npts
  real(dp),allocatable :: shape_x(:)
  real(dp),allocatable :: shape_y(:)
  
  integer :: i
  
! START

  CALL generate_circle_points(npts,shape_x,shape_y,x,y,r)
  
  do i=1,npts
    
    write(unit,8000)shape_x(i),shape_y(i)
8000 format (4E14.6)

    xmin=min(xmin,shape_x(i))
    xmax=max(xmax,shape_x(i))
    ymin=min(ymin,shape_y(i))
    ymax=max(ymax,shape_y(i))
  
  end do
  
  write(unit,*)
  write(unit,*)

  DEALLOCATE( shape_x )
  DEALLOCATE( shape_y )
   
  RETURN
  
END SUBROUTINE write_circle
!
!SUBROUTINE write_rectangle
!
! NAME
!     SUBROUTINE write_rectangle
!
! DESCRIPTION
!     write a rectangle with specified x,y centre, width, height and angle to file for plotting with gnuplot
!
!     
! COMMENTS
!     return the extent of the plotting area...
!
! HISTORY
!
!     started 23/9/2016 CJS
!     using generate_shapes.F90 20/4/2017 CJS
!
!
SUBROUTINE write_rectangle(x,y,w,h,theta,unit,xmin,xmax,ymin,ymax)
   
USE type_specifications
USE constants

IMPLICIT NONE

  real(dp),intent(IN)    :: x,y,w,h,theta         ! centre x and y coordinates, width, height and angle of rectangle
  real(dp),intent(INOUT) :: xmin,xmax,ymin,ymax   ! extent of the plotting area
  integer, intent(IN)    :: unit                  ! unit to write to  
  
! local variables  
  
  integer :: npts
  real(dp),allocatable :: shape_x(:)
  real(dp),allocatable :: shape_y(:)
  
  integer :: i
 
! START

  CALL generate_rectangle_points(npts,shape_x,shape_y,x,y,theta,w,h)
  
  do i=1,npts
  
    write(unit,8000)shape_x(i),shape_y(i)
8000 format (4E14.6)

    xmin=min(xmin,shape_x(i))
    xmax=max(xmax,shape_x(i))
    ymin=min(ymin,shape_y(i))
    ymax=max(ymax,shape_y(i))
  
  end do
  
  write(unit,*)
  write(unit,*)

  DEALLOCATE( shape_x )
  DEALLOCATE( shape_y )
   
  RETURN
  
END SUBROUTINE write_rectangle

!
! NAME
!     SUBROUTINE write_Dshape
!
! DESCRIPTION
!     write a Dshape with specified x,y centre, width1, width2, conductor separation, shell offset and angle to file for plotting with gnuplot
!
!     
! COMMENTS
!     Also return the extent of the plotting area
!
! HISTORY
!
!     started 15/11/2016 CJS
!     using generate_shapes.F90 20/4/2017 CJS
!
!
SUBROUTINE write_Dshape(x,y,w1,w2,s,r,theta,unit,xmin,xmax,ymin,ymax)
   
USE type_specifications
USE constants

IMPLICIT NONE

  real(dp),intent(IN)    :: x,y,w1,w2,s,r,theta  ! centre x and y coordinates, width1, width2, separation of wire rows, offset of D shapeand angle
  real(dp),intent(INOUT) :: xmin,xmax,ymin,ymax  ! extent of the plotting area
  integer, intent(IN)    :: unit                 ! unit to write to  
  
! local variables  

  
! local variables  
  
  integer :: npts
  real(dp),allocatable :: shape_x(:)
  real(dp),allocatable :: shape_y(:)
  
  integer :: i
 
! START

  CALL generate_Dshape_points(npts,shape_x,shape_y,x,y,w1,w2,s,r,theta)
  
  do i=1,npts
  
    write(unit,8000)shape_x(i),shape_y(i)
8000 format (4E14.6)

    xmin=min(xmin,shape_x(i))
    xmax=max(xmax,shape_x(i))
    ymin=min(ymin,shape_y(i))
    ymax=max(ymax,shape_y(i))
  
  end do
  
  write(unit,*)
  write(unit,*)

  DEALLOCATE( shape_x )
  DEALLOCATE( shape_y )
   
  RETURN
  
END SUBROUTINE write_Dshape
