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
!SUBROUTINE generate_circle_points
!SUBROUTINE generate_rectangle_points
!SUBROUTINE generate_Dshape_points
!SUBROUTINE generate_arc_points
!
!SUBROUTINE generate_circle_points
!
! NAME
!     SUBROUTINE generate_circle_points
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
!
!
SUBROUTINE generate_circle_points(npts,shape_x,shape_y,x,y,r)
   
USE type_specifications
USE constants

IMPLICIT NONE

  integer,intent(OUT)                :: npts
  real(dp),allocatable,intent(INOUT)  :: shape_x(:)
  real(dp),allocatable,intent(INOUT)  :: shape_y(:)
  real(dp),intent(IN)                :: x,y,r                 ! centre x and y coordinates and radius
  
! local variables  

  real(dp) x1,y1
  real(dp) x2,y2
  real(dp) x3,y3
  real(dp) x4,y4
  
  integer :: loop
  
! START

! write the circle as four arcs

  x1=x
  y1=y+r
  
  x2=x-r
  y2=y
  
  x3=x
  y3=y-r
  
  x4=x+r
  y4=y

  do loop=1,2  ! first pass to count the points, second pass to set the point coordinates
    
    npts=0
    
    CALL generate_arc_points(npts,shape_x,shape_y,loop,x,y,x1,y1,x2,y2)
  
    CALL generate_arc_points(npts,shape_x,shape_y,loop,x,y,x2,y2,x3,y3)
  
    CALL generate_arc_points(npts,shape_x,shape_y,loop,x,y,x3,y3,x4,y4)
  
    CALL generate_arc_points(npts,shape_x,shape_y,loop,x,y,x4,y4,x1,y1)
  
    if (loop.EQ.1) then
      ALLOCATE( shape_x(1:npts) )
      ALLOCATE( shape_y(1:npts) )
    end if
  
  end do ! next loop
   
  RETURN
  
END SUBROUTINE generate_circle_points
!
!SUBROUTINE generate_rectangle_points
!
! NAME
!     SUBROUTINE generate_rectangle_points
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
!
!
SUBROUTINE generate_rectangle_points(npts,shape_x,shape_y,x,y,theta,w,h)
   
USE type_specifications
USE constants

IMPLICIT NONE

  integer,intent(OUT)                :: npts
  real(dp),allocatable,intent(INOUT)  :: shape_x(:)
  real(dp),allocatable,intent(INOUT)  :: shape_y(:)

  real(dp),intent(IN)    :: x,y,w,h,theta         ! centre x and y coordinates, width, height and angle of rectangle
  
! local variables  

  real(dp) x1,y1
  real(dp) x2,y2
  real(dp) x3,y3
  real(dp) x4,y4
  real(dp) xt,yt
  
  integer :: loop
 
! START

! first point

  xt=w/2d0
  yt=h/2d0
  
! rotate then translate, note save the first point
  
  x1=x+xt*cos(theta)-yt*sin(theta)
  y1=y+xt*sin(theta)+yt*cos(theta)
  
! second point

  xt=-w/2d0
  yt=h/2d0
! rotate then translate 
  x2=x+xt*cos(theta)-yt*sin(theta)
  y2=y+xt*sin(theta)+yt*cos(theta)
  
! third point

  xt=-w/2d0
  yt=-h/2d0
! rotate then translate 
  x3=x+xt*cos(theta)-yt*sin(theta)
  y3=y+xt*sin(theta)+yt*cos(theta)
  
! fourth point

  xt=w/2d0
  yt=-h/2d0
! rotate then translate 
  x4=x+xt*cos(theta)-yt*sin(theta)
  y4=y+xt*sin(theta)+yt*cos(theta)
  
  do loop=1,2  ! first pass to count the points, second pass to set the point coordinates
    
    npts=0
    
    CALL generate_line_points(npts,shape_x,shape_y,loop,x1,y1,x2,y2)
  
    CALL generate_line_points(npts,shape_x,shape_y,loop,x2,y2,x3,y3)
  
    CALL generate_line_points(npts,shape_x,shape_y,loop,x3,y3,x4,y4)
  
    CALL generate_line_points(npts,shape_x,shape_y,loop,x4,y4,x1,y1)
  
    if (loop.EQ.1) then
      ALLOCATE( shape_x(1:npts) )
      ALLOCATE( shape_y(1:npts) )
    end if
  
  end do ! next loop
   
  RETURN
  
END SUBROUTINE generate_rectangle_points

!
! NAME
!     SUBROUTINE generate_Dshape_points
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
!
!
SUBROUTINE generate_Dshape_points(npts,shape_x,shape_y,x,y,w1,w2,s,r,theta)
   
USE type_specifications
USE constants

IMPLICIT NONE

  integer,intent(OUT)                :: npts
  real(dp),allocatable,intent(INOUT)  :: shape_x(:)
  real(dp),allocatable,intent(INOUT)  :: shape_y(:)

  real(dp),intent(IN)    :: x,y,w1,w2,s,r,theta  ! centre x and y coordinates, width1, width2, separation of wire rows, offset of D shapeand angle
  
! local variables  

  real(dp) x1,y1
  real(dp) x2,y2
  real(dp) x3,y3
  real(dp) x4,y4
  real(dp) x5,y5
  real(dp) x6,y6
  real(dp) x7,y7
  real(dp) x8,y8
  real(dp) x9,y9
  real(dp) x10,y10
  real(dp) x11,y11
  real(dp) x12,y12
  
  real(dp) xt,yt
  real(dp) vx,vy
  real(dp) voxl,voyl
  real(dp) norm
 
  integer :: loop
  
! START

! vector from top left conductor to bottem left conductor      
  vx=w1-w2
  vy=-2d0*s

! perpendicular vector to -x edge
  norm=sqrt(vx*vx+vy*vy)
  voxl=vy*r/norm
  voyl=-vx*r/norm

! POINT 1      ! top right
  xt=w1
  yt=s+r 
  x1=x+xt*cos(theta)-yt*sin(theta)
  y1=y+xt*sin(theta)+yt*cos(theta)
      
! POINT 2      ! top left
  xt=-w1
  yt=s+r
  x2=x+xt*cos(theta)-yt*sin(theta)
  y2=y+xt*sin(theta)+yt*cos(theta)
       
! POINT 3      ! top left centre
  xt=-w1
  yt=s
  x3=x+xt*cos(theta)-yt*sin(theta)
  y3=y+xt*sin(theta)+yt*cos(theta)
     
! POINT 4      ! top left edge
  xt=-w1+voxl
  yt=s+voyl
  x4=x+xt*cos(theta)-yt*sin(theta)
  y4=y+xt*sin(theta)+yt*cos(theta)
        
! POINT 5      ! bottom left edge
  xt=-w2+voxl
  yt=-s+voyl
  x5=x+xt*cos(theta)-yt*sin(theta)
  y5=y+xt*sin(theta)+yt*cos(theta)  
      
! POINT 6      ! bottom left centre
  xt=-w2
  yt=-s
  x6=x+xt*cos(theta)-yt*sin(theta)
  y6=y+xt*sin(theta)+yt*cos(theta)
       
! POINT 7      ! bottom left
  xt=-w2
  yt=-(s+r)
  x7=x+xt*cos(theta)-yt*sin(theta)
  y7=y+xt*sin(theta)+yt*cos(theta)
         
! POINT 8      ! bottom right
  xt=w2
  yt=-(s+r)
  x8=x+xt*cos(theta)-yt*sin(theta)
  y8=y+xt*sin(theta)+yt*cos(theta)
      
! POINT 9      ! bottom right centre
  xt=w2
  yt=-s
  x9=x+xt*cos(theta)-yt*sin(theta)
  y9=y+xt*sin(theta)+yt*cos(theta)
      
! POINT 10      ! bottom right edge
  xt=w2-voxl
  yt=-s+voyl
  x10=x+xt*cos(theta)-yt*sin(theta)
  y10=y+xt*sin(theta)+yt*cos(theta)
      
! POINT 11      ! top right edge
  xt=w1-voxl
  yt=s+voyl
  x11=x+xt*cos(theta)-yt*sin(theta)
  y11=y+xt*sin(theta)+yt*cos(theta)
       
! POINT 12      ! top right centre
  xt=w1
  yt=s
  x12=x+xt*cos(theta)-yt*sin(theta)
  y12=y+xt*sin(theta)+yt*cos(theta)
     
  do loop=1,2  ! first pass to count the points, second pass to set the point coordinates
    
    npts=0
    
    CALL generate_line_points(npts,shape_x,shape_y,loop,x1,y1,x2,y2)
    
    CALL generate_arc_points(npts,shape_x,shape_y,loop,x3,y3,x2,y2,x4,y4)

    CALL generate_line_points(npts,shape_x,shape_y,loop,x4,y4,x5,y5)
  
    CALL generate_arc_points(npts,shape_x,shape_y,loop,x6,y6,x5,y5,x7,y7)

    CALL generate_line_points(npts,shape_x,shape_y,loop,x7,y7,x8,y8)

    CALL generate_arc_points(npts,shape_x,shape_y,loop,x9,y9,x8,y8,x10,y10)

    CALL generate_line_points(npts,shape_x,shape_y,loop,x10,y10,x11,y11)

    CALL generate_arc_points(npts,shape_x,shape_y,loop,x12,y12,x11,y11,x1,y1)
   
    if (loop.EQ.1) then
      ALLOCATE( shape_x(1:npts) )
      ALLOCATE( shape_y(1:npts) )
    end if
  
  end do ! next loop

  RETURN
  
END SUBROUTINE generate_Dshape_points
!
! NAME
!     SUBROUTINE generate_arc_points
!
! DESCRIPTION
!     generate a set of points describing an arc
!     The arc is always in an anticlockwise direction.
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 20/4/2017 CJS
!
!
SUBROUTINE generate_arc_points(npts,shape_x,shape_y,loop,xc,yc,x1,y1,x2,y2)

USE type_specifications
USE constants

IMPLICIT NONE

  integer,intent(OUT)                :: npts       ! point count
  real(dp),allocatable,intent(INOUT)  :: shape_x(:) ! x coordinate value list
  real(dp),allocatable,intent(INOUT)  :: shape_y(:) ! y coordinate value list
  
  integer,intent(IN)    :: loop                  ! loop indicator. If loop=1, just count the points, if loop=2 set the coordinate values

  real(dp),intent(IN)    :: xc,yc                ! centre coordinates
  real(dp),intent(IN)    :: x1,y1                ! arc end 1 coordinates
  real(dp),intent(IN)    :: x2,y2                ! arc end 2 coordinates
   
! local variables

  real(dp) :: r           ! arc radius
  real(dp) :: x,y         ! point coordinates
  real(dp) :: tmin        ! minimum angle
  real(dp) :: tmax        ! maximum angle
  real(dp) :: t           ! angle
  real(dp) :: dt          ! angle step
  integer  :: tloop       ! theta loop variable
  
  integer,parameter :: nt=40  ! number of points in an arc
   
! START
! calculate the radius
  r=sqrt((x1-xc)**2+(y1-yc)**2)

! calculate the angle of the first point from the x axis
  tmin=atan2(y1-yc,x1-xc)

! calculate the angle of the last point from the x axis
  tmax=atan2(y2-yc,x2-xc)

! ensure that tmax is greater than tmin i.e. the arc is in an anticlockwise direction  
  if (tmax.LT.tmin) tmax=tmax+2d0*pi
  
  dt=(tmax-tmin)/dble(nt)

! first point of the arc

  npts=npts+1
  if (loop.EQ.2) then
    shape_x(npts)=x1
    shape_y(npts)=y1
  end if

! write the interior points of the arc

! loop over theta    

  do tloop=1,nt-1
  
    t=tmin+dt*dble(tloop)

    x=xc+r*cos(t)
    y=yc+r*sin(t)
  
    npts=npts+1
    if (loop.EQ.2) then
      shape_x(npts)=x
      shape_y(npts)=y
    end if
  
  end do

  npts=npts+1
  if (loop.EQ.2) then
    shape_x(npts)=x2
    shape_y(npts)=y2
  end if
  
  RETURN

END SUBROUTINE generate_arc_points
!
! NAME
!     SUBROUTINE generate_line_points
!
! DESCRIPTION
!     generate a set of points describing a line
!  
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 20/4/2017 CJS
!
!
SUBROUTINE generate_line_points(npts,shape_x,shape_y,loop,x1,y1,x2,y2)

USE type_specifications
USE constants

IMPLICIT NONE

  integer,intent(OUT)                :: npts       ! point count
  real(dp),allocatable,intent(INOUT)  :: shape_x(:) ! x coordinate value list
  real(dp),allocatable,intent(INOUT)  :: shape_y(:) ! y coordinate value list
  
  integer,intent(IN)    :: loop                  ! loop indicator. If loop=1, just count the points, if loop=2 set the coordinate values

  real(dp),intent(IN)    :: x1,y1                ! line end 1 coordinates
  real(dp),intent(IN)    :: x2,y2                ! line end 2 coordinates
   
! local variables

  real(dp) :: x,y         ! point coordinates
  real(dp) :: t           ! normalised distance along line
  integer  :: tloop       ! loop variable
  
  integer,parameter :: nt=16  ! number of points in a line
   
! START

! first point of the line

  npts=npts+1
  if (loop.EQ.2) then
    shape_x(npts)=x1
    shape_y(npts)=y1
  end if

! write the interior points of the line

! loop over theta    

  do tloop=1,nt-1
  
    t=dble(tloop)/dble(nt)

    x=x1+(x2-x1)*t
    y=y1+(y2-y1)*t
  
    npts=npts+1
    if (loop.EQ.2) then
      shape_x(npts)=x
      shape_y(npts)=y
    end if
  
  end do

  npts=npts+1
  if (loop.EQ.2) then
    shape_x(npts)=x2
    shape_y(npts)=y2
  end if
  
  RETURN

END SUBROUTINE generate_line_points
