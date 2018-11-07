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
!SUBROUTINE write_license
!SUBROUTINE path_format
!FUNCTION path_exists
!SUBROUTINE check_and_make_path
!SUBROUTINE write_long_node_list
!
! NAME
!    write_license
!
! DESCRIPTION
!     writes the license agreement note 
!
! HISTORY
!
!     started 3/05/13 CJS
!
! COMMENTS
!     
SUBROUTINE write_license()

IMPLICIT NONE 

! variables passed to subroutine
  
! local variables
  
! START

  write(*,*)''
  write(*,*)'!'
  write(*,*)'! Software to generate Spice sub-circuit models to reproduce'
  write(*,*)'! impedance functions specified as either s-domain rational functions'
  write(*,*)'! or pole-residue representations. '
  write(*,*)'!'
  write(*,*)'! Copyright (C) 2018 University of Nottingham'
  write(*,*)'!'
  write(*,*)'! NETWORK_SYNTHESIS is free software: you can redistribute it and/or modify it under the '
  write(*,*)'! terms of the GNU General Public License as published by the Free Software '
  write(*,*)'! Foundation, either version 3 of the License, or (at your option) any later '
  write(*,*)'! version.'
  write(*,*)'! '
  write(*,*)'! NETWORK_SYNTHESIS is distributed in the hope that it will be useful, but '
  write(*,*)'! WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY '
  write(*,*)'! or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License '
  write(*,*)'! for more details.'
  write(*,*)'! '
  write(*,*)'! A copy of the GNU General Public License version 3 can be found in the '
  write(*,*)'! file COPYING.txt in the root or at <http://www.gnu.org/licenses/>.'
  write(*,*)'! '
  write(*,*)'! NETWORK_SYNTHESIS uses the EISPACK library. EISPACK is subject to '
  write(*,*)'! the GNU Lesser General Public License. A copy of the GNU Lesser General Public '
  write(*,*)'! License version can be found in the file COPPYING.LESSER.txt '
  write(*,*)'! or at <http://www.gnu.org/licenses/>.'
  write(*,*)'! '
  write(*,*)'! The University of Nottingham can be contacted at: ggiemr@nottingham.ac.uk'
  write(*,*)'' 
  
  RETURN
  
END SUBROUTINE write_license
!
! NAME
!    path_format
!
! DESCRIPTION
!     check the specified path format- it should end with a /
!     if not, put one on
!
! HISTORY
!
!     started 29/1/2016 CJS
!
! COMMENTS
!     
SUBROUTINE path_format(path)

USE type_specifications

IMPLICIT NONE 

! variables passed to subroutine

character(len=filename_length),intent(INOUT) :: path
  
! local variables

integer :: length
  
! START

  length=LEN_TRIM(path)
  
! Note the different forms for the directory separator in unix and windows
! The file_separator is defined in general_module.F90 for both operating systems
! and is set appropriately using conditional compilation
  
  if (path(length:length).NE.file_separator) then
    path=trim(path)//file_separator
  end if
  
  RETURN
  
END SUBROUTINE path_format
!
! NAME
!    strip_path
!
! DESCRIPTION
!     strip the path from a name based on finding the last file separator character
!     and splitting the inpput string there
!
! HISTORY
!
!     started 24/2/2017 CJS
!
! COMMENTS
!     
!
SUBROUTINE strip_path(ipstring,path,name)

USE type_specifications

IMPLICIT NONE 

! variables passed to subroutine

character(len=filename_length),intent(IN) :: ipstring
character(len=filename_length),intent(OUT) :: path
character(len=filename_length),intent(OUT) :: name
 
! local variables

integer :: length
integer :: path_length
integer :: i
  
! START

  length=LEN_TRIM(ipstring)
  
! Note the different forms for the directory separator in unix and windows
! The file_separator is defined in general_module.F90 for both operating systems
! and is set appropriately using conditional compilation

  path_length=0
  do i=length,1,-1
   
    if (ipstring(i:i).EQ.file_separator) then
      path_length=i
      EXIT
    end if
   
  end do
  
  if (path_length.NE.0) then
    path=ipstring(1:path_length)
    name=ipstring(path_length+1:length)
  else
    path=""
    name=ipstring
  end if 
    
  RETURN
  
END SUBROUTINE strip_path
!
! NAME
!    path_exists
!
! DESCRIPTION
!     check that a given path exists
!     if not, exit with an error
!
! HISTORY
!
!     started 15/1/2016 CJS
!
! COMMENTS
!     
FUNCTION path_exists(path)

USE type_specifications

IMPLICIT NONE 

logical :: path_exists

! variables passed to subroutine

character(len=filename_length),intent(IN)  :: path
  
! local variables

integer :: ierr
  
! START

! Try to create a file in the directory. If it works then assume that the directory exists
  OPEN(unit=temp_file_unit,file=trim(path)//'temp',iostat=ierr)

! Test for success
  path_exists = (ierr == 0)

! Close and delete the temporary file
  if (ierr .EQ. 0) CLOSE(unit=temp_file_unit,status='delete')
  
  RETURN
  
END FUNCTION path_exists
!
! NAME
!    check_and_make_path
!
! DESCRIPTION
!     check that a given path exists
!     if not, make the appropriate directories
!
! HISTORY
!
!     started 15/1/2016 CJS
!
! COMMENTS
! Uses mkdir -p which can build the whole path as opposed to mkdir which can only build the bottom level directory
! mkdir -p works in ubuntu 14.04 - not sure that this is a portable solution though...
! 
SUBROUTINE check_and_make_path(path)

USE type_specifications

IMPLICIT NONE 

! variables passed to subroutine

character(len=filename_length),intent(IN)  :: path
  
! local variables

character(len=line_length)      :: command
  
! START

! check whether the path already exists, if so all is OK and we can return
  if (path_exists(path)) RETURN
  
! the path doesn't exist so we must create it
! Note there are different forms for the system command on unix and windows
! The mkdir_command is defined in general_module.F90 for both operating systems

  command=mkdir_command//trim(path)

  CALL EXECUTE_COMMAND_LINE(command)
 
  RETURN
  
END SUBROUTINE check_and_make_path
!
! NAME
!    write_long_node_list
!
! DESCRIPTION
!      write a long list of nodes such that each line doesn't exceed the specified length
!      This is used so that the maximum number of characters in a spice input file is not exceeded
!
! HISTORY
!
!     started 10/10/2016 CJS
!
! COMMENTS
! 
! 
SUBROUTINE write_long_node_list(n_nodes,node_list,max_length,unit)

USE type_specifications

IMPLICIT NONE 

! variables passed to subroutine

integer,intent(IN) :: n_nodes
integer,intent(IN) :: node_list(n_nodes)
integer,intent(IN) :: max_length
integer,intent(IN) :: unit
  
! local variables

integer :: first_node,last_node,numbers_per_line,number_width
integer :: i
  
! START

  number_width=6                                  ! unmber of characters to write an integer. Must be consistent with format in write below
  numbers_per_line=(max_length-1)/number_width    ! -1 for continuation character
  
! work out the first and last node for the first line    
  first_node=1
  last_node=min(n_nodes,numbers_per_line)
  
  do while (last_node.GE.first_node)
  
! write the line of the number list
    write(unit,'(A)',ADVANCE='NO')'+'
    write(unit,'(1000I6)')(node_list(i),i=first_node,last_node)

! work out the first and last node for the next line    
    first_node=last_node+1
    last_node=min(n_nodes,first_node+numbers_per_line-1)

  end do

  RETURN
  
END SUBROUTINE write_long_node_list
