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
!SUBROUTINE add_integer_to_string
!
! NAME
!     SUBROUTINE add_integer_to_string
!
! DESCRIPTION
!     
!  add an integer to string1 and return the result in string2
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 8/08/2012 CJS
!     revised notation 4/12/2016 CJS
!
  SUBROUTINE add_integer_to_string(string1,n,string2)

  IMPLICIT NONE

! variables passed to the subroutine
    
  character(LEN=*),intent(IN)  ::  string1   ! input string
  integer,intent(IN)           ::  n         ! input number to be added to the input string
  character(LEN=*),intent(OUT) ::  string2   ! output string
    
! local variables    

  character(LEN=20) :: stringn
  integer ::  name_len

  integer ::  n_char_n,i,string_pos
  integer ::  digit,num,n_digits
  integer ::  mask

! START

  name_len=len(string1)
  
  if (n.LT.0) then
  
    run_status='ERROR: in add_integer_to_string: n<0'
    CALL write_program_status()
    STOP 1
     
  else if (n.eq.0) then
  
! special case if n=0
    n_digits=1
    stringn(1:1)='0'
    
  else   
     
! add a positive integer to a string

    num=n
    n_digits=int(log10(dble(num)))+1  ! calculate the number of digits in the integer
  
    mask=1
  
! loop over digits, in the order units, tens, hundreds...
    do i=1,n_digits

      digit=mod(num,mask*10)/mask                           ! work out the current digit value
      string_pos=n_digits-i+1                               ! digit position in string
      stringn(string_pos:string_pos)=char(digit+ichar('0')) ! convert digit to character and add to string
      num=num-mask*digit                                    ! subtract this digit from the number
      mask=mask*10                      
    
    end do ! next digit
    
  end if
  
  string2=trim(string1)//stringn(1:n_digits)                ! append the number string to string 1
  
  RETURN
  
  END SUBROUTINE add_integer_to_string
