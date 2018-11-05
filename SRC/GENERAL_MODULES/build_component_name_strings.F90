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
! File Contents:
!SUBROUTINE build_name_with_conductor_and_end
!SUBROUTINE build_name_with_conductor_domainConductor_and_end
!SUBROUTINE build_name_with_domain_conductor_and_end
!SUBROUTINE build_name_with_domain_mode_and_end
!SUBROUTINE build_name_with_domain_and_mode
!SUBROUTINE build_name_with_domain_conductor_mode_and_end
!
! NAME
!     SUBROUTINE build_name_with_domain_conductor_and_end
!
! DESCRIPTION
!     
!  build spice component name string with conductor number and end number added
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 24/03/2016 CJS
!
!
  SUBROUTINE build_name_with_conductor_and_end(name_head,conductor,end,name)
!
  USE type_specifications
  
  IMPLICIT NONE

! variables passed to the subroutine
    
  character(LEN=spice_name_length),intent(IN)  ::  name_head
  character(LEN=spice_name_length),intent(OUT) ::  name
  integer,intent(IN)                           ::  conductor
  integer,intent(IN)                           ::  end

! local variables

  character(LEN=spice_name_length) ::  temp_name
  character(LEN=spice_name_length) ::  temp_name2
  
! START

! add conductor number
  temp_name=trim(name_head)//'_c'
  CALL add_integer_to_string(temp_name,conductor,temp_name2)
  
! add end number
  temp_name=trim(temp_name2)//'_e'
  CALL add_integer_to_string(temp_name,end,name)

  END SUBROUTINE build_name_with_conductor_and_end
!
! NAME
!     SUBROUTINE build_name_with_conductor_domainConductor_and_end
!
! DESCRIPTION
!     
!  build spice component name string with domain number, conductor number and end number added
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 24/03/2016 CJS
!
!
  SUBROUTINE build_name_with_conductor_domainConductor_and_end(name_head,conductor,domainConductor,end,name)
!
  USE type_specifications
  
  IMPLICIT NONE

! variables passed to the subroutine
    
  character(LEN=spice_name_length),intent(IN)  ::  name_head
  integer,intent(IN)                           ::  conductor
  integer,intent(IN)                           ::  domainConductor
  integer,intent(IN)                           ::  end
  character(LEN=spice_name_length),intent(OUT) ::  name

! local variables

  character(LEN=spice_name_length) ::  temp_name
  character(LEN=spice_name_length) ::  temp_name2
  
! START

! add conductor number
  temp_name=trim(name_head)//'_c'
  CALL add_integer_to_string(temp_name,conductor,temp_name2)

! add domainConductor number
  temp_name=trim(temp_name2)//'_dc'
  CALL add_integer_to_string(temp_name,domainConductor,temp_name2)
  
! add end number
  temp_name=trim(temp_name2)//'_e'
  CALL add_integer_to_string(temp_name,end,name)

  END SUBROUTINE build_name_with_conductor_domainConductor_and_end
!
! NAME
!     SUBROUTINE build_name_with_domain_conductor_and_end
!
! DESCRIPTION
!     
!  build spice component name string with domain number, conductor number and end number added
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 4/02/2016 CJS
!
!
  SUBROUTINE build_name_with_domain_conductor_and_end(name_head,domain,conductor,end,name)
!
  USE type_specifications
  
  IMPLICIT NONE

! variables passed to the subroutine
    
  character(LEN=spice_name_length),intent(IN)  ::  name_head
  integer,intent(IN)                           ::  domain
  integer,intent(IN)                           ::  conductor
  integer,intent(IN)                           ::  end
  character(LEN=spice_name_length),intent(OUT) ::  name

! local variables

  character(LEN=spice_name_length) ::  temp_name
  character(LEN=spice_name_length) ::  temp_name2
  
! START

! add domain number
  temp_name=trim(name_head)//'_d'
  CALL add_integer_to_string(temp_name,domain,temp_name2)

! add coonductor number
  temp_name=trim(temp_name2)//'_c'
  CALL add_integer_to_string(temp_name,conductor,temp_name2)
  
! add end number
  temp_name=trim(temp_name2)//'_e'
  CALL add_integer_to_string(temp_name,end,name)

  END SUBROUTINE build_name_with_domain_conductor_and_end
!
! NAME
!     SUBROUTINE build_name_with_domain_mode_and_end
!
! DESCRIPTION
!     
!  build spice component name string with domain number, mode number and end number added
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 4/02/2016 CJS
!
!
  SUBROUTINE build_name_with_domain_mode_and_end(name_head,domain,mode,end,name)
!
  USE type_specifications
  
  IMPLICIT NONE

! variables passed to the subroutine
    
  character(LEN=spice_name_length),intent(IN)  ::  name_head
  integer,intent(IN)                           ::  domain
  integer,intent(IN)                           ::  mode
  integer,intent(IN)                           ::  end
  character(LEN=spice_name_length),intent(OUT) ::  name

! local variables

  character(LEN=spice_name_length) ::  temp_name
  character(LEN=spice_name_length) ::  temp_name2
  
! START

! add domain number
  temp_name=trim(name_head)//'_d'
  CALL add_integer_to_string(temp_name,domain,temp_name2)

! add mode number
  temp_name=trim(temp_name2)//'_m'
  CALL add_integer_to_string(temp_name,mode,temp_name2)
  
! add end number
  temp_name=trim(temp_name2)//'_e'
  CALL add_integer_to_string(temp_name,end,name)

  END SUBROUTINE build_name_with_domain_mode_and_end  
!
! NAME
!     SUBROUTINE build_name_with_domain_and_mode
!
! DESCRIPTION
!     
!  build spice component name string with domain number and mode number added
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 4/02/2016 CJS
!
!
  SUBROUTINE build_name_with_domain_and_mode(name_head,domain,mode,name)
!
  USE type_specifications
  
  IMPLICIT NONE

! variables passed to the subroutine
    
  character(LEN=spice_name_length),intent(IN)  ::  name_head
  integer,intent(IN)                           ::  domain
  integer,intent(IN)                           ::  mode
  character(LEN=spice_name_length),intent(OUT) ::  name

! local variables

  character(LEN=spice_name_length) ::  temp_name
  character(LEN=spice_name_length) ::  temp_name2
  
! START

! add domain number
  temp_name=trim(name_head)//'_d'
  CALL add_integer_to_string(temp_name,domain,temp_name2)

! add mode number
  temp_name=trim(temp_name2)//'_m'
  CALL add_integer_to_string(temp_name,mode,name)
  
  END SUBROUTINE build_name_with_domain_and_mode

!
! NAME
!     SUBROUTINE build_name_with_domain_conductor_mode_and_end
!
! DESCRIPTION
!     
!  build spice component name string with domain number, conductor number,
!                                         mode number and end number added
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 4/02/2016 CJS
!
!
  SUBROUTINE build_name_with_domain_conductor_mode_and_end(name_head,domain,conductor,mode,end,name)
!
  USE type_specifications
  
  IMPLICIT NONE

! variables passed to the subroutine
    
  character(LEN=spice_name_length),intent(IN)  ::  name_head
  integer,intent(IN)                           ::  domain
  integer,intent(IN)                           ::  conductor
  integer,intent(IN)                           ::  mode
  integer,intent(IN)                           ::  end
  character(LEN=spice_name_length),intent(OUT) ::  name

! local variables

  character(LEN=spice_name_length) ::  temp_name
  character(LEN=spice_name_length) ::  temp_name2
  
! START

! add domain number
  temp_name=trim(name_head)//'_d'
  CALL add_integer_to_string(temp_name,domain,temp_name2)

! add coonductor number
  temp_name=trim(temp_name2)//'_c'
  CALL add_integer_to_string(temp_name,conductor,temp_name2)

! add mode number
  temp_name=trim(temp_name2)//'_m'
  CALL add_integer_to_string(temp_name,mode,temp_name2)
  
! add end number
  temp_name=trim(temp_name2)//'_e'
  CALL add_integer_to_string(temp_name,end,name)

  END SUBROUTINE build_name_with_domain_conductor_mode_and_end  
  
  
  
