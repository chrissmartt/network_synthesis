!
! NAME
!     RLRCR_test
!
! DESCRIPTION
!       look for a viable RLRCR branch in a given impedance/admittance function
!
! SEE ALSO
!
!
! HISTORY
!
!     started 14/09/17 CJS
!
  SUBROUTINE RLRCR_test(H_PR,type,CFtype,R1,L,R2,C,R3,found,HR,remainder_OK,remainder_zero)

USE type_specifications
USE general_module
USE constants
USE filter_module

IMPLICIT NONE
  
  type(Sfilter_PR),INTENT(IN) :: H_PR
  integer :: type
  integer :: CFtype
  real(dp):: R1,L,R2,C,R3
  logical :: found
  type(Sfilter),INTENT(INOUT) :: HR
  logical :: remainder_OK,remainder_zero

! local variables
  
  integer :: first_complex_pole,pole1,pole2,pole
  type(Sfilter_PR) :: HR_PR_local
  logical :: stable
  
  integer :: i,ii
  
  real(dp):: const_term
  logical :: complex_poles,positive_residues,positive_const_term
  logical :: positive_RLRCR
  
! function types
  logical :: conjugate_pair
  logical :: imaginary_pair
  logical :: complex_pair

!START

  write(*,*)'CALLED: RLRCR_test'

  found=.FALSE.

  first_complex_pole=H_PR%n_real_poles+1
  
! loop over complex pole pairs
  do i=1,H_PR%n_complex_pole_pairs
  
    pole1=first_complex_pole+(i-1)*2
    pole2=pole1+1
    
! check that we have a conjugate pole pair
    if ( .NOT.conjugate_pair(H_PR%poles(pole1),H_PR%poles(pole2)) ) then
      write(*,*)'ERROR in RLRCR_test: poles are not a conjugate pair'
      write(*,*)'pole 1:',H_PR%poles(pole1)
      write(*,*)'pole 2:',H_PR%poles(pole2)
      STOP
    end if
    
! check that we have a conjugate residue pair
    if ( .NOT.conjugate_pair(H_PR%residues(pole1),H_PR%residues(pole2)) ) then
      write(*,*)'ERROR in RLRCR_test: residues are not a conjugate pair'
      write(*,*)'residues 1:',H_PR%residues(pole1)
      write(*,*)'residues 2:',H_PR%residues(pole2)
      STOP
    end if

! test for whether we have a RLRCR branch here...
    
    const_term       = H_PR%R
    positive_const_term=(const_term.GT.0d0)
    complex_poles    = complex_pair(H_PR%poles(pole1),H_PR%poles(pole2))
    positive_residues= (dble(H_PR%residues(pole1)).GT.0d0)

! From the pole and residue values and the constant term,
! Calculate the values of R1, L, R2, C and R3 for the RLRCR branch

! Check whether these values are positive i.e. do we have a viable branch in terms
! of component values   
    
    if (verbose) then
      write(*,*)'Testing pole pair',i
      write(*,*)'Tests:'
      write(*,*)'positive costant term test  : ',positive_const_term,' ',const_term
      write(*,*)'complex pole test           : ',complex_poles,' ',H_PR%poles(pole1)
      write(*,*)'positive residue test       : ',positive_residues,' ',H_PR%residues(pole1)
    end if
    
! **** QUESTION: DO THE RESIDUES HAVE TO BE POSITIVE IN THIS CASE?****
! check for stable poles which are not on the imaginary s=jw axis AND positive constant term

    if (complex_poles.AND.positive_const_term) then

! Calculate the values of the 5 components of the branch R1, L, R2, C, R3 
! before checking that they are positive

! Calculate the admittance function, Y of the RLRCR branch

      positive_RLRCR=.FALSE.
      
    end if
    
    if (positive_RLRCR) then
    
      if (verbose) write(*,*)'Found possible RLRCR branch'

! this could be a viable RLRCR branch - calculate the remainder when this pole pair is removed

      CALL deallocate_Sfilter(HR)
      
! Test whether the remainder is zero
      
! build a local pole-residue filter without the test pole pair. 

! allocate the structure for the local pole-residue form function and copy the 
! required information across            
      HR_PR_local%wnorm=H_PR%wnorm
      HR_PR_local%order=H_PR%order-2
      HR_PR_local%n_real_poles=H_PR%n_real_poles
      HR_PR_local%n_complex_poles=H_PR%n_complex_poles-2
      HR_PR_local%n_complex_pole_pairs=H_PR%n_complex_pole_pairs-1
      HR_PR_local%n_real_poles=H_PR%n_real_poles

! constant term and sL term
      
      HR_PR_local%R=H_PR%R
      HR_PR_local%L=H_PR%L
      
! Test whether the remainder is zero
      if ( (HR_PR_local%order.EQ.0).AND.(HR_PR_local%R.EQ.0d0).AND.(HR_PR_local%L.EQ.0d0) ) then
        remainder_zero=.TRUE.
        GOTO 8000    
      end if

! copy any poles/ residues in the remainder
      
      ALLOCATE( HR_PR_local%complex_pole(HR_PR_local%order) )
      ALLOCATE( HR_PR_local%poles(HR_PR_local%order) )
      ALLOCATE( HR_PR_local%residues(HR_PR_local%order) )

! copy real poles      
      pole1=0
      pole2=0
      do ii=1,HR_PR_local%n_real_poles
        pole1=pole1+1
        pole2=pole2+1
        HR_PR_local%complex_pole(pole1)=.FALSE.
        HR_PR_local%poles(pole1)   =H_PR%poles(pole2)
        HR_PR_local%residues(pole1)=H_PR%residues(pole2)
      end do ! next real pole

! copy complex poles in pairs 
      do ii=1,H_PR%n_complex_pole_pairs
      
        if (ii.NE.i) then
! this is not the pole pair we are trying to remove so copy them over
          pole1=pole1+1
          pole2=pole2+1
          HR_PR_local%complex_pole(pole1)=.TRUE.
          HR_PR_local%poles(pole1)=H_PR%poles(pole2)
          HR_PR_local%residues(pole1)=H_PR%residues(pole2)
          pole1=pole1+1
          pole2=pole2+1
          HR_PR_local%complex_pole(pole1)=.TRUE.
          HR_PR_local%poles(pole1)=H_PR%poles(pole2)
          HR_PR_local%residues(pole1)=H_PR%residues(pole2)
        else
! this is the pair of poles we wish to remove so just increase the pole2 counter by 2. 
          pole2=pole2+2 
        end if
        
      end do ! next real pole

! convert to a rational function form      
      HR=Convert_filter_S_PR_to_S(HR_PR_local)

! Check the transfer funcion for stability and for whether it is positive real

      CALL check_transfer_function(HR,stable)
      
      CALL deallocate_Sfilter_PR(HR_PR_local)
      
      if (verbose) then
        if (stable) then
          write(*,*)'Remainder is stable'
        else
          write(*,*)'Remainder is unstable'        
        end if
      end if  
      
      if (stable) then
        remainder_zero=.FALSE.
        GOTO 8000
      end if
    
! Test whether the remainder is positive real    
    
    end if ! positive residues for this complex pole pair
  
  end do ! next complex pole pair
  
! we only get here if we have not found a viable RLRCR branch

  remainder_OK=.FALSE.
  found=.FALSE.
  remainder_zero=.FALSE.
  
  RETURN
  
8000 CONTINUE
! jump here if we have found a viable RLRCR branch

  remainder_OK=.TRUE.
  found=.TRUE.
  CALL deallocate_Sfilter_PR(HR_PR_local)
  pole1=first_complex_pole+(i-1)*2
  pole2=pole1+1
  if (type.EQ.type_impedance) then
    CFtype=series_RLRCR
    C=0d0
    L=0d0
    R1=0d0
    R2=0d0
    R3=0d0
  else
    CFtype=shunt_RLRCR
    L=0d0
    C= 0d0
    R1=0d0
    R2=0d0
    R3=0d0
  end if
  if (verbose) then
    write(*,*)'FOUND VIABLE RLRCR BRANCH'
    write(*,*)'R1=',R1
    write(*,*)'L =',L
    write(*,*)'R2=',R2
    write(*,*)'C =',C
    write(*,*)'R3=',R3
    write(*,*)'remainder_OK   :',remainder_OK
    write(*,*)'remainder_zero :',remainder_zero
  end if
  
  RETURN

  END SUBROUTINE RLRCR_test
