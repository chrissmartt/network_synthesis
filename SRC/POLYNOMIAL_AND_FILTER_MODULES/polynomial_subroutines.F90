!
! NAME
!      deallocate_poly
!      
! DESCRIPTION
!      deallocate the polynomial structure coefficients
!     
!
! SEE ALSO
!
!
! HISTORY
!
!     started 31/8/17 CJS
!
  SUBROUTINE deallocate_poly(H)

USE type_specifications
USE general_module
USE constants
  
  type(polynomial),intent(INOUT)  :: H

! local variables  

!START
       
  if(allocated( H%coeff )) DEALLOCATE( H%coeff )
   
  RETURN
  END SUBROUTINE deallocate_poly
!      
! DESCRIPTION
!      deallocate the polynomial structure coefficients
!     
!
! SEE ALSO
!
!
! HISTORY
!
!     started 31/8/17 CJS
!
  SUBROUTINE deallocate_complex_poly(H)

USE type_specifications
USE general_module
USE constants
  
  type(complex_polynomial),intent(INOUT)  :: H

! local variables  

!START
       
  if(allocated( H%coeff )) DEALLOCATE( H%coeff )
   
  RETURN
  END SUBROUTINE deallocate_complex_poly
!
! NAME
!      get_min_order_poly
!      
! DESCRIPTION
!      return the minimum order polynomial i.e. remove any high order zero coefficients
!     
!
! SEE ALSO
!
!
! HISTORY
!
!     started 31/8/17 CJS
!
  SUBROUTINE get_min_order_poly(H)

USE type_specifications
USE general_module
USE constants
  
  type(polynomial),intent(INOUT)  :: H

! local variables  

  integer :: on
  type(polynomial)  :: Hlocal
  real(dp)          :: max_coeff
  integer :: i

!START

! copy H
  Hlocal=H

! Get the maximum coefficient value
  max_coeff=0d0
  do i=0,Hlocal%order
    max_coeff=max (abs(Hlocal%coeff(i)),max_coeff)
  end do
       
! check the true order of the function
  on=0
  
  if (max_coeff.NE.0d0) then
    do i=0,Hlocal%order
      if (abs(Hlocal%coeff(i)/max_coeff).GT.zero_test_small) on=i
    end do
  end if
  
  if (on.NE.Hlocal%order) then
  
    DEALLOCATE( H%coeff )
    H%order=on
    ALLOCATE( H%coeff(0:on) )
    do i=0,on
      H%coeff(i)=Hlocal%coeff(i)
    end do
  
  end if
       
  CALL deallocate_poly(Hlocal)
   
  RETURN
  END SUBROUTINE get_min_order_poly
!
! NAME
!     write_polyrat_local
!
! DESCRIPTION
!       write Laplace domain polynomial rational function
!
! SEE ALSO
!
!
! HISTORY
!
!     started 23/01/09 CJS
!
  SUBROUTINE write_polyrat_local(a,b)


USE type_specifications
USE general_module
USE constants

IMPLICIT NONE
  
  type(Polynomial),intent(IN):: a,b
  
  integer i

!START

! write to screen
  write(*,*)'Laplace domain filter'
  write(*,*)'a order=',a%order
  write(*,*)'b order=',b%order
  write(*,*)''
  
  write(*,'(A)',advance='NO')'F(s) = '
  
  do i=a%order,1,-1
    write(*,8000,advance='NO')a%coeff(i),' s^',i,' + '
  end do
  write(*,8010)a%coeff(0)
  
  write(*,'(A)',advance='NO')'       '
  do i=a%order,1,-1
    write(*,'(A)',advance='NO')'------------'
  end do
  write(*,'(A)')'-----'
  
  write(*,'(A)',advance='NO')'       '
  do i=b%order,1,-1
    write(*,8000,advance='NO')b%coeff(i),' s^',i,' + '
  end do
  write(*,8010)b%coeff(0)
      
8000 format(F5.2,A3,I1,A3)
8010 format(F5.2)
  
  RETURN
  END SUBROUTINE write_polyrat_local
!
! NAME
!     write_poly_local
!
! DESCRIPTION
!       write Laplace domain polynomial
!
! SEE ALSO
!
!
! HISTORY
!
!     started 23/01/09 CJS
!
SUBROUTINE write_poly_local(a)

USE type_specifications
USE general_module
USE constants
  
IMPLICIT NONE

  type(Polynomial),intent(IN):: a
  
  integer i

!START

  if (.NOT.allocated(a%coeff)) then
    write(*,*)'Polynomial not allocated'
    RETURN
  end if
  
! write to screen
  write(*,*)'Polynomial'
  write(*,*)'a order=',a%order
  write(*,*)''
  
  write(*,'(A)',advance='NO')'F(s) = '
  
  do i=a%order,1,-1
    write(*,8000,advance='NO')a%coeff(i),' s^',i,' + '
  end do
  write(*,8010)a%coeff(0)
      
8000 format(F6.2,A3,I1,A3)
8010 format(F6.2)
  
RETURN
END SUBROUTINE write_poly_local
!
! NAME
!     write_complex_poly_local2
!
! DESCRIPTION
!       write complex polynomial
!
! SEE ALSO
!
!
! HISTORY
!
!     started 23/01/09 CJS
!
SUBROUTINE write_complex_poly_local2(f,x,a)

USE type_specifications
USE general_module
USE constants
  
IMPLICIT NONE

  character,intent(IN) :: f
  character,intent(IN) :: x
  type(complex_Polynomial),intent(IN):: a
  
  integer i

!START

  if (.NOT.allocated(a%coeff)) then
    write(*,*)'Polynomial not allocated'
    RETURN
  end if
  
! write to screen
  write(*,*)'Polynomial'
  write(*,*)'a order=',a%order
  write(*,*)''
  
  write(*,'(A,A,A,A4)',advance='NO')f,'(',x,') = '
    
  do i=a%order,1,-1
    if (i.NE.1) then
      write(*,8000,advance='NO')'(',real(a%coeff(i)),',',aimag(a%coeff(i)),')',' ',x,'^',i,' + '
    else
      write(*,8005,advance='NO')'(',real(a%coeff(i)),',',aimag(a%coeff(i)),')',' ',x,' + '
    end if
  end do
  write(*,8010)'(',real(a%coeff(0)),',',aimag(a%coeff(0)),')'
      
8000 format(A1,F6.2,A1,F6.2,A1,A,A,A,I1,A3)
8005 format(A1,F6.2,A1,F6.2,A1,A,A,A3)
8010 format(A1,F6.2,A1,F6.2,A1)
  
RETURN
END SUBROUTINE write_complex_poly_local2
!
! NAME
!     write_poly_local2
!
! DESCRIPTION
!       write Laplace domain polynomial. Different format to write_poly_local
!
! SEE ALSO
!
!
! HISTORY
!
!     started 6/09/2017 CJS
!
SUBROUTINE write_poly_local2(f,x,a)

USE type_specifications
USE general_module
USE constants
  
IMPLICIT NONE

  character,intent(IN) :: f
  character,intent(IN) :: x
  type(Polynomial),intent(IN):: a
  
  integer i

!START

  if (.NOT.allocated(a%coeff)) then
    write(*,*)'Polynomial not allocated'
    RETURN
  end if

! write to screen
  
  write(*,'(A,A,A,A4)',advance='NO')f,'(',x,') = '
  
  do i=a%order,1,-1
    if (i.NE.1) then
      write(*,8000,advance='NO')a%coeff(i),' ',x,'^',i,' + '
    else
      write(*,8005,advance='NO')a%coeff(i),' ',x,' + '
    end if
  end do
  write(*,8010)a%coeff(0)
      
8000 format(F6.2,A,A,A,I1,A3)
8005 format(F6.2,A,A,A3)
8010 format(F6.2)
  
RETURN
END SUBROUTINE write_poly_local2
!
! NAME
!     write_poly_local3
!
! DESCRIPTION
!       write Laplace domain polynomial
!
! SEE ALSO
!
!
! HISTORY
!
!     started 5/10/2017 CJS
!
SUBROUTINE write_poly_local3(a)

USE type_specifications
USE general_module
USE constants
  
IMPLICIT NONE

  type(Polynomial),intent(IN):: a
  
  integer i

!START

  if (.NOT.allocated(a%coeff)) then
    write(*,*)'Polynomial not allocated'
    RETURN
  end if
  
! write to screen
  write(*,*)'Polynomial'
  write(*,*)'a order=',a%order
  write(*,*)''
  
  write(*,'(A)',advance='NO')'F(s) = '
  
  do i=a%order,1,-1
    write(*,8000,advance='NO')a%coeff(i),' s^',i,' + '
  end do
  write(*,8010)a%coeff(0)
      
8000 format(ES10.3,A3,I1,A3)
8010 format(ES10.3)
  
RETURN
END SUBROUTINE write_poly_local3
  
  
!
! NAME
!     divide_poly
!
! DESCRIPTION
!     divide polynomials A(x)/B(x)
!     return the result, C(x) and remainder, R(x) where B(x)=A(x)C(x)+R(x)
!
! SEE ALSO
!
!
! HISTORY
!
!     started 04/03/09 CJS
!
  
SUBROUTINE divide_poly(Ain,Bin,Q,R,local_verbose)

USE type_specifications
USE general_module
USE constants

IMPLICIT NONE

! variables passed to subroutine  
  type(polynomial),intent(IN)      :: Ain
  type(polynomial),intent(IN)      :: Bin
  type(polynomial),intent(INOUT)   :: Q
  type(polynomial),intent(INOUT)   :: R
  logical :: local_verbose

! local variables

  type(polynomial)   :: A,B,T,Q2,R2,K
  integer :: torder

! START

  if (local_verbose) then
    write(*,*)'CALLED: divide_poly'
    CALL write_poly_local2('A','x',Ain)
    write(*,*)'B(s)='
    CALL write_poly_local2('B','x',Bin)
  end if

! deallocate result and remainder polynomials if required
  if (allocated(Q%coeff)) deallocate (Q%coeff)
  if (allocated(R%coeff)) deallocate (R%coeff)
  
! clean up the input polynomials - remove any zero high order coefficients
  A=Ain
  B=Bin
  
  CALL get_min_order_poly(A)
  CALL get_min_order_poly(B)
    
! check for zero denominator
  if ( polynomial_is_zero(B) ) then
    write(*,*)'ERROR in divide_poly: denominator =0'
    STOP  
  end if
      
! work out the order of the result and the remainder (assuming the general case)
  
  Q=0d0
  
  R=A
  
  do while ( (R%order.GE.B%order).AND.(.NOT.polynomial_is_zero(R)) )

    if (local_verbose) then
      write(*,*)'Calculate new term'
      CALL write_poly_local2('R','x',R)
      CALL write_poly_local2('B','x',B)
    end if

! get the contribution to Q  
    Torder=R%order-B%order
    T=allocate_polynomial(Torder)
    T%coeff(Torder)=R%coeff(R%order)/B%coeff(B%order)
    
    if (local_verbose) then
      write(*,*)'T='
      CALL write_poly_local2('T','x',T)
    end if
    
! add the contribution to Q. Q=Q+T
    Q2=Q+T
    Q=Q2
    CALL deallocate_poly(Q2)
    
    if (local_verbose) then
      write(*,*)'Q='
      CALL write_poly_local2('Q','x',Q)
    end if
    
! Calculate the next remainder. R=R-T*B
    K=T*B
    R2=R-K
    R=R2
    CALL get_min_order_poly(R)
    CALL deallocate_poly(K)
    CALL deallocate_poly(R2)
    if (local_verbose) then
      write(*,*)'R='
      CALL write_poly_local2('R','x',R)
      
      do torder=0,R%order
        write(*,*)R%coeff(torder)
      end do
    end if
    
  end do
  
  CALL deallocate_poly(A)
  CALL deallocate_poly(B)
  
  if (local_verbose) then
    write(*,*)'_____________________________________________________________'
    write(*,*)' '
    write(*,*)'FINISHED'
    write(*,*)' '
    CALL write_poly_local2('Q','x',Q)
    CALL write_poly_local2('R','x',R)
    write(*,*)'_____________________________________________________________'
  end if
  
  RETURN
  
END SUBROUTINE divide_poly
  
!
! NAME
!      remove_even_multiple_zeros
!      
! DESCRIPTION
!      look for the zeros of A(x) and remove any zeros with even multiplicity
!     
!
! SEE ALSO
!
!
! HISTORY
!
!     started 21/9/17 CJS
!
  SUBROUTINE remove_even_multiple_zeros(A,odd_multiplicity_flag,local_verbose)

USE type_specifications
USE general_module
USE constants
  
IMPLICIT NONE
  
  type(polynomial),intent(INOUT)  :: A
  logical,intent(OUT) :: odd_multiplicity_flag
  logical,intent(IN)  :: local_verbose

! local variables  

  integer :: i,loop
  integer      :: poly_order
  integer      :: final_poly_order
  complex(dp),allocatable  :: roots(:)
  complex(dp),allocatable  :: final_roots(:)
  integer,allocatable      :: root_check(:)
  type(complex_polynomial)  :: term
  type(complex_polynomial)  :: Alocal
  type(complex_polynomial)  :: Alocal2
  type(polynomial)  :: Ap
 
  logical :: real_root
  
  real(dp) ::gain
  
  integer :: root1,root2
  integer :: multiplicity
  
  integer,parameter :: root_under_test=1
  integer,parameter :: root_removed=2
  
! This a small number used in testing for double roots etc. It reflects the accuracy with which the 
! roots are found.

  real(dp),parameter :: eps=1D-6
  integer,parameter  :: n_newton=10
  complex(dp)        :: new_root

!START
  
  if (local_verbose) then
    write(*,*)'CALLED: remove_even_multiple_zeros'
    write(*,*)'Input polynomial:'
    CALL write_poly_local2('A','x',A)  
    write(*,*)'Full precision:'
    do i=0,A%order
      write(*,*)'coeff ',i,' ',A%coeff(i)
    end do
  end if
    
  gain=A%coeff(A%order)
  odd_multiplicity_flag=.FALSE.
  
! calculate the zeros of A(x)

  poly_order=A%order
  
  ALLOCATE( roots(1:poly_order) )
  
  ALLOCATE( root_check(1:poly_order) )
  root_check(1:poly_order)=0
  
  ALLOCATE( final_roots(1:poly_order) )
  final_roots(1:poly_order)=(0d0,0d0)
  
  CALL findroots(A,roots,poly_order)
  
  if (local_verbose) then
  
    write(*,*)'Gain=',gain
    do i=1,poly_order
      write(*,*)'zero ',i,' =',roots(i),evaluate_polynomial(A,roots(i))
    end do
        
  end if
  
! Apply Newton's method to 'polish' the roots

! Calculate the differential polynomial Ap

  Ap=allocate_polynomial(a%order-1)

  do i=0,ap%order
    ap%coeff(i)=a%coeff(i+1)*dble(i+1)
  end do

  do i=1,poly_order
    do loop=1,n_newton
      new_root=roots(i)-evaluate_polynomial(A,roots(i))/evaluate_polynomial(Ap,roots(i))
      roots(i)=new_root
    end do
  end do
  
  deallocate( AP%coeff )
  
  if (local_verbose) then
    write(*,*)'After root polishing:'
    write(*,*)'Gain=',gain
    do i=1,poly_order
      write(*,*)'zero ',i,' =',roots(i),evaluate_polynomial(A,roots(i))
    end do
        
  end if
  
! look for multiple real zeros
  
  do root1=1,poly_order-1
  
    real_root=( abs(aimag(roots(root1))).LT.eps )
    
    if (local_verbose) then
      write(*,*)'Checking root ',root1,' real root:',real_root
    end if
  
    if (real_root.AND.( root_check(root1).NE.root_removed) ) then
! this root has not been tested

      root_check(root1)=root_under_test
      multiplicity=1
      do root2=root1+1,poly_order
    
        if ( (root_check(root2).NE.root_removed).AND.(abs(roots(root1)-roots(root2)).LT.eps) ) then
! this is a repeated root so increase the multiplicity and flag the root as a root_under_test
          multiplicity=multiplicity+1
          root_check(root2)=root_under_test
          if (local_verbose) write(*,*)'Found repeated root:',root2
         
        end if
  
      end do ! next root 2
      
      if (local_verbose) write(*,*)'Root multiplicity is ',  multiplicity  
        
      if (multiplicity.GT.1)then
! we have multiple roots and some will need to be removed
        if (mod(multiplicity,2).EQ.1) then
! we must keep one root i.e. root1
          if (local_verbose) write(*,*)'we must keep one root i.e. root1=',root1
          odd_multiplicity_flag=.TRUE.
          root_check(root1)=0
        else
! we must remove the initial root
          if (local_verbose) write(*,*)'Removing repeated root:',root1
          root_check(root1)=root_removed
        end if
        
! remove pairs of roots
        do root2=root1,poly_order
      
          if ( (root_check(root2).EQ.root_under_test) ) then
! this is a repeated root so remove it
            root_check(root2)=root_removed
            if (local_verbose) then
              write(*,*)'Removing repeated root:',root2
            end if
          end if
  
        end do
        
      else if (multiplicity.EQ.1)then
      
        odd_multiplicity_flag=.TRUE.

      end if ! multiplicity.GT.1
            
    end if ! this is a real root which has not been removed i.e. a root to test
    
  end do ! next test root
  
! work out the final polynomial order
! and assemble the polynomial with roots removed.

  Alocal=allocate_complex_polynomial(0)
  Alocal%coeff(0)=1d0
  
  term=allocate_complex_polynomial(1)
  term%coeff(1)=1d0
  
  if (local_verbose) write(*,*)'Assembling the final polynomial'
      
  if (local_verbose) then
    write(*,*)'Initial polynomial:'
    CALL write_complex_poly_local2('Alocal','x',Alocal)  
  end if
  
  final_poly_order=0
  do root1=1,poly_order
  
    if (root_check(root1).NE.root_removed) then
    
      if (local_verbose) write(*,*)'Adding root number',root1,' :',roots(root1)
      
      final_poly_order=final_poly_order+1
      term%coeff(0)=-roots(root1)
      
      Alocal2=Alocal*term
      CALL deallocate_complex_poly(Alocal)
      Alocal=Alocal2
      
      if (local_verbose) then
        write(*,*)'Cumulative polynomial:'
        CALL write_complex_poly_local2('Alocal','x',Alocal)  
      end if

    else
    
      if (local_verbose) write(*,*)'Not including root number',root1,' :',roots(root1)
    
    end if
      
    if (local_verbose) then
      write(*,*)'Cumulative polynomial:'
      CALL write_complex_poly_local2('Alocal','x',Alocal)  
    end if
    
  end do
  
! copy the final polynomial to A

  A=allocate_polynomial(Alocal%order)
  
  A%coeff(0:A%order)=dble(Alocal%coeff(0:A%order)*gain)
  
  if (local_verbose) then
    write(*,*)'Polynomial with even multiplicity repeated roots removed:'
    CALL write_poly_local2('A','x',A)  
  end if
            
  CALL deallocate_complex_poly(term)
  CALL deallocate_complex_poly(Alocal)
  CALL deallocate_complex_poly(Alocal2)
  DEALLOCATE( roots )
  DEALLOCATE( root_check )
  DEALLOCATE( final_roots )
   
  RETURN
  END SUBROUTINE remove_even_multiple_zeros
!
! NAME
!      check_for_multiple_roots
!      
! DESCRIPTION
!      look for multiple roots of of A(x) 
!     
!
! SEE ALSO
!
!
! HISTORY
!
!     started 25/9/17 CJS
!
  SUBROUTINE check_for_multiple_roots(A,multiple_roots,local_verbose)

USE type_specifications
USE general_module
USE constants
  
IMPLICIT NONE
  
  type(polynomial),intent(INOUT)  :: A
  logical,intent(OUT) :: multiple_roots
  logical,intent(IN)  :: local_verbose

! local variables  

  integer :: i
  integer      :: poly_order
  complex(dp),allocatable  :: roots(:)
  
  integer :: root1,root2
  
!START
  
  if (local_verbose) then
    write(*,*)'CALLED: check for multiple roots'
    write(*,*)'Input polynomial:'
    CALL write_poly_local2('A','x',A)  
    write(*,*)'Full precision:'
    do i=0,A%order
      write(*,*)'coeff ',i,' ',A%coeff(i)
    end do
  end if
        
! calculate the zeros of A(x)

  poly_order=A%order
  
  ALLOCATE( roots(1:poly_order) )
  
  CALL findroots(A,roots,poly_order)
  
  if (local_verbose) then
  
    do i=1,poly_order
      write(*,*)'zero ',i,' =',roots(i),evaluate_polynomial(A,roots(i))
    end do
        
  end if
  
! look for multiple roots
  multiple_roots=.FALSE.
  
  do root1=1,poly_order-1
  
    do root2=root1+1,poly_order
    
        if ( (abs(roots(root1)-roots(root2)).LT.zero_test_small) ) then
        
! this is a repeated root set the mutliple roots flag
          multiple_roots=.TRUE.
          if (local_verbose) write(*,*)'Found repeated root:',root2
         
        end if
  
      end do ! next root 2
         
  end do ! next test root
 
   
  RETURN
  END SUBROUTINE check_for_multiple_roots

