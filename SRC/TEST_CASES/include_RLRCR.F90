! Set the numbers of each component type and the number of intermediate filters 
! required for the impedance calculation

  filename='test_RLRCR'

  NR=3
  NL=1
  NC=1
  NZ=1
  NY=5

! Allocate the R,L and C arrays plus Y and Z filter arrays
  write(*,*)'Allocate the R,L and C arrays plus Y and Z filter arrays'
  allocate( R(1:NR) )
  allocate( L(1:NL) )
  allocate( C(1:NC) )
  allocate( Z(1:NZ) )
  allocate( Y(1:NY) )
  
! Set the component values
  write(*,*)'Set the component values'
  R(1)=1d0
  L(1)=3d0
  R(2)=2d0
  C(1)=5d0
  R(3)=7d0
  
! Set the initial branch impedance/ admittances
  write(*,*)'Set the initial branch impedance/ admittances'

  Y(1)=allocate_Sfilter(0,0)
  Y(1)%a%coeff(0)=1d0
  
  Y(1)%b%coeff(0)=R(1)

  Y(2)=allocate_Sfilter(0,1)
  Y(2)%a%coeff(0)=1d0
  
  Y(2)%b%coeff(0)=R(2)
  Y(2)%b%coeff(1)=L(1)
  
  Y(3)=allocate_Sfilter(1,1)
  Y(3)%a%coeff(0)=0d0
  Y(3)%a%coeff(1)=C(1)
  
  Y(3)%b%coeff(0)=1D0
  Y(3)%b%coeff(1)=C(1)*R(3)
  
! Calculate the impedance filter function from the individual branch 
! impedance/admittance matrices
  write(*,*)'Calculate the impedance filter function'
  
  Y(4)=Y(1)+Y(2)
  Y(5)=Y(4)+Y(3)
  Z(1)=reciprocal_Sfilter(Y(5))
  
  H=Z(1)
