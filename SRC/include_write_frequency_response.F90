! evaluate the frequency response of the continued fraction function

  open(unit=30,file='Continued_fraction.fout')

  wstep=(wmax-wmin)/(nw-1)
  
  do i=1,nw
  
    s=(0d0,1d0)*(wmin+wstep*(i-1))
      
    last_CF_term=(0d0,0d0)
    last_type=CFtype(n_branches)
    
    do loop=n_branches,1,-1    ! note evaluate the continued fraction from the bottom up
       
      R=CFterm(loop,1)
      L=CFterm(loop,2)
      C=CFterm(loop,3)
      
      type=CFtype(loop)
      
      if (type.EQ.series_RLC) then  
       
        CF_term=s*L/(s*s*L*C+s*L/R+1D0)
    
      else if (type.EQ.series_LC) then   
       
        CF_term=s*L/(s*s*L*C+1d0)
        
      else if (type.EQ.series_RC) then   
       
        CF_term=R/(s*C*R+1d0)
    
      else if (type.EQ.series_RL) then   
       
        CF_term=s*L/(s*L/R+1D0)
    
      else if (type.EQ.series_C) then   
       
        CF_term=1d0/(s*C)
    
      else if (type.EQ.series_L) then   
       
        CF_term=s*L
    
      else if (type.EQ.series_R) then   
       
        CF_term=R    
   
! admittance terms
    
      else if (type.EQ.shunt_RLC) then   
        
        CF_term=s*C/(s*s*L*C+s*R*C+1D0)
      
      else if (type.EQ.shunt_LC) then   
       
        CF_term=s*C/(s*s*L*C+1D0)
    
      else if (type.EQ.shunt_RC) then   
       
        CF_term=s*C/(s*C*R+1D0)
    
      else if (type.EQ.shunt_RL) then   
       
        CF_term=1d0/(s*L+R)
    
      else if (type.EQ.shunt_C) then   
       
        CF_term=(s*C)
    
      else if (type.EQ.shunt_L) then   
       
        CF_term=1d0/(s*L)
    
      else if (type.EQ.shunt_R) then   
    
        CF_term=1d0/R

      end if 
      
      if (loop.NE.max_order) then 
        
        if (last_type*type.LT.0) then
! we switch from impedance to admittance or vice versa in the network
          CF_term=CF_term+1d0/last_CF_term
        else
          CF_term=CF_term+last_CF_term         
        end if
        
      end if
      
      last_CF_term=CF_term
      last_type=type
   
    end do
    
    if (last_type.GT.0d0) then
      H_CF=last_CF_term
    else
      H_CF=1d0/last_CF_term
    end if
    write(30,8030)wnorm_save*(wmin+wstep*(i-1))/6.28318530718,real(H_CF),aimag(H_CF)
  
  end do
  
  close(unit=30)
