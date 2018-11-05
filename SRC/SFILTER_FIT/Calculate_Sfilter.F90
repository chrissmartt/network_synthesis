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
! File contents:
! SUBROUTINE Calculate_Sfilter
!
! NAME
!     Calculate_Sfilter
!
! DESCRIPTION
!
! Fit rational function model to complex frequency domain input data over a range of model orders
! If the model order is positive on input then the process calculates the filter function of that model order only
! If the model order is negative on input then the 'optimium' model  is returned
! between model order 0 and abs(model_order)
!
! Two types of model are considered - s general model in which no constraints are placed on the filter function
! and a model for the propagation correction in which the constraint that the value of the filter function at
! zero frequency is equal to 1.0
!     
! COMMENTS
!     The stopping criterion for calculation of the optimum model is an inflection in the function log_error(model_order)
!
! HISTORY
!
!     started 28/04/2016 CJS 
!     10/05/2016 CJS use automatic order criterion if the specified order is negative, otherwise
!                    return the order model specified.
!
!
SUBROUTINE Calculate_Sfilter(Hp,f,n_frequencies,Hfilter,model_order,relative_border,fit_type)

USE type_specifications
USE general_module
USE constants
USE filter_module
USE maths

IMPLICIT NONE

! variables passed to subroutine

integer,intent(IN)                 :: n_frequencies           ! number of frequencies at which we have data points
complex(dp),intent(IN)             :: Hp(1:n_frequencies)     ! list of complex function values at each frequency
real(dp),intent(IN)                :: f(1:n_frequencies)      ! list of frequencies at which we have data
type(Sfilter),intent(OUT)          :: Hfilter                 ! Output 'best fit' rational (filter) function
integer,intent(IN)                 :: model_order             ! Numerator model order to calculate
integer,intent(IN)                 :: relative_border         ! Denominator model order to calculate relative to numerator order, border=aorder+relative_border
integer,intent(IN)                 :: fit_type                ! Type of model fit. fit_type=0 has no restrictions on f=0 value
                                                              !                    fit_type=1 imposes Hfilter(f=0)=1.0
! local variables
  
  integer  :: order
  integer  :: aorder
  integer  :: border
  integer  :: max_order
  
  real(dp) :: MSE
  real(dp) :: test
  
  real(dp) :: ftest(0:abs(model_order))
  real(dp) :: MSE_save(0:abs(model_order))
  type(Sfilter)  :: Hfilter_save(0:abs(model_order))
  integer  :: min_MSE_order
  real(dp) :: min_MSE

! START

  if (verbose) write(*,*)'CALLED: Calculate_Sfilter'
  
  if (model_order.LT.0) then

! choose the best model fit up to the abs value of the model order specified  

    min_MSE_order=0
    min_MSE=1D30    ! set minimum mean square error to a very large number initially
    
    do order=0,abs(model_order)

! calculate the numerator and denominator model orders 
      aorder=order
      border=aorder+relative_border

! Calculate the best fit model for this numerator and denominator model order using the Weiner Hopf method 
      CALL Weiner_Hopf(Hp,f,n_frequencies,Hfilter_save(order),aorder,border,fit_type) 
      
      max_order=order ! this is the maximum order filter calculated so far

! calculate the mean square error (MSE) for this model order      
      CALL filter_MSE(Hp,f,n_frequencies,Hfilter_save(order),MSE_save(order))

! save the minimum MSE and the corresponding model order      
      if (MSE_save(order).LT.min_MSE) then
        min_MSE=MSE_save(order)
        min_MSE_order=order
      end if

! evaluate the function which is used to determine the 'optimum' model order      
      ftest(order)=log10(MSE_save(order))

! evaluate the test criterion for  the optimum model order. the function test is an
! estimatre of the second derivative of the mean square error as a funtion of model order      
      if(aorder.GE.2) then
        test=ftest(order)-2d0*ftest(order-1)+ftest(order-2)
      else
        test=-1d0
      end if
  
      if (verbose) write(*,*)'A order=',aorder,'B order=',border,' Mean Square Error=',MSE_save(order),' test',test
    
      if(test.GT.0d0) then
! the rate of decrease of the mean square error is slowing so we are getting to the diminishing returns stage...
! therefore stop here and choose the best model order so far.
        EXIT   
      end if
    
    end do ! next model order to calculate

! Return the best filter model calculated 
   Hfilter=Hfilter_save(min_MSE_order)

! deallocate all the filters used    
    do order=0,max_order
      CALL deallocate_Sfilter(Hfilter_save(order))
    end do
  
  else
  
! evaluate the model order specified and use that

    aorder=model_order
    border=aorder+relative_border
 
! Calculate the best fit model for this numerator and denominator model order using the Weiner Hopf method 
    CALL Weiner_Hopf(Hp,f,n_frequencies,Hfilter,aorder,border,fit_type) 
  
! calculate the mean square error (MSE) for this model order      
    CALL filter_MSE(Hp,f,n_frequencies,Hfilter,MSE)
      
    if (verbose) write(*,*)'A order=',aorder,'B order=',border,' Mean Square Error=',MSE  
  
  end if
  
  if (verbose) write(*,*)'FINISHED: Calculate_Sfilter'
                   
  RETURN

END SUBROUTINE Calculate_Sfilter
