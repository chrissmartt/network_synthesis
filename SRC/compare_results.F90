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
!
! MODULE compare_results_data
!   CONTAINS 
!   SUBROUTINE read_compare_results_input_data_file
!   SUBROUTINE deallocate_compare_results_data
!   PROGRAM compare_results
!
! DESCRIPTION
!
! Software to read two datasets from files and compare them to produce a numerical measure of the difference
!
! The input data is two columns,   x and f(x) thus the data could be either time domain or frequency domain (magnitude) data
! the columns for x and f(x) must be specified.
!
! The measure of the difference is calculated as  (integral |f1(x)-f2(x)|dx) /(integral dx) where the 
! integration is over the range of x for which the data sets overlap.
!     
! COMMENTS
!
!
! HISTORY
!
!     started 24/11/2015 CJS 
!
! _________________________________________________________________
!
! MODULE NAME
!    compare_results_data
!
! DESCRIPTION
!
! module containing the following:
! compare_results input data structure 
! subroutine to allocate and read this data from a file
! subroutine to deallocate the data structure
! program to calculate the comparison measure.
!     
! COMMENTS
!
!
! HISTORY
!
!     started 24/11/2015 CJS 
!
!
MODULE compare_results_data

  USE type_specifications

  IMPLICIT NONE

! Structure to hold the data sets for comparison
  TYPE comparison_dataset
  
    integer                 :: n_samples    ! number of data samples
    real(dp),allocatable    :: x(:)         ! list of x values
    real(dp),allocatable    :: fx(:)        ! list of function of x values
    real(dp)                :: xmin,xmax    ! range of x for which data is available

  END TYPE comparison_dataset
  
  CONTAINS

! SUBROUTINE NAME
!    read_compare_results_input_data_file
!
! DESCRIPTION
!     
! read data from a file into x and fx arrays within the camparison_dataset structure
! return the complete camparison_dataset structure: x and fx arrays and also the min and max x values and number of data samples 
!
! COMMENTS
! If the expected data columns cannot be read in then the line is ignored and the process moves on to the next line
! This allows for comment lines in the data files.
!
! HISTORY
!
!     started 20/11/2015 CJS 
!
!
  SUBROUTINE read_compare_results_input_data_file(dataset)

  USE type_specifications
  USE general_module

  IMPLICIT NONE

! variables passed to the subroutine

    type(comparison_dataset),intent(OUT)    :: dataset

! local variables

    character(len=filename_length)    :: filename
    character(len=line_length)        :: line
    
    logical        :: file_exists
    
    integer        :: xcol,fxcol,maxcol    ! columns for x, f(x) and the maximum column number ot be read
    
    real(dp),allocatable    :: r_in(:)     ! array of reals for to read the columns of data into

! temporary local variables    
    integer        :: read_loop
    integer        :: sample
    integer        :: comment
    integer        :: i
   
    integer        :: ierr     ! integer to return error codes from file reads

! START

    write(*,*)'Enter the filename for frequency domain data'

    read(*,'(A)')filename

    inquire(file=trim(filename),exist=file_exists)
    if (.NOT.file_exists) then
      write(*,*)'Error in compare_results. Cannot find the file:',trim(filename)
      run_status='ERROR: in read_compare_results_input_data_file, Cannot find the file:'//trim(filename)
      CALL write_program_status()
      STOP 1
    end if 
   
! open and read the file
  
    OPEN(unit=local_file_unit,file=filename)

    write(*,*)'opened file:',trim(filename)
    
    write(*,*)'Enter the column for x data'
    read(*,*)xcol
    
    write(*,*)'Enter the column for f(x) data'
    read(*,*)fxcol
    
    maxcol=max(xcol,fxcol)
    ALLOCATE( r_in(1:maxcol) )

! read the file twice, first to count the number of valid samples then allocate data, second to set the frequency and data values
    do read_loop=1,2
    
      sample=0
      comment=0
      
      do

! read the next line of the file into a string   
        read(local_file_unit,'(A)',IOSTAT=ierr)line
        if (ierr.LT.0) EXIT    ! end of file
    
! try to read a list of numbers from the string
        read(line,*,IOSTAT=ierr)(r_in(i),i=1,maxcol)
    
        if (ierr.eq.0) then 
! this looks like a valid line of data so increase the sample number and set the frequency and value (on second read) 
          sample=sample+1
          if (read_loop.EQ.2) then
      
            dataset%x(sample) =r_in(xcol)
            dataset%fx(sample)=r_in(fxcol)
! work out the frequency range of the data
            dataset%xmin=min( dataset%xmin,dataset%x(sample) )
            dataset%xmax=max( dataset%xmax,dataset%x(sample) )
        
          end if  ! read_loop=2
          
        else ! not a valid data line i.e. a comment
        
          comment=comment+1
          
        end if ! comment line
    
      end do ! read the next line from the file
      
      if (read_loop.EQ.1) then
      
        if (sample.LE.1) then
          write(*,*)'Error in compare_results. No data found in the file:',trim(filename)
          run_status='ERROR. No data found in the file:'//trim(filename)
          CALL write_program_status()
          STOP 1
        end if 
      
        dataset%n_samples=sample
        ALLOCATE( dataset%x(1:dataset%n_samples) )
        ALLOCATE( dataset%fx(1:dataset%n_samples) )
        dataset%xmin=1D30
        dataset%xmax=-1D30
        REWIND(unit=local_file_unit)
    
      end if ! read_loop.EQ.1
    
    end do ! next read_loop
    
    DEALLOCATE( r_in )
    
! close the file
  
    CLOSE(unit=local_file_unit)
    
    write(*,*)'Number of samples read =',dataset%n_samples
    write(*,*)'Number of comment lines=',comment

  END SUBROUTINE read_compare_results_input_data_file
!
! ________________________________________________________________
!
!  
! SUBROUTINE NAME
!    deallocate_compare_results_data
!
! DESCRIPTION
!     
! deallocate the allocatable arrays in the compare_results_data data structure
!
! COMMENTS
! 
!
! HISTORY
!
!     started 20/11/2015 CJS 
!
!
  SUBROUTINE deallocate_compare_results_data(dataset)

  USE type_specifications

  IMPLICIT NONE

! variables passed to the subroutine

    type(comparison_dataset),intent(INOUT)    :: dataset

! START

    if (ALLOCATED( dataset%x ))  DEALLOCATE ( dataset%x )
    if (ALLOCATED( dataset%fx )) DEALLOCATE ( dataset%fx )

  END SUBROUTINE deallocate_compare_results_data

END MODULE compare_results_data
!
! ________________________________________________________________
!
!  
!
! PROGRAM NAME
!    compare_results
!
! DESCRIPTION
!    Calculate the comparison measure of the difference between two datasets.
!    PROCESS:
!    1. Read the two datasets
!    2. Calculate the x range for which the datasets overlap
!    3. Calculate the measure of the difference between the datasets as (integral |f1(x)-f2(x)|dx) /(integral dx) 
!       where the integration is over the range of x for which the data sets overlap.
!    4. Write the result to screen and file
!     
! COMMENTS
!    The integration is based on the x values in dataset 1, the function value from dataset 2 is calculated by linear 
!    interpolation hence the datasets need not have common x values.
!
!    The end points of the integration range needs special treatment maybe to make the process more rigorous.
!    
!
! HISTORY
!
!     started 20/11/2015 CJS based on compare_results.F90 from the GGI_TLM IELF code (www.github.com/ggiemr/GGI_TLM)
!
!
PROGRAM compare_results

USE type_specifications
USE compare_results_data
USE general_module

IMPLICIT NONE

! local variables
  
! frequency domain data to be read from input files
  type(comparison_dataset)    :: dataset1
  type(comparison_dataset)    :: dataset2
  
! filename for the comaprison result
  character(len=filename_length)    :: opfilename

! output quantities
  real(dp)    :: compare_results_value
  
! range of x for which the datasets overlap i.e. the integration range
  real(dp)    :: xmin_overlap,xmax_overlap
  
  integer    :: sample_min,sample_max,last_sample,sample
  
  real(dp)    :: fx_1,fx_2      ! function1 and function 2 values at the current x value 
  real(dp)    :: x_1,x_2,x      ! x values for function 2 which bracket the current x value - used for interpolation
          
  real(dp)    :: error   ! magnitude of the difference between the function values at the current x value
 
  real(dp)    :: local_compare_results_value
  integer     :: i
  
  integer    :: n_samples1,n_samples2  ! number of samples in dataset1 and dataset 2

! START

  program_name='compare_results'
  run_status='Started'
  CALL write_program_status()
  
  CALL read_version()
    
  CALL write_license()
  
  write(*,*)'compare_results Analysis'
  
  write(*,*)' '
  
! STAGE 1. read the two datasets to be compared
  CALL read_compare_results_input_data_file(dataset1)
  CALL read_compare_results_input_data_file(dataset2)
  
! STAGE 2. Calculate the overlapping x range 
  write(*,*)'Minimum x in dataset 1=',dataset1%xmin
  write(*,*)'Minimum x in dataset 2=',dataset2%xmin
  write(*,*)'Maximum x in dataset 1=',dataset1%xmax
  write(*,*)'Maximum x in dataset 2=',dataset2%xmax
  
  xmin_overlap=max(dataset1%xmin,dataset2%xmin)
  xmax_overlap=min(dataset1%xmax,dataset2%xmax)
  
  write(*,*)'The x range over which the datasets overlap is'
  write(*,*)'xmin:',xmin_overlap,' Hz, xmax:',xmax_overlap,' Hz'
  
  if ((xmax_overlap-xmin_overlap).LE.0D0) then
    run_status='ERROR: the x ranges of the datasets do not overlap'
    CALL write_program_status()
    STOP 1
  end if
    
  n_samples1=dataset1%n_samples
  n_samples2=dataset2%n_samples
  
! loop over samples of dataset1 to work out sample_min and sample_max for local_compare_results_value

  sample_min=0
  sample_max=0
  
  do sample=2,n_samples1-1
  
    if ( ( dataset1%x(sample-1).LE.xmin_overlap ).AND.( sample_min.eq.0 ) ) sample_min=sample
    if ( ( dataset1%x(sample+1).GE.xmax_overlap ).AND.( sample_max.eq.0 ) ) sample_max=sample
    
  end do

! ensure that the sample range is within the array bounds  
  if (sample_min.eq.0) sample_min=2
  if (sample_max.eq.0) sample_max=n_samples1-1
  
  write(*,*)'n_samples1 =',n_samples1
  write(*,*)'Sample_min=',sample_min,' xmin=',dataset1%x(sample_min)
  write(*,*)'Sample_max=',sample_max,' xmax=',dataset1%x(sample_max)

! STAGE 3. calculate the difference measure over the full overlap range       
       
! loop over samples in frequency range and do local_compare_results_value calculation

  local_compare_results_value=0d0

  last_sample=1

  do sample=sample_min,sample_max
  
    fx_1=dataset1%fx(sample)
    x=dataset1%x(sample)
    
! find the frequencies which lie either side of f1 and interpolate to give the value from file 2
    do i=last_sample,n_samples2-1
    
     x_1=dataset2%x(i)
     x_2=dataset2%x(i+1)
     
      if ( (x_1.le.x).AND.    &
           (x_2.gt.x)  ) then

! the required x value is bracketed by x_1 and x_2 so calculate fx_2 using linear interpolation then exit the loop       
        fx_2=dataset2%fx(i)+( (x-x_1)/(x_2-x_1) )*(dataset2%fx(i+1)-dataset2%fx(i))
        
        last_sample=i
        
        GOTO 1000
    
      end if
      
    end do  ! next sample of file 2

    write(*,*)'Sample not found in file 2, x=',x
    write(*,*)'First x=',dataset2%x(last_sample)
    write(*,*)'Last x=',dataset2%x(n_samples2-1)

    run_status='ERROR: Sample not found in file 2'
    CALL write_program_status()
    STOP 1

1000  CONTINUE

! we now have fx_1 and fx_2 at the current x value

! error is the magnitude of the difference between the values         
    error=abs(fx_1-fx_2)
    
! The contribution to the difference measure is (1/2)*f(x_i)*((x_i+1)-x(i-1))/(total range of x)
! i.e. an integration over x using a piecewise linear integration scheme
    
    local_compare_results_value=local_compare_results_value+    &
      error*( ((dataset1%x(sample+1))-(dataset1%x(sample-1)) )/2d0)     &
            /( (dataset1%x(sample_max+1))-(dataset1%x(sample_min-1)) )

  end do ! next sample in the integration over x

! STAGE 4. write the result to screen and file
              
  write(*,8000)'compare_results value:',local_compare_results_value
  write(*,*)'compare_results value:',local_compare_results_value

  write(*,*)'Enter the filename for the compare_results data'
  read(*,*)opfilename
  
  OPEN(unit=local_file_unit,file=opfilename)
  
  write(local_file_unit,8000)' ',local_compare_results_value
  
8000 format(A,F14.8)
  
! deallocate the data  
  CALL deallocate_compare_results_data(dataset1)
  CALL deallocate_compare_results_data(dataset2)

  run_status='Finished_Correctly'
  CALL write_program_status()

END PROGRAM compare_results
