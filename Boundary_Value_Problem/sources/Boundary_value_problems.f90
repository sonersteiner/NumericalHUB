!******************************************************************************
!
! Author: Juan A Hernandez (juanantonio.hernandez@upm.es) & Pablo Lopez Negro & Javier Escoto 
!******************************************************************************
module Boundary_value_problems 

use Boundary_value_problems1D
use Boundary_value_problems2D
use Boundary_value_problems3D

implicit none  

private

public :: Boundary_Value_Problem
public :: linear2D, linear1D, linear3D 

 
 
 interface Boundary_Value_Problem
      module procedure Boundary_Value_Problem1D, &
                       Boundary_Value_Problem2D, & 
                       Boundary_Value_Problem2D_system
 end interface

 
 
contains 
    
 
    
end module
