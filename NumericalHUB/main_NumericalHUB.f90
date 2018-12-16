!***************************************************
!* Book:  How to learn Applied maths
!***************************************************    
program main_NumericalHUB

       
       use API_Example_Systems_of_Equations
       use API_Example_Lagrange_Interpolation
       use API_Example_Chebyshev_Interpolation
       use API_Example_Cauchy_Problem
       use API_Example_Finite_Differences
       use API_Example_Boundary_Value_Problem
       
       use API_Example_Initial_Value_Boundary_Problem
       use API_Example_IVBP_and_BVP
       
       use my_examples    
          
       implicit none 
       
   
   call my_Kepler_orbit    
        
   call Test_preliminaries 
   call Systems_of_Equations_examples 
   call Lagrange_Interpolation_examples 
   call Chebyshev_Interpolation_examples
   
   call Cauchy_problem_examples
   call Finite_difference_examples
  
   call BVP_examples
   
   call IVBP_examples 
   
   call Test_IVBP_and_BVP
 
 
 
end program  

