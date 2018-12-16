
module my_examples

       use dislin 
       use Cauchy_Problem
       implicit none 

       real, parameter :: PI = 4 * atan(1d0) 
       
contains  
    

subroutine my_Kepler_orbit 

    integer, parameter :: N= 200000, M = 4
    real :: U(0:N, M), Time (0:N)
   
    real :: t0 = 0, tf = 2*PI*365*10 ! 10 years 
    integer :: i
  
    Time = [ (t0 + (tf -t0 )*i/N, i=0, N ) ]
    
    U(0,:) = [ 1, 0, 0, 1 ]
    
    call Cauchy_ProblemS( Time_Domain = Time ,                      & 
                          Differential_operator = Kepler,           & 
                          Scheme = Cash_Karp , Solution = U )
      
  
    CALL PAGE (2000, 2000)
    call scrmod("reverse")
    call metafl("xwin")
    call qplot( U(:,1), U(:,2), N+1) 
    
    call Cauchy_ProblemS( Time_Domain = Time ,                      & 
                          Differential_operator = Kepler,           & 
                          Scheme = Euler , Solution = U )
    
    call qplot( U(:,1), U(:,2), N+1) 
   
contains  
!-------------------------------------------------  
function Kepler(U, t) result(F) 
  real :: U(:), t
  real :: F ( size(U) ) 
    
     real :: r(2), drdt(2)   
     
     r = U(1:2);  drdt = U(3:4)
     
     F = [ drdt, - r/norm2(r)**3] 
         
end function  

end subroutine



end module  

