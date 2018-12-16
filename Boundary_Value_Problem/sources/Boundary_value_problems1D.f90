!******************************************************************************
!
! Author: Juan A Hernandez (juanantonio.hernandez@upm.es) & Javier Escoto 
!******************************************************************************
module Boundary_value_problems1D 

use Linear_systems
use Non_Linear_Systems 
use Finite_differences

use Dependencies
use Linearity
 
implicit none  

private
public :: Boundary_Value_Problem1D
public :: linear1D
   
abstract interface  

       real function DifferentialOperator1D(x, u, ux, uxx) 
                        real, intent(in) :: x, u, ux, uxx 
       end function  

       real function      BC1D(x, u, ux) 
           real, intent(in) :: x, u, ux 
       end function  

       subroutine  NonLinearSolver(Function, x)
          use Jacobian_module
          procedure(FunctionRN_RN) :: Function
          real, intent(inout) :: x(:)
       end subroutine  


 end interface

logical :: linear1D = .true. 
logical :: dU(2) ! matrix of dependencies( order )       
 
  
contains 
    
!**************************************************************************************************************
! Boundary value problem 
!
!       Differential_operator(x, u, ux, uxx) = 0. Linear and non linear 
!       Boundary_conditions(x, u, ux) at x=0 and x=L 
!
!**************************************************************************************************************
subroutine Boundary_Value_Problem1D( x_nodes, Order, Differential_operator, & 
                                     Boundary_conditions, Solution, Solver )
       
     real, intent(in) :: x_nodes(0:)
     integer, intent(in) :: Order
     procedure (DifferentialOperator1D) :: Differential_operator
     procedure (BC1D) ::  Boundary_conditions
     real, intent(inout) :: Solution(0:) 
     procedure (NonLinearSolver), optional:: Solver
     
     
     dU = Dependencies_BVP_1D( Differential_operator ) 
     
     linear1D = Linearity_BVP_1D( Differential_operator ) 
     
     if (linear1D) then 
         
       call Linear_Boundary_Value_Problem1D( x_nodes, Order,                &
                                             Differential_operator,         &
                                             Boundary_conditions, Solution)
         
     else 
         
       call Non_Linear_Boundary_Value_Problem1D( x_nodes, Order,            & 
                                                 Differential_operator,     & 
                                                 Boundary_conditions,       & 
                                                 Solver, Solution)
         
     end if 
   

end subroutine 


!**************************************************************************************************************
! Linear Boundary Value Problem 1D
!
!       Differential_operator(x, u, ux, uxx)
!       Boundary_conditions(x, u, ux) at x=0 and x=L 
!
!**************************************************************************************************************
subroutine Linear_Boundary_Value_Problem1D( x_nodes, Order, Differential_operator, &  
                                            Boundary_conditions, Solution)
        
     real, intent(in) :: x_nodes(0:)
     integer, intent(in) :: Order
     procedure (DifferentialOperator1D) :: Differential_operator
     procedure (BC1D) ::  Boundary_conditions
     real, intent(out) :: Solution(0:) 
     

!  *** auxiliary variables
       integer ::  i, N 
       
       real, allocatable ::  bi(:), F(:), U(:), Ux(:), Uxx(:),  &
                             Difference_operator(:,:)       

!  *** Integration domain 
       N = size(x_nodes) - 1 
       allocate( bi(0:N), F(0:N), U(0:N), Ux(0:N), Uxx(0:N),    & 
                 Difference_operator(0:N, 0:N) ) 
          
!  *** independent term  A U = b  ( U = inverse(A) b )         
       U = 0 
       call Difference_equation1DL(x_nodes, U, Ux, Uxx, bi) 
       
       
       
!  *** Delta kronecker to calculate the difference operator        
       do i=0, N
          U(0:N) = 0 
          U(i) = 1.0 
          
          call Difference_equation1DL(x_nodes, U, Ux, Uxx, F) 
          Difference_operator(0:N, i) = F - bi
       enddo 
  
!  *** solve the linear system of equations  
       call LU_Factorization(Difference_operator)

       Solution = Solve_LU( Difference_operator, -bi ) 
      
      deallocate( U, Ux, Uxx, F, Difference_operator, bi ) 

contains 
!-----------------------------------------------------------------
subroutine Difference_equation1DL(x, W, Wx, Wxx, F) 
           real, intent(in) :: x(0:), W(0:)
           real, intent(out) :: Wx(0:), Wxx(0:), F(0:)  

    integer :: i 
    real :: D, C 
    
    
        if (dU(1)) call Derivative( "x", 1, W, Wx) 
        if (dU(2)) call Derivative( "x", 2, W, Wxx) 
        
!  ***  boundary conditions
        do i=0, N, N 
             D = Differential_operator( x(i), W(i), Wx(i), Wxx(i) ) 
             C = Boundary_conditions( x(i), W(i), Wx(i) )
             if (C == FREE_BOUNDARY_CONDITION) then 
                 F(i) = D 
             else 
                 F(i) = C 
             end if 
         end do   
           
!  *** inner grid points         
        do i=1, N-1
            F(i) = Differential_operator( x(i), W(i), Wx(i), Wxx(i) ) 
        enddo 
                 
end subroutine 

end subroutine 

!**************************************************************************************************************
! Non linear BVP 
!
!       Differential_operator(x, y, u, ux, uy, uxx, uyy, uxy)
!       Boundary_conditions(x, u, ux, uy)  at boundary [a, b]x[c,d] 
!
!**************************************************************************************************************
subroutine Non_Linear_Boundary_Value_Problem1D( x_nodes, Order,             &
              Differential_operator,  Boundary_conditions, Solver, Solution)

     real, intent(in) :: x_nodes(0:)
     integer, intent(in) :: Order
     procedure (DifferentialOperator1D) :: Differential_operator
     procedure (BC1D) ::  Boundary_conditions
     procedure (NonLinearSolver), optional:: Solver
     real, intent(inout) :: Solution(0:)


!  *** variables specification
       integer ::   N
     

!  *** Integration domain
       N = size(x_nodes) - 1
  
!  *** Non linear solver
       if (present(Solver))  then
                   call Solver(System_BVP, Solution)
       else
                   call Newton(System_BVP, Solution)
       end if


    contains
!-----------------------------------------------------------------------
 function System_BVP(U) result(F)
        real, intent (in) :: U(:)
        real :: F(size(U)) 
        

            call Difference_equation1DNL(x_nodes, U, F)

end function 
!-----------------------------------------------------------------------
subroutine Difference_equation1DNL(x, W, F)
           real, intent(in) :: x(0:), W(0:)
           real, intent(out) :: F(0:)

    real :: Wx(0:N), Wxx(0:N)
    integer :: i
    real :: C, D 

     if (dU(1)) call Derivative( "x", 1, W, Wx)
     if (dU(2)) call Derivative( "x", 2, W, Wxx)

 !  ***  boundary conditions
         do i=0, N, N 
             D = Differential_operator( x(i), W(i), Wx(i), Wxx(i) ) 
             C = Boundary_conditions( x(i), W(i), Wx(i) )
             if (C == FREE_BOUNDARY_CONDITION) then 
                 F(i) = D 
             else 
                 F(i) = C 
             end if 
         end do   
           
!  ***   inner grid points         
         do i=1, N-1
            F(i) = Differential_operator( x(i), W(i), Wx(i), Wxx(i) ) 
         enddo     
     

end subroutine


end subroutine



    
    
    
end module
