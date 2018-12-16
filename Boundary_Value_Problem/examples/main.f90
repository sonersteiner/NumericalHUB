module examples

use Boundary_value_problems
implicit none

contains

!********************************************************************************************
!*
!********************************************************************************************
subroutine Test_BVP1D

      
       
       integer, parameter :: N = 30 
       real :: x(0:N), U(0:N) 
       real :: x0 = 0 , xf = 1
       integer :: i
       real :: pi = 4 * atan(1.0)

        x = [ (x0 + (xf-x0)*i/N, i=0, N) ]

       call Linear_Boundary_Value_Problem(  x_nodes = x, Order = 4, Differential_operator = L, Boundary_conditions = BCs, Solution = U )

       write(*,*) " maxval, minval U =", maxval(U), minval(U) 
       
       call scrmod('revers') 
       
       call metafl('xwin')
       call disini() 
       CALL TEXMOD('ON')
       CALL TITLIN (" $ y_{xx} + \pi y = 0  $ ", 2)
       CALL TITLIN (" $ y(0) = 1, \quad y_{x}(1) = - \pi $  ", 4)
       call qplot(x, U, N+1)

contains 
!---------------------------------------------------
real function L(x, y, yx, yxx) 
           real, intent(in) :: x, y, yx, yxx 
    
           L = yxx + pi*pi*y
           
end function 
!-----------------------------------------------------
real function BCs(x, y, yx) 
           real, intent(in) :: x, y, yx 
           

   if (x==x0) then
                      BCs = y -1
   elseif (x==xf) then
                      BCs = yx + pi
   else 
        write(*,*) " Error BCs x=", x  
        write(*,*) " a, b=", x0, xf
        stop  
   endif 
           
                 
end function  

end subroutine 

 
   

!********************************************************************************************
!*
!********************************************************************************************
subroutine Test_BVP2D

            
       integer, parameter :: Nx = 15, Ny = 15
       real :: x(0:Nx), y(0:Ny), U(0:Nx, 0:Ny) 
      
       integer :: i, j
       real :: a=-1, b=1
       real :: pi = 4 * atan(1.0) 

    
       x = [ (a + (b-a)*i/Nx, i=0, Nx) ]
       y = [ (a + (b-a)*j/Ny, j=0, Ny) ]
       
       
      call Linear_Boundary_Value_Problem( x_nodes = x, y_nodes = y, Order = 10, Differential_operator = L, Boundary_conditions = BCs, Solution = U )
                               
       
       call scrmod('revers') 
       call disini() 
       CALL TEXMOD('ON')
       CALL TITLIN (" $ u_{xx} + u_{yy} =  \pi^2  \sin( \pi x)  \sin( \pi y ) $  ", 2)
       CALL TITLIN (" $  u=0 \quad $ boundary $ [-1,1] x [-1,1] $ ", 4)
       call qplcon( U, Nx+1, Ny+1, 20);
   
contains 
!-----------------------------------------------------------------------
real function L(x, y, u, ux, uy, uxx, uyy, uxy) 
           real, intent(in) :: x, y, u, ux, uy, uxx, uyy, uxy 
        
        real :: pi = 4 * atan(1.0)
           
           L = uxx + uyy - pi**2 * sin(pi*x) * sin(pi*y) 
                 
end function 
!-----------------------------------------------------------------------
real function BCs(x, y, u, ux, uy) 
           real, intent(in) :: x, y, u, ux, uy 
           

   if (x==a) then
                      BCs = u 
   elseif (x==b) then
                      BCs = u 
   elseif (y==a) then
                      BCs = u 
   elseif (y==b) then
                      BCs = u                    
   else 
        write(*,*) " Error BCs x=", x  
        write(*,*) " a, b=", a, b
        stop  
   endif 
           
                 
end function  

end subroutine 






!********************************************************************************************
!*
!********************************************************************************************

subroutine Test_non_linear_BVP1D

      
       integer, parameter :: N = 30
       real :: x(0:N), U(0:N), err(0:N)
       real :: x0=1d0 , xf=2d0
       integer :: i

       x = [ (x0 + (xf-x0)*i/N, i=0, N) ]
       U = 0


       call Non_Linear_Boundary_Value_Problem( x_nodes = x, Order = 4, Differential_operator = L, Boundary_conditions = BCs, Solution = U) 
       
       write(*,*) " maxval, minval U =", maxval(U), minval(U)

       do i=0,N
        err(i)=U(i)-1/(x(i)+1)  ! exact solution
       end do

       call scrmod('revers')
       call disini() 
       CALL TEXMOD('ON')
       CALL TITLIN (" $ y_{xx} - y^3 -y y_x  =  0 $  ", 2)
       CALL TITLIN (" $  y(1)=0.5, \quad y(2) = 1/3  $ ", 4)
       call qplot(x, U, N+1)
       
       
       call disini() 
       CALL TEXMOD('ON')
       CALL TITLIN (" $ Error = y(x) - y_{exact} $  ", 2)
       call qplot(x, err, N+1)




contains
!-----------------------------------------------------------------------
real function L(x, y, yx, yxx)
           real, intent(in) :: x, y, yx, yxx

           real, parameter :: lambda=1

           L = yxx - y**3 - y * yx

end function

!-----------------------------------------------------------------------
real function BCs(x, y, yx)
           real, intent(in) :: x, y, yx


   if (x==x0) then
                      BCs = y - 0.5
   elseif (x==xf) then
                      BCs = y - 1/3.
   else
        write(*,*) " Error NLBCs x=", x
        write(*,*) " a, b=", x0, xf
        read(*,*)
   endif


end function

end subroutine



!********************************************************************************************
!*
!********************************************************************************************
subroutine Test_non_linear_BVP2D

   

       integer, parameter :: Nx = 10, Ny = 10
       real :: x(0:Nx), y(0:Ny), U(0:Nx, 0:Ny)

       integer :: i, j
       real :: a=0, b=1
       real :: pi = 4 * atan(1.0)


       x = [ (a + (b-a)*i/Nx, i=0, Nx) ]
       y = [ (a + (b-a)*j/Ny, j=0, Ny) ]
       U = 1 
     

      call Non_Linear_Boundary_Value_Problem( x_nodes = x, y_nodes = y, Order = 6,  & 
                                              Differential_operator = L, Boundary_conditions = BCs, Solution = U ) 
                          
   
       call scrmod('revers')
       call disini() 
       CALL TEXMOD('ON')
       CALL TITLIN (" $ ( u_{xx} + u_{yy} ) u = 0  $ ", 2)
       CALL TITLIN (" $ u(0,y) = 0, \quad u(1,y) = y, \quad u(x,0) = 0, \quad u(x,1) = x $  ", 4)
       call qplcon( U, Nx+1, Ny+1, 20);


contains
!-----------------------------------------------------------------------
real function L(x, y, u, ux, uy, uxx, uyy, uxy)
           real, intent(in) :: x, y, u, ux, uy, uxx, uyy, uxy

           L = ( uxx + uyy) * u

end function
!-----------------------------------------------------------------------
real function BCs(x, y, u, ux, uy)
           real, intent(in) :: x, y, u, ux, uy


   if (x==a) then
                      BCs = u
   elseif (x==b) then
                      BCs = u-y
   elseif (y==a) then
                      BCs = u
   elseif (y==b) then
                      BCs = u-x
   else
        write(*,*) " Error BCs x=", x
        write(*,*) " a, b=", a, b
        stop 
   endif


end function

end subroutine

end module 



!*************************************
program Bounday_value_examples 

   
    use examples
    implicit none 
 
    
    Write(*,*) "   Test linear_BVP_1D" 
    call Test_BVP1D
    
    Write(*,*) "   Test linear_BVP_2D" 
    call Test_BVP2D
    
    Write(*,*) " Test non_linear_BVP1D"
    call Test_non_linear_BVP1D
  
    Write(*,*) " Test non_linear_BVP2D"
    call Test_non_linear_BVP2D
    


end program 










 

