module BVP_Validation
    
    use Boundary_value_problems
    use Finite_differences
       
implicit none

abstract interface  

       real function DifferentialOperator1D(x, u, ux, uxx) 
                        real, intent(in) :: x, u, ux, uxx 
       end function  
       
       real function DifferentialOperator2D(x, y, u, ux, uy, uxx, uyy, uxy) 
                        real, intent(in) :: x, y, u, ux, uy, uxx, uyy, uxy
       end function  

        function DifferentialOperator2D_system(x, y, u, ux, uy, uxx, uyy, uxy) 
                        real, intent(in) :: x, y, u(:), ux(:), uy(:), uxx(:), uyy(:), uxy(:)
                        real :: DifferentialOperator2D_system(size(u))
        end function
        
       real function      BC1D(x, u, ux) 
           real, intent(in) :: x, u, ux 
       end function
       
       real function      BC2D(x, y, u, ux, uy) 
           real, intent(in) :: x, y, u, ux, uy 
       end function  
       
       function BC2D_system(x, y, u, ux, uy) 
           real, intent(in) :: x, y, u(:), ux(:), uy(:) 
           real :: BC2D_system(size(u)) 
       end function  
       
        function Analytical_Solution1D(x)
            real, intent(in) :: x(:)
            real :: Analytical_Solution1D( size(x) )
        end function
        
        function Analytical_Solution2D(x, y)
            real, intent(in) :: x(:), y(:)
            real :: Analytical_Solution2D( size(x), size(y) )
        end function
        
        real function DifferentialOperator1DIVBP( x, t, u, ux, uxx) 
            real, intent(in) ::  x, t, u, ux, uxx
        end function
        
        real function DifferentialOperator2DIVBP(  x, y, t, U, Ux, Uy, Uxx, Uyy, Uxy ) 
            real, intent(in) :: x, y, t, U, Ux, Uy, Uxx, Uyy, Uxy 
        end function
        
        function IC1D(x)
            real, intent(in) :: x(:)
            real :: IC1D(size(x))
        end function
        
        function IC2D(x,y)
            real, intent(in) :: x(:), y(:)
            real ::  IC2D(size(x), size(y))
        end function
        
            
        
        real function BC1D_IVBP(x, t, u, ux)
            real, intent(in) :: x, t, u, ux
        end function
        real function BC2D_IVBP( x, y, t, U, Ux, Uy )
            real, intent(in) ::  x, y, t, U, Ux, Uy 
        end function
        
end interface

    contains

    
    
    
    
    
    
    
    
    
    
    
    
    
!*****     LATEX - 90 (60)    *****    
subroutine BVP_1D_Validation( Differential_operator, Boundary_conditions, Order, Spatial_Domain, Analytical_Solution, log_e, log_N)

     procedure (DifferentialOperator1D) :: Differential_operator
     procedure (BC1D) ::  Boundary_conditions
     procedure (Analytical_Solution1D) :: Analytical_Solution
     integer, intent(in) :: Order
     real, intent(in) :: Spatial_Domain(:)
     real, allocatable, intent(out) :: log_e(:), log_N(:)
     
     real :: x0 , xf, error
     real, allocatable :: x_nodes(:),  U(:), U_an(:)
     
     integer :: i, j, N
    
    integer :: N0 = 30, Nf = 200
    integer :: gap = 10
     
     allocate( log_e(0:(Nf-N0)/gap), log_N(0:(Nf-N0)/gap) )
     
     x0 = Spatial_Domain(1)
     xf = Spatial_Domain(size(Spatial_Domain))

     

     do i = 0, (Nf-N0)/gap
         
         N = N0 + gap*i
        
        allocate ( x_nodes(0:N), U(0:N), U_an(0:N) )
        
        x_nodes(0) = x0; x_nodes(N) = xf

        call Grid_Initialization( grid_spacing = "uniform", direction = "x", q = Order, nodes = x_nodes )

        U_an = Analytical_Solution(x_nodes)

        
        call Boundary_Value_Problem( x_nodes = x_nodes, Order = Order, &                  
                    Differential_operator = Differential_operator,     & 
                    Boundary_conditions = Boundary_conditions,         & 
                    Solution = U )
        

        error = abs( norm2(U_an - U) )
        
        log_e(i) = log10( error )
        log_N(i) = log10( real(N) ) 
        
        deallocate( x_nodes, U, U_an)

     end do
     
end subroutine

!*****     LATEX - 114     *****  
subroutine BVP_2D_Validation( Differential_operator, Boundary_conditions, Order, Spatial_Domain_x, Spatial_Domain_y, Analytical_Solution, log_e, log_N)

     procedure (DifferentialOperator2D) :: Differential_operator
     procedure (BC2D) ::  Boundary_conditions
     procedure (Analytical_Solution2D) :: Analytical_Solution
     integer, intent(in) :: Order
     real, intent(in) :: Spatial_Domain_x(:), Spatial_Domain_y(:)
     real, allocatable, intent(out) :: log_e(:), log_N(:)
     
     real :: x0, y0, xf, yf, error
     real, allocatable :: x_nodes(:), y_nodes(:),  U(:,:), U_an(:,:)
     
     integer :: i, j, N
    
     integer :: N0 = 20, Nf = 50
     integer :: gap = 5
     
     allocate( log_e(0:(Nf-N0)/gap), log_N(0:(Nf-N0)/gap) )
     
     x0 = Spatial_Domain_x(1)
     xf = Spatial_Domain_x( size(Spatial_Domain_x) )
     
     y0 = Spatial_Domain_y(1)
     yf = Spatial_Domain_y( size(Spatial_Domain_y) )

     

     do i = 0, (Nf-N0)/gap
         
         N = N0 + gap*i
        
        allocate ( x_nodes(0:N), y_nodes(0:N), U(0:N, 0:N), U_an(0:N, 0:N) )
        
        x_nodes(0) = x0; x_nodes(N) = xf
        y_nodes(0) = y0; y_nodes(N) = yf

        call Grid_Initialization( grid_spacing = "uniform", direction = "x", q = Order, nodes = x_nodes )
        call Grid_Initialization( grid_spacing = "uniform", direction = "y", q = Order, nodes = y_nodes )

        U_an = Analytical_Solution(x_nodes, y_nodes)

        call Boundary_Value_Problem( x_nodes = x_nodes, y_nodes = y_nodes, Order = Order,   & 
                                     Differential_operator = Differential_operator,         &
                                     Boundary_conditions = Boundary_conditions,             &
                                     Solution = U ) 
        
        

        error = abs( norm2(U_an - U) )
        log_e(i) = log10( error ) 
        log_N(i) = log10 ( real(N) ) 
        
        write(*,*) 'N = ', N, 'error = ', error

        
        deallocate( x_nodes, y_nodes, U, U_an)

     end do  
end subroutine
















!*****     LATEX - 190     *****  
subroutine BVP_2D_System_Validation( Differential_operator, Boundary_conditions, Order, Nvariables, Spatial_Domain_x, Spatial_Domain_y, Analytical_Solution, log_e, log_N)
     
     procedure (DifferentialOperator2D_system) :: Differential_operator
     procedure (BC2D_system) ::  Boundary_conditions
     procedure (Analytical_Solution2D) :: Analytical_Solution
     integer, intent(in) :: Order, Nvariables
     real, intent(in) :: Spatial_Domain_x(:), Spatial_Domain_y(:)
     real, allocatable, intent(out) :: log_e(:), log_N(:)
     
      real :: x0, y0, xf, yf, error
     real, allocatable :: x_nodes(:), y_nodes(:),  U(:,:,:), U_an(:,:)
     
     integer :: i, j, N
    
     integer :: N0 = 20, Nf = 40
     integer :: gap = 5
     
     allocate( log_e(0:(Nf-N0)/gap), log_N(0:(Nf-N0)/gap) )
     
     x0 = Spatial_Domain_x(1)
     xf = Spatial_Domain_x( size(Spatial_Domain_x) )
     
     y0 = Spatial_Domain_y(1)
     yf = Spatial_Domain_y( size(Spatial_Domain_y) )
     
     
     do i = 0, (Nf-N0)/gap
         
         N = N0 + gap*i
        
        allocate ( x_nodes(0:N), y_nodes(0:N), U(0:N, 0:N, Nvariables), U_an(0:N, 0:N) )
        
        x_nodes(0) = x0; x_nodes(N) = xf
        y_nodes(0) = y0; y_nodes(N) = yf

        call Grid_Initialization( grid_spacing = "uniform", direction = "x", q = Order, nodes = x_nodes )
        call Grid_Initialization( grid_spacing = "uniform", direction = "y", q = Order, nodes = y_nodes )

        U_an = Analytical_Solution(x_nodes, y_nodes)
        
            call Boundary_Value_Problem( x_nodes = x_nodes, y_nodes = y_nodes, Order = Order,   & 
                                 N_variables = Nvariables,                          &
                                 Differential_operator = Differential_operator,     & 
                                 Boundary_conditions = Boundary_conditions,         &
                                 Solution = U )
        
        error = abs( norm2(U_an - U(:,:,1)) )
        log_e(i) = log10( error ) 
        log_N(i) = log10( real(N) ) 
        
        write(*,*) 'N = ', N, 'error = ', error

        
        deallocate( x_nodes, y_nodes, U, U_an)

     end do  
end subroutine



subroutine crea_name(name_fichero_, numero, extension, name_creado)

     character,intent(in)    :: name_fichero_*(*)
     character,intent(in)    :: extension*(*)
     integer, intent(in)     :: numero
     character, intent(out)  :: name_creado*(*)

     integer     :: numero_caracteres
     character   :: indice_texto*2

     if (numero == 0) then
         name_creado = trim(name_fichero_)//extension
     else
         call intcha(numero, numero_caracteres, indice_texto)
         indice_texto = adjustl(indice_texto)
         name_creado = trim(name_fichero_)//&
                         indice_texto(1:numero_caracteres)//extension
     endif


endsubroutine


!******************************************
!
!           VALIDATION
!
!*******************************************

!LINEAR 1D
subroutine Save_validation_1D_linear(x, y, N, q)
     real, intent(in) :: x(:), y(:)
    integer, intent(in) :: N, q

    integer :: i
    character(len=100) ::  name, names
    character(len=4) :: extension = '.txt'
    
    name = 'BVP1D_Validation_'
    
    call crea_name(trim(name), q, extension, names)
    
    open(62, file = names)
    
        do i = 1, N
            write (62,*) x(i), y(i)
        end do
        
    close (62)
    
    write(*,*) 'Saved as: ', names


end subroutine




subroutine Save_validation_2D(x, y, N, q)
    real, intent(in) :: x(:), y(:)
    integer, intent(in) :: N, q

    integer :: i
    character(len=100) ::  name, names
    character(len=4) :: extension = '.txt'
    
    name = 'BVP2D_Validation_'
    
    
    call crea_name(trim(name), q, extension, names)
    
    open(62, file = names)
    
        do i = 1, N
            write (62,*) x(i), y(i)
        end do
        
    close (62)
    
    write(*,*) 'Saved as: ', names


end subroutine


!LINEAR 2D System
subroutine Save_validation_2D_System(x, y, N, q)
    real, intent(in) :: x(:), y(:)
    integer, intent(in) :: N, q

    integer :: i
    character(len=100) ::  name, names
    character(len=4) :: extension = '.txt'
    
    name = 'BVP2DS_Validation_'
       
    call crea_name(trim(name), q, extension, names)
    
    open(62, file = names)
    
        do i = 1, N
            write (62,*) x(i), y(i)
        end do
        
    close (62)
    
    write(*,*) 'Saved as: ', names

end subroutine


end module 




