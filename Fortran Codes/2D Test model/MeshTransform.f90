!!! This file contains all grid generation operations including:
!!! 1) 1D GLL and Fourier Points
!!! 2) 2D GLL-GLL, GLL-Fourier, Fourier-Fourier points
!!! 3) Coordinate Mapping functions like Numerical Jacobian Matrix calculations
subroutine OneDGrid(x_ini,x_end,x_grid,N)
	integer i
	! the number of grid points
	integer, intent(in):: N
	!the start and end coordinates of mesh
	real*8, intent(in):: x_ini,x_end
	!the grid coordinate array
	real*8, intent(inout),dimension(N):: x_grid
	! the length of domain and the uniform distance between points
	real*8 domLength, dx
	
	! the domain length is calculated as
	domLength=x_end -x_ini
	! let's calculate the uniform distance between points
	dx=domLength/(N-1)
	! first point of grid starts with x_ini
	x_grid(1)=x_ini
	
	do i=1,(N-1)
		x_grid(i+1)=x_grid(i)+dx
		!print*, x_grid(i+1)
	end do
		
end subroutine OneDGrid