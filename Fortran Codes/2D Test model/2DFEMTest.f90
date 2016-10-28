Program FEM2DTest
	!this program is written to test 2-D deformable elements
	!it is an initial model to test effects of plasma on thruster structures
	integer i,j,k
	! the start and the end point of the grid 
	real*8 x_ini, x_end
	! the grid points array
	real*8, allocatable:: x_grid(:)
	! the number of points in grid 
	integer N
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                                                     !
	!     2D Spectral Element Model Simulation Parameter  !
	!                                                     !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	N=9
	x_ini=0.1
	x_end=0.9
	
	allocate(x_grid(N))
	! let's generate a sample grid  
	call	OneDGrid(x_ini,x_end,x_grid,N)
	
	print*, "Sample Grid Points"
	print*, x_grid
	
End Program FEM2DTest	