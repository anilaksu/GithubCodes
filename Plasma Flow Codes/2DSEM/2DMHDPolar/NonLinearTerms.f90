!! this routine generates 1-D and 2-D nonlinear convective terms 
subroutine OneDNterms(UUx,U,Ux,N)
	!! This function returns nonlinear convective terms in 1D
	integer i,j,k 
	! the number of grid points
	integer, intent(in):: N
	!the differentiation matrix
	real*8, intent(inout),dimension(N):: U,Ux,UUx
	! the derivative of UX
	real*8, allocatable::dUdx(:)
	!! Differentiation Matrix 
	real*8,allocatable:: D1(:,:)
	! 1D GLL points and corresponding weight
	allocate(dUdx(N))
	! the allocation of first derivative matrix 
	allocate(D1(N,N))
	! let's first generate differentiation matrix
	call FirstDiff(D1,N)
	call matvect(D1,UX,dUdx,N)
	do i=1,N
		UUx(i)=U(i)*dUdx(i)
	end do
	
	!call OneDGLLFilter(UF,U,x_grid,w,N)
	
end subroutine OneDNterms

subroutine TimeStepSize1D(dt,U,x_grid,N)
	!! This function returns step size under CFL condition
	integer i
	! the number of grid points
	integer, intent(in):: N
	!the velocity and the grid
	real*8, intent(inout),dimension(N):: U,x_grid
	! the time step
	real*8, intent(inout):: dt
	! the distance between each grid
	real*8, allocatable::dx(:)
	! the maximum values of dx and U
	real*8 dx_max,U_max
	
	allocate(dx(N-1))
	! let's calculate dx
	do i=1,N-1
		dx(i)=abs(x_grid(i+1)-x_grid(i))
	end do
	!print*,U
	! the maximum values 
	dx_max=MAXVAL(dx)
	U_max=MAXVAL(abs(U))
	!print*,dx_max,U_max
	dt=0.5*dx_max/U_max
	
end subroutine TimeStepSize1D
