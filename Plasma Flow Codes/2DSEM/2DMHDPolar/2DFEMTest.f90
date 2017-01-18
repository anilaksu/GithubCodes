
Program FEM2DTest
	!this program is written to test 2-D deformable elements
	!it is an initial model to test effects of plasma on thruster structures
	integer i,j,k,jstart,jend
	! the start and the end point of the grid 
	real*8 x_ini, x_end
	! the grid points array and corresponding weight
	real*8, allocatable:: x_grid(:),w(:),x_total(:)
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                                                     !
	!     Independent Variables  (Time, Coordinate)       !
	!                                                     !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! the test grid in coordinate transformation
	real*8, allocatable:: x_trans(:,:),x_surf(:,:),x_test(:,:),x_2Dgrid(:,:)
	! the time array 
	real*8, allocatable:: SimTime(:)
	! time step size
	real*8 dt
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                                                     !
	!               Dependent Variables                   !
	!                                                     !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! the velocity field
	real*8, allocatable::u(:,:),v(:,:)
	! the auxilary velocity field
	real*8, allocatable::ua(:),va(:)
	! the convective terms
	real*8, allocatable::uux(:,:),vuy(:,:),uvx(:,:),vvy(:,:)
	! the filtered convective terms
	real*8, allocatable::fuux(:,:),fvuy(:,:),fuvx(:,:),fvvy(:,:)
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                                                     !
	!         The dimensions of dependent and 			  !
	!             independent variables                   !
	!                                                     !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! the number of points in grid and the order of legendre polynomial
	integer N,Nl,Ngrid
	! the number of subgrids and the number of time steps
	integer Nsub,Ntime
	! the length of each subsomain
	real*8, allocatable:: delx(:)
	! the number of points in each subgrid in x and y direction
	integer Nx,Ny
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                                                     !
	!                Material Properties 			      !
	!                                                     !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! the kinematic viscosity and dynamic viscosity
	real*8 xnu, xmu
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                                                     !
	!                Numerical Operators			      !
	!                                                     !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! 1-D Discontinuous SEM Laplacian Operator and Boundary Condition Matrix
	real*8, allocatable:: LapDis1D(:,:),SEMBCMatrix1D(:,:)
	! the identity matrix 
	real*8, allocatable:: ID(:,:),ID2D(:,:)
	! 2-D Laplacian Operator and Boundary Condition Matrix
	real*8, allocatable:: Lap2D(:,:),BCMatrix2D(:,:)
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                                                     !
	!        2-D Grid and related parameters			  !
	!                                                     !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! the grid points
	real*8, allocatable:: x_grid2D(:,:),wx(:),wy(:)
	
	
	! the number of points on the surface
	integer, allocatable::Nsurf(:)
	! legendre polynomials and associated q
	real*8, allocatable:: Ln(:),Lpn(:),q(:),qp(:)
	! Test matrix 
	real*8, allocatable:: A(:,:),Ainv(:,:)
	! Boundary condition, Governing Equation and System Matrix and inverse of the system matrix
	real*8, allocatable:: BCMatrix(:,:),GEMatrix(:,:),SysMatrix(:,:),InvSysMatrix(:,:)
	! Test Function
	real*8, allocatable:: F(:),dF(:),FB(:),LegenData(:)
	! the right hand side of the equation
	real*8, allocatable:: RHS(:)
	!complex*8, allocatable:: F(:),FB(:)
	complex,allocatable:: FourF(:)
	integer ssign
	! the differentiation matrix
	real*8, allocatable:: D1(:,:)
	! the laplacian matrix 
	real*8, allocatable:: LaplaceMatrix(:,:)
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                                                     !
	!                  Useful Constants                   !
	!                                                     !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! pi number
	real*8 pi
	! for complex exponential representation
	complex phi
	pi=4.*datan(1.d0)
	phi=CMPLX(0,2.*pi)
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                                                     !
	!     2D Spectral Element Model Simulation Parameter  !
	!                                                     !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	!! test grid parameters
	!Nsurf=40
	Ngrid=20
		
	N=32
	Nl=2
	x_ini=-1.
	x_end=1.
	
	!! let's allocate location arrays
	allocate(x_grid(Ngrid))
	allocate(x_2Dgrid(Ngrid*Ngrid,2))
	allocate(LaplaceMatrix(Ngrid*Ngrid,Ngrid*Ngrid))
	allocate(w(Ngrid))
	allocate(Ln(N))
	allocate(Lpn(N))
	allocate(q(N))
	allocate(qp(N)) 
	
	!! test matrix 
	!allocate(A(N,N))
	!allocate(Ainv(N,N))
	!! test function
	!allocate(F(N))
	allocate(F(Ngrid*Ngrid))
	allocate(dF(Ngrid*Ngrid))
	!allocate(FourF(N))
	!allocate(LegenData(N+1))
	
	!! test grid allocation
	!allocate(x_test(Ngrid,2))
	!allocate(x_trans(Ngrid,2))
	!allocate(x_surf(Nsurf,2))
	
	!! first differentiation matrix 
	allocate(D1(Ngrid,Ngrid))
	

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                                                     !
	!             Sample Coordinate Mapping               !
	!                                                     !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	
	!do i=1,N
		!f(i)=x_grid(i)**2.
	!end do 
	!call OneDIntegration(int_f,f,w,N)
	!print*,"Sum of weights and integral of the function"
	!print*,sum(w),int_f
	
	
	!deallocate(A)
	!allocate(A(Ngrid,Nsurf))
	!call GLLPoints(x_test(:,2),w,Ngrid)write(139,*) i-1,ULeg(i)
	!do i=1,Nsurf
	!	x_surf(i,1)=2.*dcos(pi*(i-1)/Nsurf)
	!	x_surf(i,2)=2.*dsin(pi*(i-1)/Nsurf)
	!end do
	
	!call OneDTransform(x_trans,x_surf,x_test,Nsurf,Ngrid)
	
	!! test grid allocation
	
	
	!deallocate(x_test)
	!deallocate(x_trans)

	
	!allocate(Nsurf(4))
	!Nsurf=5
	
	!allocate(x_test(Ngrid*Ngrid,2))
	!allocate(x_trans(Ngrid*Ngrid,2))
	!allocate(x_surf(sum(Nsurf),2))
	
	!do i=1,4
	
	!	if(i==1)then
	 !   	jstart=1
	!	else
	!		jstart=sum(Nsurf(1:(i-1)))+1
	!	end if
		
	!	do j=jstart,sum(Nsurf(1:i))
	!		if(i==1) then
	!			x_surf(j,1)=0.5
	!			x_surf(j,2)=-0.5+0.25*(j-jstart)
	!		else if (i==2) then
	!				x_surf(j,1)=0.5-0.25*(j-jstart)
	!				x_surf(j,2)=0.5
	!			 else if (i==3) then
	!			 		x_surf(j,1)=-0.5
	!			 		x_surf(j,2)=0.5-0.25*(j-jstart)
	!			    	 x_surf(j,1)=-0.5+0.25*(j-jstart)
	!				 	 x_surf(j,2)=-0.5			 
	!			 end if
	!	end do
	!end do
	
	!do i=1,Nsurf(1)
	!	print*,x_surf(i,:)
	!end do
	!call TwoDMeshTransform(x_trans,x_surf,x_test,Nsurf,Ngrid)
	
	!open(128,file='1DMap.dat',status='unknown')
	!open(129,file='OMap.dat',status='unknown')
	!do i=1,Ngrid
	!		write(128,*) x_trans(i,1),x_trans(i,1)
	!		print*, x_test(i,:)
	!end do 
	!do i=1,Nsurf
	!	write(129,*) x_surf(i,1),x_surf(i,1)
	!end do 
	
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                                                     !
	!        1D Laplace Equation Solution on GLL          !
	!                                                     !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	!deallocate(F)
	!deallocate(dF)
	! let's allocate all the related matrices
	!allocate(BCMatrix(Ngrid,Ngrid))
	!allocate(GEMatrix(Ngrid,Ngrid))
	!allocate(SysMatrix(Ngrid,Ngrid))
	!allocate(InvSysMatrix(Ngrid,Ngrid))
	!allocate(SysMatrix(Ngrid-2,Ngrid-2))
	!allocate(InvSysMatrix(Ngrid-2,Ngrid-2))
	!allocate(FB(Ngrid))
	!allocate(F(Ngrid))
	!allocate(dF(Ngrid))
	!allocate(RHS(Ngrid))

	!do i=1,Ngrid
		!FB(i)=dsin(2*pi*x_grid(i))
		!FB(i)=dexp(-1.*x_grid(i)**2.)
		!FB(i)=x_grid(i)
		!if(i==1. .or. i==Ngrid) then
			!FB(i)=0.
		!else
			!FB(i)=1.
			!FB(i)=x_grid(i)
		!end if
		
	!end do
	
	!call OneDBC(BCMatrix,0,0,Ngrid)
	!call Lagder(GEMatrix,Ngrid)
	!call OneDBCApply(BCMatrix,GEMatrix,SysMatrix,Ngrid)
	!call invertSVD(Ngrid,SysMatrix,InvSysMatrix)
	!call OneDRHS(RHS,FB,Ngrid)
	!print*,RHS
	!F=0.
	!call matvect(InvSysMatrix,RHS,F,Ngrid)
	
	!print*,F
	
	!print*,GEMatrix(1,:)
	!print*,GEMatrix(2,:)
	!print*,GEMatrix(3,:)
	
	!open(135,file='1DAnalytical.dat',status='unknown')
	!open(136,file='1DNumerical.dat',status='unknown')
	!do i=1,Ngrid
	!	write(135,*) x_grid(i),FB(i)/(-4.*pi*pi)
	!	write(136,*) x_grid(i),F(i)
	!end do
	
	! let's generate the differentiation matrix
	!call FirstDiff(D1,Ngrid)
	! let's differentiate F
	!call matvect(D1,FB,dF,Ngrid)
	
	!open(131,file='Ffunction.dat',status='unknown')
	!open(132,file='dFdx.dat',status='unknown')
	!do i=1,Ngrid
		!print*,D1(i,:)
	!	write(131,*) x_grid(i),FB(i)
	!	write(132,*) x_grid(i),dF(i)
	!end do 
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                                                     !
	!          1-D Burger's Equation Solutiom             !
	!                                                     !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	! the number time steps
	Ntime=2000
	! the number of subdomains
	Nsub=2
	! the number of grid points within each domain
	Ngrid=11

	! each domain has same amount of points
	allocate(x_total(Nsub*Ngrid))
	deallocate(x_grid)
	deallocate(w)
	allocate(x_grid(Ngrid))
	allocate(w(Ngrid))
	allocate(delx(Nsub))
	allocate(SimTime(Ntime))
	! the velocity field
	allocate(u(Ntime,Nsub*Ngrid))
	allocate(ua(Nsub*Ngrid))
	! the convective terms
	allocate(uux(Ntime,Nsub*Ngrid))
	allocate(fuux(Ntime,Nsub*Ngrid))
	! numerical operators
	allocate(LapDis1D(Nsub*Ngrid,Nsub*Ngrid))
	allocate(SEMBCMatrix1D(Nsub*Ngrid,Nsub*Ngrid))
	allocate(SysMatrix(Nsub*Ngrid,Nsub*Ngrid))
	allocate(InvSysMatrix(Nsub*Ngrid,Nsub*Ngrid))
	allocate(ID(Nsub*Ngrid,Nsub*Ngrid))
	! let's generate GLL points
	call GLLPoints(x_grid,w,Ngrid)

	
	
	! let's fill the domain which from -1 to 1
	do i=1,Nsub
		do j=1,Ngrid
		 	x_total((i-1)*Ngrid+j)=1.*(i-2)+0.5*(x_grid(j)+1.)
		end do
	end do
	! the time step size
	dt=0.2
	! the domain lengths
	delx(1)=1.
	delx(2)=1.
	! the initial condition
	do i=1,Nsub*Ngrid
		u(1,i)=dsin(pi*x_total(i))
	end do
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                                                     !
	!             The Material Properties                 !
	!                                                     !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! the kinematic viscosity
	xnu=10E-1
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                                                     !
	!             Diffusive Operator in SEM               !
	!                                                     !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	call OneDDisLaplace(LapDis1D,delx,Nsub,Ngrid)
	call OneDSEMBC(SEMBCMatrix1D,Nsub,Ngrid)
	call OneDSEMBCApply(SEMBCMatrix1D,LapDis1D,SysMatrix,Nsub,Ngrid)
	call IdMatrix(ID,Nsub*Ngrid)
	!! it is given as this form as the system integrates viscous part implicitly
	SysMatrix=ID-dt*xnu*SysMatrix
	call invertSVD(Ngrid*Nsub,SysMatrix,InvSysMatrix)
	
	open(140,file='SysMatrix.dat',status='unknown')
	do i=1,Ngrid*Nsub
		write(140,*) SysMatrix(i,:)
	end do
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                                                     !
	!        Time Integration of The Equation             !
	!                                                     !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	SimTime(1)=0.
	! let's time integrate the equation
	do i=1,Ntime-1
		! let's calculate the nonlinear terms for each grid 
		do j=1,Nsub
			jstart=(j-1)*Ngrid+1
			jend=j*Ngrid
			call OneDNterms(uux(i,jstart:jend),u(i,jstart:jend),u(i,jstart:jend),Ngrid)
			call OneDGLLFilter(fuux(i,jstart:jend),uux(i,jstart:jend),x_grid,w,Ngrid)
		end do
		! let's find the time step
		!call TimeStepSize1D(dt,u(i,1:Ngrid),x_grid,Ngrid)
		! let's integrate the velocity in time
		!u(i+1,:)=u(i,:)+dt*uux(i,:)		
		! the integration of convective part 
		ua=u(i,:)+dt*fuux(i,:)
		! the integration of diffusive part
		call matvect(InvSysMatrix,ua,u(i+1,:),Ngrid*Nsub)
		!u(i+1,:)=ua
		SimTime(i+1)=SimTime(i)+dt
		!print*,SimTime(i+1),dt
	end do 
	
	
	open(137,file='UVelocity20.dat',status='unknown')
	open(138,file='UVelocity40.dat',status='unknown')
	open(139,file='UVelocity60.dat',status='unknown')
	do i=1,Ngrid*Nsub
		write(137,*) x_total(i),u(20,i)
		write(138,*) x_total(i),u(40,i)
		write(139,*) x_total(i),u(100,i)
	end do
	
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                                                     !
	!    2D Poisson Equation Solution on GLL-GLL Grid     !
	!                                                     !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	Nx=20
	Ny=20
	deallocate(F)
	deallocate(dF)
	! let's allocate all the related matrices
	allocate(BCMatrix2D(Nx*Ny,Nx*Ny))
	allocate(Lap2D(Nx*Ny,Nx*Ny))
	deallocate(SysMatrix)
	deallocate(InvSysMatrix)
	allocate(SysMatrix(Nx*Ny,Nx*Ny))
	allocate(InvSysMatrix(Nx*Ny,Nx*Ny))
	
	allocate(F(Nx*Ny))
	allocate(RHS(Nx*Ny))
	allocate(x_grid2D(Nx*Ny,2))
	allocate(wx(Nx))
	allocate(wy(Ny))
	
	
	do i=1,Ny
		do j=1,Nx
			if(j==1) then
				F((i-1)*Nx+j)=1.
			else
				F(i)=0.
			end if		
		end do
	end do
	
	call TwoDBC(BCMatrix2D,0,0,0,0,Nx,Ny)
	call Laplace(Lap2D,Nx,Ny)
	call TwoDBCApply(BCMatrix2D,Lap2D,SysMatrix,Nx,Ny)
	call invertSVD(Nx*Ny,SysMatrix,InvSysMatrix)
	call TwoDRHS(RHS,F,Nx,Ny)
	call matvect(InvSysMatrix,RHS,F,Nx*Ny)
	! let's generate the grid
	call TwoDMapXY(x_grid2D,wx,wy,Nx,Ny)

	open(140,file='2DLaplace.dat',status='unknown')
	open(141,file='2DNumerical.dat',status='unknown')
	do i=1,Nx*Ny
	!	write(140,*) BCMatrix2D(i,:)
		write(140,*) SysMatrix(i,:)
		write(141,*) x_grid2D(i,:),F(i)
	end do
		
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                                                     !
	!    2D Heat Equation Solution on GLL-GLL Grid     !
	!                                                     !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	!Nx=20
	!Ny=20
	!deallocate(F)
	!deallocate(dF)
	! the velocity field
	!deallocate(u)
	!deallocate(ua)
	! let's allocate all the related matrices
	!allocate(BCMatrix2D(Nx*Ny,Nx*Ny))
	!allocate(Lap2D(Nx*Ny,Nx*Ny))
	!deallocate(SysMatrix)
	!deallocate(InvSysMatrix)
	!allocate(SysMatrix(Nx*Ny,Nx*Ny))
	!allocate(InvSysMatrix(Nx*Ny,Nx*Ny))
	! the velocity field
	!allocate(u(Ntime,Nx*Ny))
	!allocate(ua(Nx*Ny))
	!allocate(F(Nx*Ny))
	!allocate(RHS(Nx*Ny))
	!allocate(x_grid2D(Nx*Ny,2))
	!allocate(wx(Nx))
	!allocate(wy(Ny))
	!allocate(ID2D(Nx*Ny,Nx*Ny))
	
	!do i=1,Nx*Ny
	!	if(i<=Nx) then
	!		u(1,i)=1.
	!	else
	!		u(1,i)=0.
	!	end if		
	!end do
	
	!call TwoDIDMatrix(ID2D,Nx,Ny)
	!call TwoDBC(BCMatrix2D,0,0,0,1,Nx,Ny)
	!call Laplace(Lap2D,Nx,Ny)
	!SysMatrix=ID2D-dt*xnu*Lap2D
	!call TwoDBCApply(BCMatrix2D,SysMatrix,Lap2D,Nx,Ny)
	!call invertSVD(Nx*Ny,Lap2D,InvSysMatrix)
	
	! let's generate the grid
	!call TwoDMapXY(x_grid2D,wx,wy,Nx,Ny)


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                                                     !
	!        Time Integration of The Equation             !
	!                                                     !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!SimTime(1)=0.
	! let's time integrate the equation
	!do i=1,Ntime-1
		! the integration of diffusive part
	!	call matvect(InvSysMatrix,u(i,:),u(i+1,:),Nx*Ny)
		!u(i+1,:)=ua
	!	SimTime(i+1)=SimTime(i)+dt
		!print*,SimTime(i+1),dt
	!end do 
	
	!open(140,file='2DLaplace.dat',status='unknown')
	!open(141,file='2DNumerical.dat',status='unknown')
	!do i=1,Nx*Ny
	!	write(140,*) BCMatrix2D(i,:)
	!	write(140,*) InvSysMatrix(i,:)
	!	write(141,*) x_grid2D(i,:),u(Ntime,i)
	!end do
			
End Program FEM2DTest	