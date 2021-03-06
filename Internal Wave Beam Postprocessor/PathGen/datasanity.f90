Program DataSanity 

 
  integer i
! the size of the character file 
  integer ::M
! the size of the coordinate Pressure and Velocity Array
  integer :: nxr,nzr
! the character where the file names kept
  CHARACTER(20),allocatable:: pfiles(:),vfiles(:)
! pressure velocity 
  real*8,allocatable::Press(:),u(:),v(:)
! the coordinate data
  real*8,allocatable::x1(:,:),z1(:,:)

M=128

nxr=512

nzr=491

allocate(pfiles(M))

allocate(vfiles(M))

allocate(Press(nxr*nzr))

allocate(u(nxr*nzr))

allocate(v(nxr*nzr))

allocate(x1(nxr,nzr))

allocate(z1(nxr,nzr))

!print*,"okey"
 call pread(pfiles,M)

 call vread(vfiles,M)

 call DataRead(pfiles(1),vfiles(1),Press,u,v,nxr,nzr,x1,z1)

do i=1,nzr
print*,Press(i),u(i),v(i)
end do 

contains

subroutine pread(pfiles,M)
 integer i
! the size of the character file 
  integer,intent(in)::M
! the character where the file names kept
  CHARACTER(20),dimension(M),intent(inout):: pfiles

! here I generate File names and open it as well



open(15,file='depot.dat',status='unknown')

do i=1,M
if (i<10) then 
write(15,210) 'TPress_',i,'.dat'
210 FORMAT(A12,I1,A4)
goto 300
end if 

if (i<100) then 
write(15,220) 'TPress_',i,'.dat'
220 FORMAT(A12,I2,A4)
goto 300
end if 

if (i<1000) then 
write(15,230) 'TPress_',i,'.dat'
230 FORMAT(A12,I3,A4)
goto 300
end if 

300 continue
end do

 close(15)
open(16,file='depot.dat',status='old')
do i=1,M
read(16,*) pfiles(i)
end do

 close(16)

end subroutine pread 

subroutine vread(pfiles,M)
 integer i
! the size of the character file 
  integer,intent(in)::M
! the character where the file names kept
  CHARACTER(20),dimension(M),intent(inout):: pfiles

! here I generate File names and open it as well



open(15,file='depot.dat',status='unknown')

do i=1,M
if (i<10) then 
write(15,210) 'TVelocity_',i,'.dat'
210 FORMAT(A12,I1,A4)
goto 300
end if 

if (i<100) then 
write(15,220) 'TVelocity_',i,'.dat'
220 FORMAT(A12,I2,A4)
goto 300
end if 

if (i<1000) then 
write(15,230) 'TVelocity_',i,'.dat'
230 FORMAT(A12,I3,A4)
goto 300
end if 

300 continue
end do

 close(15)
open(16,file='depot.dat',status='old')
do i=1,M
read(16,*) pfiles(i)
end do

 close(16)

end subroutine vread 

subroutine DataRead(pfile,vfile,Press,u,v,nxr,nzr,x1,z1)
  integer i,j
! the size of the data in x and y direction
  integer,intent(in):: nxr,nzr
! Pressure and Velocity Data
  real*8,dimension(nxr*nzr),intent(inout)::Press,u,v
! the coordinate data
  real*8 ,dimension(nxr,nzr),intent(inout)::x1,z1

! the character where the file names kept
  CHARACTER(20),intent(inout):: pfile ,vfile





open(20,file=pfile,status='unknown')
  
 open(21,file=vfile,status='unknown') 

  do j=1,nzr

  do i=1,nxr

  read(20,*) x1(i,j),z1(i,j),Press((j-1)*nxr+i)
  
  read(21,*) x1(i,j),z1(i,j),u((j-1)*nxr+i),v((j-1)*nxr+i)
  
  end do 

  end do




end subroutine DataRead


end Program DataSanity 
