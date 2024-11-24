Subroutine growth

USE CLUSTER
USE PARAMETERS

Implicit none

Integer :: i
Real*8    :: xcentr_sum,ycentr_sum,zcentr_sum,ddist,ddist_tmp
Real*8    :: gsphere_radius,theta,phi,mod_v
Real*8    :: fattore_conversione
Real*8, parameter :: gradius=6.d0


if(ipas .eq. growth_steps)then
	!This is the first call to growth:
  open(23,file='growth.out',status='unknown')
  write(23,*)'***********************'
  write(23,*)'Growth at ipas = ',ipas
else
  open(23,file='growth.out',status='old',position='append')
  write(23,*)'***********************'
  write(23,*)'Growth at ipas = ',ipas
endif

natom = natom+1
write(23,*)'Growth time: an atom is added' 
write(23,*) 'natom is now = ',natom

!Choice of the atom type:
call ranmar(zzz)

spec(natom)=elem(nelem)
itype(natom)=nelem
do i=1,nelem
   if((zzz.ge.sum_p1(i-1)).and.(zzz.lt.sum_p1(i))) then
      spec(natom)=elem(i)
      itype(natom)=i
   endif
enddo
write(23,*)'The new atom is a ',spec(natom),' atom.'

!Setting forces acting on natom to be = 0
fx(natom)=0.
fy(natom)=0.
fz(natom)=0.
dfx(natom)=0.
dfy(natom)=0.
dfz(natom)=0.
!Setting velocities of natom to be = 0
vx(natom)=0.
vy(natom)=0.
vz(natom)=0.


!Looking for the center of mass of the cluster
xcentr_sum=0.d0
ycentr_sum=0.d0
zcentr_sum=0.d0
do i=1,natom-1
   xcentr_sum=xcentr_sum+x(i)
   ycentr_sum=ycentr_sum+y(i)
   zcentr_sum=zcentr_sum+z(i)
enddo
xcentr_sum=xcentr_sum/dfloat(natom-1)
ycentr_sum=ycentr_sum/dfloat(natom-1)
zcentr_sum=zcentr_sum/dfloat(natom-1)
write(123,*)'Cluster com coordinates: ',xcentr_sum, ycentr_sum, zcentr_sum

!How far from the com is the farthest atom of the cluster?
ddist=0.d0
do i=1,natom-1
   ddist_tmp=(x(i)-xcentr_sum)**2+(y(i)-ycentr_sum)**2+(z(i)-zcentr_sum)**2
   if (ddist_tmp.gt.ddist) ddist=ddist_tmp
enddo
ddist=dsqrt(ddist)
write(23,*)'The largest distance between the com and the cluster atoms is: ',ddist

!Setting the radius of the growth sphere
gsphere_radius=ddist+gradius 
write(23,*)'The radius of the deposition sphere is thus = ',ddist,'+',gradius,'=',gsphere_radius
!Extracting a random position on the sphere of radius gsphere_radius
call ranmar(zzz)
theta=pi*zzz
call ranmar(zzz)
phi=2*pi*zzz
write(123,*)'Random point on the sphere has r, theta, phi = ', gsphere_radius, theta, phi
!Assigning coordinates to the new atom:
x(natom)=xcentr_sum+gsphere_radius*dsin(theta)*dsin(phi)
y(natom)=ycentr_sum+gsphere_radius*dsin(theta)*dcos(phi)
z(natom)=zcentr_sum+gsphere_radius*dcos(theta)

!Let's assign a velocity to the newcomer:


!!TO EXTRACT VEL FROM A GAUSSIAN DISTRIB
fattore_conversione=dsqrt(evsujoule/uasukg)/angsum
call gauss   
vx(natom)=g1*fattore_conversione*dsqrt((cbol*temperatura)/mass(itype(natom)))
vy(natom)=g2*fattore_conversione*dsqrt((cbol*temperatura)/mass(itype(natom)))
vz(natom)=g3*fattore_conversione*dsqrt((cbol*temperatura)/mass(itype(natom)))
mod_v=dsqrt(vx(natom)**2+vy(natom)**2+vz(natom)**2)


!Calculating the x,y,z components of the velocity vector pointing to the center of mass:
vx(natom)=-mod_v*dsin(theta)*dsin(phi)
vy(natom)=-mod_v*dsin(theta)*dcos(phi)
vz(natom)=-mod_v*dcos(theta)

write(23,*)'The new atom has a velocity = ',mod_v,' A/s'

close(23)


End subroutine growth
