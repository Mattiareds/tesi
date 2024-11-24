Subroutine growth

USE CLUSTER
USE PARAMETERS

Implicit none

Integer :: i
Real    :: xcentr_sum, ycentr_sum, zcentr_sum, ddist, ddist_tmp
Real    :: gsphere_radius, theta, phi, new_atom_mass, mod_v
Real    :: fattore_conversione
Real, parameter :: gradius=6.


if(ipas .eq. nterm+growth_steps)then
	!This is the first call to growth:
	open(123,file='growth.out',status='unknown')
	write(123,*)'***********************'
	write(123,*)'Growth at ipas = ',ipas
else
	open(123,file='growth.out',status='old',position='append')
	write(123,*)'***********************'
	write(123,*)'Growth at ipas = ',ipas

endif

natom = natom + 1
write(123,*)'It"s growth time: let"s add an atom. Natom is now = ',natom

!Choice of the atom type:
call ranmar(zzz)
if(zzz .lt. growth_p1)then
   elem(natom)=elem1
   itype(natom)=1
   new_atom_mass=mass(1)
else
   elem(natom)=elem2
   itype(natom)=2
   new_atom_mass=mass(2)
endif
write(123,*)'The new atom is a ',elem(natom),' atom.'

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
xcentr_sum=0.
ycentr_sum=0.
zcentr_sum=0.
do i=1,natom-1
   xcentr_sum=xcentr_sum+x(i)+u(i)
   ycentr_sum=ycentr_sum+y(i)+v(i)
   zcentr_sum=zcentr_sum+z(i)+w(i)
enddo
xcentr_sum=xcentr_sum/dfloat(natom-1)
ycentr_sum=ycentr_sum/dfloat(natom-1)
zcentr_sum=zcentr_sum/dfloat(natom-1)
write(123,*)'Cluster com coordinates: ',xcentr_sum, ycentr_sum, zcentr_sum

!How far from the com is the farthest atom of the cluster?
ddist=0.
do i=1,natom-1
   ddist_tmp=((x(i)+u(i))-xcentr_sum)**2+((y(i)+v(i))-ycentr_sum)**2+((z(i)+w(i))-zcentr_sum)**2
   if(ddist_tmp .gt. ddist)ddist=ddist_tmp
enddo
ddist=sqrt(ddist)
write(123,*)'The largest distance between the com and the cluster atoms is: ',ddist

!Setting the radius of the growth sphere
gsphere_radius=ddist+gradius 
write(123,*)'The radius of the deposition sphere is thus = ',ddist,'+',gradius,'=',gsphere_radius
!Extracting a random position on the sphere of radius gsphere_radius
call ranmar(zzz)
theta=pi*zzz
call ranmar(zzz)
phi=2*pi*zzz
write(123,*)'Random point on the sphere has r, theta, phi = ', gsphere_radius, theta, phi
!Assigning coordinates to the new atom:
x(natom)=xcentr_sum+gsphere_radius*sin(theta)*sin(phi)
y(natom)=ycentr_sum+gsphere_radius*sin(theta)*cos(phi)
z(natom)=zcentr_sum+gsphere_radius*cos(theta)
u(natom)=0.
v(natom)=0.
w(natom)=0.


!Let's assign a velocity to the newcomer:


!************************************
!!TO ASSIGN ALWAYS THE SAME VEL:
!!kB = Boltzmann constant = 1.3806E-23
!!1 amu = 1.660538 E-27 kg
!!mod_v = sqrt(3*kB*growth_T/m)
!! 1 A/s = 10^10 m/s
!mod_v=10E10*sqrt(3*1.3806*growth_temp*(10**4)/new_atom_mass)


!************************************
!!TO EXTRACT VEL FROM A GAUSSIAN DISTRIB
fattore_conversione=sqrt(evsujoule/uasukg)/angsum
call gauss   
vx(natom)=g1*fattore_conversione*sqrt((cbol*temperatura)/mass(itype(natom)))
vy(natom)=g2*fattore_conversione*sqrt((cbol*temperatura)/mass(itype(natom)))
vz(natom)=g3*fattore_conversione*sqrt((cbol*temperatura)/mass(itype(natom)))
mod_v=sqrt(vx(natom)**2+vy(natom)**2+vz(natom)**2)


!Calculating the x,y,z components of the velocity vector pointing to the center of mass:
vx(natom)=-mod_v*sin(theta)*sin(phi)
vy(natom)=-mod_v*sin(theta)*cos(phi)
vz(natom)=-mod_v*cos(theta)

write(123,*)'The new atom has a velocity = ',mod_v,' A/s'
!write(*,*)'For comparison, the atom natom-1 has a velocity = ',sqrt(vx(natom-1)**2+vy(natom-1)**2+vz(natom-1)**2),'A/s'
!write(*,*)'For comparison, the atom natom-2 has a velocity = ',sqrt(vx(natom-2)**2+vy(natom-2)**2+vz(natom-2)**2),'A/s'
!write(*,*)'For comparison, the atom natom-3 has a velocity = ',sqrt(vx(natom-3)**2+vy(natom-3)**2+vz(natom-3)**2),'A/s'
!write(*,*)'For comparison, the atom natom-4 has a velocity = ',sqrt(vx(natom-4)**2+vy(natom-4)**2+vz(natom-4)**2),'A/s'
!write(*,*)'For comparison, the atom natom-5 has a velocity = ',sqrt(vx(natom-5)**2+vy(natom-5)**2+vz(natom-5)**2),'A/s'
!write(*,*)'For comparison, the atom natom-6 has a velocity = ',sqrt(vx(natom-6)**2+vy(natom-6)**2+vz(natom-6)**2),'A/s'
!write(*,*)'For comparison, the atom natom-7 has a velocity = ',sqrt(vx(natom-7)**2+vy(natom-7)**2+vz(natom-7)**2),'A/s'
!write(*,*)'For comparison, the atom natom-8 has a velocity = ',sqrt(vx(natom-8)**2+vy(natom-8)**2+vz(natom-8)**2),'A/s'

close(123)


End subroutine growth