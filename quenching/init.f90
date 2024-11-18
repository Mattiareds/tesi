subroutine init

!This subroutine initializes velocities to start the MD run

USE PARAMETERS  !uso il modulo di definizione dei parametri
USE CLUSTER     !uso il modulo dove defin1:natomco variabili e parametri cluster

Implicit None

!Real, Intrinsic :: dble 
! variabili locali
Integer :: is
Real*8 :: ecinetica
Real*8 :: vxsum,vysum,vzsum
Real*8 :: vxsum2,vysum2,vzsum2
Real*8 :: kchevorrei,kcheho,costante

ecinetica=0.d0
!simulation starts at temperature tinit
!write(*,*)'INIT, tinit=',tinit
!temperatura=tinit

vx(1:natom)=0.d0
vy(1:natom)=0.d0
vz(1:natom)=0.d0
dx(1:natom)=0.d0                                                       
dy(1:natom)=0.d0                                                       
dz(1:natom)=0.d0                                                       
fx(1:natom)=0.d0                                                       
fy(1:natom)=0.d0                                                       
fz(1:natom)=0.d0                                                      
dfx(1:natom)=0.d0                                                       
dfy(1:natom)=0.d0                                                       
dfz(1:natom)=0.d0 

!random extraction of velocities
do is=nfix+1,natom
   call ranmar(zzz)
   vx(is)=dble(2.d0*zzz-1.d0) 
   call ranmar(zzz)
   vy(is)=dble(2.d0*zzz-1.d0)  
   call ranmar(zzz)
   vz(is)=dble(2.d0*zzz-1.d0)
enddo

!drift velocity
   vxsum=0.d0
   vysum=0.d0
   vzsum=0.d0
   do is=1,natom
      vxsum=vxsum+vx(is)
      vysum=vysum+vy(is)
      vzsum=vzsum+vz(is)
   enddo
   vxsum=vxsum/dfloat(natom)
   vysum=vysum/dfloat(natom)
   vzsum=vzsum/dfloat(natom)
   
   write(*,*)'in INIT the initial drift is: '
   write(*,*)'vxsum,    vysum,    vzsum [m/s]'
   write(*,'(3f10.6)')vxsum,vysum,vzsum
   write(*,*)
   
   !removing the drift velocity:
   do is=nfix+1,natom     
      vx(is)=vx(is)-vxsum
      vy(is)=vy(is)-vysum
      vz(is)=vz(is)-vzsum
   enddo
   
vxsum2=0.d0
vysum2=0.d0
vzsum2=0.d0

!kinetic energy from ideal temperature T distribution:
kchevorrei=1.5d0*natom*cbol*evsujoule*temperatura
!write(*,*)'kchevorrei',kchevorrei

!kinetic energy from random velocities (units of random v. are m/s)
do is=nfix+1,natom
   vxsum2=vxsum2+vx(is)*vx(is)*mass(itype(is))*uasukg
   vysum2=vysum2+vy(is)*vy(is)*mass(itype(is))*uasukg
   vzsum2=vzsum2+vz(is)*vz(is)*mass(itype(is))*uasukg
enddo

kcheho=0.5d0*(vxsum2+vysum2+vzsum2)
!scaling factor from random velocities to velocity distribution@temperatura
costante=sqrt(kchevorrei/kcheho)

!new velocities [m/s]
do is=nfix+1,natom
   vx(is)=vx(is)*costante
   vy(is)=vy(is)*costante
   vz(is)=vz(is)*costante
enddo

!new velocities [A/s]
!write(*,*)'INIT'
do is=nfix+1,natom
   vx(is)=vx(is)/angsum
   vy(is)=vy(is)/angsum
   vz(is)=vz(is)/angsum
   !vx(is)=0.
   !vy(is)=0.
   !vz(is)=0.
enddo

vxsum=0.d0
vysum=0.d0
vzsum=0.d0
do is=nfix+1,natom
   vxsum=vxsum+vx(is)
   vysum=vysum+vy(is)
   vzsum=vzsum+vz(is)
enddo

End Subroutine init
