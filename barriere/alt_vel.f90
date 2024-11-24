sUBROUTINE VEL

USE PARAMETERS  !uso il modulo di definizione dei parametri
USE CLUSTER     !uso il modulo dove defin1:natomco variabili e parametri cluster

Implicit None
Real(8) :: dlaran
Real*8 :: qx(nmax),qy(nmax),qz(nmax)
Real*8 :: qxx(nmax),qyy(nmax),qzz(nmax)
Real*8 :: pcmx,pcmy,pcmz
Real*8 :: versx,versy,versz,versmod,velmod,vel2
Real*8 :: poguy,vargau
Real*8 :: vdfa
Real*8 :: vdfae,fattore_conversione
Integer :: is
Integer :: ir(4)
Integer :: i

!la frequenza di termostatazione e' fissata in module_parameters.f90
!vdfa e' un numero 0<vdfa<1 che serve per far intervenire il termostato
vdfa=vnu*tstep

ecinet=0.d0
poguy=0.d0
pcmx=0.d0
pcmy=0.d0
pcmz=0.d0




pot=0.d0
do i=nfix+1,natom 
   pot=pot+fx(i)*vx(i)+fy(i)*vy(i)+fz(i)*vz(i) 
enddo

do i=nfix+1,natom  
   fattore_conversione=evsujoule/(uasukg*angsum*angsum)
   !le df sono quelle vecchie, le f sono state appena calcolate da force
   !le forze sono espresse in eV/A, vanno convertite in J/m
   !a questo punto il secondo addendo va convertito in da m/s a A/s
   if (pot.lt.0) then
      qx(i)=0.d0
   else
      qx(i)=vx(i)+fattore_conversione*0.5d0*(fx(i)+dfx(i))*tstep/mass(itype(i)) 
   endif
   if (pot.lt.0) then
      qy(i)=0.d0
   else
      qy(i)=vy(i)+fattore_conversione*0.5d0*(fy(i)+dfy(i))*tstep/mass(itype(i)) 
   endif
   if (pot.lt.0) then
      qz(i)=0.d0
   else
      qz(i)=vz(i)+fattore_conversione*0.5d0*(fz(i)+dfz(i))*tstep/mass(itype(i)) 
   endif
   vx(i)=qx(i) 
   vy(i)=qy(i)
   vz(i)=qz(i)      
   !in A*A/s*s   
   vel2=vx(i)*vx(i)+vy(i)*vy(i)+vz(i)*vz(i)
   !l-energia cinetica la esprimo in eV
   ecinet=ecinet+(0.5d0*mass(itype(i))*uasukg*vel2*angsum*angsum)/evsujoule
enddo 




!temp e' la temperatura calcolata sulla base delle velocita':
temp=2.d0*ecinet/(3.d0*(natom-nfix)*cbol) 

!tutte le enrgie a questo punto sono espresse in eV:
etot=ener+ecinet 
!write(*,*)etot
!open(102,file='e2.dat',status='unknown',position='append')
!write(102,'(i10,3f16.8)')ipas,ecinet,etot,temp
!close(102)
!write(*,*)'In the routine vel.f90 natom=',natom
!write(*,*)'In the routine vel.f90 the last atom has velocity: ',sqrt(vel2)

!if (ipas .eq. 301) stop

END SUBROUTINE VEL
        
        
        
