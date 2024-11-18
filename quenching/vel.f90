sUBROUTINE VEL

USE PARAMETERS  !uso il modulo di definizione dei parametri
USE CLUSTER     !uso il modulo dove defin1:natomco variabili e parametri cluster

Implicit None

Real*8 :: qx(nmax),qy(nmax),qz(nmax)
Real*8 :: pcmx,pcmy,pcmz
Real*8 :: vdfa,vel2
Real*8 :: vdfae,fattore_conversione
Integer :: last_atom_2bthermalized
Integer :: i

!la frequenza di termostatazione e' fissata in module_parameters.f90
!vdfa e' un numero 0<vdfa<1 che serve per far intervenire il termostato

vdfa=vnu*tstep 

ecinet=0.d0
pcmx=0.d0
pcmy=0.d0
pcmz=0.d0
!fattore_conversione=evsujoule/(uasukg*angsum*angsum)
!write (*,*) fattore_conversione
!stop

   !le df sono quelle vecchie, le f sono state appena calcolate da force
   !le forze sono espresse in eV/A, vanno convertite in J/m
   !a questo punto il secondo addendo va convertito in da m/s a A/s

if (choice_sim.eq.'qq') then
   fattore_conversione=evsujoule/(uasukg*angsum*angsum)
   do i=nfix+1,natom 
      if (tempo.lt.1000*timelag) then
         qx(i)=vx(i)+fattore_conversione*0.5d0*(fx(i)+dfx(i))*tstep/mass(itype(i)) 
         qy(i)=vy(i)+fattore_conversione*0.5d0*(fy(i)+dfy(i))*tstep/mass(itype(i)) 
         qz(i)=vz(i)+fattore_conversione*0.5d0*(fz(i)+dfz(i))*tstep/mass(itype(i))   
      else 
         if (fx(i)*vx(i).lt.0) then
            qx(i)=0.d0
         else
            qx(i)=vx(i)+fattore_conversione*0.5d0*(fx(i)+dfx(i))*tstep/mass(itype(i)) 
         endif
         if (fy(i)*vy(i).lt.0) then
            qy(i)=0.d0
         else
            qy(i)=vy(i)+fattore_conversione*0.5d0*(fy(i)+dfy(i))*tstep/mass(itype(i)) 
         endif
         if (fz(i)*vz(i).lt.0) then
            qz(i)=0.d0
         else
            qz(i)=vz(i)+fattore_conversione*0.5d0*(fz(i)+dfz(i))*tstep/mass(itype(i)) 
         endif
      endif    
      if (nxfix.eq.1) then
         qx(natom)=0.d0       
      endif
      vx(i)=qx(i) 
      vy(i)=qy(i)
      vz(i)=qz(i)      
      !in A*A/s*s   
      vel2=vx(i)*vx(i)+vy(i)*vy(i)+vz(i)*vz(i)
      !l-energia cinetica la esprimo in eV
      ecinet=ecinet+(0.5d0*mass(itype(i))*uasukg*vel2*angsum*angsum)/evsujoule
  enddo 

else IF ((choice_sim.eq.'ev').or.(choice_sim.eq.'hf').or.((choice_sim.eq.'gr').and.(ipas.le.growth_steps)))   THEN
!write(*,*)'at step ',ipas,'velocities are assigned as usual'
        last_atom_2bthermalized=natom
	!quello che segue lo faccio solo in termalizzazione:
        do i=1,last_atom_2bthermalized !termalizzo gli atomi con vicini      
           if (thermostat .eq. 'si') then
              call ranmar(zzz)
              vdfae=zzz
           if (vdfae .ge. vdfa) then   
                 fattore_conversione=evsujoule/(uasukg*angsum*angsum)
              !le df sono quelle vecchie, le f sono state appena calcolate da force
              !le forze sono espresse in eV/A, vanno convertite in J/m
              !a questo punto il secondo addendo va convertito in da m/s a A/s
                 qx(i)=vx(i)+fattore_conversione*0.5d0*(fx(i)+dfx(i))*tstep/mass(itype(i)) 
                 qy(i)=vy(i)+fattore_conversione*0.5d0*(fy(i)+dfy(i))*tstep/mass(itype(i))     
                 qz(i)=vz(i)+fattore_conversione*0.5d0*(fz(i)+dfz(i))*tstep/mass(itype(i))     
                 vx(i)=qx(i) 
                 vy(i)=qy(i)
                 vz(i)=qz(i)
              else ! altrimenti interviene il termostato a riscalare le velocita':
                 fattore_conversione=sqrt(evsujoule/uasukg)/angsum
                 call gauss   
                 vx(i)=g1*fattore_conversione*sqrt((cbol*temperatura)/mass(itype(i)))
                 vy(i)=g2*fattore_conversione*sqrt((cbol*temperatura)/mass(itype(i)))
                 vz(i)=g3*fattore_conversione*sqrt((cbol*temperatura)/mass(itype(i)))
              endif  !sulla chiamata a GAUSS
              !se thermostat=no e ho superato il numero di passi di termalizzazione, non termostato:
           elseif (thermostat .eq. 'no') then
              fattore_conversione=evsujoule/(uasukg*angsum*angsum)
              !le df sono quelle vecchie, le f sono state appena calcolate da force
              !le forze sono espresse in eV/A, vanno convertite in J/m
              !a questo punto il secondo addendo va convertito in da m/s a A/s
              qx(i)=vx(i)+fattore_conversione*0.5d0*(fx(i)+dfx(i))*tstep/mass(itype(i)) 
              qy(i)=vy(i)+fattore_conversione*0.5d0*(fy(i)+dfy(i))*tstep/mass(itype(i))     
              qz(i)=vz(i)+fattore_conversione*0.5d0*(fz(i)+dfz(i))*tstep/mass(itype(i))     
              vx(i)=qx(i) 
              vy(i)=qy(i)
              vz(i)=qz(i)      
           endif
	   !in A*A/s*s
           vel2=vx(i)*vx(i)+vy(i)*vy(i)+vz(i)*vz(i)
	   !l'energia cinetica in eV
           ecinet=ecinet+(0.5d0*mass(itype(i))*uasukg*vel2*angsum*angsum)/evsujoule
        enddo   !su i atomi con vicini

ELSE IF ((choice_sim.eq.'gr') .and. (ipas.ge.(growth_steps+1))) THEN
       last_atom_2bthermalized=natom-1
	!quello che segue lo faccio solo in termalizzazione:
       do i=1,last_atom_2bthermalized !termalizzo gli atomi con vicini      
           !estraggo con una certa frequenza il termostato nei due casi in cui:
           !1. thermostat=si
           !2. thermostat=no ma sono ancora in termalizzazione
           if (thermostat .eq. 'si')then
              call ranmar(zzz)
              vdfae=zzz
              if (vdfae .ge. vdfa) then   
                 fattore_conversione=evsujoule/(uasukg*angsum*angsum)
                 !le df sono quelle vecchie, le f sono state appena calcolate da force
                 !le forze sono espresse in eV/A, vanno convertite in J/m
                 !a questo punto il secondo addendo va convertito in da m/s a A/s
                 qx(i)=vx(i)+fattore_conversione*0.5d0*(fx(i)+dfx(i))*tstep/mass(itype(i)) 
                 qy(i)=vy(i)+fattore_conversione*0.5d0*(fy(i)+dfy(i))*tstep/mass(itype(i))     
                 qz(i)=vz(i)+fattore_conversione*0.5d0*(fz(i)+dfz(i))*tstep/mass(itype(i))     
                 vx(i)=qx(i) 
                 vy(i)=qy(i)
                 vz(i)=qz(i)
              else ! altrimenti interviene il termostato a riscalare le velocita':
                 fattore_conversione=sqrt(evsujoule/uasukg)/angsum
                 call gauss   
                 vx(i)=g1*fattore_conversione*sqrt((cbol*temperatura)/mass(itype(i)))
                 vy(i)=g2*fattore_conversione*sqrt((cbol*temperatura)/mass(itype(i)))
                 vz(i)=g3*fattore_conversione*sqrt((cbol*temperatura)/mass(itype(i)))
              endif  !sulla chiamata a GAUSS
              !se thermostat=no e ho superato il numero di passi di termalizzazione, non termostato:
           elseif(thermostat .eq. 'no')then
              fattore_conversione=evsujoule/(uasukg*angsum*angsum)
              !le df sono quelle vecchie, le f sono state appena calcolate da force
              !le forze sono espresse in eV/A, vanno convertite in J/m
              !a questo punto il secondo addendo va convertito in da m/s a A/s
              qx(i)=vx(i)+fattore_conversione*0.5d0*(fx(i)+dfx(i))*tstep/mass(itype(i)) 
              qy(i)=vy(i)+fattore_conversione*0.5d0*(fy(i)+dfy(i))*tstep/mass(itype(i))     
              qz(i)=vz(i)+fattore_conversione*0.5d0*(fz(i)+dfz(i))*tstep/mass(itype(i))     
              vx(i)=qx(i) 
              vy(i)=qy(i)
              vz(i)=qz(i)      
           endif
	   !in A*A/s*s   
           vel2=vx(i)*vx(i)+vy(i)*vy(i)+vz(i)*vz(i)
	   !l'energia cinetica in eV
           ecinet=ecinet+(0.5d0*mass(itype(i))*uasukg*vel2*angsum*angsum)/evsujoule   
       enddo   !su i atomi con vicini
       !the atom that has just come in from the source is given a velocity with no thermalization:
       fattore_conversione=evsujoule/(uasukg*angsum*angsum)
       !le df sono quelle vecchie, le f sono state appena calcolate da force
       !le forze sono espresse in eV/A, vanno convertite in J/m
       !a questo punto il secondo addendo va convertito in da m/s a A/s
       qx(natom)=vx(natom)+fattore_conversione*0.5d0*(fx(natom)+dfx(natom))*tstep/mass(itype(natom)) 
       qy(natom)=vy(natom)+fattore_conversione*0.5d0*(fy(natom)+dfy(natom))*tstep/mass(itype(natom))     
       qz(natom)=vz(natom)+fattore_conversione*0.5d0*(fz(natom)+dfz(natom))*tstep/mass(itype(natom))     
       vx(natom)=qx(natom) 
       vy(natom)=qy(natom)
       vz(natom)=qz(natom)
       !in A*A/s*s   
       vel2=vx(natom)*vx(natom)+vy(natom)*vy(natom)+vz(natom)*vz(natom)
       !write(*,*)'In vel.f90, the velocity of the adatom is: ',sqrt(vel2)
       !l-energia cinetica la esprimo in eV
       ecinet=ecinet+(0.5d0*mass(itype(i))*uasukg*vel2*angsum*angsum)/evsujoule
ENDIF
  
!temp e' la temperatura calcolata sulla base delle velocita':
temp=2.d0*ecinet/(3.d0*(natom-nfix)*cbol) 

!tutte le energie a questo punto sono espresse in eV:
etot=ener+ecinet 

!write (55,*) temp
!do i=1,natom
!   write (55,*) i, vx(i),vy(i),vz(i)
!enddo
!stop

END SUBROUTINE VEL
        
        
        
