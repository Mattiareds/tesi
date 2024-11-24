SUBROUTINE THERMA

  USE PARAMETERS
  USE CLUSTER     !uso il modulo dove definiscovariabili e parametri cluster
  
  Implicit None
  
  Real*8 :: fattore_conversione
  Integer :: is
  
  fattore_conversione=evsujoule/(angsum*angsum*uasukg)
  
! solo gli atomi da nfix+1 a natom si muovono  
if (choice_sim.eq.'qq') then

  do is=nfix+1,natom   !evoluzione da Velocity-Verlet 
     if (tempo.lt.1000*timelag) then 
        dx(is)=fattore_conversione*0.5d0*tstep*tstep*fx(is)/mass(itype(is))+vx(is)*tstep
        dy(is)=fattore_conversione*0.5d0*tstep*tstep*fy(is)/mass(itype(is))+vy(is)*tstep
        dz(is)=fattore_conversione*0.5d0*tstep*tstep*fz(is)/mass(itype(is))+vz(is)*tstep       
     else
     !(la forza e' espressa in eV/A e va convertita in J/m.
     ! fatto questo, il primo addendo va convertito in A.
     ! le velocita' sono espresse in A/m e vanno bene.) 
        if (fx(is)*vx(is).lt.0) then
           dx(is)=0.d0
        else            
           dx(is)=fattore_conversione*0.5d0*tstep*tstep*fx(is)/mass(itype(is))+vx(is)*tstep
        endif
        if (fy(is)*vy(is).lt.0) then
           dy(is)=0.d0
        else            
           dy(is)=fattore_conversione*0.5d0*tstep*tstep*fy(is)/mass(itype(is))+vy(is)*tstep
        endif
        if (fz(is)*vz(is).lt.0) then
           dz(is)=0.d0
        else            
           dz(is)=fattore_conversione*0.5d0*tstep*tstep*fz(is)/mass(itype(is))+vz(is)*tstep
        endif
     endif  
     if (nxfix.eq.1) then
        dx(natom)=0.d0       
     endif
     x(is)=x(is)+dx(is)
     y(is)=y(is)+dy(is)
     z(is)=z(is)+dz(is)     
   enddo
   
else

  do is=1,natom   !evoluzione da Velocity-Verlet 
     !(la forza e' espressa in eV/A e va convertita in J/m.
     ! fatto questo, il primo addendo va convertito in A.
     ! le velocita' sono espresse in A/m e vanno bene.)            
     dx(is)=fattore_conversione*00.5d0*tstep*tstep*fx(is)/mass(itype(is))+vx(is)*tstep
     dy(is)=fattore_conversione*00.5d0*tstep*tstep*fy(is)/mass(itype(is))+vy(is)*tstep
     dz(is)=fattore_conversione*00.5d0*tstep*tstep*fz(is)/mass(itype(is))+vz(is)*tstep     
     x(is)=x(is)+dx(is)
     y(is)=y(is)+dy(is)
     z(is)=z(is)+dz(is)         
  enddo
  
endif  

END SUBROUTINE THERMA
