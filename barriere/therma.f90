SUBROUTINE THERMA

  USE PARAMETERS
  USE CLUSTER     !uso il modulo dove definiscovariabili e parametri cluster
  
  Implicit None
  
  Real*8 :: fattore_conversione
  Real    :: xcentr_sum, ycentr_sum, zcentr_sum, ddist
  Integer :: is,i
  
  fattore_conversione=evsujoule/(angsum*angsum*uasukg)
  
! solo gli atomi da nfix+1 a natom si muovono  
  do is=nfix+1,natom   !evoluzione da Velocity-Verlet 
     if (tempo.lt.1000*timelag) then 
        du(is)=fattore_conversione*0.5d0*tstep*tstep*fx(is)/mass(itype(is))+vx(is)*tstep
        dv(is)=fattore_conversione*0.5d0*tstep*tstep*fy(is)/mass(itype(is))+vy(is)*tstep
        dw(is)=fattore_conversione*0.5d0*tstep*tstep*fz(is)/mass(itype(is))+vz(is)*tstep       
     else
     !(la forza e' espressa in eV/A e va convertita in J/m.
     ! fatto questo, il primo addendo va convertito in A.
     ! le velocita' sono espresse in A/m e vanno bene.) 
        if (fx(is)*vx(is).lt.0) then
           du(is)=0.d0
        else            
           du(is)=fattore_conversione*0.5d0*tstep*tstep*fx(is)/mass(itype(is))+vx(is)*tstep
        endif
        if (fy(is)*vy(is).lt.0) then
           dv(is)=0.d0
        else            
           dv(is)=fattore_conversione*0.5d0*tstep*tstep*fy(is)/mass(itype(is))+vy(is)*tstep
        endif
        if (fz(is)*vz(is).lt.0) then
           dw(is)=0.d0
        else            
           dw(is)=fattore_conversione*0.5d0*tstep*tstep*fz(is)/mass(itype(is))+vz(is)*tstep
        endif
     endif  
     if (nxfix.eq.1) then
        du(natom)=0.d0       
     endif
     u(is)=u(is)+du(is)
     v(is)=v(is)+dv(is)
     w(is)=w(is)+dw(is)     
     !write(*,*)'du',is,vx(is)*tstep
  enddo
  
  !se boundary_conditions=si, periodicizzo la scatola e salvo le coordinate assolute: 
!  if(boundary_cond .eq. 'si')then
     
     !box_edge e' sempre positivo (calcolato sulla base di una distanza)     
!     do is=1,natom
!        if ( (x(is)+u(is)) .gt. x_boxedge/2. )        u(is)=u(is)-x_boxedge
!        if ( (x(is)+u(is)) .lt. (x_boxedge*(-1)/2.) ) u(is)=u(is)+x_boxedge   
!        if ( (y(is)+v(is)) .gt. y_boxedge/2.)        v(is)=v(is)-y_boxedge
!        if ( (y(is)+v(is)) .lt. (y_boxedge*(-1)/2.) ) v(is)=v(is)+y_boxedge   
!        IF(zboundary_cond .eq. 'si')THEN
!          if ( (z(is)+w(is)) .gt. z_boxedge/2.)        w(is)=w(is)-z_boxedge
!          if ( (z(is)+w(is)) .lt. (z_boxedge*(-1)/2.) ) w(is)=w(is)+z_boxedge   
!        ENDIF
!     enddo
     
!  endif
  

END SUBROUTINE THERMA
