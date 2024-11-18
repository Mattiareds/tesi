subroutine leggi_input

  USE PARAMETERS  
  USE CLUSTER 
  Implicit None
  
  ! variabili locali
  Integer :: ind,i,j,inderr
  Real*8 :: x0sum,y0sum,z0sum
  Real*8 :: rate_temp,summax
  Character(20):: param_file
  Character(20):: cluster_file

!*************************************************************
!READING PARAMETERS FOR MOLECULAR DYNAMICS RUN
!*************************************************************

  cluster_file='cluster.in'
  write (*,*) "START READING FILE settings.in" 
  write (*,*)
  open(14,file='settings.in',status='old')
! general settings
  read(14,*)
  read(14,*) potential
  read(14,*) param_file
!*************************************************************
!READING POTENTIAL PARAMETERS
!*************************************************************
write (*,*) "START READING FILE ", param_file
! it is read now to know nelem, which is necessary for growth probabilities
  if (potential.eq.'l_j') then
     write (*,*)
     write (*,*) "Lennard-Jones potential"
     write (*,*)
     call leggi_pot_lj(param_file)
  else if (potential.eq.'mrs') then
     write (*,*)
     write (*,*) "Morse potential"   
    write (*,*)
    call leggi_pot_morse(param_file)
  else if (potential.eq.'gpt') then
     write (*,*)
     write (*,*) "Morse potential"   
     write (*,*)
     call leggi_pot_gupta(param_file)
  else
    write (*,*)
    write (*,*) 'potential is unknown'
    stop
  endif 
  write (*,*)
!stop
!*************************************************************
!  close(1)
!  open(1,file='settings.in',status='old')
! simulation settings
!  read(14,*)
!  read(14,*) 
!  read(14,*) 
  read(14,*) cluster_file
  read(14,*) tstep
  read(14,'(a2)') thermostat
  read(14,*) vnu
  read(14,*) rshell
  read(14,'(a2)') cutoff
  read(14,*) nfix
  read(14,'(a2)') choice_sim
  read(14,*) irand1,irand2
! for quenching simulations 
  if (choice_sim.eq.'qq') then
     read(14,*)
     read(14,*) que_temp
     read(14,*) timefinqq
     read(14,*) timelag
     read(14,*) nxfix
  else
     read(14,*)
     read(14,*)
     read(14,*) 
     read(14,*) 
     read(14,*) 
  endif
! for growth simulations    
  if (choice_sim.eq.'gr') then
     read(14,*)
     read(14,*) growth_rate
     read(14,*) growth_temp
     read(14,*) ndepo
     read(14,*) (growth_p1(i),i=1,nelem)
  else
     read(14,*)
     read(14,*)
     read(14,*) 
     read(14,*) 
     read(14,*) 
  endif
! for evolution simulations 
  if (choice_sim.eq.'ev') then
     read(14,*)
     read(14,*) evo_temp
     read(14,*) timefin
     read(14,*) nsections
  else
     read(14,*)
     read(14,*)
     read(14,*) 
     read(14,*) 
  endif
! for heating/freezing simulations 
  if (choice_sim.eq.'hf') then
     read(14,*)
     read(14,*) tinit
     read(14,*) tfin
     read(14,*) rate_temp
  else
     read(14,*)
     read(14,*)
     read(14,*) 
     read(14,*) 
  endif
! for output 
  read(14,*)  
  read(14,*) nsections
  read(14,*) photo_rate
  
  close(14)

  write(*,*)  
    
  write(*,*) 'potential = ', potential
  write(*,*) 'param_file = ', param_file
  write(*,*) 'cluster_file = ', cluster_file
  write(*,*) 'tstep = ', tstep
  write(*,*) 'thermostat = ', thermostat
  write(*,*) 'vnu = ', vnu
  write(*,*) 'rshell = ', rshell
  write(*,*) 'cutoff = ', cutoff
  write(*,*) 'nfix = ', nfix
  write(*,*) 'choice_sim = ', choice_sim
  write(*,*) 'irand1,irand2 = ', irand1,irand2
  if (choice_sim.eq.'qq') then 
!    for quenching simulations 
     write(*,*) 'que_temp = ', que_temp
     write(*,*) 'timefinqq = ', timefinqq
     write(*,*) 'timelag = ', timelag
     write(*,*) 'nxfix = ', nxfix
  else if (choice_sim.eq.'gr') then
     write(*,*) 'growth_rate =', growth_rate
     write(*,*) 'growth_p1 =', (growth_p1(i),i=1,nelem)
     write(*,*) 'growth_temp =', growth_temp
     write(*,*) 'ndepo =', ndepo 
  else if (choice_sim.eq.'ev') then
     write(*,*) 'evo_temp = ', evo_temp
     write(*,*) 'timefin = ', timefin
     write(*,*) 'nsections = ', nsections
  else if (choice_sim.eq.'hf') then
     write(*,*) 'tinit = ', tinit
     write(*,*) 'tfin = ', tfin
     write(*,*) 'rate_temp =', rate_temp
  endif
  
! for output 
  write(*,*) 'nsections = ', nsections
  write(*,*) 'photo_rate = ', photo_rate

  write(*,*)  
  
  write (*,*) "END READING FILE settings.in"
  write (*,*)
  if (nxfix.gt.1) then
     write (*,*) 'nxfix sbagliato'
     stop
  endif
  !numeri casuali
  call rmarin(irand1,irand2)
  !poi quando serve un numero casuale: call ranmar(zzz)

  
  !calcolo i numeri di passi totali e parziali, fisso delta_temp
  !e stabilisco la temperatura iniziali
  delta_temp=0.d0
  if (choice_sim.eq.'qq') then
     temperatura=que_temp
     npas=int(timefinqq*1.d-9/tstep,8)+1
     evo_steps=npas/nsections
     write (*,*) "evo_steps = ", evo_steps
     write (*,*) "npas = ", npas
 else if (choice_sim.eq.'hf') then
     temperatura=tinit
     aumento_temp=int(abs(1.d0/(tstep*rate_temp*1.d9)),8)
     npas=int(abs(tinit-tfin)*aumento_temp,8)
     write (*,*) "aumento_temp = ", aumento_temp
     if (tfin.gt.tinit) then
        delta_temp=1.d0
     else if (tfin.lt.tinit) then
        delta_temp=-1.d0
     else
        write (*,*) 'error: temperature is constant in hf simulation'
        stop
     endif
     if ((rate_temp*(tinit-tfin)).gt.0.d0) then
        write (*,*) 'error: no concordance in rate_temp, tfin-tinit"'
        stop
     endif    
  else if (choice_sim.eq.'gr') then
     sum_p1(0)=0.d0
     do i=1,nelem
        sum_p1(i)=sum_p1(i-1)+growth_p1(i)
     enddo
     summax=sum_p1(nelem)
     do i=1,nelem
        sum_p1(i)=sum_p1(i)/summax
        growth_p1(i)=growth_p1(i)/summax
        write (*,*) i, elem(i),growth_p1(i),sum_p1(i)       
      enddo
!     stop
     temperatura=growth_temp
     growth_steps=int(1.d0/(tstep*growth_rate*1.d9),8)+1
     npas=growth_steps*(ndepo+1)
     write (*,*) "growth_steps = ", growth_steps     
  else if (choice_sim.eq.'ev') then
     temperatura=evo_temp
     npas=int(timefin*1.d-9/tstep,8)+1
     evo_steps=npas/nsections
     write (*,*) "evo_steps = ", evo_steps
     write (*,*) "npas = ", npas
  else
     write (*,*) 'unknown type of simulation'
     stop
  endif
  write (*,*) "npas = ", npas
  
 !stabilisce ogni quanti passi deve prendere una foto
  scrivo=int(1.d-9/(photo_rate*tstep),8)
  write (*,*) "scrivo = ", scrivo
  
!*************************************************************
!READING INITIAL STRUCTURE
!*************************************************************
  write(*,*)
  write(*,*) "START READING FILE ", cluster_file 
  write(*,*)
     open(25,file=cluster_file,status='unknown')
     read (25,*) natom
     if (natom.gt.nmax) then
        write(*,*) 'The input cluster is too big: natom > nmax.'
        write(*,*)'Program stops here.'
        stop
     endif
     if (natom.le.nfix) then
        write(*,*) 'no atom can move: nfix>=natom'
        write(*,*)'Program stops here.'
        stop
     endif     
     read(25,*)      
     ! reading input cluster: atomic positions and atomic types
     do i=1,nelem
        ntype(i)=0
     enddo
     do ind=1,natom
        read(25,*) spec(ind),x0(ind),y0(ind),z0(ind)
        inderr=0
        do i=1,nelem
!           write (*,*) ind,i,'',spec(ind),' ',elem(i)
           if (spec(ind).eq.elem(i)) then 
              itype(ind)=i
              ntype(i)=ntype(i)+1
              inderr=1
           endif        
        enddo
          if(inderr.eq.0) then
             write (*,*) "wrong element",ind,spec(ind)
             stop
          endif
     enddo
     
!     if (boundary_cond .eq. 'si') then
!        read(25,*) x_boxedge
!        read(25,*) y_boxedge
!        read(25,*) z_boxedge
!     endif
     close(25)
!     ntype1=natom-ntype2
     do i=1,nelem
        write(*,*) cluster_file,' contains',ntype(i),'atoms of species no.',i, elem(i)
     enddo
     
  write(*,*)
  write(*,*) "END READING FILE ", cluster_file 
  write(*,*)


!     stop
!System is centered in (0.,0.,0.)
     x0sum=0.d0
     y0sum=0.d0
     z0sum=0.d0
     do ind=1,natom
        x0sum=x0sum+x0(ind)
        y0sum=y0sum+y0(ind)
        z0sum=z0sum+z0(ind)
     enddo
     do ind=1,natom
        x(ind)=x0(ind)-x0sum/dfloat(natom)
        y(ind)=y0(ind)-y0sum/dfloat(natom)
        z(ind)=z0(ind)-z0sum/dfloat(natom)
     enddo     
  


  !cutoff_end and cutoff_start are set to very large values if cutoff is not required
  if(cutoff .eq. 'no')then
     do i=1,nelem
        do j=1,nelem
           cutoff_start(i,j)=2000.
           cutoff_end(i,j)=2000.
        enddo
     enddo
  endif
! calcolo del massimo cutoff_end
  cutoff_end_max=0.d0
  do i=1,nelem
     do j=1,nelem
     if(cutoff_end(i,j).gt.cutoff_end_max) cutoff_end_max=cutoff_end(i,j)
!     write(*,*) "cutoff_end", i, cutoff_end(i) 
     enddo
!     write(*,*) "cutoff_end_max = ", cutoff_end_max 
  enddo
  if ((potential.ne.'gpt').and.(cutoff.eq.'si')) then
     write (*,*) 'Please put no in cutoff and restart. Cutoff is allowed only for Gupta potential'
     stop 
  endif


!stop
End Subroutine leggi_input
