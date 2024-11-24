program MolDin
    
  USE CLUSTER
  USE PARAMETERS
  
  implicit none
  
  integer :: i
  real*8 :: dmax,secd,firstd
  real*8 :: partialx(nmax),partialy(nmax),partialz(nmax)
  Real*8 :: delta,sum_tempmedia,sum_enermedia

  partialx(:)=0.
  partialy(:)=0.  
  partialz(:)=0.
  
  write (*,*) "chiamo leggi"
  call leggi_input
!  stop
  !to check that the reading of the initial coordinates was ok:
  open(12,file='out.xyz',status='unknown')
  write(12,*) natom
  write(12,*) (elem(i),i=1,nelem)
  do i=1,natom
     write(12,'(a2,1x,3f16.6)')spec(i),x(i),y(i),z(i)
  enddo
  close(12)
   
  !initialization of velocities @ T=tinit
  call init
  
  !neighbor list update
  call bigvoi
  
  !neighbor list
  call voisin
  
  !energy and forces
  if (potential.eq.'l_j') then
     call force_lj
  else if (potential.eq.'mrs') then
     call force_morse
  else if (potential.eq.'gpt') then
     call force_gupta
!     write (83,*) ipas,ener
  else
     write (*,*) 'wrong potential'
     stop
  endif      


  !position update
  call therma

  write(*,*) 'INITIALIZATION COMPLETED'
  write (*,*) "ener=", ener

  !headings of medie.out  
  open (4100, file='medie.out', status='unknown')
       write(4100,*)'#uni, ipas, time, temp, average temper, average ener'     
  close(4100)  


  !INIZIO LOOP SUI TEMPI
  
  sum_tempmedia=0.d0
  sum_enermedia=0.d0
  write(*,*)'Total simulation steps:',npas
  do ipas=1,npas
     tempo=ipas*tstep*1.d12        !time units are ps          
     !temperature changes in freezing/melting simulations:     
     !check if neighbor list needs to be updated by bigvoi:
     firstd=0.d0
     secd=0.d0
     do i=1,natom
        partialx(i)=partialx(i)+dx(i)
        partialy(i)=partialy(i)+dy(i)
        partialz(i)=partialz(i)+dz(i)
        delta=partialx(i)**2+partialy(i)**2+partialz(i)**2
        if ((delta.gt.secd).and.(delta.le.firstd)) then
           secd=delta
        endif
        if (delta.gt.firstd) then
           secd=firstd
           firstd=delta
        endif
     enddo
     
     dmax=sqrt(firstd)+sqrt(secd)
     if (dmax.gt.rshell) then !update is needed
        do i=1,natom
           partialx(i)=0.d0
           partialy(i)=0.d0
           partialz(i)=0.d0
        enddo

        call bigvoi
     endif
     
     !proseguono i passi necessari all'uso del velocity verlet
     call voisin

     if (potential.eq.'l_j') then
        call force_lj
     else if (potential.eq.'mrs') then
        call force_morse
     else if (potential.eq.'gpt') then
        call force_gupta
 !       write (93,*) ipas,ener
     endif
     
!    PROPAGATION OF TRAJECTORIES     
     call vel
     !updating atomic positions
     call therma

     sum_tempmedia=sum_tempmedia+temp
     sum_enermedia=sum_enermedia+etot

     if ((mod(ipas,scrivo).eq.0).or.(ipas.eq.1)) then
        call write_photo
     endif


!    WRITING OF AVERAGES IN QUENCHING OR EVOLUTION SIMULATIONS
     
     if ((choice_sim.eq.'qq').or.(choice_sim.eq.'ev')) then
!        write (*,*) 'sono qui',evo_steps
        if (mod(ipas,evo_steps).eq.0) then
           tempmedia=sum_tempmedia/dfloat(evo_steps)
           enermedia=sum_enermedia/dfloat(evo_steps)
           call write_avg
           if (ipas.lt.npas) then
              sum_tempmedia=0.d0
              sum_enermedia=0.d0
           endif
        endif
     endif 



     !CHANGE OF TEMPERATURE IN HEATING/FREEZING AND WRITING OF AVERAGES
     
     if (choice_sim.eq.'hf') then
        if (mod(ipas,aumento_temp).eq.0) then
           tempmedia=sum_tempmedia/dfloat(aumento_temp)
           enermedia=sum_enermedia/dfloat(aumento_temp)
           call write_avg
           if (ipas.lt.npas) then
              temperatura=temperatura+delta_temp 
              sum_tempmedia=0.d0
              sum_enermedia=0.d0
           endif
        endif
     endif 
               
     !DEPOSITION OF A NEW ATOM GROWTH SIMULATIONS AND WRITING OF AVERAGES
     
     if (choice_sim.eq.'gr') then
        if (mod(ipas,growth_steps).eq.0) then
           tempmedia=sum_tempmedia/dfloat(growth_steps)
           enermedia=sum_enermedia/dfloat(growth_steps)
           call write_avg
           if (ipas.lt.npas) then
              call growth !the number of atoms is increased there
              sum_tempmedia=0.d0
              sum_enermedia=0.d0
           endif
        endif
     endif  

     
 
  enddo

  
end program MolDin

