program MolDin
  
  ! Molecular Dynamics code
  ! last changed: 25 January 2010
  ! author: Giulia Rossi
  
  USE CLUSTER
  USE PARAMETERS
  
  implicit none
  
  integer :: i
  real(8) :: rvec,dmax,secd,firstd
  real*8 :: partialx(nmax),partialy(nmax),partialz(nmax)
  Real*8 :: delta,sum_tempmedia,sum_enermedia
  Character*7 :: gruppo
  Integer :: ordine_gruppo
  Logical :: minim

  partialx(:)=0.
  partialy(:)=0.  
  partialz(:)=0.
  
  write (*,*) "chiamo leggi"
  call leggi_input

  !to check that the reading of cluster.in was ok:
  open(12,file='out.xyz',status='unknown')
  write(12,*)natom
  write(12,'(a2,1x,a2)')elem1,elem2
  do i=1,natom
     write(12,'(a2,1x,3f16.6)')elem(i),x(i),y(i),z(i)
  enddo
  close(12)
   
  !initialization of velocities @ T=tinit
  call init
  
  !neighbor list update
  call bigvoi
  
  !neighbor list
  call voisin
  
  !energy and forces
  call force_rgl

!  if(oxide_substrate .eq. 'si')then
!    minim=.false.
!    call force_met_mgo(minim)
!  endif

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
        partialx(i)=partialx(i)+du(i)
        partialy(i)=partialy(i)+dv(i)
        partialz(i)=partialz(i)+dw(i)
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
     
     call force_rgl
!     if(oxide_substrate .eq. 'si')then
!       minim=.false.
!       call force_met_mgo(minim)
!     endif
      
     !write(*,*)'before calling vel'
     call vel
     !updating atomic positions
     call therma

     sum_tempmedia=sum_tempmedia+temp
     sum_enermedia=sum_enermedia+etot

     if ((mod(ipas,scrivo).eq.0).or.(ipas.eq.1)) then
        call write_photo
     endif


     ! WRITING OF AVERAGES IN QUENCHING SIMULATIONS
     
     if (choice_sim.eq.'qq') then
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


     
 
  enddo

  
end program MolDin

