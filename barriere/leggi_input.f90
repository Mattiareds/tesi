subroutine leggi_input

  USE PARAMETERS  
  USE CLUSTER 
  USE STRUCTURES 
  Implicit None
  
  ! variabili locali
  Integer :: ind,i,j,k
  Real*8 :: x0sum,y0sum,z0sum
  Real*8 :: distanza, distanza_max
  Real*8 :: rate_temp,volume_scatola
  Real*8 :: ar,br,cr,ab,bb,cb,dik0 
  Character(12) :: sis
  Character(20):: param_file
  Character(2) :: metal(2)

!*************************************************************
!READING PARAMETERS FOR MOLECULAR DYNAMICS RUN
!*************************************************************

  write (*,*) "entro in leggi" 
  open(1,file='leggi.in',status='old')
! general settings
  read(1,*)
  read(1,*)irand1
  read(1,*)irand2
  read(1,*)potential
  read(1,*)param_file
  read(1,*)tstep
  read(1,'(a2)')thermostat
  read(1,*)rshell
  read(1,'(a2)')cutoff
  read(1,'(a2)')choice_sim
  read(1,*)photo_rate
  read(1,*)nfix
  read(1,*)nxfix
! for evolution simulations 
  read(1,*)
  read(1,*)evo_temp
  read(1,*)timefin
  read(1,*)timelag
  read(1,*)nsections
  close(1)
  write (*,*) "finito leggi.in"
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
     temperatura=evo_temp
     npas=idint(timefin*1.d-9/tstep)+1
     evo_steps=npas/nsections
     write (*,*) "evo_steps = ", evo_steps
  else
     write (*,*) 'unknown type of simulation'
     stop
  endif
  write (*,*) "npas = ", npas
  
 !stabilisce ogni quanti passi deve prendere una foto
  scrivo=idint(1.d-9/(photo_rate*tstep))
  write (*,*) "scrivo = ", scrivo
  
  
!*************************************************************
!READING INPUT STRUCTURE
!*************************************************************
     
     open(25,file='cluster.in',status='unknown')
     read (25,*) natom
     if (natom.gt.nmax) then
        write (*,*) 'The input cluster is too big: natom > nmax.'
        write(*,*)'Program stops here.'
        stop
     endif
     if (natom.le.nfix) then
        write (*,*) 'no atom can move: nfix>=natom'
        write(*,*)'Program stops here.'
        stop
     endif     
     read (25,*) elem1,elem2  
     sis='monometallic'
     if(elem1.ne.elem2) sis='_bimetallic_'
     write(*,'(a8,1x,a13,1x,a2,1x,a2)') 'System =',sis,elem1,elem2
     
     ! reading input cluster: atomic positions and atomic types
     itype(1:natom)=1
     ntype1=0
     ntype2=0
     do ind=1,natom
        read(25,*) elem(ind),x0(ind),y0(ind),z0(ind)
        if(elem(ind).eq.elem2) then 
           itype(ind)=2
           ntype2=ntype2+1
        endif
        if((elem(ind).ne.elem1).and.(elem(ind).ne.elem2)) then
           write(*,*)'WARNING: the species declared in the heading of cluster.in are not correct.'
           write(*,*)'Program stops while reading atom ',ind,' of cluster.in'
           stop
        endif
     enddo
!     if (boundary_cond .eq. 'si') then
!        read(25,*) x_boxedge
!        read(25,*) y_boxedge
!        read(25,*) z_boxedge
!     endif
     close(25)
     ntype1=natom-ntype2
     write(*,*)'cluster.in contains',ntype1,'atoms of species no.1'
     write(*,*)'cluster.in contains',ntype2,'atoms of species no.2'
     
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
  
!*************************************************************
!READING POTENTIAL PARAMETERS
!*************************************************************
  
  open(19,file=param_file,status='old')
  read(19,*)
  read(19,*)
  read(19,*)metal(1), metal(2)
            !check consistency with metals declared in seed.in
            if ((metal(1) .ne. elem1) .or. (metal(2) .ne. elem2))then
               write(*,*)'####################################################'
               write(*,*)'ERROR:'
               write(*,*)'Metals declared in cluster.in'
               write(*,*)'MUST BE THE SAME AND IN THE SAME ORDER' 
               write(*,*)'as in the parametrization file.'
               write(*,*)'Program stops here.'
               write(*,*)'####################################################'
               stop
            endif
            !end check consistency
  read(19,*)
  read(19,*)
  read(19,*)p(1),p(2),p(3)
  read(19,*)q(1),q(2),q(3)
  read(19,*)a(1),a(2),a(3)
  read(19,*)qsi(1),qsi(2),qsi(3)
  read(19,*)
  read(19,*)
  read(19,*)ecoh(1),ecoh(2)
  read(19,*)rat(1),rat(2)
  read(19,*)mass(1),mass(2)
  read(19,*)
  read(19,*)
  read(19,*)cutoff_start(1),cutoff_end(1)
  read(19,*)cutoff_start(2),cutoff_end(2)
  read(19,*)cutoff_start(3),cutoff_end(3)
!  read(19,*)cutoff_start,cutoff_end
  write (*,*) 'END READING FILE ', param_file
  close(19)

  !Unit conversions

  arete(1)=rat(1)*dsqrt(8.d0)
  arete(2)=rat(2)*dsqrt(8.d0)
  arete(3)=(arete(1)+arete(2))/2.0d0

  !width for the box to apply boundary conditions, in arete units
  x_boxedge_a=x_boxedge/arete(1)
  y_boxedge_a=y_boxedge/arete(1)
  z_boxedge_a=z_boxedge/arete(1)

  !nn are the nearest neighbours distances in A
  nn(1)=arete(1)/dsqrt(2.d0)
  nn(2)=arete(2)/dsqrt(2.d0)
  nn(3)=arete(3)/dsqrt(2.d0)

  !dist are the nearest neighbours distances in arete(1) units
  dist(1)=1.d0/dsqrt(2.d0)
  dist(2)=nn(2)/arete(1)
  dist(3)=nn(3)/arete(1)

  !cutoff_end and cutoff_start are set to very large values if cutoff is not required
  if(cutoff .eq. 'no')then
     do i=1,3
        cutoff_start(i)=2000.
        cutoff_end(i)=2000.
     enddo
  endif
! calcolo del massimo cutoff_end
  cutoff_end_max=0.d0
  do i=1,3
     if(cutoff_end(i).gt.cutoff_end_max) cutoff_end_max=cutoff_end(i)
     write(*,*) "cutoff_end", i, cutoff_end(i) 
  enddo
     write(*,*) "cutoff_end_max = ", cutoff_end_max 

  !Cutoff parameters a5,a4,a3,x5,x4,x3
  do i=1,3 
        dik0=nn(i)
        ar=-a(i)*dexp(-p(i)*(cutoff_start(i)/dik0-1.d0))/(cutoff_end(i)-cutoff_start(i))**3 
        br=-(p(i)/dik0)*a(i)*dexp(-p(i)*(cutoff_start(i)/dik0-1.d0))/(cutoff_end(i)-cutoff_start(i))**2
        cr=-((p(i)/dik0)**2) &
             *a(i)*dexp(-p(i)*(cutoff_start(i)/dik0-1.d0))/(cutoff_end(i)-cutoff_start(i))
        ab=-qsi(i)*dexp(-q(i)*(cutoff_start(i)/dik0-1.d0))/(cutoff_end(i)-cutoff_start(i))**3
        bb=-(q(i)/dik0)*qsi(i)*dexp(-q(i)*(cutoff_start(i)/dik0-1.d0))/(cutoff_end(i)-cutoff_start(i))**2 
        cb=-((q(i)/dik0)**2) &
             *qsi(i)*dexp(-q(i)*(cutoff_start(i)/dik0-1.d0))/(cutoff_end(i)-cutoff_start(i))
        x5(i)=(12.d0*ab-6.d0*bb+cb)/(2.d0*(cutoff_end(i)-cutoff_start(i))**2)
        x4(i)=(15.d0*ab-7.d0*bb+cb)/(cutoff_end(i)-cutoff_start(i))
        x3(i)=(20.d0*ab-8.d0*bb+cb)/2.d0
        a5(i)=(12.d0*ar-6.d0*br+cr)/(2.d0*(cutoff_end(i)-cutoff_start(i))**2)
        a4(i)=(15.d0*ar-7.d0*br+cr)/(cutoff_end(i)-cutoff_start(i))
        a3(i)=(20.d0*ar-8.d0*br+cr)/2.d0
     enddo

End Subroutine leggi_input
