MODULE PARAMETERS
    
  Implicit None
  
  !COSTANTI
  Integer, parameter :: nelemmax=10 
  Integer, parameter :: nmax=3000 
  Integer, parameter :: nvmax=3000
  
!  Real*8, Parameter :: evsujoule=1.60219d-19
!  Real*8, Parameter :: angsum=1.d-10
!  Real*8, Parameter :: uasukg=1.d0/6.023d26
  
!  Real*8, Parameter :: cbol=8.62d-05           !in eV/K
 Real*8, Parameter :: evsujoule=1.602176634d-19
 Real*8, Parameter :: angsum=1.d-10
 Real*8, Parameter :: uasukg=1.660538921d-27
  
 Real*8, Parameter :: cbol=8.6173303d-05           !in eV/K
  Real*8, Parameter :: pi=3.141592653589793d0  !pi greco
  
  !assegnati dalla lettura di leggi.in:
  Integer :: irand1,irand2
  Integer :: nfix,nxfix
  Character(3) :: potential
  Real*8  :: tstep 
  Real*8  :: vnu       !frequenza di termostatazione!

  Integer(kind=8) :: npas
  Integer(kind=8) :: ipas
  Integer(kind=8) :: growth_steps,evo_steps
  Integer(kind=8) :: scrivo,ndepo,nsections,nevo_ave
  Integer(kind=8) :: aumento_temp
  Real*8    :: growth_temp,evo_temp,que_temp
  Real*8    :: growth_rate,photo_rate
  Real*8    :: growth_p1(nelemmax),sum_p1(0:nelemmax)
  Real*8  :: delta_temp !incremento della temperatura
  Real*8  :: tinit,tfin,timefin,timefinqq,timelag  !temperatura iniziale
  Character(2) :: choice_sim
  Character(2) :: thermostat
  Character(2) :: cutoff
 

  Real*8  :: zzz !!seme numero casuale
  
  ! PER IL CALCOLO INTERAZIONI:
  !Lennard-Jones
  Real*8 :: epsi(nelemmax,nelemmax),quattroepsi(nelemmax,nelemmax)
  Real*8 :: rmin(nelemmax,nelemmax),sigma(nelemmax,nelemmax),sigma6(nelemmax,nelemmax)
  Real*8 :: sigma12(nelemmax,nelemmax),duesigma12(nelemmax,nelemmax)
  ! Morse
  Real*8 :: u0(nelemmax,nelemmax),r0(nelemmax,nelemmax),alpha(nelemmax,nelemmax)
  !Gupta
  Real*8 :: p(nelemmax,nelemmax),q(nelemmax,nelemmax),a(nelemmax,nelemmax),qsi(nelemmax,nelemmax)
  Real*8 :: rat(nelemmax,nelemmax),nn(nelemmax,nelemmax)
  Real*8 :: a5(nelemmax,nelemmax),a4(nelemmax,nelemmax),a3(nelemmax,nelemmax)
  Real*8 :: x5(nelemmax,nelemmax),x4(nelemmax,nelemmax),x3(nelemmax,nelemmax)
  !for all potentials  
  Real*8 :: mass(nelemmax)
  Real*8 :: cutoff_start(nelemmax,nelemmax),cutoff_end(nelemmax,nelemmax)
  Real*8 :: cutoff_end_max
  
  ! QUANTI VICINI HANNO GLI ATOMI?
  Real*8  :: rshell
  Integer :: nv4(nmax),iv4(nvmax,nmax)
  Integer :: nvois(nmax),ivois(nvmax,nmax)
  
  ! COORDINATE, VELOCITA', FORZE:
  Real*8 :: temperatura,ener,etot,temp,ecinet,tempo
  Real*8 :: ener_min
  Real*8 :: dx(nmax),dy(nmax),dz(nmax)
  Real*8 :: fx(nmax),fy(nmax),fz(nmax)             
  Real*8 :: dfx(nmax),dfy(nmax),dfz(nmax)
  Real*8 :: vx(nmax),vy(nmax),vz(nmax)
  Real*8 :: g1,g2,g3
  ! PER IL QUENCHING
  Real*8 :: pot    
  ! PER L'OUTPUT:
  Real*8 :: tempmedia,enermedia

END MODULE PARAMETERS
