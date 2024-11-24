MODULE PARAMETERS
  
  !modulo per costanti e evoluzione
  !data ultima modifica: 4 gennaio 2005
  !autore: Giulia Rossi
  
  Implicit None
  
  !COSTANTI
  Integer, parameter :: nmax=4000 
  Integer, parameter :: nvmax=150 
  
  Real*8, Parameter :: vnu=5.d11       !frequenza di termostatazione!
  
  Real*8, Parameter :: evsujoule=1.60219d-19
  Real*8, Parameter :: angsum=1.d-10
  Real*8, Parameter :: uasukg=1.d0/6.023d26
  
  Real*8, Parameter :: cbol=8.62d-05           !in eV/K
  Real*8, Parameter :: pi=3.141592653589793d0  !pi greco
  Real*8, Parameter :: rac2=1.414213562373095d0
  Real*8, Parameter :: rac3=1.732050807568877d0
  Real*8, Parameter :: rac5=2.236067977499790d0
  Real*8, Parameter :: rac8=2.828427124746190d0
  
  !assegnati dalla lettura di leggi.in:
  Integer :: irand1,irand2
  Integer :: nfix,nxfix
  Character(3) :: potential
  Real*8  :: tstep 
  Integer(kind=8) :: npas
  Integer :: nterm
  Integer :: growth_steps,evo_steps
  Real*8    :: growth_temp,evo_temp
  Real*8    :: growth_rate,photo_rate
  Real*8    :: growth_p1
  Integer :: scrivo,ndepo,nsections,nevo_ave
  Integer :: aumento_temp
  Real*8  :: tinit,tfin,timefin,timelag  !temperatura iniziale
  Character(2) :: choice_sim
  Character(2) :: thermostat
  Real*8  :: delta_temp !incremento della temperatura
  Character(2) :: boundary_cond,zboundary_cond
  Character(2) :: cutoff
 

  Real*8  :: zzz !!seme numero casuale
  Integer :: ipas
  
  !BOX DIMENSIONS
  Real*8 :: x_boxedge,y_boxedge,z_boxedge
  Real*8 :: x_boxedge_a,y_boxedge_a,z_boxedge_a !in arete units

  ! PER CALCOLARE IL POTENZIALE DEVO SAPERE CHE:
  Real*8 :: dmas(2),ecoh(2),p(3),q(3),a(3),qsi(3),arete(3)
  Real*8 :: rat(2),mass(2)
  Real*8 :: cutoff_start(3),cutoff_end(3)
  Real*8 :: cutoff_end_max
  Real*8 :: dist(3),dc1,dc2,dc3,dc4,dc5,dc33,dc55
  Real*8 :: a5(3),a4(3),a3(3),x5(3),x4(3),x3(3)
  
  ! QUANTI VICINI HANNO GLI ATOMI?
  Real*8  :: rshell
  Integer :: nv4(nmax),iv4(nvmax,nmax)
  Integer :: nvois(nmax),ivois(nvmax,nmax)
  
  ! COORDINATE, VELOCITA', FORZE:
  Real*8 :: temperatura,ener,etot,temp,ecinet,tempo
  Real*8 :: ener_dyn,ener_min
  Real*8 :: u(nmax),v(nmax),w(nmax)                           
  Real*8 :: du(nmax),dv(nmax),dw(nmax)
  Real*8 :: fx(nmax),fy(nmax),fz(nmax)             
  Real*8 :: dfx(nmax),dfy(nmax),dfz(nmax)
  Real*8 :: vx(nmax),vy(nmax),vz(nmax)
  Real*8 :: g1,g2,g3
  ! PER IL QUENCHING
  Real*8 :: pot    
  ! PER L'OUTPUT:
  Real*8 :: tempmedia,enermedia

END MODULE PARAMETERS
