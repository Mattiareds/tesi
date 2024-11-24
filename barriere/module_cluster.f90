MODULE CLUSTER

USE PARAMETERS

! modulo che contiene le caratteristiche del cluster: numero di paricelle, loro coordinate, loro specie chimica
!data ultima modifica: 20 dicenbre 2004
!autore: Giulia Rossi

integer :: natom ! il numero di particelle
integer :: itype(nmax) ! il numero della specie chimica (1 oppure 2)
integer :: ntype1,ntype2 ! il numero di atomi della specie1, della specie2
integer :: imet1, imet2 ! i numeri delle due specie (da 1 a 9)
real(8) :: x0(nmax),y0(nmax),z0(nmax) ! le coordinate del cluster in input
real(8) :: x(nmax),y(nmax),z(nmax) ! le coordinate del cluster durante l'evoluzione
double precision :: xat(nmax),yat(nmax),zat(nmax)
character(2) :: elem1,elem2,elem(nmax) ! i simboli chimici

END MODULE CLUSTER
