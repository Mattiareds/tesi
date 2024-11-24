MODULE CLUSTER

USE PARAMETERS

! modulo che contiene le caratteristiche del cluster: numero di paricelle, loro coordinate, loro specie chimica

integer :: nelem! il numero di specie chimiche, letto due volte
integer :: natom ! il numero di particelle
integer :: itype(nmax) ! il numero della specie chimica (1,...,nelem)
integer :: ntype(nelemmax)! il numero di atomi delle varie specie
real(8) :: x0(nmax),y0(nmax),z0(nmax) ! le coordinate del cluster in input
real(8) :: x(nmax),y(nmax),z(nmax) ! le coordinate del cluster durante l'evoluzione
double precision :: xat(nmax),yat(nmax),zat(nmax)
character(2) :: spec(nmax) ! specie degli atomi
character(2) :: elem(nelemmax) ! simboli chimici delle specie
END MODULE CLUSTER
