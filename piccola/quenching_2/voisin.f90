subroutine voisin

USE PARAMETERS  !uso il modulo di definizione dei parametri
USE CLUSTER     !uso il modulo dove defin1:natomco variabili e parametri cluster

Implicit None
Integer :: i,jv,j,ik,jk
Real*8 :: dij2,xij,yij,zij
Real*8 :: cutoff_end2(nelem,nelem)
 do i=1,nelem
    do j=1,nelem
       cutoff_end2(i,j)=cutoff_end(i,j)*cutoff_end(i,j)! !voisin.f90 works with Angstrom coordinates
!       write (*,*) i,j, cutoff_end(i,j)
    enddo
 enddo
! stop
 
do i=1,natom
   nvois(i)=0
enddo

do i=1,natom-1                                                   
   do jv=1,nv4(i)                                                   
      j=iv4(jv,i)
      ik=itype(i)
      jk=itype(j)                         
      xij=x(j)-x(i)
      yij=y(j)-y(i)
      zij=z(j)-z(i)       
      dij2=xij*xij+yij*yij+zij*zij
      !neighbor list update for force_lj.f90
      if (dij2.lt.cutoff_end2(ik,jk)) then
         nvois(i)=nvois(i)+1
         nvois(j)=nvois(j)+1
         ivois(nvois(i),i)=j
         ivois(nvois(j),j)=i
      endif      
   enddo      !fine ciclo su jv
enddo       !fine ciclo su i

!do i=1,natom
!   write (*,*) i, nvois(i)
!enddo


!stop
EndSubroutine voisin
