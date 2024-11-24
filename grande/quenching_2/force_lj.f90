subroutine force_lj
use iso_c_binding!in realta` e` uno statement per fortran 2003 ma se compila non mi lamento
USE PARAMETERS  !uso il modulo di definizione dei parametri
USE CLUSTER     !uso il modulo dove defin1:natomco variabili e parametri cluster

Implicit None
Integer :: i,j,k
Integer :: itypi,itypk
Real*8 :: dik2,dik6,dik12,xik,yik,zik
Real*8 :: adden,addf
Real*8 :: energia(natom)

ener=0.d0                                                        

do i=1,natom
   dfx(i)=fx(i)
   dfy(i)=fy(i)
   dfz(i)=fz(i)
   fx(i)=0.d0
   fy(i)=0.d0
   fz(i)=0.d0
   energia(i)=0.d0
enddo

do i=1,natom
   !se l-atomo i non ha vicini, pongo la forza = 0
   do j=1,nvois(i) 
      k=ivois(j,i)  
!      write (*,*) "prima di k>i, i,k=",i,k
      if (k.gt.i) then
         itypi=itype(i)
         itypk=itype(k)
         xik=x(k)-x(i)
         yik=y(k)-y(i)
         zik=z(k)-z(i) 
         dik2=xik*xik+yik*yik+zik*zik
         dik6=dik2**3
         dik12=dik6*dik6
!     energia     
         adden=quattroepsi(itypi,itypk)*(sigma12(itypi,itypk)/dik12-sigma6(itypi,itypk)/dik6)
         energia(i)=energia(i)+adden
         energia(k)=energia(k)+adden
!     forze
         addf=6.d0*quattroepsi(itypi,itypk)*(duesigma12(itypi,itypk)/(dik12*dik2)-sigma6(itypi,itypk)/(dik6*dik2))   
         fx(i)=fx(i)-addf*xik
         fy(i)=fy(i)-addf*yik
         fz(i)=fz(i)-addf*zik
         fx(k)=fx(k)+addf*xik
         fy(k)=fy(k)+addf*yik
         fz(k)=fz(k)+addf*zik      
      endif !su k>i   
!      write (89,*) spec(i),spec(k),i,k,itype(i),itype(k)
   enddo !su j
   ener=ener+energia(i)/2.d0
!   write (90,*) i, energia(i)/2
!   write (90+i,*) i, energia(i)/2,fx(i),fy(i),fz(i)

enddo !su i

!stop

EndSubroutine force_lj
