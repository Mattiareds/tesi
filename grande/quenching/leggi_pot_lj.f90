subroutine leggi_pot_lj(pfile)

  USE PARAMETERS  
  USE CLUSTER 
  Implicit None
  
  ! variabili locali
  Integer :: i,j
  Character(20):: pfile
  Character(2) :: sym
  Character(2) :: pippo(nelemmax)

  
!*************************************************************
!READING POTENTIAL PARAMETERS
!*************************************************************
  write (*,*) "START READING FILE ", pfile 
  write (*,*)
 
   cutoff='no'
   open(19,file=pfile,status='old') 
     read(19,*) nelem
     read(19,*)
     read(19,*) (elem(i),i=1,nelem)
     do i=1,nelem
        pippo(i)=elem(i) 
     enddo
     write(*,*) "System " 
     write(*,*) (pippo(i),i=1,nelem) 
     write(*,*)
     read(19,*)
     read(19,*) (mass(i),i=1,nelem)
     write(*,*) "Masses "
     write (*,*) '             ', (pippo(j),j=1,nelem)
     write(*,*) ' ', (mass(i),i=1,nelem)
     write (*,*)
     read(19,*)
     read(19,*)
     read(19,*)
     do i=1,nelem
        read(19,*) sym,(epsi(i,j),j=1,i)           
     enddo   
     do i=2,nelem
        do j=1,i-1
           epsi(j,i)=epsi(i,j)
        enddo
     enddo
     write (*,*) 'epsi(i,j)'
     write (*,*) '             ', (pippo(j),j=1,nelem)
     do i=1,nelem
        write(*,*) elem(i), (epsi(i,j),j=1,nelem)           
     enddo
     write (*,*)
     read(19,*)
     read(19,*)
     read(19,*)
     do i=1,nelem
        read(19,*) sym,(rmin(i,j),j=1,i)           
     enddo   
     do i=2,nelem
        do j=1,i-1
           rmin(j,i)=rmin(i,j)
        enddo
     enddo
     write(*,*) 'rmin(i,j)'
     write(*,*) '             ', (pippo(j),j=1,nelem)
     do i=1,nelem
        write(*,*) elem(i), (rmin(i,j),j=1,nelem)           
     enddo       
     write (*,*)     
     write (*,*) 'END READING FILE ', pfile
  close(19)

  !Calculation of useful constants
     do i=1,nelem
        do j=1,nelem
        sigma(i,j)=rmin(i,j)/2.d0**(1.d0/6.d0)
        sigma6(i,j)=sigma(i,j)**6
        sigma12(i,j)=sigma6(i,j)**2
        duesigma12(i,j)=2.d0*sigma12(i,j)
        quattroepsi(i,j)=4.d0*epsi(i,j)
     enddo
   enddo

write (*,*)
!stop
End Subroutine leggi_pot_lj
