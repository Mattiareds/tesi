subroutine leggi_pot_morse(pfile)

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
     read (19,*)
     read (19,*) (mass(i),i=1,nelem)
     write (*,*) "Masses " 
     write (*,*) '             ', (pippo(j),j=1,nelem)
     write(*,*) ' ', (mass(i),i=1,nelem)
     write(*,*)
     read(19,*)
     read(19,*)
     read(19,*)
     do i=1,nelem
        read(19,*) sym,(u0(i,j),j=1,i)           
     enddo   
     do i=2,nelem
        do j=1,i-1
           u0(j,i)=u0(i,j)
        enddo
     enddo
     write (*,*) 'u0(i,j)'
     write (*,*) '             ', (pippo(j),j=1,nelem)
     do i=1,nelem
        write(*,*)  elem(i), (u0(i,j),j=1,nelem)           
     enddo
     write(*,*)
     read(19,*)
     read(19,*)
     read(19,*)
     do i=1,nelem
        read(19,*) sym,(r0(i,j),j=1,i)           
     enddo   
     do i=2,nelem
        do j=1,i-1
           r0(j,i)=r0(i,j)
        enddo
     enddo
     write(*,*) 'r0(i,j)'
     write(*,*) '             ', (pippo(j),j=1,nelem)
     do i=1,nelem
        write(*,*) elem(i), (r0(i,j),j=1,nelem)           
     enddo       
     write(*,*)
     read(19,*)
     read(19,*)
     read(19,*)
     do i=1,nelem
        read(19,*) sym,(alpha(i,j),j=1,i)           
     enddo   
     do i=2,nelem
        do j=1,i-1
           alpha(j,i)=alpha(i,j)
        enddo
     enddo
     write (*,*) 'alpha(i,j)'
     write (*,*) '             ', (pippo(j),j=1,nelem)
     do i=1,nelem
        write(*,*) elem(i), (alpha(i,j),j=1,nelem)           
     enddo  
     write (*,*)     
     write (*,*) 'END READING FILE ', pfile
  close(19)


write (*,*)
!stop
End Subroutine leggi_pot_morse
