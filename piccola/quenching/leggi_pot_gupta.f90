subroutine leggi_pot_gupta(pfile)

  USE PARAMETERS  
  USE CLUSTER 
  Implicit None
  
  ! variabili locali
  Integer :: i,j
  Character(20):: pfile
  Character(2) :: sym
  Character(2) :: pippo(nelemmax)
  Real*8 :: ar,br,cr,ab,bb,cb,dik0
  
  
 
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
        read(19,*) sym,(p(i,j),j=1,i)           
     enddo   
     do i=2,nelem
        do j=1,i-1
           p(j,i)=p(i,j)
        enddo
     enddo
     write (*,*) 'p(i,j)'
     write (*,*) '             ', (pippo(j),j=1,nelem)
     do i=1,nelem
        write(*,*)  elem(i), (p(i,j),j=1,nelem)           
     enddo
     write(*,*)
     read(19,*)
     read(19,*)
     read(19,*)
     do i=1,nelem
        read(19,*) sym,(q(i,j),j=1,i)           
     enddo   
     do i=2,nelem
        do j=1,i-1
           q(j,i)=q(i,j)
        enddo
     enddo
     write(*,*) 'q(i,j)'
     write(*,*) '             ', (pippo(j),j=1,nelem)
     do i=1,nelem
        write(*,*) elem(i), (q(i,j),j=1,nelem)           
     enddo       
     write(*,*)
     read(19,*)
     read(19,*)
     read(19,*)
     do i=1,nelem
        read(19,*) sym,(a(i,j),j=1,i)           
     enddo   
     do i=2,nelem
        do j=1,i-1
           a(j,i)=a(i,j)
        enddo
     enddo
     write (*,*) 'a(i,j)'
     write (*,*) '             ', (pippo(j),j=1,nelem)
     do i=1,nelem
        write(*,*) elem(i), (a(i,j),j=1,nelem)           
     enddo  
     write(*,*)
     read(19,*)
     read(19,*)
     read(19,*)
     do i=1,nelem
        read(19,*) sym,(qsi(i,j),j=1,i)           
     enddo   
     do i=2,nelem
        do j=1,i-1
           qsi(j,i)=qsi(i,j)
        enddo
     enddo
     write (*,*) 'qsi(i,j)'
     write (*,*) '             ', (pippo(j),j=1,nelem)
     do i=1,nelem
        write(*,*) elem(i), (qsi(i,j),j=1,nelem)           
     enddo  
     write(*,*)
     read(19,*)
     read(19,*)
     read(19,*)
     do i=1,nelem
        read(19,*) sym,(rat(i,j),j=1,i)           
     enddo   
     do i=2,nelem
        do j=1,i-1
           rat(j,i)=rat(i,j)
        enddo
     enddo
     write (*,*) 'rat(i,j)'
     write (*,*) '             ', (pippo(j),j=1,nelem)
     do i=1,nelem
        write(*,*) elem(i), (rat(i,j),j=1,nelem)           
     enddo  
     write(*,*)
     read(19,*)
     read(19,*)
     read(19,*)
     do i=1,nelem
        read(19,*) sym,(cutoff_start(i,j),j=1,i)           
     enddo   
     do i=2,nelem
        do j=1,i-1
           cutoff_start(j,i)=cutoff_start(i,j)
        enddo
     enddo
     write (*,*) 'cutoff_start(i,j)'
     write (*,*) '             ', (pippo(j),j=1,nelem)
     do i=1,nelem
        write(*,*) elem(i), (cutoff_start(i,j),j=1,nelem)           
     enddo  
     write(*,*)
     read(19,*)
     read(19,*)
     read(19,*)
     do i=1,nelem
        read(19,*) sym,(cutoff_end(i,j),j=1,i)           
     enddo   
     do i=2,nelem
        do j=1,i-1
           cutoff_end(j,i)=cutoff_end(i,j)
        enddo
     enddo
     write (*,*) 'cutoff_end(i,j)'
     write (*,*) '             ', (pippo(j),j=1,nelem)
     do i=1,nelem
        write(*,*) elem(i), (cutoff_end(i,j),j=1,nelem)           
     enddo  
     write (*,*)     
     write (*,*) 'END READING FILE ', pfile
  close(19)
!stop
  !Calculation of useful constants
do i=1,nelem
   do j=1,nelem
      nn(i,j)=2.d0*rat(i,j)
      
      dik0=nn(i,j)
        ar=-a(i,j)*dexp(-p(i,j)*(cutoff_start(i,j)/dik0-1.d0))/(cutoff_end(i,j)-cutoff_start(i,j))**3 
        br=-(p(i,j)/dik0)*a(i,j)*dexp(-p(i,j)*(cutoff_start(i,j)/dik0-1.d0))/(cutoff_end(i,j)-cutoff_start(i,j))**2
        cr=-((p(i,j)/dik0)**2) &
             *a(i,j)*dexp(-p(i,j)*(cutoff_start(i,j)/dik0-1.d0))/(cutoff_end(i,j)-cutoff_start(i,j))
        ab=-qsi(i,j)*dexp(-q(i,j)*(cutoff_start(i,j)/dik0-1.d0))/(cutoff_end(i,j)-cutoff_start(i,j))**3
        bb=-(q(i,j)/dik0)*qsi(i,j)*dexp(-q(i,j)*(cutoff_start(i,j)/dik0-1.d0))/(cutoff_end(i,j)-cutoff_start(i,j))**2 
        cb=-((q(i,j)/dik0)**2) &
             *qsi(i,j)*dexp(-q(i,j)*(cutoff_start(i,j)/dik0-1.d0))/(cutoff_end(i,j)-cutoff_start(i,j))
        x5(i,j)=(12.d0*ab-6.d0*bb+cb)/(2.d0*(cutoff_end(i,j)-cutoff_start(i,j))**2)
        x4(i,j)=(15.d0*ab-7.d0*bb+cb)/(cutoff_end(i,j)-cutoff_start(i,j))
        x3(i,j)=(20.d0*ab-8.d0*bb+cb)/2.d0
        a5(i,j)=(12.d0*ar-6.d0*br+cr)/(2.d0*(cutoff_end(i,j)-cutoff_start(i,j))**2)
        a4(i,j)=(15.d0*ar-7.d0*br+cr)/(cutoff_end(i,j)-cutoff_start(i,j))
        a3(i,j)=(20.d0*ar-8.d0*br+cr)/2.d0

   enddo
enddo
   
write (*,*)
!stop
End Subroutine leggi_pot_gupta
