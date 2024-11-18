subroutine bigvoi
  USE PARAMETERS  !uso il modulo di definizione dei parametri
  USE CLUSTER     !uso il modulo dove defin1:natomco variabili e parametri cluster

  Implicit None
  Integer :: i,j
  Real*8 :: dij2,xij,yij,zij
  Real*8 :: bigvoi_cage

  nv4(1:natom)=0.d0
  iv4(:,:)=0
  bigvoi_cage=(cutoff_end_max+rshell)**2
!  write (*,*) 'bigvoi_cage=', bigvoi_cage
!  stop
    
  do i=1,natom-1                                              
     do j=i+1,natom                                                  
        xij=x(j)-x(i)                                          
        yij=y(j)-y(i)                                           
        zij=z(j)-z(i)             
        dij2=xij*xij+yij*yij+zij*zij
     
      
      if (dij2.lt.bigvoi_cage) then
         nv4(i)=nv4(i)+1 
         iv4(nv4(i),i)=j
      endif
      
!     
   enddo
enddo

EndSubroutine bigvoi
