Subroutine bigvoi

  USE PARAMETERS  !uso il modulo di definizione dei parametri
  USE CLUSTER     !uso il modulo dove defin1:natomco variabili e parametri cluster

  Implicit None
  Integer :: i,j
  Real*8 :: dij2,xij,yij,zij
  Real*8 :: bigvoi_cage

  nv4(1:natom)=0.d0
  iv4(:,:)=0
  !bigvoi_cage=(arete(1)*cutoff_end+rshell)**2 !bigvoi works with Angstrom coordinates
  bigvoi_cage=(cutoff_end_max+rshell)**2
    
  do i=1,natom-1                                              
     do j=i+1,natom                                                  
        xij=x(j)+u(j)-x(i)-u(i)                                           
        yij=y(j)+v(j)-y(i)-v(i)                                           
        zij=z(j)+w(j)-z(i)-w(i)
!        xij=distanzex(i,j)
!        yij=distanzey(i,j)
!        zij=distanzez(i,j)        

        !periodic boundary conditions
!        if(boundary_cond .eq. 'si')then
!           if (abs(xij) .gt. x_boxedge/2.) then
!              if (xij .gt. 0.)then
!                 xij=xij-x_boxedge
!              else
!                 xij=xij+x_boxedge
!              endif
!         endif
!         if (abs(yij) .gt. y_boxedge/2.) then
!            if (yij .gt. 0.)then
!               yij=yij-y_boxedge
!            else
!               yij=yij+y_boxedge
!            endif
!         endif
!         if(zboundary_cond .eq. 'si')then
!           if (abs(zij) .gt. z_boxedge/2.) then
!              if (zij .gt. 0.)then
!                 zij=zij-z_boxedge
!              else
!                 zij=zij+z_boxedge
!              endif
!           endif
!         endif
!         distanzex(i,j)=xij
!         distanzey(i,j)=yij
!         distanzez(i,j)=zij

!      endif
      

      dij2=xij*xij+yij*yij+zij*zij
     
      
      if (dij2.lt.bigvoi_cage) then
         nv4(i)=nv4(i)+1 
         iv4(nv4(i),i)=j
      endif
      
!       if((j .eq. natom) .and. (ipas .gt. nterm+growth_steps))then
      !write(*,*)'In bigvoi the distance between the last atom and ',i,'is: ',sqrt(dij2)
      !write(*,*)'In bigvoi nv4(natom)=',nv4(natom)
!      endif
      
      
   enddo
   
!      if (nv4(i).ge.nmax) then                                        
!      write (*,'(i4,1x,f6.3,a18)'),i,nv4(i),'too many neighbours in bigvoi.f90'
!      write(*,*)'Program stops here.'
!      stop                                                          
!   endif
   
enddo

!check:
!write(*,*)'bigvoi:'
!do i=1,natom
!   write(*,*)i,'ha ', nv4(i),' vicini'
!   write(*,*)'lista:'
!   do j=1,nv4(i)
!      write(*,*)iv4(j,i)
!   enddo
!enddo
!stop


!making the distance matrix symmetric:
!do i=1,natom-1                                                   
!   do j=i+1,natom        
!      distanzex(j,i)=-distanzex(i,j)   
!      distanzey(j,i)=-distanzey(i,j)   
!      distanzez(j,i)=-distanzez(i,j)   
!   enddo
!enddo

EndSubroutine bigvoi
