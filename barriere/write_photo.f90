subroutine write_photo
  
  USE PARAMETERS
  USE CLUSTER
  
  Implicit None
  
  Integer :: i,temp_intera
  Integer,save :: uni
  Character(18) :: nomefile !output file name p.#####t####
  
  temp_intera=int(temperatura)

  !output only after the thermalization steps, or at ipas=1
  if (ipas .eq. 1) then 
     uni=0
  endif

  if (mod(ipas,scrivo) .eq. 0) then
     uni=uni+1 !logic unit for output files p.#####t####
  endif
    
  !building up the name of the output file p.#####t####  
  nomefile(1:2)='p.'
  nomefile(8:8)='t'
  nomefile(13:13)='n'

  if (mod(ipas,scrivo) .eq. 0) then
     nomefile(3:7)='00000'
     write(nomefile(3:7),'(i5.5)')uni
     nomefile(9:12)='0000'
     write(nomefile(9:12),'(i4.4)')temp_intera
     nomefile(14:18)='00000'
     write(nomefile(14:18),'(i5.5)')natom     
     open (uni, file=nomefile, status='unknown')
     write(uni,*)natom
     write(uni,'(a2,1x,a2,2x,f16.6,f12.4)')elem1, elem2, ener, temperatura      
     do i=1,natom
        write(uni,'(a2,1x,f18.12,1x,f18.12,1x,f18.12)')elem(i),x(i)+u(i)&
	&,y(i)+v(i),z(i)+w(i)
     enddo
     close(uni)
     call flush(uni)     
  endif
  
!  if(ipas .eq. npas)then
!     uni=0
!  endif

end subroutine write_photo
