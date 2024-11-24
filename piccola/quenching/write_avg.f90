subroutine write_avg
  
  USE PARAMETERS
  USE CLUSTER
  
  Implicit None

  tempo=ipas*tstep*1.d12
  open(4100, file='medie.out', status='unknown', position='append')
  write(4100,'(i6,2x,i12,1x,f18.4,4x,f12.6,3x,f18.6,3x,f18.10)') natom,ipas,tempo,&
  &temperatura,tempmedia,enermedia   
  close(4100)
  

end subroutine write_avg
