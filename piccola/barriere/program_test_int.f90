program test_int
  integer :: i = 42
  real*8 :: aa=12366616618818188.d0
  complex :: z = (-3.7, 1.0)
  print *, int(i)
  print *, int(aa), int(aa,8)
end program
