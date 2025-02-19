program loop_loop
implicit none

integer imin, imax
real i, root

imin = 1
imax = 10

do i = imin, imax
  root = sqrt(i)
  write(*,*) "Sqrt of", i, ":", root
enddo

write(*,*) "loop loop!"

end