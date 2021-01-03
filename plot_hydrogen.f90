program plot
  implicit none
  double precision, external :: e_loc

  double precision :: x(50), energy, dx, r(3), a(6)
  integer :: i, j

  a = (/ 0.1d0, 0.2d0, 0.5d0, 1.d0, 1.5d0, 2.d0 /)

  dx = 10.d0/(size(x)-1)
  do i=1,size(x)
     x(i) = -5.d0 + (i-1)*dx
  end do

  r(:) = 0.d0

  do j=1,size(a)
     print *, '# a=', a(j)
     do i=1,size(x)
        r(1) = x(i)
        energy = e_loc( a(j), r )
        print *, x(i), energy
     end do
     print *, ''
     print *, ''
  end do

end program plot
