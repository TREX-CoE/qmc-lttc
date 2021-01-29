program energy_hydrogen
  implicit none
  double precision, external :: e_loc, psi
  double precision :: x(50), w, delta, energy, dx, r(3), a(6), norm
  integer :: i, k, l, j

  a = (/ 0.1d0, 0.2d0, 0.5d0, 1.d0, 1.5d0, 2.d0 /)

  dx = 10.d0/(size(x)-1)
  do i=1,size(x)
     x(i) = -5.d0 + (i-1)*dx
  end do

  do j=1,size(a)

     ! TODO

     print *, 'a = ', a(j), '    E = ', energy
  end do

end program energy_hydrogen

program energy_hydrogen
  implicit none
  double precision, external :: e_loc, psi
  double precision :: x(50), w, delta, energy, dx, r(3), a(6), norm
  integer          :: i, k, l, j

  a = (/ 0.1d0, 0.2d0, 0.5d0, 1.d0, 1.5d0, 2.d0 /)

  dx = 10.d0/(size(x)-1)
  do i=1,size(x)
     x(i) = -5.d0 + (i-1)*dx
  end do

  delta = dx**3

  r(:) = 0.d0

  do j=1,size(a)
     energy = 0.d0
     norm   = 0.d0
     
     do i=1,size(x)
        r(1) = x(i)

        do k=1,size(x)
           r(2) = x(k)

           do l=1,size(x)
              r(3) = x(l)

              w = psi(a(j),r)
              w = w * w * delta

              energy = energy + w * e_loc(a(j), r)
              norm   = norm   + w 
           end do

        end do

     end do

     energy = energy / norm
     print *, 'a = ', a(j), '    E = ', energy
  end do

end program energy_hydrogen
