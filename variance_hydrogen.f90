program variance_hydrogen
  implicit none

  double precision :: x(50), w, delta, energy, energy2
  double precision :: dx, r(3), a(6), norm, e_tmp, s2
  integer          :: i, k, l, j

  double precision, external :: e_loc, psi

  a = (/ 0.1d0, 0.2d0, 0.5d0, 1.d0, 1.5d0, 2.d0 /)

  dx = 10.d0/(size(x)-1)
  do i=1,size(x)
     x(i) = -5.d0 + (i-1)*dx
  end do

  delta = dx**3

  r(:) = 0.d0

  do j=1,size(a)
     energy  = 0.d0
     energy2 = 0.d0
     norm    = 0.d0

     do i=1,size(x)
        r(1) = x(i)

        do k=1,size(x)
           r(2) = x(k)

           do l=1,size(x)
              r(3) = x(l)

              w = psi(a(j),r)
              w = w * w * delta

              e_tmp = e_loc(a(j), r)

              energy  = energy  + w * e_tmp
              energy2 = energy2 + w * e_tmp * e_tmp
              norm   = norm     + w 
           end do

        end do

     end do

     energy  = energy  / norm
     energy2 = energy2 / norm

     s2 = energy2 - energy*energy

     print *, 'a = ', a(j), ' E = ', energy, ' s2 = ', s2
  end do

end program variance_hydrogen
