double precision function potential(r)
  implicit none
  double precision, intent(in) :: r(3)
  potential = -1.d0 / dsqrt( r(1)*r(1) + r(2)*r(2) + r(3)*r(3) )
end function potential

double precision function psi(a, r)
  implicit none
  double precision, intent(in) :: a, r(3)
  psi = dexp(-a * dsqrt( r(1)*r(1) + r(2)*r(2) + r(3)*r(3) ))
end function psi

double precision function kinetic(a,r)
  implicit none
  double precision, intent(in) :: a, r(3)
  kinetic = -0.5d0 * (a*a - (2.d0*a) / &
       dsqrt( r(1)*r(1) + r(2)*r(2) + r(3)*r(3) ) ) 
end function kinetic

double precision function e_loc(a,r)
  implicit none
  double precision, intent(in) :: a, r(3)
  double precision, external   :: kinetic, potential
  e_loc = kinetic(a,r) + potential(r)
end function e_loc

subroutine drift(a,r,b)
  implicit none
  double precision, intent(in)  :: a, r(3)
  double precision, intent(out) :: b(3)
  double precision :: ar_inv
  ar_inv = -a / dsqrt(r(1)*r(1) + r(2)*r(2) + r(3)*r(3))
  b(:) = r(:) * ar_inv
end subroutine drift
