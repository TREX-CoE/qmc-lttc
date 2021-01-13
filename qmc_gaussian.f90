double precision function gaussian(r)
  implicit none
  double precision, intent(in) :: r(3)
  double precision, parameter :: norm_gauss = 1.d0/(2.d0*dacos(-1.d0))**(1.5d0)
  gaussian = norm_gauss * dexp( -0.5d0 * (r(1)*r(1) + r(2)*r(2) + r(3)*r(3) ))
end function gaussian


subroutine gaussian_montecarlo(a,nmax,energy)
  implicit none
  double precision, intent(in)  :: a
  integer*8       , intent(in)  :: nmax 
  double precision, intent(out) :: energy

  integer*8 :: istep

  double precision :: norm, r(3), w

  double precision, external :: e_loc, psi, gaussian

  energy = 0.d0
  norm   = 0.d0
  do istep = 1,nmax
     call random_gauss(r,3)
     w = psi(a,r) 
     w = w*w / gaussian(r)
     norm = norm + w
     energy = energy + w * e_loc(a,r)
  end do
  energy = energy / norm
end subroutine gaussian_montecarlo

program qmc
  implicit none
  double precision, parameter :: a = 0.9
  integer*8       , parameter :: nmax = 100000
  integer         , parameter :: nruns = 30

  integer :: irun
  double precision :: X(nruns)
  double precision :: ave, err

  do irun=1,nruns
     call gaussian_montecarlo(a,nmax,X(irun))
  enddo
  call ave_error(X,nruns,ave,err)
  print *, 'E = ', ave, '+/-', err
end program qmc
