subroutine uniform_montecarlo(a,nmax,energy)
  implicit none
  double precision, intent(in)  :: a
  integer*8       , intent(in)  :: nmax 
  double precision, intent(out) :: energy

  integer*8        :: istep
  double precision :: norm, r(3), w

  double precision, external :: e_loc, psi

  energy = 0.d0
  norm   = 0.d0

  do istep = 1,nmax

     call random_number(r)
     r(:) = -5.d0 + 10.d0*r(:)

     w = psi(a,r)
     w = w*w

     energy = energy + w * e_loc(a,r)
     norm   = norm   + w

  end do

  energy = energy / norm

end subroutine uniform_montecarlo

program qmc
  implicit none
  double precision, parameter :: a     = 0.9
  integer*8       , parameter :: nmax  = 100000
  integer         , parameter :: nruns = 30

  integer          :: irun
  double precision :: X(nruns)
  double precision :: ave, err

  do irun=1,nruns
     call uniform_montecarlo(a, nmax, X(irun))
  enddo

  call ave_error(X, nruns, ave, err)

  print *, 'E = ', ave, '+/-', err
end program qmc
