subroutine metropolis_montecarlo(a,nmax,tau,energy,accep)
  implicit none
  double precision, intent(in)  :: a
  integer*8       , intent(in)  :: nmax 
  double precision, intent(in)  :: tau
  double precision, intent(out) :: energy
  double precision, intent(out) :: accep

  integer*8 :: istep

  double precision :: norm, r_old(3), r_new(3), psi_old, psi_new
  double precision :: v, ratio, n_accep
  double precision, external :: e_loc, psi, gaussian

  energy = 0.d0
  norm   = 0.d0
  n_accep = 0.d0
  call random_number(r_old)
  r_old(:) = tau * (2.d0*r_old(:) - 1.d0)
  psi_old = psi(a,r_old)
  do istep = 1,nmax
     call random_number(r_new)
     r_new(:) = r_old(:) + tau * (2.d0*r_new(:) - 1.d0)
     psi_new = psi(a,r_new)
     ratio = (psi_new / psi_old)**2
     call random_number(v)
     if (v < ratio) then
        r_old(:) = r_new(:)
        psi_old = psi_new
        n_accep = n_accep + 1.d0
     endif
     norm = norm + 1.d0
     energy = energy + e_loc(a,r_old)
  end do
  energy = energy / norm
  accep  = n_accep / norm
end subroutine metropolis_montecarlo

program qmc
  implicit none
  double precision, parameter :: a = 0.9d0
  double precision, parameter :: tau = 1.3d0
  integer*8       , parameter :: nmax = 100000
  integer         , parameter :: nruns = 30

  integer :: irun
  double precision :: X(nruns), Y(nruns)
  double precision :: ave, err

  do irun=1,nruns
     call metropolis_montecarlo(a,nmax,tau,X(irun),Y(irun))
  enddo
  call ave_error(X,nruns,ave,err)
  print *, 'E = ', ave, '+/-', err
  call ave_error(Y,nruns,ave,err)
  print *, 'A = ', ave, '+/-', err
end program qmc
