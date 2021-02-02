subroutine metropolis_montecarlo(a,nmax,dt,energy,accep)
  implicit none
  double precision, intent(in)  :: a
  integer*8       , intent(in)  :: nmax 
  double precision, intent(in)  :: dt
  double precision, intent(out) :: energy
  double precision, intent(out) :: accep

  double precision :: r_old(3), r_new(3), psi_old, psi_new
  double precision :: v, ratio
  integer*8        :: n_accep
  integer*8        :: istep

  double precision, external :: e_loc, psi, gaussian

  energy  = 0.d0
  n_accep = 0_8

  call random_number(r_old)
  r_old(:) = dt * (2.d0*r_old(:) - 1.d0)
  psi_old = psi(a,r_old)

  do istep = 1,nmax
     energy = energy + e_loc(a,r_old)

     call random_number(r_new)
     r_new(:) = r_old(:) + dt*(2.d0*r_new(:) - 1.d0)

     psi_new = psi(a,r_new)

     ratio = (psi_new / psi_old)**2
     call random_number(v)

     if (v <= ratio) then

        n_accep = n_accep + 1_8

        r_old(:) = r_new(:)
        psi_old = psi_new

     endif

  end do

  energy = energy / dble(nmax)
  accep  = dble(n_accep) / dble(nmax)

end subroutine metropolis_montecarlo

program qmc
  implicit none
  double precision, parameter :: a = 1.2d0
  double precision, parameter :: dt = 1.0d0
  integer*8       , parameter :: nmax = 100000
  integer         , parameter :: nruns = 30

  integer          :: irun
  double precision :: X(nruns), Y(nruns)
  double precision :: ave, err

  do irun=1,nruns
     call metropolis_montecarlo(a,nmax,dt,X(irun),Y(irun))
  enddo

  call ave_error(X,nruns,ave,err)
  print *, 'E = ', ave, '+/-', err

  call ave_error(Y,nruns,ave,err)
  print *, 'A = ', ave, '+/-', err

end program qmc
