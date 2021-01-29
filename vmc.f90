subroutine variational_montecarlo(a,dt,nmax,energy)
  implicit none
  double precision, intent(in)  :: a, dt
  integer*8       , intent(in)  :: nmax 
  double precision, intent(out) :: energy

  integer*8 :: istep
  double precision :: norm, r_old(3), r_new(3), d_old(3), sq_dt, chi(3)
  double precision, external :: e_loc

  sq_dt = dsqrt(dt)
  
  ! Initialization
  energy = 0.d0
  norm   = 0.d0
  call random_gauss(r_old,3)

  do istep = 1,nmax
     call drift(a,r_old,d_old)
     call random_gauss(chi,3)
     r_new(:) = r_old(:) + dt * d_old(:) + chi(:)*sq_dt
     norm = norm + 1.d0
     energy = energy + e_loc(a,r_new)
     r_old(:) = r_new(:)
  end do
  energy = energy / norm
end subroutine variational_montecarlo

program qmc
  implicit none
  double precision, parameter :: a = 0.9
  double precision, parameter :: dt = 0.2
  integer*8       , parameter :: nmax = 100000
  integer         , parameter :: nruns = 30

  integer :: irun
  double precision :: X(nruns)
  double precision :: ave, err

  do irun=1,nruns
     call variational_montecarlo(a,dt,nmax,X(irun))
  enddo
  call ave_error(X,nruns,ave,err)
  print *, 'E = ', ave, '+/-', err
end program qmc
