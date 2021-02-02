subroutine variational_montecarlo(a,dt,nmax,energy,accep)
  implicit none
  double precision, intent(in)  :: a, dt
  integer*8       , intent(in)  :: nmax 
  double precision, intent(out) :: energy, accep

  integer*8        :: istep
  integer*8        :: n_accep
  double precision :: sq_dt, chi(3), d2_old, prod, u
  double precision :: psi_old, psi_new, d2_new, argexpo, q
  double precision :: r_old(3), r_new(3)
  double precision :: d_old(3), d_new(3)

  double precision, external :: e_loc, psi

  sq_dt = dsqrt(dt)

  ! Initialization
  energy  = 0.d0
  n_accep = 0_8

  call random_gauss(r_old,3)

  call drift(a,r_old,d_old)
  d2_old  = d_old(1)*d_old(1) + &
            d_old(2)*d_old(2) + &
            d_old(3)*d_old(3)

  psi_old = psi(a,r_old)

  do istep = 1,nmax
     energy = energy + e_loc(a,r_old)

     call random_gauss(chi,3)
     r_new(:) = r_old(:) + dt*d_old(:) + chi(:)*sq_dt

     call drift(a,r_new,d_new)
     d2_new = d_new(1)*d_new(1) + &
              d_new(2)*d_new(2) + &
              d_new(3)*d_new(3)

     psi_new = psi(a,r_new)

     ! Metropolis
     prod = (d_new(1) + d_old(1))*(r_new(1) - r_old(1)) + &
            (d_new(2) + d_old(2))*(r_new(2) - r_old(2)) + &
            (d_new(3) + d_old(3))*(r_new(3) - r_old(3))

     argexpo = 0.5d0 * (d2_new - d2_old)*dt + prod

     q = psi_new / psi_old
     q = dexp(-argexpo) * q*q

     call random_number(u)

     if (u <= q) then

        n_accep = n_accep + 1_8

        r_old(:) = r_new(:)
        d_old(:) = d_new(:)
        d2_old   = d2_new
        psi_old  = psi_new

     end if

  end do

  energy = energy / dble(nmax)
  accep  = dble(n_accep) / dble(nmax)

end subroutine variational_montecarlo

program qmc
  implicit none
  double precision, parameter :: a     = 1.2d0
  double precision, parameter :: dt    = 1.0d0
  integer*8       , parameter :: nmax  = 100000
  integer         , parameter :: nruns = 30

  integer          :: irun
  double precision :: X(nruns), accep(nruns)
  double precision :: ave, err

  do irun=1,nruns
     call variational_montecarlo(a,dt,nmax,X(irun),accep(irun))
  enddo

  call ave_error(X,nruns,ave,err)
  print *, 'E = ', ave, '+/-', err

  call ave_error(accep,nruns,ave,err)
  print *, 'A = ', ave, '+/-', err

end program qmc
