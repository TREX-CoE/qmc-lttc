subroutine variational_montecarlo(a,tau,nmax,energy,accep_rate)
  implicit none
  double precision, intent(in)  :: a, tau
  integer*8       , intent(in)  :: nmax 
  double precision, intent(out) :: energy, accep_rate

  integer*8 :: istep
  double precision :: sq_tau, chi(3), d2_old, prod, u
  double precision :: psi_old, psi_new, d2_new, argexpo, q
  double precision :: r_old(3), r_new(3)
  double precision :: d_old(3), d_new(3)
  double precision, external :: e_loc, psi

  sq_tau = dsqrt(tau)

  ! Initialization
  energy = 0.d0
  accep_rate = 0.d0
  call random_gauss(r_old,3)
  call drift(a,r_old,d_old)
  d2_old = d_old(1)*d_old(1) + d_old(2)*d_old(2) + d_old(3)*d_old(3)
  psi_old = psi(a,r_old)

  do istep = 1,nmax
     call random_gauss(chi,3)
     r_new(:) = r_old(:) + tau * d_old(:) + chi(:)*sq_tau
     call drift(a,r_new,d_new)
     d2_new = d_new(1)*d_new(1) + d_new(2)*d_new(2) + d_new(3)*d_new(3)
     psi_new = psi(a,r_new)
     ! Metropolis
     prod = (d_new(1) + d_old(1))*(r_new(1) - r_old(1)) + &
            (d_new(2) + d_old(2))*(r_new(2) - r_old(2)) + &
            (d_new(3) + d_old(3))*(r_new(3) - r_old(3))
     argexpo = 0.5d0 * (d2_new - d2_old)*tau + prod
     q = psi_new / psi_old
     q = dexp(-argexpo) * q*q
     call random_number(u)
     if (u<q) then
        accep_rate = accep_rate + 1.d0
        r_old(:) = r_new(:)
        d_old(:) = d_new(:)
        d2_old = d2_new
        psi_old = psi_new
     end if
     energy = energy + e_loc(a,r_old)
  end do
  energy = energy / dble(nmax)
  accep_rate = dble(accep_rate) / dble(nmax)
end subroutine variational_montecarlo

program qmc
  implicit none
  double precision, parameter :: a = 0.9
  double precision, parameter :: tau = 1.0
  integer*8       , parameter :: nmax = 100000
  integer         , parameter :: nruns = 30

  integer :: irun
  double precision :: X(nruns), accep(nruns)
  double precision :: ave, err

  do irun=1,nruns
     call variational_montecarlo(a,tau,nmax,X(irun),accep(irun))
  enddo
  call ave_error(X,nruns,ave,err)
  print *, 'E = ', ave, '+/-', err
  call ave_error(accep,nruns,ave,err)
  print *, 'A = ', ave, '+/-', err
end program qmc
