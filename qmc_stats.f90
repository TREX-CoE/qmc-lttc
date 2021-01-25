subroutine ave_error(x,n,ave,err)
  implicit none
  integer, intent(in)           :: n 
  double precision, intent(in)  :: x(n) 
  double precision, intent(out) :: ave, err
  ! TODO
end subroutine ave_error

subroutine ave_error(x,n,ave,err)
  implicit none
  integer, intent(in)           :: n 
  double precision, intent(in)  :: x(n) 
  double precision, intent(out) :: ave, err
  double precision :: variance
  if (n < 1) then
     stop 'n<1 in ave_error'
  else if (n == 1) then
     ave = x(1)
     err = 0.d0
  else
     ave = sum(x(:)) / dble(n)
     variance = sum( (x(:) - ave)**2 ) / dble(n-1)
     err = dsqrt(variance/dble(n))
  endif
end subroutine ave_error

subroutine random_gauss(z,n)
  implicit none
  integer, intent(in) :: n
  double precision, intent(out) :: z(n)
  double precision :: u(n+1)
  double precision, parameter :: two_pi = 2.d0*dacos(-1.d0)
  integer :: i

  call random_number(u)
  if (iand(n,1) == 0) then
     ! n is even
     do i=1,n,2
        z(i)   = dsqrt(-2.d0*dlog(u(i))) 
        z(i+1) = z(i) * dsin( two_pi*u(i+1) )
        z(i)   = z(i) * dcos( two_pi*u(i+1) )
     end do
  else
     ! n is odd
     do i=1,n-1,2
        z(i)   = dsqrt(-2.d0*dlog(u(i))) 
        z(i+1) = z(i) * dsin( two_pi*u(i+1) )
        z(i)   = z(i) * dcos( two_pi*u(i+1) )
     end do
     z(n)   = dsqrt(-2.d0*dlog(u(n))) 
     z(n)   = z(n) * dcos( two_pi*u(n+1) )
  end if
end subroutine random_gauss
