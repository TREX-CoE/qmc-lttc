subroutine ave_error(x,n,ave,err)
  implicit none
  integer, intent(in)           :: n 
  double precision, intent(in)  :: x(n) 
  double precision, intent(out) :: ave, err
  double precision :: variance
  if (n == 1) then
     ave = x(1)
     err = 0.d0
  else
     ave = sum(x(:)) / dble(n)
     variance = sum( (x(:) - ave)**2 ) / dble(n-1)
     err = dsqrt(variance/dble(n))
  endif
end subroutine ave_error
