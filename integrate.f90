module integrate
 implicit none

 public :: integrate_trap,integrate_trap_log

 private

 abstract interface
  pure real function func(x)
   real, intent(in) :: x
  end function func
 end interface

 interface integrate_trap
  module procedure integrate_trap_func,integrate_trap_array
 end interface

 interface integrate_trap_log
  module procedure integrate_trap_log,integrate_trap_log_func
 end interface

contains
!-------------------------------------------------------------------
! helper routine to integrate a function using the trapezoidal rule
!-------------------------------------------------------------------
pure real function integrate_trap_func(f,xmin,xmax) result(g)
 real, intent(in) :: xmin,xmax
 procedure(func) :: f
 real :: fx,fprev,dx,x
 integer, parameter :: npts = 4096
 integer :: i

 g = 0.
 dx = (xmax-xmin)/(npts)
 fprev = f(xmin)
 do i=2,npts
    x = xmin + i*dx
    fx = f(x)
    g = g + 0.5*dx*(fx + fprev)
    fprev = fx
 enddo

end function integrate_trap_func

!-------------------------------------------------------------------
! helper routine to integrate a function using the trapezoidal rule
!-------------------------------------------------------------------
pure real function integrate_trap_array(n,x,f) result(g)
 integer, intent(in) :: n
 real, intent(in) :: x(n),f(n)
 integer :: i

 g = 0.
 do i=2,n
    g = g + 0.5*(x(i)-x(i-1))*(f(i) + f(i-1))
 enddo

end function integrate_trap_array

!--------------------------------------------------------
!+
!  Integrate function on evenly spaced logarithmic grid
!  i.e. \int f(x) dx = \int x f(x) d(ln x)
!  using trapezoidal rule
!+
!--------------------------------------------------------
real function integrate_trap_log(n,x,f) result(fint)
 integer, intent(in) :: n
 real, intent(in) :: x(:),f(:)
 real :: dlogx
 integer :: i

 fint = 0.
 if (n < 2) return
 do i=2,n
    dlogx = log(x(i)/x(i-1))
    fint = fint + 0.5*(f(i)*x(i) + f(i-1)*x(i-1))*dlogx
 enddo

end function integrate_trap_log

!--------------------------------------------------------
!+
!  Integrate function on evenly spaced logarithmic grid
!  with a function call
!+
!--------------------------------------------------------
real function integrate_trap_log_func(n,x,f) result(fint)
 integer, intent(in) :: n
 real, intent(in) :: x(:)
 procedure(func) :: f
 real :: dlogx
 integer :: i

 fint = 0.
 if (n < 2) return
 do i=2,n
    dlogx = log(x(i)/x(i-1))
    fint = fint + 0.5*(f(x(i))*x(i) + f(x(i-1))*x(i-1))*dlogx
 enddo

end function integrate_trap_log_func

end module integrate
