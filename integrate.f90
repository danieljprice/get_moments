module integrate
 implicit none

 public :: integrate_trap

 private

 abstract interface
  pure real function func(x)
   real, intent(in) :: x
  end function func
 end interface

 interface integrate_trap
  module procedure integrate_trap_func,integrate_trap_array
 end interface

contains
!-------------------------------------------------------------------
! helper routine to integrate a function using the trapezoidal rule
!-------------------------------------------------------------------
pure real function integrate_trap_func(f,xmin,xmax) result(g)
 real, intent(in) :: xmin,xmax
 procedure(func) :: f
 real :: fx,fprev,dx,x
 integer, parameter :: npts = 128
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

end module integrate
