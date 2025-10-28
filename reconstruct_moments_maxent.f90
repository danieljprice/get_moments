module reconstruct_maxent
 implicit none
 real, parameter :: pi = 4.*atan(1.)
 integer, parameter :: rk = kind(0.d0)

 public :: reconstruct,reconstruct_moments_maxent,integrand

 logical, private :: use_log_grid = .false.
 logical, private :: log_grid_is_uniform = .false.
 private

contains
!
! reconstruct a function from a limited number
! of specified moments mu_n, corresponding to
!
!  mu_n  = \int f(x) x^n dx
!
! Arguments:
!  IN:
!    mu(:) : array of moments of arbitrary length
!    x(:)  : grid of positions on which to evaluate function
!  OUT:
!    f(size(x)) : reconstructed function evaluated
!                 on grid positions
!    lambsol(size(mu)) : parameters used in max entropy reconstruction
!    ierr : error code (>0 if something went wrong)
!
subroutine reconstruct(mu,x,f,lambsol,ierr,lambguess)
 real, intent(in) :: mu(:)
 real, intent(in) :: x(:)
 real, intent(out) :: f(size(x))
 real, intent(out) :: lambsol(size(mu))
 integer, intent(out) :: ierr
 real, intent(in), optional :: lambguess(size(mu))
 real :: lsum(size(mu))

 ! in principle can choose method here
 if (present(lambguess)) then
    call reconstruct_moments_maxent(mu,x,f,lambsol,lsum,ierr,lambguess)
 else
    call reconstruct_moments_maxent(mu,x,f,lambsol,lsum,ierr)
 endif

end subroutine reconstruct

!
! maximum entropy reconstruction of function given moments
!
subroutine reconstruct_moments_maxent(mu,x,f,lambsol,lsum,ierr,lambguess,use_log,weights)
 use integrate,  only:check_log_grid_is_uniform
 use fsolve_mod, only:fsolve
 real, intent(in) :: mu(:)
 real, intent(in) :: x(:)
 real, intent(out) :: f(size(x))
 real, intent(out) :: lambsol(size(mu))
 real, intent(out) :: lsum(size(mu))
 integer, intent(out) :: ierr
 real, intent(in), optional :: lambguess(size(mu))
 real, intent(in), optional :: weights(size(x))
 logical, intent(in), optional :: use_log
 integer :: n_moments
 real, parameter :: tol = 1.e-10
 real, allocatable :: logx(:)
 logical :: use_gauss_legendre = .true.


 lambsol = 0.
 ierr = 0
 f = 0.
 n_moments = size(mu)

 ! initial guesses for Lagrange multipliers
 if (present(lambguess)) then
    ! use guesses passed as arguments (e.g. from previous attempt)
    lambsol = lambguess
    !print*,' using lambguess = ',lambguess
 else
    ! use predefined guess
    lambsol(1) = -log(sqrt(2.*pi))
 endif

 !
 ! if we are using a log grid, precalculate the grid
 ! of logx (this is an optimisation)
 !
 use_log_grid = .false.
 if (present(use_log)) use_log_grid = use_log
 if (use_log_grid) then
    log_grid_is_uniform = check_log_grid_is_uniform(x)
    allocate(logx(size(x)))
    logx = log(x)
 endif
 !
 ! if pre-calculated weights are sent, use Gauss-Legendre quadrature
 !
 use_gauss_legendre = present(weights)

 call fsolve(residual,n_moments,lambsol,lsum,tol,ierr)
 f = integrand(x,lambsol,n_moments)

contains
!
!  residual of the moment approximation function
!
!    Calculates the residual of the moment approximation function.
!
!    Parameters:
!        lamb (array): an array of Lagrange constants used to approximate the distribution
!        x (array):
!        k (integer): order of the moment
!        mu (array): an array of the known moments needed to approximate the distribution function
!
!    Returns:
!        rhs: the integrated right hand side of the moment approximation function
!
subroutine residual(n,lamb,l_sum)
 use integrate, only:integrate_trap,integrate_trap_log,integrate_trap_log_uniform,integrate_gauss_legendre
 integer, intent(in) :: n
 real ( kind = rk ), intent(in)  :: lamb(n) ! guess for  solution
 real ( kind = rk ), intent(out) :: l_sum(n) ! function evaluated for given lamb
 integer :: k

 if (use_gauss_legendre) then
    do k=1,n
       l_sum(k) = integrate_gauss_legendre(size(x),weights,integrand(x,lamb,k-1,logx=logx)) - mu(k)
    enddo
 else
    if (use_log_grid) then
       if (log_grid_is_uniform) then
          do k=1,n
             l_sum(k) = integrate_trap_log_uniform(size(x),x,integrand(x,lamb,k-1,logx=logx)) - mu(k)
          enddo
       else
          do k=1,n
             l_sum(k) = integrate_trap_log(size(x),x,integrand(x,lamb,k-1,logx=logx)) - mu(k)
          enddo
       endif
    else
       do k=1,n
          l_sum(k) = integrate_trap(size(x),x,integrand(x,lamb,k-1)) - mu(k)
       enddo
    endif
 endif
 !print*,' moments are: ',l_sum(:) + mu(:)
 !print*,' residuals are: ',l_sum(:)

end subroutine residual

end subroutine reconstruct_moments_maxent

!
! compute the integrand (bit inside the integral) of the kth moment
!
!  IN:
!    x   : grid on which to evaluate the integrand
!    lamb: array of Lagrange multipliers
!    k   : order of moment to calculate (e.g. 0, 1, 2...)
!
function integrand(x,lamb,k,logx) result(y)
 real( kind = rk ), intent(in) :: x(:)
 real( kind = rk ), intent(in) :: lamb(:)
 real( kind = rk ) :: y(size(x))  ! result
 real( kind = rk ), intent(in), optional :: logx(:)
 integer, intent(in) :: k
 real( kind = rk ) :: xi(size(lamb))
 integer :: i,j
 real, parameter :: e100 = exp(-100.)
 real( kind = rk ) :: term

 do i=1,size(x)
    if (use_log_grid) then
       if (present(logx)) then
          do j=1,size(lamb)
             xi(j) = logx(i)**(j-1)
          enddo
       else
          do j=1,size(lamb)
             xi(j) = log(x(i))**(j-1)
          enddo
       endif
    else
       do j=1,size(lamb)
          xi(j) = x(i)**(j-1)
       enddo
    endif
    ! function is exp( lam_0 + lam_1 x + lam_2 x^2 + ...)
    ! or if log_grid is set exp( lam_0 + lam_1 log(x) + lam_2 log(x)^2 + ...)
    term = dot_product(lamb,xi)
    if (term < -100.) then ! avoid wasted effort on exp(-100) etc
       y(i) = x(i)**(k/3.)*e100
    else
       y(i) = x(i)**(k/3.)*exp(term)
    endif

 enddo

end function integrand

end module reconstruct_maxent
