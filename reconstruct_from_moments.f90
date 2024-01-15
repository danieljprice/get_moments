module reconstruct_from_moments
 implicit none
 real, parameter :: pi = 4.*atan(1.)
 integer, parameter :: rk = kind(0.d0)

 public :: reconstruct,reconstruct_maxent,fsolve_error,integrand
 public :: reconstruct_gamma_dist

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
    call reconstruct_maxent(mu,x,f,lambsol,lsum,ierr,lambguess)
 else
    call reconstruct_maxent(mu,x,f,lambsol,lsum,ierr)
 endif

end subroutine reconstruct

!
! maximum entropy reconstruction of function given moments
!
subroutine reconstruct_maxent(mu,x,f,lambsol,lsum,ierr,lambguess,use_log,weights)
 use integrate, only:check_log_grid_is_uniform
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
    !print*,' lsum analytically = ',l_sum(:)
    do k=1,n
       l_sum(k) = integrate_gauss_legendre(size(x),weights,integrand(x,lamb,k-1,logx=logx)) - mu(k)
    enddo
    !print*,' lsum numerically = ',l_sum(:)
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

end subroutine reconstruct_maxent

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
 real( kind = rk ) :: xi(size(lamb)), term
 integer :: i,j
 real, parameter :: e100 = exp(-100.)
 real :: prefac,theta,d_on_p,p,beta,d

 beta = lamb(1)   ! overall normalisation factor
 theta = lamb(2)  ! scale parameter theta on wikipedia
 d_on_p = lamb(3)      ! shape parameter k on wikipedia
 p = lamb(4)      ! shape parameter
 prefac = beta * p / (theta * Gamma(d_on_p))
 d = d_on_p * p
 !prefac = beta / (theta * Gamma(d/p) * (1./(p*theta))**(1./p))
 !print*, 1./Gamma(d/p) !prefac,beta,a,d,p
 !read*
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
    !if (.not.present(logx)) stop 'logx is not present'

    ! GENERALIZED GAMMA DISTRIBUTION (standard Gamma distribution if p=1)
    term = -(x(i)/theta)**p
    if (term < -100.) then
       y(i) = x(i)**(k/3.)*e100*prefac !prefac*(x(i)/theta)**(d-1)*e100
    else
       y(i) = x(i)**(k/3.)*prefac*(x(i)/theta)**(d-1)*exp(term)
    endif
    !x(i)**(lamb(2)-1)*exp(-lamb(3)*x(i))*prefac
    !y(i) = x(i)**(k/3.) * exp(max(dot_product(lamb,xi),-100.))
 enddo

end function integrand

!
! print the error message corresponding to the error code
!
function fsolve_error(ierr) result(string)
 integer, intent(in) :: ierr
 character(len=62) :: string

 select case(ierr)
 case(0)
     string = 'improper input parameters'
 case(1)
     string = 'relative error between X and solution is at most tol'
 case(2)
     string = 'number of calls to residual function exceeded 200*(N+1)'
 case(3)
     string = 'TOL is too small. No further improvement in solution possible'
 case(4)
     string = 'the iteration is not making good progress'
 case default
    string = ''
 end select

end function fsolve_error

!
! maximum entropy reconstruction of function given moments
!
subroutine reconstruct_gamma_dist(mu,lambsol,lsum,ierr,lambguess)
 real, intent(in) :: mu(:)
 real, intent(out) :: lambsol(size(mu))
 real, intent(out) :: lsum(size(mu))
 integer, intent(out) :: ierr
 real, intent(in), optional :: lambguess(size(mu))
 integer :: n_moments,k
 real, parameter :: tol = 1.e-6
 real :: d_on_p,theta,p

 lambsol = 0.
 ierr = 0
 n_moments = 2 !size(mu)

 ! initial guesses for Lagrange multipliers
 if (present(lambguess)) then
    lambsol = lambguess
 else
    lambsol(1) = 1.7
 endif

 call fsolve(residual_fit_gamma,n_moments,lambsol,lsum,tol,ierr)

 d_on_p = lambsol(1)
 p = lambsol(2)
 theta = (mu(2)/mu(1) * Gamma(d_on_p) / Gamma(d_on_p + 1./(3.*p)))**3
 !print*,'INPUT  moments = ',mu
 !print*,'SOLVED moments = ',(mu(1)*theta**(k/3.)*Gamma(d_on_p + k/(3.*p))/Gamma(d_on_p),k=0,3)

! call residual_fit_gamma(n_moments,lambsol,lsum)

contains

subroutine residual_fit_gamma(n,lamb,l_sum)
 integer, intent(in) :: n
 real ( kind = rk ), intent(in)  :: lamb(n) ! guess for  solution
 real ( kind = rk ), intent(out) :: l_sum(n) ! function evaluated for given lamb
 real ( kind = rk ) :: scalefac, p, theta, d_on_p
 integer :: k

 d_on_p = lamb(1)
 p = lamb(2)
 theta = (mu(2)/mu(1) * Gamma(d_on_p) / Gamma(d_on_p + 1./(3.*p)))**3

 ! solve for moments 3 and 4 (k=2,3)
 l_sum(1) = theta**(2./3.) * Gamma(d_on_p + 2./(3.*p)) / Gamma(d_on_p) * mu(1)/mu(3) - 1.! /mu(1)
 l_sum(2) = theta**(3./3.) * Gamma(d_on_p + 3./(3.*p)) / Gamma(d_on_p) * mu(1)/mu(4) - 1.! /mu(1)

 !print*,'d/p,p = ',d_on_p,p,' want moments = ',mu(3:4), ' got ',&
 !       (theta**(k/3.) * Gamma(d_on_p + k/(3.*p)) / Gamma(d_on_p) * mu(1),k=2,3),' err = ',l_sum(1:2)
 !do k=1,4
 !   l_sum(k) = scalefac * theta ** (1.+(k-1)/3.) / lamb(2) * Gamma(d_on_p + (k-1)/3.) - mu(k)
! enddo

end subroutine residual_fit_gamma

end subroutine reconstruct_gamma_dist


end module reconstruct_from_moments
