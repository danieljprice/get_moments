module reconstruct_from_moments
 implicit none
 real, parameter :: pi = 4.*atan(1.)
 integer, parameter :: rk = kind(0.d0)

 public :: reconstruct,reconstruct_maxent,fsolve_error,integrand

 logical, private :: use_log_grid = .false.
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

 ! in principle can choose method here
 if (present(lambguess)) then
    call reconstruct_maxent(mu,x,f,lambsol,ierr,lambguess)
 else
    call reconstruct_maxent(mu,x,f,lambsol,ierr,lambguess)
 endif

end subroutine reconstruct

!
! maximum entropy reconstruction of function given moments
!
subroutine reconstruct_maxent(mu,x,f,lambsol,ierr,lambguess,use_log)
 real, intent(in) :: mu(:)
 real, intent(in) :: x(:)
 real, intent(out) :: f(size(x))
 real, intent(out) :: lambsol(size(mu))
 integer, intent(out) :: ierr
 real, intent(in), optional :: lambguess(size(mu))
 logical, intent(in), optional :: use_log
 integer :: n_moments
 real, parameter :: tol = 1.e-6
 real :: lsum(size(mu))

 lambsol = 0.
 ierr = 0
 f = 0.
 n_moments = size(mu)

 ! initial guesses for Lagrange multipliers
 if (present(lambguess)) then
    ! use guesses passed as arguments (e.g. from previous attempt)
    lambsol = lambguess
 else
    ! use predefined guess
    lambsol(1) = -log(sqrt(2.*pi))
 endif
 use_log_grid = .false.
 if (present(use_log)) use_log_grid = use_log

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
 use integrate, only:integrate_trap,integrate_trap_log
 integer, intent(in) :: n
 real ( kind = rk ), intent(in)  :: lamb(n) ! guess for  solution
 real ( kind = rk ), intent(out) :: l_sum(n) ! function evaluated for given lamb

 integer :: k

 if (use_log_grid) then
    do k=1,n
       l_sum(k) = integrate_trap_log(size(x),x,integrand(x,lamb,k-1)) - mu(k)
    enddo
 else
    do k=1,n
       l_sum(k) = integrate_trap(size(x),x,integrand(x,lamb,k-1)) - mu(k)
    enddo
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
function integrand(x,lamb,k) result(y)
 real( kind = rk ), intent(in) :: x(:)
 real( kind = rk ), intent(in) :: lamb(:)
 real( kind = rk ) :: y(size(x))  ! result
 integer, intent(in) :: k
 real( kind = rk ) :: xi(size(lamb))
 integer :: i,j

 do i=1,size(x)
    if (use_log_grid) then
       do j=1,size(lamb)
          xi(j) = log(x(i))**(j-1)
       enddo
    else
       do j=1,size(lamb)
          xi(j) = x(i)**(j-1)
       enddo
    endif
    ! function is exp( lam_0 + lam_1 x + lam_2 x^2 + ...)
    ! or if log_grid is set exp( lam_0 + lam_1 log(x) + lam_2 log(x)^2 + ...)
    y(i) = x(i)**k * exp(dot_product(lamb,xi))
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

end module reconstruct_from_moments
