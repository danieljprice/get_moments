program test_moments
 use reconstruct_maxent, only:reconstruct_moments_maxent,integrand
 use integrate,          only:integrate_trap,integrate_trap_log
 use fsolve_mod,         only:fsolve_error
 implicit none
 integer, parameter :: nx = 1024
 integer :: i,k,ierr
 real :: xmin,xmax,dx,m,c
 real :: x(nx),f(nx),ftmp(nx)
 real :: mu(4),mu0,mu1,mu2,mu3,lambsol(4),lambguess(4),lsum(4)
 logical :: log_grid

 xmin=0.1
 xmax=3.
 log_grid = .true.
 if (log_grid) then
    dx = log10(xmax/xmin)/real(nx-1)
    do i=1,nx
       x(i) = log10(xmin) + (i-1)*dx
    enddo
    x = 10.**x
 else
    dx = (xmax-xmin)/nx
    do i=1,nx
       x(i) = xmin + (i-0.5)*dx
    enddo
 endif

 m = 5.
 c = 3.145
 k = 0

 if (log_grid) then
    k=0; mu0=integrate_trap_log(nx,x,y)
    k=1; mu1=integrate_trap_log(nx,x,y)
    k=2; mu2=integrate_trap_log(nx,x,y)
    k=3; mu3=integrate_trap_log(nx,x,y)
 else
    k=0; mu0=integrate_trap(y,xmin,xmax)
    k=1; mu1=integrate_trap(y,xmin,xmax)
    k=2; mu2=integrate_trap(y,xmin,xmax)
    k=3; mu3=integrate_trap(y,xmin,xmax)
 endif
 k=0
 print*,'moments are: ',mu0,mu1,mu2,mu3

 mu=[mu0,mu1,mu2,mu3]

 !lambguess = [-7.97995483,-36.31574354,151.04895024,-118.41603063]
 lambguess = [0.,-m,0.,0.]
 !print*,' using lambguess = ',lambguess

 call reconstruct_moments_maxent(mu,x,f,lambsol,lsum,ierr,use_log=log_grid,lambguess=lambguess)
 if (ierr /= 1) print*,' INFO: '//fsolve_error(ierr)

 !lambsol = lambguess
 print*,' got lambsol = ',lambsol
 f = integrand(x, lambsol, 0)

 do k=0,size(mu)-1
    ftmp = f*x**k
    print*,' moment ',k,' is',integrate_trap(nx,x,ftmp) !,sum(integrand(x,lambguess,k))
 enddo

 ! print results
 open(unit=1,file='recon.out',status='replace',action='write')
 do i=1,nx
    write(1,*) x(i),f(i)
 enddo
 close(1)

 k = 0
 ! print exact
 open(unit=2,file='exact.out',status='replace',action='write')
 do i=1,nx
    write(2,*) x(i),y(x(i))
 enddo
 close(2)

contains
!
!  internal function used for testing the moment reconstruction
!
pure real function y(x)
 real, intent(in) :: x
 real, parameter :: logxmean = log(10.**(0.2)),log_x_sigma = 0.1

 !return ((x+0.1)**(-m) + c)*x**k
 !y = (x**(-m) + 0.*c)*x**k
 !y = (m*x**2 + c)*x**k
 y = exp(-(log(x) - logxmean)**2/(2.*log_x_sigma**2))*x**k
 !y = exp(-0.5*(x-0.5)**2/(0.1**2))*x**k

end function y

end program test_moments
