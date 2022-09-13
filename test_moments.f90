program test_moments
 use reconstruct_from_moments, only:reconstruct_maxent,&
                                    integrand,fsolve_error
 use integrate,                only:integrate_trap
 integer, parameter :: nx = 100
 integer :: i,k
 real :: xmin,xmax,dx
 real :: x(nx),f(nx),ftmp(nx)
 real :: mu(4),mu0,mu1,mu2,mu3,lambsol(4),lambguess(4)

 xmin=0.
 xmax=2.
 dx = (xmax-xmin)/nx
 do i=1,nx
    x(i) = xmin + (i-0.5)*dx
 enddo

 m = 5.
 c = 3.145
 k = 0

 k=0; mu0=integrate_trap(y,xmin,xmax)
 k=1; mu1=integrate_trap(y,xmin,xmax)
 k=2; mu2=integrate_trap(y,xmin,xmax)
 k=3; mu3=integrate_trap(y,xmin,xmax)
 k=0
 print*,'moments are: ',mu0,mu1,mu2,mu3

 mu=[mu0,mu1,mu2,mu3]

 lambguess = [-7.97995483,-36.31574354,151.04895024,-118.41603063]
 !print*,' using lambguess = ',lambguess

 call reconstruct_maxent(mu,x,f,lambsol,ierr)
 if (ierr /= 1) print*,' INFO: '//fsolve_error(ierr)

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
 ! print exact
 open(unit=2,file='exact.out',status='replace',action='write')
 do i=1,nx
    write(2,*) x(i),y(x(i))
 enddo
 close(2)

contains
!
!  internal function used for testing
!
pure real function y(x)
 real, intent(in) :: x

 !return ((x+0.1)**(-m) + c)*x**k
 !return (x**(-m) + c)*x**k
 !return (m*x**2 + c)*x**k
 y = exp(-0.5*(x**2-0.5)**2/(0.1**2))*x**k

end function y

end program test_moments
