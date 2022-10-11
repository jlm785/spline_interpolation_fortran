program test_grid_interp

! example of test program to check the spline interpolation.
! using gnuplot with
! gnuplot> pl 'fort.100', 'fort.102'
! will allow to check the interpolation

! Change to any other function, grid, etc...

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

  integer        ::  idum
  real(REAL64)   ::  small
  real(REAL64), external   ::  ran2

  integer        ::  nin
  real(REAL64)   ::  xin(101), fin(101)
  integer        ::  i
  integer        ::  ng
  real(REAL64)   ::  xgmin, xgmax, xdif
  real(REAL64)   ::  xg(1000)

  real(REAL64)   ::  a, b
  integer        ::  l
  real(REAL64)   ::  alfa, rmax

  real(REAL64)   ::  fexact

  integer        ::  isx, ierr
  real(REAL64)   ::  yp(101)
  real(REAL64)   ::  ypp(101)
  real(REAL64)   ::  w(101,3)

  real(REAL64)   ::  fgs(1000)
  real(REAL64)   ::  ypg(1000)
  real(REAL64)   ::  yppg(1000)

  real(REAL64), parameter    ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64


  l = 0

  rmax = 20.0_REAL64
!  alfa = 0.25
  alfa = 0.04

  nin = 101
  idum = 12345
  small = 0.001


  b = 0.1
  a = rmax / (exp(b*(nin-1))-UM)
  do i = 1,nin
    xin(i) = a*(exp(b*(i-1))-UM)
    fin(i) = xin(i)**l * exp(-alfa*xin(i)*xin(i) )
!    fin(i) = xin(i)**l * exp(-alfa*xin(i) )
  enddo
  do i = 1,nin
    write(100,'(2f16.6)') xin(i),fin(i)
  enddo
!  do i = 1,nin
!    xin(i) = 0.5*(exp(0.1*i)-1.0d0)
!    fin(i) = sin(xin(i)) + small*ran2(idum)
!    write(100,'(2f16.6)') xin(i),fin(i)
!  enddo

  xgmin = 0.0
  xgmax = 20.0

  ng = 501
  xdif = (xgmax - xgmin)/(ng-1)
  do i = 2,ng
    xg(i) = xgmin + (i-1)*xdif
  enddo

! Spline interpolation

  isx = 0
  call splift (xin,fin,yp,ypp,nin,w,ierr,isx,ZERO,ZERO,ZERO,ZERO)

  write(6,*) '  splift  ierr = ',ierr

  call splint (xin,fin,ypp,nin,xg,fgs,ypg,yppg,ng,ierr)

  write(6,*) '  splint  ierr = ',ierr

  do i = 1,ng
    fexact = xg(i)**l * exp(-alfa*xg(i)*xg(i) / rmax*rmax)
    write(102,'(2f16.6,f20.10)') xg(i),fgs(i),fgs(i) - fexact
  enddo

  stop

end program test_grid_interp

