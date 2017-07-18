program diffv

!calculates the dipole autocorrelation function
!then fourier transforms it to get the Infrared spectra
!
! to compile:
! gfortran -O3 Dauto.f90 four1.f -o dauto.x
! Input is a file with all dipoles at all time-steps
!
! Option to limit number of points outputted / smooth spectrum , added by D. Elton 3/27/14
!

implicit none

integer, parameter :: dp=kind(1.0d0)

real(dp), dimension(:,:,:), allocatable :: Da ! saves ell the dipoles for all the times
integer :: trun, tcor,na,i,j,k,t0,tt0,tt0max,t,nmol,istep,ia,l,ix
integer :: n,m,count,ja, tread, kk,nsoft, ierr
real(dp) :: timestep, omega,IR,magn, normal, avgMag, MaxFreqOut, length, vol, PointsAvailable
character*180 :: fileinp,fileout
character*2, dimension(:), allocatable :: sym
complex, dimension(:), allocatable :: aux1,aux2,vcross,ACF
real(dp), parameter :: Cspeed=3.00d10 ! cm/s
real(dp), parameter :: Kb=1.38d-23 !J/K
real(dp), parameter :: Temp=300 !K
real(dp), parameter :: hbar=1.0545718*10d-34 ! J*s
real(dp), parameter :: fs2s=1d-15 ! 1fs in s
real(dp), parameter :: Debye2SI = 3.33564d-30
real(dp), parameter :: eps0 = 8.854187d-12 !SI units
real(dp) :: pi, A0,A1,A2,A3, window
integer :: NumAvgOver, NumPointsOut


NumPointsOut = 250
open(unit=1, file="inputV.dat", form="formatted", status="unknown")

read(1,*) fileinp  ! name of the file with all the dipoles
read(1,*) fileout  ! name of output file
read(1,*) timestep ! times step in fs
read(1,*) nmol     ! number of molecules, if 1 then assumes that it uses the total dipole M
read(1,*) length   ! length of box size (assumed a cubic box) (Ang)
read(1,*) MaxFreqOut !maximum frequency to plot (cm^-1)

close(1)

open(unit=1, file=fileinp, form="formatted", status="unknown")

tread=0
ierr =0
!!count number of timesteps
do while (ierr .eq. 0)
	do j=1,nmol
		read(1,*,IOSTAT=ierr)
	enddo
	tread = tread  + 1
	! if (mod(i,100) .eq. 1) write(*,*) i
enddo
write(*,*) "File contains ", tread, " points"

!find closest power of 2 less than tread
trun = 2**(    floor( dlog(  dble(tread) )/dlog(2d0)  )  + 1	 )

rewind(1)

allocate(Da(3,nmol,0:2*trun-1))
allocate(aux1(0:2*trun-1)) !Auxiliary array, for fourier transform
allocate(aux2(0:2*trun-1)) ! second Aux arrays
allocate(Vcross(0:2*trun-1)) ! Stores the cross correlation in reciprocal space
allocate(ACF(0:2*trun-1))  ! stores  the autocorrelation Fucntion (real space)


!read in data
do i = 1, tread-1
  do j=1,nmol
    read(1,*)  (Da(ix,j,i), ix=1,3) ! format of the file: x, y z of each dipole in 3 columns
  enddo
  if (mod(i,100) .eq. 1) write(*,*) i
enddo


tread = trun

close(1)





ACF=0.0d0
count=0
normal = 0.0
Do i=1,nmol
write(*,*) "doing mol.", i
  do j=i,nmol
   do l=1,3
	! Call to the direct FFT
         aux1=cmplx(Da(l,i,:))
         aux2=cmplx(Da(l,j,:))
         n=size(aux1)
          call four1(aux1,n,-1)
          call four1(aux2,n,-1)
	  	  !Build the cross correlation in fourier space
          Vcross=aux1*conjg(aux2)
		  ! Inverse transform CrossV
          call four1(Vcross,n,1)
          Vcross=Vcross/real(size(Vcross))
          !Store the autocorrelation function
          ACF=ACF+Vcross
    enddo
  count=count+1
  enddo
enddo

!!!Save the correlation function to file
!open(unit=10, file="corr_function.dat", form="formatted", status="unknown")
!do i = 1, tread
!	write(10,*) i*timestep, real(ACF(i))/real(ACF(1))
!enddo
! close(10)

! Fourier transform the ACF
aux1=cmplx(ACF)
call four1(aux1,n,-1)
n=size(aux1)!necessary
aux1=aux1/real(n)


!--------------------- Save the IR spectrum in the file -----------------
open(unit=11, file=fileout, form="formatted", status="unknown")
write(*,*) "n = ",  n
write(*,*) "MaxFreqOut = ",  MaxFreqOut


PointsAvailable = MaxFreqOut/( 1d0/(timestep*fs2s*n*Cspeed) )
write(*,*) "points available = ", PointsAvailable

!if (PointsAvailable .lt. NumPointsOut) NumPointsOut = PointsAvailable

NumAvgOver =PointsAvailable/NumPointsOut
if (NumAvgOver == 0) NumAvgOver = 1
write(*,*) "Averaging over", NumAvgOver

vol = (length**3)*1d-30

do t = 0, NumPointsOut-1

  avgMag = 0
  do i = 1, NumAvgOver
	omega=( t*NumAvgOver+i )/(timestep*n*fs2s) !get freq in 1/s (Hz

	avgMag = avgMag + (omega**2)*real(aux1(t*NumAvgOver+i))

! 	!avgMag = avgMag + omega*tanh(hbar*omega/(Kb*2.0d0*Temp))*real(aux1(t*NumAvgOver+i))
	!avgMag = avgMag + omega*(1d0 - dexp(-(hbar*omega/(Kb*Temp))))*real(aux1(t*NumAvgOver+i))
  enddo
  avgMag = avgMag/real(NumAvgOver)


  IR = (2d0*(3.14159d0)*(Debye2SI**2)*avgMag)/(3d0*vol*2.99d8*kb*Temp*eps0)

  ! Use the prefactor with harmonic quantum correction (See Ramirez paper)
  !IR = (2d0*(3.14159d0)*(Debye2SI**2)*avgMag)/(3d0*vol*2.99d8*hbar)

  !IR = 3.14159d0*(Debye2SI**2)*avgMag/(3d0*vol*hbar*2.99d8)

  IR = IR*fs2s ! convert units of dt from Fourier transform into seconds

  IR = IR*.01 !convert 1/m to 1/cm

  IR = IR*2   !fudge factor -D. Elton )

  omega=( floor((t+.5)*numAvgOver) )/(timestep*n*fs2s) !get central freq in 1/s (Hz
  omega = omega/Cspeed  	! convert frequency to cm-1

  WRITE(11,*)  omega, IR
enddo

close(11)
end
