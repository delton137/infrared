!------------------------------------------------------------------------
!--- compute infrared spectrum and write out smoothed spectra ----------
!------------------------------------------------------------------------
! format of the dipole moment file : x, y z of each dipole in 3 columns
!
! 2013-2016 Daniel C. Elton
!------------------------------------------------------------------------
program infrared
 use math
 implicit none
 Integer  :: NumPointsOut = 350 !maximum number of points to use
 real, dimension(:,:), allocatable :: dip_moms
 character(len=400)     :: filein, fileout
 integer :: trun, tcor, ix, i, iH, t, n,  tread,  Nhyd, PointsAvailable, temp
 double precision :: omega, IR, magn, avgMag, MinFreqOut, vol, timestep, length, MaxFreqOut
 double precision, dimension(:), allocatable :: ACF, output, DFT, allfreqs
 double precision, dimension(:), allocatable :: spectra_smoothed, allfreqs_smoothed
 complex, dimension(:), allocatable :: aux1
 integer :: NumAvgOver, ierr, j, nmol, lun=30, BlockSize
 double precision, parameter :: Cspeed=3.00d10 ! cm/s
 double precision, parameter :: hbarSI=1.0545718*10d-34  ! J*s
 double precision, parameter :: fs2s=1d-15 ! 1fs in s
 double precision, parameter :: Debye2SI = 3.33564d-30
 double precision, parameter :: kb = 0.0019872041d0 !Boltzmann constant in kcal/mole/K

 open(unit=1, file="infrared.inp", form="formatted", status="unknown")
 read(1,*) filein  ! name of the file with all the dipoles
 read(1,*) fileout  ! name of output file
 read(1,*) timestep ! times step in fs
 read(1,*) nmol     ! number of molecules, if 1 then assumes that it uses the total dipole M
 read(1,*) length   ! length of box size (assumed a cubic box) (Ang)
 read(1,*) MaxFreqOut !maximum frequency to plot (cm^-1)
 close(1)

 open(unit=1, file=trim(filein), form="formatted", status="unknown")

 !-------------- count number of timesteps ---------------------
 ierr=0
 tread=0
 do while (ierr .eq. 0)
	read(1,*,IOSTAT=ierr)
	tread = tread  + 1
 enddo
 write(*,*) "File contains ", tread, " points"

 !----- find closest power of 2 less than tread
 trun = 2**(    floor( dlog(  dble(tread) )/dlog(2d0)  )  + 1	 )

 rewind(1)

 allocate(dip_moms(3,tread))
 vol = (length**3)*1d-30

 !------------ read in data
 do i = 1, tread-1
	read(1,*)  (dip_moms(ix, i), ix=1,3) ! format of the file: x, y z of each dipole in 3 columns
	if (mod(i,1) .eq. 1) write(*,*) i
 enddo

 allocate(output(tread))! Stores the cross correlation in reciprocal space
 allocate(ACF(tread))   ! stores the autocorrelation Function (real space)
 allocate(allfreqs(tread))
 allocate(DFT(tread))   ! stores the autocorrelation Function (real space)

 ACF=0
 !find correlation function
 do ix = 1, 3
	call calc_corr_function(dip_moms(ix, 1:tread), output)
	ACF = ACF + output
 enddo
 ACF = ACF/ACF(1)

 !! Save the correlation function to file--------------------
  open(40, file="out_"//trim(fileout)//"_vel_vel_corr_function.dat")
  do i = 1, tread
   	write(40,*) i*timestep, real(ACF(i))
  enddo
  close(40)
 !-----------------------------------------------------------

 call calc_DFT(ACF, DFT, allfreqs, timestep, size(ACF))

 allfreqs = allfreqs/fs2s !convert to Hz

 !apply quantum harmonic correction
 do i = 1, PointsAvailable
	DFT(i) = allfreqs(i)*tanh(hbarSI*allfreqs(i)/(Kb*2.0d0*Temp))*DFT(i)
 enddo

 !Use the prefactor with harmonic quantum correction (See Ramirez paper)
 DFT = (2d0*3.14159d0*(Debye2SI**2)*DFT)/(3d0*vol*2.99d8*hbarSI*Cspeed)

 allfreqs = allfreqs/Cspeed !convert to cm^-1

 if (size(allfreqs) .gt. 1) MinFreqOut =  allfreqs(2)

 PointsAvailable = floor(MaxFreqOut/MinFreqOut)

 if (PointsAvailable .lt. NumPointsOut) NumPointsOut = PointsAvailable

 allocate(allfreqs_smoothed(NumPointsOut))
 allocate(spectra_smoothed(NumPointsOut))

 write(*,*) size(allfreqs), NumPointsOut

 BlockSize = floor(real(size(allfreqs))/N)

 do i = 1, N
	allfreqs_smoothed(i) = sum(allfreqs((i-1)*BlockSize+1:i*BlockSize) ) /BlockSize
 enddo

  BlockSize = floor(real(size(DFT))/N)

 do i = 1, N
	spectra_smoothed(i) = sum(DFT((i-1)*BlockSize+1:i*BlockSize) ) /BlockSize
 enddo


 !! Save the DOS to file--------------------
 open(lun, file="out_"//trim(fileout)//"_IR.dat")
 do i = 1, size(allfreqs) !NumPointsOut
 	write(lun,'(2g12.4)')  allfreqs(i), DFT(i)
!	write(lun,'(2g12.4)')  allfreqs_smoothed(i), spectra_smoothed(i)
 enddo
 close(lun)
 !!----------------------------------------

end program infrared
