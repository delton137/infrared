
DEBUGOPS=""

debug: DEBUGOPS=--debug --backtrace
debug: infrared.x
	
infrared.x:  four1.f   math.f90  infrared.f90
	gfortran -O3 four1.f math.f90 infrared.f90  -o infrared.x $(DEBUGOPS)
	
clean:
	rm *.o *.x
