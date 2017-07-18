Dauto.x: Dauto.f90 four1.f 
	gfortran -O3 Dauto.f90 four1.f -o Dauto.x --debug --bounds-check --backtrace
clean:
	rm *.o

pristine: 
	rm *.x
