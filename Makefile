#
# Makefile for REAFFOX-AlOH version
# 
# Byoungseon Jeon
# School of Engineering and Applied Sciences
# Harvard University
# May 2011
#
.SUFFIXES: .o .f90
#
F90 = gfortran
F90 = ifort
LD  = xild
OMP = -openmp
FLAGS = -xHOST -O3 -ipo -no-prec-div #-heap-arrays
FLAG2 = ${FLAGS} ${OMP}

#FLAGS = -g  -Wall  -fbounds-check   -O3
#FLAGS =  -march=native -ffast-math -funroll-loops -O3 
#OMP = -fopenmp
#FLAGS = -march=native -ffast-math -funroll-loops -O3 -g -fbounds-check
#FLAG2 =  ${FLAGS} ${OMP}
OBJ = datafmt.o main.o force.o vverlet.o eem.o util.o  cellsort.o
TARGET = AlOH_Efield_OMP
${TARGET}:${OBJ}
	${F90} -o ${TARGET} ${OMP} ${OBJ} ${LIB} 

.f90.o:
	${F90} ${FLAGS} -c $< ${INC}
main.o:main.f90
	${F90} ${FLAG2} -c main.f90
force.o:force.f90
	${F90} ${FLAG2} -c force.f90
eem.o:eem.f90
	${F90} ${FLAG2} -c eem.f90

clean:
	rm -rf *.mod *.o *.f90~ core ${TARGET}
