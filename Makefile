
CFLAGS	        =  
FFLAGS	        =
CPPFLAGS        =
FPPFLAGS        =

PETSC_DIR = /home/peggyhuang/Software/petsc-3.6.3
PETSC_ARCH=arch-linux2-c-debug

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

all: 	main.o nrutil.o boun_func.o init_cond.o init_setting.o matrices_update.o user_param.o interface_marker.o fluid.o chkopts
	mpicc -o NonlinSolid main.o nrutil.o boun_func.o init_cond.o init_setting.o matrices_update.o user_param.o interface_marker.o fluid.o  ${PETSC_LIB}
	${RM} main.o nrutil.o boun_func.o init_cond.o init_setting.o matrices_update.o user_param.o interface_marker.o fluid.o 



