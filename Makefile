
CFLAGS	        = 
FFLAGS	        =
CPPFLAGS        =
FPPFLAGS        =
CLEANFILES = PETSCtest.o

PETSC_DIR = /home/peggyhuang/Software/petsc-3.5.3
PETSC_ARCH=arch-linux2-c-debug

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

all: 	main.o nrutil.o boun_func.o init_cond.o init_setting.o matrices_update.o user_param.o chkopts
	mpicc -o NonlinSolid main.o nrutil.o boun_func.o init_cond.o init_setting.o matrices_update.o user_param.o ${PETSC_LIB}
	${RM} main.o nrutil.o boun_func.o init_cond.o init_setting.o matrices_update.o user_param.o



