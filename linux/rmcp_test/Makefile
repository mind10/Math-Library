LIB = ../../lib
SRC = .
CC = g++

CFLAGS = -O2
INC = -I../../include -I../../vender/linux/include
LIB = -L../../lib -L../../vender/linux/lib

SRC_FILES = $(SRC)/dcspline_test.cpp \
			$(SRC)/dgqmatrix_test.cpp \
			$(SRC)/dmatrix_test.cpp \
			$(SRC)/dquaternion_test.cpp \
			$(SRC)/duquaternion_test.cpp \
			$(SRC)/drobotics_test.cpp \
			$(SRC)/dse3_test.cpp \
			$(SRC)/dso3_test.cpp \
			$(SRC)/dsqoptimize_test.cpp \
			$(SRC)/dvector_test.cpp \
			$(SRC)/main.cpp \
			$(SRC)/utility.cpp

LIB_FILES = -lm -lrmcp -lmkl_ia32 -lmkl_lapack -lguide -lmkl -lmkl_solver

all : lib test

lib :
	cd ../rmcp_lib; $(MAKE)

test : $(SRC_FILES)
	$(CC) -o test $(CFLAGS) $(INC) $(LIB) $(SRC_FILES) $(LIB_FILES)

clean :
	rm -f test *.o

