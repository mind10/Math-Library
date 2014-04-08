LIB = ../../lib
SRC = ../../src
CC = g++

CFLAGS = -O2
INC = -I../../include -I../../vender/linux/include
LIBS = 

RMCP_LIB = librmcp.a

SRC_FILES = $(SRC)/dcspline.cpp \
			$(SRC)/dgqmatrix.cpp \
			$(SRC)/dinertia.cpp \
			$(SRC)/djacobian.cpp \
			$(SRC)/dmatrix.cpp \
			$(SRC)/donlp2.cpp \
			$(SRC)/donlp2_sub.cpp \
			$(SRC)/donmatrix.cpp \
			$(SRC)/dormatrix.cpp \
			$(SRC)/dpolynomial.cpp \
			$(SRC)/dquaternion.cpp \
			$(SRC)/drobotics.cpp \
			$(SRC)/dse3.cpp \
			$(SRC)/dso3.cpp \
			$(SRC)/dsqoptimize.cpp \
			$(SRC)/dssmatrix.cpp \
			$(SRC)/duquaternion.cpp \
			$(SRC)/dvector.cpp \
			$(SRC)/newx.cpp \
			$(SRC)/qsort_wkey.cpp \
			$(SRC)/user_eval.cpp

OBJ_FILES = dcspline.o \
			dgqmatrix.o \
			dinertia.o \
			djacobian.o \
			dmatrix.o \
			donlp2.o \
			donlp2_sub.o \
			donmatrix.o \
			dormatrix.o \
			dpolynomial.o \
			dquaternion.o \
			drobotics.o \
			dse3.o \
			dso3.o \
			dsqoptimize.o \
			dssmatrix.o \
			duquaternion.o \
			dvector.o \
			newx.o \
			qsort_wkey.o \
			user_eval.o


all : lib

lib : $(SRC_FILES)
	$(CC) -c $(CFLAGS) $(INC) $(SRC_FILES)
	$(AR) rcv $(RMCP_LIB) $(OBJ_FILES)
	ranlib $(RMCP_LIB)
	mv $(RMCP_LIB) $(LIB)

clean :
	rm -f $(OBJ_FILES)

distclean :
	rm -f $(OBJ_FILES) $(LIB)/$(RMCP_LIB)

