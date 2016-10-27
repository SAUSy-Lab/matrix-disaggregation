# Makefile - Written by Wim R. Cardoen
TARGET=ODdis
CC=gcc
DEBUG=   
OBJ=info.o input.o init.o main.o loop.o dump.o
CFLAGS=
LDFLAGS= -lm
OPT= -O2 

$(TARGET):$(OBJ)
	${CC} -o $(TARGET) ${DEBUG} $(OBJ) ${LDFLAGS} 

%.o: %.c
	${CC} -c ${CFLAGS} ${DEBUG} ${OPT} $*.c

clean:
	rm -rf $(OBJ) 

veryclean:
	make clean
	rm -rf ${TARGET}
