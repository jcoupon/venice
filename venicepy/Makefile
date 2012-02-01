
#------------------------------------------------#
#Makefile for venice. Edit the following lines 
#to change the current configuration:

CC	= gcc
GSL	= /usr

CFLAGS	= -I$(GSL)/local/include -Wall -Wuninitialized -O3 -fPIC
LDFLAGS	= -L$(GSL)/local/lib -lgsl -lgslcblas  -lm
EXEC	= venice
SRC	= main.c
OBJ	= $(SRC:.c=.o)

all: $(EXEC)

$(EXEC) : $(OBJ)
	$(CC) -o $@ $^ $(LDFLAGS)

#main.o: main.h
#	$(CC) $(CFLAGS) -c main.c

%.o: %.c
	$(CC) -o $@ -c $< $(CFLAGS)

.PHONY: clean mrproper

clean:
	rm -rf *.o $(EXEC).tar

mrproper: clean
	rm -rf $(EXEC)

tar: 
	tar cvf $(EXEC).tar Makefile main.c main.h README


swig: mask.cpp mask.i main.o
	g++ $(CFLAGS) -c mask.cpp -o mask.o
	swig -python -c++ mask.i
	g++ $(CFLAGS) -c mask_wrap.cxx -I/usr/include/python2.7 
	g++ $(CFLAGS) -shared  mask_wrap.o mask.o main.o -o _venicepy.so  $(LDFLAGS)
