
#------------------------------------------------#
#Makefile for venice. Edit the following lines 
#to change the current configuration:

CC	= icc
CFLAGS	= -I$(HOME)/local/include #-Wall -Wuninitialized -O3 -fPIC
LDFLAGS	= -L$(HOME)/local/lib -lgsl -lgslcblas  -lm -lcfitsio
#LDFLAGS	= -L$(HOME)/local/lib -lgsl -lgslcblas  -lm 	
EXEC	= venice
SRC	= main.c
#SRC	= main_nofits.c
OBJ	= $(SRC:.c=.o)

all: $(EXEC)

$(EXEC) : $(OBJ)
	$(CC) -o $@ $^ $(LDFLAGS)

main.o: main.h

%.o: %.c
	$(CC) -o $@ -c $< $(CFLAGS)

.PHONY: clean mrproper

clean:
	rm -rf *.o $(EXEC).tar

mrproper: clean
	rm -rf $(EXEC)

tar: 
	tar cvf $(EXEC).tar Makefile *.c *.h README
