#------------------------------------------------#
#Makefile for venice                             #	
#------------------------------------------------#

CC	= icc
FITS    = yes

ifeq ($(FITS),yes)
   SRC	   = main.c fits.c
   LDFLAGS = -L$(HOME)/local/lib -lgsl -lgslcblas  -lm -lcfitsio
else
   SRC  = main.c withoutFits.c
   LDFLAGS = -L$(HOME)/local/lib -lgsl -lgslcblas  -lm
endif

CFLAGS	= -I$(HOME)/local/include #-Wall -Wuninitialized -O3 -fPIC
EXEC	= venice
OBJ	= $(SRC:.c=.o)

all: $(EXEC)

$(EXEC) : $(OBJ)
	$(CC) -o $@ $^ $(LDFLAGS)

%.o: %.c
	$(CC) -o $@ -c $< $(CFLAGS)

.PHONY: clean mrproper

clean:
	rm -rf *.o $(EXEC).tar

mrproper: clean
	rm -rf $(EXEC)

tar: 
	tar cvf $(EXEC).tar Makefile *.c *.h README
