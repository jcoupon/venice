# Makefile for venice
SHELL := /bin/bash

# Options
RM = rm -f
EXEC = bin/venice
CC = gcc
CFLAGS  = -fPIC -Wall -Wextra -O3 #-g
LDFLAGS =
PREFIX_GSL = # default: /usr/local
PREFIX_CFITSIO = # default: /usr/local
# Where python is installed
CFLAGS_PYTHON  = -I/usr/include/python2.6
LDFLAGS_PYTHON = -ldl -lpython2.6


# read locations for gsl and CFITSIO
ifneq ($(PREFIX_GSL), )
	GSL = $(PREFIX_GSL)
else
	GSL = /usr/local
endif

ifneq ($(PREFIX_CFITSIO), )
	CFITSIO = $(PREFIX_CFITSIO)
else
	CFITSIO = /usr/local
endif

# source files
SRCS    = utils.c fits.c init.c main.c
OBJS    = $(SRCS:.c=.o)

# extra headers
CFLAGS += -Iinclude -I$(CFITSIO)/include  -I$(GSL)/include
LFLAGS += -lm  -lcfitsio -lgsl -lgslcblas -L$(CFITSIO)/lib -L$(GSL)/lib
LDFLAGS +=  -Wl,-rpath,$(FFTW)/lib -Wl,-rpath,$(GSL)/lib

.PHONY: all
all: $(EXEC)

vpath %.h include
vpath %.c src

$(EXEC):  $(OBJS)
	$(CC)  -o $@ $^ $(CFLAGS) $(LFLAGS) $(LDFLAGS)

%.o:  %.c
	$(CC) -c -o $@ $< $(CFLAGS)

%.h:

.PHONY: clean
clean:
	-${RM} ${OBJS}

tar_test:
	tar czvf test_venice.tgz test
	mv test_venice.tgz $(HOME)/gdrive/public

# TODO: the code below is probalably broken
python:
	swig -python main.i
	$(CC) -c $(SRC) main_wrap.c $(CFLAGS) $(CFLAGS_PYTHON)
	$(CC) -bundle $(SRC:.c=.o) main_wrap.o -o _$(EXEC)py.so $(LDFLAGS) $(LDFLAGS_PYTHON)
