# Makefile for venice

# Where cfitsio library is installed
CFITSIO =  # /usr/cfitsio

# Where GSL library is installed
GSL =  # /usr/local

CC = gcc

# Where python is installed
CFLAGS_PYTHON  = -I/usr/include/python2.6
LDFLAGS_PYTHON = -ldl -lpython2.6

CFLAGS      = -Iinclude # -use-asm # for old icc versions
LDFLAGS     =  -lgsl -lgslcblas -lm -lcfitsio
RM          = rm -f
EXEC        = bin/venice
# SRC         = main.c
# OBJ         = $(SRC:.c=.o)

# source files
SRCS    = utils.c fits.c init.c main.c
OBJS    = $(SRCS:.c=.o)

# Headers for libraries
ifneq ($(CFITSIO), )
	CFLAGS     +=  -I$(CFITSIO)/include
	LDFLAGS    +=  -L$(CFITSIO)/lib
endif

ifneq ($(GSL), )
	CFLAGS     +=  -I$(GSL)/include
	LDFLAGS    +=  -L$(GSL)/lib
endif

.PHONY: all
all: $(EXEC)

vpath %.h include
vpath %.c src

$(EXEC):  $(OBJS)
	$(CC)  -o $@ $^ $(CFLAGS) $(LDFLAGS)

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
