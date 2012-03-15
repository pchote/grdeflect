CFLAGS = -g -Wall
CC = gcc
FF = gfortran
OSX_LIBS = -I/sw/lib/pgplot -L/sw/lib/pgplot -L/usr/X11R6/lib -lX11 -L/sw/lib -laquaterm -Wl,-framework -Wl,Foundation  -lpng
OSX_INCL = -I/sw/lib/pgplot
LIBS =  -lm -lcpgplot -lpgplot -lcurses

grdeflect: lightDeflection.o
	${FF} ${CFLAGS} -o $@ lightDeflection.o  ${OSX_LIBS} ${LIBS}

.c.o:
	${CC} ${CFLAGS} ${OSX_INCL} -c $<

install:
	cp grdeflect /usr/bin/
	chmod -R 755 /usr/bin/grdeflect

uninstall:
	rm /usr/bin/grdeflect
