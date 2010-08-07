# Please enter here the locations for CloogLib include and libraries if they
# aren't the default values (/usr/lib and /usr/include).
CLOOG_INC = $(HOME)/progs/linux/include
CLOOG_LIB = $(HOME)/progs/linux/lib

CC = gcc
LDLIBS= -lcloog
CFLAGS= -I $(CLOOG_INC) -L $(CLOOG_LIB)


example: example.c
	@echo "             /*-----------------------------------------------*"
	@echo "              *                Making example                 *"
	@echo "              *-----------------------------------------------*/"
	$(CC) example.c -o example $(CFLAGS) $(LDLIBS)

clean:
	@echo "             /*-----------------------------------------------*"
	@echo "              *               Cleaning example                *"
	@echo "              *-----------------------------------------------*/"
	-rm -f example example.exe core
