PROG=tree64
CPPFLAGS=-Wall -Wextra -O2 -march=native -g

$(PROG): $(PROG).o

clean:
	rm *.o $(PROG)
