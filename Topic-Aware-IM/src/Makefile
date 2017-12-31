CC=g++
CFLAGS=-O3 -Wall -Wunused-result
# CFLAGS=-g3 -Wall
LDFLAGS=-lm 
TAGS=tests experiments offline
target_tests=tic.o tim.o tests.o
target_offline=tic.o tim.o offline.o
target_experiments=tic.o tim.o experiments.o

all:$(TAGS)

experiments: $(target_experiments)
		$(CC) $(CFLAGS) $(target_experiments) -o experiments $(LDFLAGS)
experiments.o: experiments.cpp
		$(CC) $(CFLAGS) -o experiments.o -c experiments.cpp

tests: $(target_tests)
		$(CC) $(CFLAGS) $(target_tests) -o tests $(LDFLAGS)
tests.o: tests.cpp
		$(CC) $(CFLAGS) -o tests.o -c tests.cpp

offline: $(target_offline)
		$(CC) $(CFLAGS) $(target_offline) -o offline $(LDFLAGS)
offline.o: offline.cpp
		$(CC) $(CFLAGS) -o offline.o -c offline.cpp

tim.o: tim.cpp
		$(CC) $(CFLAGS) -o tim.o -c tim.cpp
tic.o: tic.cpp
		$(CC) $(CFLAGS) -o tic.o -c tic.cpp
clean:
	rm -rf *.o $(TAGS)
