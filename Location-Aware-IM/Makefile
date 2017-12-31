CC= g++
CFLAGS= -O3 -Wall
# CFLAGS= -Wall -O0 -fno-inline -g
LDFLAGS= -lm
TAGS=mip
target_mip= Quadtree.o heap.o Cache.o limit.o general_cascade.o graph.o pmia_upgrade.o mip.o main.o 

all:$(TAGS)

mip: $(target_mip)
		$(CC) $(CFLAGS) $(target_mip) -o mip $(LDFLAGS)
main.o: main.cpp
		$(CC) $(CFLAGS) -o main.o -c main.cpp
mip.o: mip.cpp
		$(CC) $(CFLAGS) -o mip.o -c mip.cpp
pmia_upgrade.o: pmia_upgrade.cpp
		$(CC) $(CFLAGS) -o pmia_upgrade.o -c pmia_upgrade.cpp
graph.o: graph.cpp        
		$(CC) $(CFLAGS) -o graph.o -c graph.cpp
general_cascade.o: general_cascade.cpp        
		$(CC) $(CFLAGS) -o general_cascade.o -c general_cascade.cpp
limit.o: limit.cpp        
		$(CC) $(CFLAGS) -o limit.o -c limit.cpp
Cache.o: Cache.cpp        
		$(CC) $(CFLAGS) -o Cache.o -c Cache.cpp
heap.o: heap.cpp        
		$(CC) $(CFLAGS) -o heap.o -c heap.cpp
Quadtree.o: Quadtree.cpp        
		$(CC) $(CFLAGS) -o Quadtree.o -c Quadtree.cpp
clean:
	rm -rf *.o $(TAGS)
