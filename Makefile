CC=gcc
CFLAGS=-I. -lstdc++ -lm -O3 -std=c++17 -fopenmp 
DEPS = consts.h lin_alg.h nbd_object.h nbd_sys.h octtree.h 
OBJ = main.o lin_alg.o nbd_object.o nbd_sys.o octtree.o

%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

main: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)

clean:
	rm -f *.o
