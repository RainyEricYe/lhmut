CC= g++ -Wall -O3

INCLUDE= -I /home/yerui/miniconda2/include -I /home/yerui/src/alglib-3.14.0/src
LIBS= -L /home/yerui/miniconda2/lib -L.
OBJS= main.o likelihood.o option.o

lhmut: lhmut.cc option.h $(OBJS)
	$(CC) $< -o $@ $(OBJS) /home/yerui/src/alglib-3.14.0/src/*.o $(INCLUDE) $(LIBS)

main.o: main.cc main.h option.h
	$(CC) -c $< -o $@

likelihood.o: likelihood.cc likelihood.h main.h option.h
	$(CC) -c $< $(INCLUDE) -o $@

option.o: option.cc option.h
	$(CC) -c $< -o $@

.PHONY: clean

clean:
	-rm lhmut $(OBJS)
