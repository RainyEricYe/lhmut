CC= g++ -Wall -O3

INCLUDE= -I /home/yerui/miniconda2/include -I /home/yerui/src/alglib-3.14.0/src
LIBS= -L /home/yerui/miniconda2/lib -L.
OBJS= main.o likelihood.o

lhmut: lhmut.cc option.hpp $(OBJS)
	$(CC) $< $(OBJS) $(INCLUDE) $(LIBS) -o $@ /home/yerui/src/alglib-3.14.0/src/*.o

main.o: main.cc main.h option.hpp
	$(CC) -c $< -o $@

likelihood.o: likelihood.cc likelihood.h main.h option.hpp
	$(CC) -c $< $(INCLUDE) $(INCLUDE) -o $@

.PHONY: clean

clean:
	-rm lhmut $(OBJS)
