CC= g++ -Wall -O3

#INCLUDE= -I /home/yerui/src/SeqLib/ -I /home/yerui/src/SeqLib/htslib/ /home/yerui/src/SeqLib/bin/libseqlib.a /home/yerui/src/SeqLib/bin/libhts.a /home/yerui/src/SeqLib/bin/libbwa.a /home/yerui/src/SeqLib/bin/libfml.a
#LIBS= -llzma -lbz2 -L. -lz

INCLUDE= -I /home/yerui/miniconda2/include -I /home/yerui/src/alglib-3.14.0/src
LIBS= -L /home/yerui/miniconda2/lib -L.
OBJS= main.o likelihood.o

lhmut: lhmut.cc $(OBJS)
	$(CC) $< $(OBJS) $(INCLUDE) $(LIBS) -o $@ /home/yerui/src/alglib-3.14.0/src/*.o

main.o: main.cc main.h
	$(CC) -c $< -o $@

likelihood.o: likelihood.cc likelihood.h main.h
	$(CC) -c $< $(INCLUDE) $(INCLUDE) -o $@

.PHONY: clean

clean:
	-rm lhmut $(OBJS)
