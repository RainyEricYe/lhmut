CC= g++ -Wall

#INCLUDE= -I /home/yerui/src/SeqLib/ -I /home/yerui/src/SeqLib/htslib/ /home/yerui/src/SeqLib/bin/libseqlib.a /home/yerui/src/SeqLib/bin/libhts.a /home/yerui/src/SeqLib/bin/libbwa.a /home/yerui/src/SeqLib/bin/libfml.a
#LIBS= -llzma -lbz2 -L. -lz

OBJS= main.o compass_search.o likelihood.o

lhmut: lhmut.cc $(OBJS)
	$(CC) $< $(OBJS) -o $@

main.o: main.cc main.h
	$(CC) -c $< -o $@

compass_search.o:

likelihood.o: likelihood.cc likelihood.h compass_search.* main.*
	$(CC) -c $< -o $@

.PHONY: clean

clean:
	-rm lhmut $(OBJS)
