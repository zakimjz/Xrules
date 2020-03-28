CC = g++
CFLAGS = -O3
HDRS = xminer.h timetrack.h calcdb.h eqclass.h stats.h idlist.h hashtable.h
OBJS = calcdb.o eqclass.o stats.o idlist.o xminer.o hashtable.o
LIBS = 
TARGET = xminer

default: $(TARGET)

clean:
	rm -rf *~ *.o $(TARGET)

xminer: $(OBJS) $(HFILES)
	$(CC) $(CFLAGS) -o xminer $(OBJS) $(LIBS)
#	strip xminer

.SUFFIXES: .o .cpp

.cpp.o:
	$(CC) $(CFLAGS) -c $<


# dependencies
# $(OBJS): $(HFILES)
xminer.o: $(HDRS)
calcdb.o: calcdb.h xminer.h idlist.h
eqclass.o: eqclass.h xminer.h idlist.h calcdb.h idlist.h
stats.o: stats.h
idlist.o: idlist.h
hashtable.o: hashtable.h xminer.h eqclass.h
