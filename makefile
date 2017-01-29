OBJECTST = General.o Sequence.o DataSet.o Model.o CPTParameterisation.o Train.o
HEADERST = General.h Sequence.h DataSet.h Model.h CPTParameterisation.h

CC	= g++
CFLAGS  = -g

Train :	$(OBJECTST) 
	$(CC) $(CFLAGS) -o $@ $(OBJECTST) $(LLIBS)

Train.o : Train.cpp $(HEADERST)
	$(CC) -c $(CFLAGS) $*.cpp

General.o : General.cpp General.h
	$(CC) -c $(CFLAGS) $*.cpp

Sequence.o : Sequence.cpp Sequence.h General.o
	$(CC) -c $(CFLAGS) $*.cpp

DataSet.o : DataSet.cpp DataSet.h Sequence.h General.h Sequence.o
	$(CC) -c $(CFLAGS) $*.cpp

CPTParameterisation.o : CPTParameterisation.cpp CPTParameterisation.h General.o
	$(CC) -c $(CFLAGS) $*.cpp

Model.o : Model.cpp Model.h Sequence.o CPTParameterisation.o
	$(CC) -c $(CFLAGS) $*.cpp

clean:
	rm *.o # Train
