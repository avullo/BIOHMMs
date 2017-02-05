OBJECTST = General.o Alphabet.o Instance.o DataSet.o Model.o CPTParameterisation.o Train.o
HEADERST = General.h Alphabet.o Instance.h DataSet.h Model.h CPTParameterisation.h

CC	= g++
CFLAGS  = -g

Train :	$(OBJECTST) 
	$(CC) $(CFLAGS) -o $@ $(OBJECTST) $(LLIBS)
	mv $@ examples/secondary_structure

Train.o : Train.cpp $(HEADERST)
	$(CC) -c $(CFLAGS) $*.cpp

General.o : General.cpp General.h
	$(CC) -c $(CFLAGS) $*.cpp

Alphabet.o : Alphabet.cpp Alphabet.h
	$(CC) -c $(CFLAGS) $*.cpp

Instance.o : Instance.cpp Instance.h General.o
	$(CC) -c $(CFLAGS) $*.cpp

DataSet.o : DataSet.cpp DataSet.h General.o Alphabet.o Instance.o
	$(CC) -c $(CFLAGS) $*.cpp

CPTParameterisation.o : CPTParameterisation.cpp CPTParameterisation.h General.o
	$(CC) -c $(CFLAGS) $*.cpp

Model.o : Model.cpp Model.h Instance.o CPTParameterisation.o
	$(CC) -c $(CFLAGS) $*.cpp

clean:
	rm *.o examples/secondary_structure/Train
