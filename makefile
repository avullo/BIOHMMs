OBJECTST = General.o Sequence.o DataSet.o # Model.o BIOHMM.o Train.o
HEADERST = General.h Sequence.h DataSet.h # Model.h BIOHMM.h Sequence.h 

CC	= g++
CFLAGS  = -O3

# Train :	$(OBJECTST) 
# 	$(CC) $(CFLAGS) -o $@ $(OBJECTST) $(LLIBS)

# Train.o : Train.cpp $(HEADERST)
# 	$(CC) -c $(CFLAGS) $*.cpp

General.o : General.cpp General.h
	$(CC) -c $(CFLAGS) $*.cpp

Sequence.o : Sequence.cpp Sequence.h General.h
	$(CC) -c $(CFLAGS) $*.cpp

DataSet.o : DataSet.cpp DataSet.h Sequence.h General.h
	$(CC) -c $(CFLAGS) $*.cpp

# Model.o : Model.cpp $(HEADERST)
# 	$(CC) -c $(CFLAGS) $*.cpp

# BIOHMM.o : BIOHMM.cpp $(HEADERST)
# 	$(CC) -c $(CFLAGS) $*.cpp


# clean:
# 	rm *.o Train
