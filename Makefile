CC=g++-8


CPPFLAGS  += `root-config --cflags --libs`
LDLIBS    += -lboost_system -lboost_filesystem 
LDFLAGS   += $(shell root-config --ldflags)
OBJFILES += caenTextFile_to_Root.o rootSupport.o
TARGET += caenTextFile_to_Root


all: caenTextFile_to_Root

caenTextFile_to_Root.o: caenTextFile_to_Root.cpp
	$(CC) -Wall -c caenTextFile_to_Root.cpp -o caenTextFile_to_Root.o $(CPPFLAGS) $(LDLIBS) 

fitFunctions.o: fitFunctions.cpp
	$(CC) -Wall -c fitFunctions.cpp -o fitFunctions.o $(CPPFLAGS) $(LDLIBS)

caenTextFile_to_Root: $(OBJFILES)
	$(CC) -o caenTextFile_to_Root $(OBJFILES) $(CPPFLAGS) $(LDLIBS)

clean:
	rm ./*.o
