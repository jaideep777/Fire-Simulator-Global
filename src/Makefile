# makefile for fire model

TARGET = fire
LIBPATH = -L/usr/local/netcdf-cxx-legacy/lib -L/home/jaideep/codes/libgsm/lib
INCPATH = -I/usr/local/netcdf-c-4.3.2/include -I/usr/local/netcdf-cxx-legacy/include
INCPATH += -I/home/jaideep/codes/libgsm/include
LDFLAGS =  
CPPFLAGS = -g

CPPFILES = $(wildcard *.cpp) # globals.cpp init.cpp io.h main.cpp
OBJFILES = $(CPPFILES:.cpp=.o)

LIBS = -lgsm -lnetcdf_c++


all: $(TARGET)

$(TARGET): $(OBJFILES)
	g++ -o $(LDFLAGS) $(TARGET) $(OBJFILES) $(LIBPATH) $(LIBS)

.cpp.o:
	g++ -c $(CPPFLAGS) $(INCPATH) $< -o $@ 

clean:
	rm -f *.o $(TARGET)

re: clean all


