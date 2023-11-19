#Gnu makefile for gnu+Linux system

CC=g++

IDIR =include
ODIR=obj
SRCDIR=src
OUTDIR=bin
OUTNAME=InsertName

CXXFLAGS=-std=c++17 -Wall -O2 -Wextra -Wpedantic -Wdouble-promotion -I$(IDIR)

LIBS = -lSDL2 -lSDL2_image -lSDL2_ttf -lSDL2_mixer -lGLEW -lGL -ljsoncpp



#Yes, there are better ways of doing makefiles but this is a rather small project
_DEPS = SFXwrap.hpp IO.hpp load_shader.hpp get_assets.hpp
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))


_OBJ = main.o SFXwrap.o load_shader.o get_assets.o IO_graphics.o IO_audio.o IO_input_devices.o IO.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

default: $(OUTDIR)/$(OUTNAME)


#Create object files (deprecated, that is stupidly slow)
#$(ODIR)/%.o: $(SRCDIR)/%.cpp $(DEPS)
#	$(CC) -c -o $@ $< $(CXXFLAGS)

#The main file depends on EVERYTHING
$(ODIR)/main.o: $(SRCDIR)/main.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CXXFLAGS)


$(ODIR)/IO_graphics.o: $(SRCDIR)/IO_graphics.cpp  $(IDIR)/IO.hpp $(IDIR)/load_shader.hpp
	$(CC) -c -o $@ $< $(CXXFLAGS)

$(ODIR)/IO_audio.o: $(SRCDIR)/IO_audio.cpp  $(IDIR)/IO.hpp $(IDIR)/SFXwrap.hpp
	$(CC) -c -o $@ $< $(CXXFLAGS)

$(ODIR)/IO_input_devices.o: $(SRCDIR)/IO_input_devices.cpp  $(IDIR)/IO.hpp
	$(CC) -c -o $@ $< $(CXXFLAGS)

$(ODIR)/SFXwrap.o: $(SRCDIR)/SFXwrap.cpp $(IDIR)/SFXwrap.hpp
	$(CC) -c -o $@ $< $(CXXFLAGS)

$(ODIR)/load_shader.o: $(SRCDIR)/load_shader.cpp $(IDIR)/load_shader.hpp $(IDIR)/IO.hpp
	$(CC) -c -o $@ $< $(CXXFLAGS)

$(ODIR)/get_assets.o: $(SRCDIR)/get_assets.cpp $(IDIR)/get_assets.hpp
	$(CC) -c -o $@ $< $(CXXFLAGS)

$(ODIR)/IO.o: $(SRCDIR)/IO.cpp  $(IDIR)/IO.hpp
	$(CC) -c -o $@ $< $(CXXFLAGS)

#Compile the final program
$(OUTDIR)/$(OUTNAME):	$(OBJ)
	$(CC) -o $@ $^ $(CXXFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o
	rm -f $(OUTDIR)/$(OUTNAME)
