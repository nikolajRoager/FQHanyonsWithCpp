#I don't like literal tabs
.RECIPEPREFIX := ;

CC = g++

IDIR =include
ODIR=obj
SRCDIR=src
OUTDIR=bin

CXXFLAGS=-Wall -I$(IDIR)

LIBS = -std=c++17 -O2 -larmadillo  -lpthread

default:  log.txt


#Calculate some static potentials at different positions, these two are the two individual steps we need to rotate clockwise or counterclockwise

log.txt : $(OUTDIR)/static_potential.exe
;/usr/bin/time -f "%P %M" $(OUTDIR)/static_potential.exe 6 6 3 n0663.tsv raw663.tsv > log.txt
;/usr/bin/time -f "%P %M" $(OUTDIR)/static_potential.exe 6 6 4 n0664.tsv raw664.tsv >> log.txt
;/usr/bin/time -f "%P %M" $(OUTDIR)/static_potential.exe 6 6 5 n0665.tsv raw665.tsv >> log.txt
;/usr/bin/time -f "%P %M" $(OUTDIR)/static_potential.exe 6 6 6 n0666.tsv raw666.tsv >> log.txt
;/usr/bin/time -f "%P %M" $(OUTDIR)/static_potential.exe 7 6 3 n0673.tsv raw673.tsv >> log.txt
;/usr/bin/time -f "%P %M" $(OUTDIR)/static_potential.exe 7 6 4 n0674.tsv raw674.tsv >> log.txt
;/usr/bin/time -f "%P %M" $(OUTDIR)/static_potential.exe 7 6 5 n0675.tsv raw675.tsv >> log.txt
;/usr/bin/time -f "%P %M" $(OUTDIR)/static_potential.exe 7 6 6 n0676.tsv raw676.tsv >> log.txt
;/usr/bin/time -f "%P %M" $(OUTDIR)/static_potential.exe 8 6 3 n0683.tsv raw683.tsv >> log.txt
;/usr/bin/time -f "%P %M" $(OUTDIR)/static_potential.exe 8 6 4 n0684.tsv raw684.tsv >> log.txt
;/usr/bin/time -f "%P %M" $(OUTDIR)/static_potential.exe 8 6 5 n0685.tsv raw685.tsv >> log.txt
;/usr/bin/time -f "%P %M" $(OUTDIR)/static_potential.exe 8 6 6 n0686.tsv raw686.tsv >> log.txt
;/usr/bin/time -f "%P %M" $(OUTDIR)/static_potential.exe 9 6 3 n0693.tsv raw693.tsv >> log.txt
;/usr/bin/time -f "%P %M" $(OUTDIR)/static_potential.exe 9 6 4 n0694.tsv raw694.tsv >> log.txt
;/usr/bin/time -f "%P %M" $(OUTDIR)/static_potential.exe 9 6 5 n0695.tsv raw695.tsv >> log.txt
#;/usr/bin/time -f "%P %M" $(OUTDIR)/static_potential.exe 9 6 6 n0696.tsv >> log.txt


$(OUTDIR)/static_potential.exe: $(ODIR)/static_potential.o $(ODIR)/get_state.o $(ODIR)/generate_states.o $(ODIR)/approx.o $(ODIR)/minimize.o
;$(CC) -o $@ $^ $(CXXFLAGS) $(LIBS)


$(ODIR)/static_potential.o: $(SRCDIR)/static_potential.cpp #$(IDIR)/generate_states.hpp $(IDIR)/approx.hpp  $(IDIR)/minimize.hpp
;$(CC) -c -o $@ $< $(CXXFLAGS) $(LIBS)

$(ODIR)/generate_states.o: $(SRCDIR)/generate_states.cpp $(IDIR)/generate_states.hpp  $(IDIR)/approx.hpp
;$(CC) -c -o $@ $< $(CXXFLAGS) $(LIBS)

$(ODIR)/get_state.o: $(SRCDIR)/get_state.cpp $(IDIR)/get_state.hpp  $(IDIR)/approx.hpp
;$(CC) -c -o $@ $< $(CXXFLAGS) $(LIBS)

$(ODIR)/approx.o: $(SRCDIR)/approx.cpp $(IDIR)/approx.hpp
;$(CC) -c -o $@ $< $(CXXFLAGS) $(LIBS)

$(ODIR)/minimize.o: $(SRCDIR)/minimize.cpp $(IDIR)/minimize.hpp
;$(CC) -c -o $@ $< $(CXXFLAGS) $(LIBS)

clean:
;rm -f $(ODIR)/*.o
;rm -f $(OUTDIR)/*.exe
;rm -f log*.txt
;rm -f q*.tsv
;rm -f  *.pdf
