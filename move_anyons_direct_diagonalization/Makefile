#I don't like literal tabs
.RECIPEPREFIX := ;

CC = g++

IDIR =include
ODIR=obj
SRCDIR=src
OUTDIR=bin

CXXFLAGS=-Wall -I$(IDIR)

LIBS = -std=c++17 -O2 -larmadillo  -lpthread

default:  q.pdf


q.pdf : q.tsv #q_4.tsv
;gnuplot plot.gpi
;notify-send 'Computation finished'

#Calculate some static potentials at different positions, these two are the two individual steps we need to rotate clockwise or counterclockwise
q.tsv: $(OUTDIR)/move_potential.exe pos0.tsv pos1.tsv pos2.tsv pos3.tsv pos4.tsv pos5.tsv pos6.tsv pos7.tsv pos8.tsv pos9.tsv pos10.tsv pos11.tsv pos12.tsv Makefile
;$(OUTDIR)/move_potential.exe 6 6 q.tsv q_all.tsv 4 pos0.tsv pos1.tsv pos2.tsv pos3.tsv pos4.tsv pos5.tsv pos6.tsv pos7.tsv pos8.tsv pos9.tsv pos10.tsv pos11.tsv pos12.tsv 1> log.txt 2> F.tsv



q_crazy.tsv: $(OUTDIR)/move_potential.exe pos1.tsv pos11.tsv  Makefile
;$(OUTDIR)/move_potential.exe 6 6 q_crazy.tsv q_all_crazy.tsv 4 pos1.tsv pos7.tsv pos1.tsv 1> log_crazy.txt 2> F_crazy.tsv

q_4.tsv: $(OUTDIR)/move_potential.exe pos0_4.tsv pos1_4.tsv pos2_4.tsv pos3_4.tsv pos4_4.tsv pos5_4.tsv pos6_4.tsv pos7_4.tsv pos8_4.tsv pos9_4.tsv pos10_4.tsv pos11_4.tsv pos12_4.tsv Makefile
;$(OUTDIR)/move_potential.exe 6 6 q_4.tsv  q_all_4.tsv 5 pos0_4.tsv pos1_4.tsv pos2_4.tsv pos3_4.tsv pos4_4.tsv pos5_4.tsv pos6_4.tsv pos7_4.tsv pos8_4.tsv pos9_4.tsv pos10_4.tsv pos11_4.tsv pos12_4.tsv 1> log.txt 2> F_4.tsv


$(OUTDIR)/move_potential.exe: $(ODIR)/move_potential.o $(ODIR)/get_state.o $(ODIR)/generate_states.o $(ODIR)/approx.o $(ODIR)/minimize.o
;$(CC) -o $@ $^ $(CXXFLAGS) $(LIBS)


$(ODIR)/move_potential.o: $(SRCDIR)/move_potential.cpp #$(IDIR)/generate_states.hpp $(IDIR)/approx.hpp  $(IDIR)/minimize.hpp
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
