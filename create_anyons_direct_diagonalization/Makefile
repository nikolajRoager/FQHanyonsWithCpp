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


q.pdf : q0.tsv q1.tsv q2.tsv q3.tsv q4.tsv q5.tsv q6.tsv q7.tsv q8.tsv q9.tsv q10.tsv q11.tsv q12.tsv

#Calculate some static potentials at different positions, these two are the two individual steps we need to rotate clockwise or counterclockwise

vacuum.tsv: $(OUTDIR)/static_potential.exe
;$(OUTDIR)/static_potential.exe 6 6 pos0.tsv vacuum.tsv &

q0.tsv: $(OUTDIR)/static_potential.exe pos0.tsv
;$(OUTDIR)/static_potential.exe 6 6 pos0.tsv V0.tsv q0.tsv    >log0.txt  &

q1.tsv: $(OUTDIR)/static_potential.exe pos1.tsv
;$(OUTDIR)/static_potential.exe 6 6 pos1.tsv V1.tsv q1.tsv    >log1.txt  &

q2.tsv: $(OUTDIR)/static_potential.exe pos2.tsv
;$(OUTDIR)/static_potential.exe 6 6 pos2.tsv V2.tsv q2.tsv    >log2.txt  &

q3.tsv: $(OUTDIR)/static_potential.exe pos3.tsv
;$(OUTDIR)/static_potential.exe 6 6 pos3.tsv V3.tsv q3.tsv    >log3.txt  &

q4.tsv: $(OUTDIR)/static_potential.exe pos4.tsv
;$(OUTDIR)/static_potential.exe 6 6 pos4.tsv V4.tsv q4.tsv    >log4.txt  &

q5.tsv: $(OUTDIR)/static_potential.exe pos5.tsv
;$(OUTDIR)/static_potential.exe 6 6 pos5.tsv V5.tsv q5.tsv    >log5.txt  &

q6.tsv: $(OUTDIR)/static_potential.exe pos6.tsv
;$(OUTDIR)/static_potential.exe 6 6 pos6.tsv V6.tsv q6.tsv    >log6.txt  &

q7.tsv: $(OUTDIR)/static_potential.exe pos7.tsv
;$(OUTDIR)/static_potential.exe 6 6 pos7.tsv V7.tsv q7.tsv    >log7.txt  &

q8.tsv: $(OUTDIR)/static_potential.exe pos8.tsv
;$(OUTDIR)/static_potential.exe 6 6 pos8.tsv V8.tsv q8.tsv   >log8.txt  &

q9.tsv: $(OUTDIR)/static_potential.exe pos9.tsv
;$(OUTDIR)/static_potential.exe 6 6 pos9.tsv V9.tsv q9.tsv   >log9.txt &

q10.tsv: $(OUTDIR)/static_potential.exe pos10.tsv
;$(OUTDIR)/static_potential.exe 6 6 pos10.tsv V10.tsv q10.tsv  >log10.txt  &

q11.tsv: $(OUTDIR)/static_potential.exe pos11.tsv
;$(OUTDIR)/static_potential.exe 6 6 pos11.tsv V11.tsv q11.tsv  >log11.txt &

q12.tsv: $(OUTDIR)/static_potential.exe pos12.tsv
;$(OUTDIR)/static_potential.exe 6 6 pos12.tsv V12.tsv q12.tsv  >log12.txt &

q13.tsv: $(OUTDIR)/static_potential.exe pos13.tsv
;$(OUTDIR)/static_potential.exe 6 6 pos13.tsv V13.tsv q13.tsv  >log13.txt &

q14.tsv: $(OUTDIR)/static_potential.exe pos14.tsv
;$(OUTDIR)/static_potential.exe 6 6 pos14.tsv V14.tsv q14.tsv >log14.txt &


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
