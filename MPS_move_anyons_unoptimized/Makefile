#I don't like literal tabs
.RECIPEPREFIX := ;

CC = g++

IDIR =include
ODIR=obj
SRCDIR=src
OUTDIR=bin

CXXFLAGS=-Wall -I$(IDIR)

LIBS = -std=c++17 -O2 -larmadillo  -lpthread

default: MoveBiggest.txt
;gnuplot plot.gpi

MoveBiggest.txt :$(OUTDIR)/mps_demonstration.exe
;$(OUTDIR)/mps_demonstration.exe 9 6 5 24 32 q9k.tsv\
 pos_k0_9.tsv\
 pos_k1_9.tsv\
 pos_k2_9.tsv\
 pos_k3_9.tsv\
 pos_k4_9.tsv\
 pos_k5_9.tsv\
 pos_k6_9.tsv\
 pos_k7_9.tsv\
 pos_k8_9.tsv\
 pos_k9_9.tsv\
 pos_k10_9.tsv\
 pos_k11_9.tsv\
 pos_k12_9.tsv\
 pos_k13_9.tsv\
 pos_k14_9.tsv\
 pos_k15_9.tsv\
 pos_k16_9.tsv\
 pos_k17_9.tsv\
 pos_k18_9.tsv\
 pos_k19_9.tsv\
 pos_k20_9.tsv\
 pos_k21_9.tsv\
 pos_k22_9.tsv\
 pos_k23_9.tsv\
 pos_k24_9.tsv\
 pos_k25_9.tsv\
 pos_k26_9.tsv\
 pos_k27_9.tsv\
 pos_k28_9.tsv\
 pos_k29_9.tsv\
 pos_k30_9.tsv\
> MoveLog.txt



$(OUTDIR)/mps_demonstration.exe: $(ODIR)/mps.o $(ODIR)/mpo.o $(ODIR)/main.o $(ODIR)/approx.o
;$(CC) -o $@ $^ $(CXXFLAGS) $(LIBS)



$(ODIR)/main.o: $(SRCDIR)/main.cpp $(IDIR)/mps.hpp $(IDIR)/mpo.hpp $(IDIR)/approx.hpp
;$(CC) -c -o $@ $< $(CXXFLAGS) $(LIBS)

$(ODIR)/mpo.o: $(SRCDIR)/mpo.cpp $(IDIR)/mpo.hpp  $(IDIR)/mps.hpp
;$(CC) -c -o $@ $< $(CXXFLAGS) $(LIBS)

$(ODIR)/mps.o: $(SRCDIR)/mps.cpp $(IDIR)/mps.hpp
;$(CC) -c -o $@ $< $(CXXFLAGS) $(LIBS)

$(ODIR)/approx.o: $(SRCDIR)/approx.cpp $(IDIR)/approx.hpp
;$(CC) -c -o $@ $< $(CXXFLAGS) $(LIBS)


clean:
;rm -f $(ODIR)/*.o
;rm -f *.pdf
;rm -f $(OUTDIR)/*.exe
