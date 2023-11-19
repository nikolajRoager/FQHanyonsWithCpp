#I don't like literal tabs
.RECIPEPREFIX := ;

CC = g++

IDIR =include
ODIR=obj
SRCDIR=src
OUTDIR=bin

CXXFLAGS=-Wall -I$(IDIR)

LIBS = -std=c++17 -O2 -larmadillo  -lpthread

default:  log.txt logBig.txt logBigger.txt logBiggest.txt

log.txt :$(OUTDIR)/mps_demonstration.exe
;$(OUTDIR)/mps_demonstration.exe 6 6 3 2 16 pos0.tsv q02.tsv raw664.tsv > log.txt
;$(OUTDIR)/mps_demonstration.exe 6 6 3 3 16 pos0.tsv q03.tsv raw664.tsv >> log.txt
;$(OUTDIR)/mps_demonstration.exe 6 6 3 4 16 pos0.tsv q04.tsv raw664.tsv >> log.txt
;$(OUTDIR)/mps_demonstration.exe 6 6 3 5 16 pos0.tsv q05.tsv raw664.tsv >> log.txt
;$(OUTDIR)/mps_demonstration.exe 6 6 3 6 16 pos0.tsv q06.tsv raw664.tsv >> log.txt
;$(OUTDIR)/mps_demonstration.exe 6 6 3 7 16 pos0.tsv q07.tsv raw664.tsv >> log.txt
;$(OUTDIR)/mps_demonstration.exe 6 6 3 8 16 pos0.tsv q08.tsv raw664.tsv >> log.txt
;$(OUTDIR)/mps_demonstration.exe 6 6 3 12 16 pos0.tsv q012.tsv raw664.tsv >> log.txt
;$(OUTDIR)/mps_demonstration.exe 6 6 3 16 16 pos0.tsv q016.tsv raw664.tsv >> log.txt


thread_speed_test/logSpeed.txt :$(OUTDIR)/mps_demonstration.exe
;$(OUTDIR)/mps_demonstration.exe 6 6 5 8  32 pos0.tsv thread_speed_test/q6_08.tsv raw666.tsv  > thread_speed_test/logSpeed.txt
;$(OUTDIR)/mps_demonstration.exe 6 6 5 16 32 pos0.tsv thread_speed_test/q6_016.tsv raw666.tsv >>  thread_speed_test/logSpeed.txt
;$(OUTDIR)/mps_demonstration.exe 6 6 5 24 32 pos0.tsv thread_speed_test/q6_024.tsv raw666.tsv >>  thread_speed_test/logSpeed.txt


logBig.txt :$(OUTDIR)/mps_demonstration.exe
;$(OUTDIR)/mps_demonstration.exe 6 6 5 2  32 pos0.tsv q6_02.tsv raw666.tsv  > logBig.txt
;$(OUTDIR)/mps_demonstration.exe 6 6 5 3  32 pos0.tsv q6_03.tsv raw666.tsv  >> logBig.txt
;$(OUTDIR)/mps_demonstration.exe 6 6 5 4  32 pos0.tsv q6_04.tsv raw666.tsv  >> logBig.txt
;$(OUTDIR)/mps_demonstration.exe 6 6 5 5  32 pos0.tsv q6_05.tsv raw666.tsv  >> logBig.txt
;$(OUTDIR)/mps_demonstration.exe 6 6 5 6  32 pos0.tsv q6_06.tsv raw666.tsv  >> logBig.txt
;$(OUTDIR)/mps_demonstration.exe 6 6 5 7  32 pos0.tsv q6_07.tsv raw666.tsv  >> logBig.txt
;$(OUTDIR)/mps_demonstration.exe 6 6 5 8  32 pos0.tsv q6_08.tsv raw666.tsv  >> logBig.txt
;$(OUTDIR)/mps_demonstration.exe 6 6 5 10 32 pos0.tsv q6_010.tsv raw666.tsv >> logBig.txt
;$(OUTDIR)/mps_demonstration.exe 6 6 5 12 32 pos0.tsv q6_012.tsv raw666.tsv >> logBig.txt
;$(OUTDIR)/mps_demonstration.exe 6 6 5 14 32 pos0.tsv q6_014.tsv raw666.tsv >> logBig.txt
;$(OUTDIR)/mps_demonstration.exe 6 6 5 16 32 pos0.tsv q6_016.tsv raw666.tsv >> logBig.txt
;$(OUTDIR)/mps_demonstration.exe 6 6 5 24 32 pos0.tsv q6_024.tsv raw666.tsv >> logBig.txt


logBigger.txt :$(OUTDIR)/mps_demonstration.exe
;$(OUTDIR)/mps_demonstration.exe 8 6 5 4  32 pos0_8.tsv q8_04.tsv  raw686.tsv > logBigger.txt
;$(OUTDIR)/mps_demonstration.exe 8 6 5 6  32 pos0_8.tsv q8_06.tsv  raw686.tsv >> logBigger.txt
;$(OUTDIR)/mps_demonstration.exe 8 6 5 8  32 pos0_8.tsv q8_08.tsv  raw686.tsv >> logBigger.txt
;$(OUTDIR)/mps_demonstration.exe 8 6 5 10 32 pos0_8.tsv q8_010.tsv raw686.tsv >> logBigger.txt
;$(OUTDIR)/mps_demonstration.exe 8 6 5 12 32 pos0_8.tsv q8_012.tsv raw686.tsv >> logBigger.txt
;$(OUTDIR)/mps_demonstration.exe 8 6 5 14 32 pos0_8.tsv q8_014.tsv raw686.tsv >> logBigger.txt
;$(OUTDIR)/mps_demonstration.exe 8 6 5 16 32 pos0_8.tsv q8_016.tsv raw686.tsv >> logBigger.txt
;$(OUTDIR)/mps_demonstration.exe 8 6 5 20 32 pos0_8.tsv q8_020.tsv raw686.tsv >> logBigger.txt
;$(OUTDIR)/mps_demonstration.exe 8 6 5 24 32 pos0_8.tsv q8_024.tsv raw686.tsv >> logBigger.txt



logBiggest.txt :$(OUTDIR)/mps_demonstration.exe
;$(OUTDIR)/mps_demonstration.exe 9 6 5 8  32 pos0_9.tsv q9_08.tsv > logBiggest.txt
;$(OUTDIR)/mps_demonstration.exe 9 6 5 16 32 pos0_9.tsv q9_016.tsv >> logBiggest.txt
;$(OUTDIR)/mps_demonstration.exe 9 6 5 24 32 pos0_9.tsv q9_024.tsv >> logBiggest.txt



log24.txt :$(OUTDIR)/mps_demonstration.exe
;$(OUTDIR)/mps_demonstration.exe 9 6 2 24  32 pos0_9.tsv q0693.tsv raw693.tsv >> log24.txt
;$(OUTDIR)/mps_demonstration.exe 9 6 3 24  32 pos0_9.tsv q0694.tsv raw694.tsv >> log24.txt
;$(OUTDIR)/mps_demonstration.exe 9 6 4 24  32 pos0_9.tsv q0695.tsv raw695.tsv >> log24.txt



NEWlogBiggest.txt :$(OUTDIR)/mps_demonstration.exe
;$(OUTDIR)/mps_demonstration.exe 9 6 5 12 32 pos0_9.tsv q9_012.tsv > NEWlogBiggest.txt
;$(OUTDIR)/mps_demonstration.exe 9 6 5 18 32 pos0_9.tsv q9_018.tsv >> NEWlogBiggest.txt
;$(OUTDIR)/mps_demonstration.exe 9 6 5 20 32 pos0_9.tsv q9_020.tsv >> NEWlogBiggest.txt
;$(OUTDIR)/mps_demonstration.exe 9 6 6 16 32 pos0_9.tsv q9_6_016.tsv >  NEWlogBiggerest.txt
;$(OUTDIR)/mps_demonstration.exe 9 6 7 16 32 pos0_9.tsv q9_7_016.tsv > NEWlogBIGgerest.txt
;$(OUTDIR)/mps_demonstration.exe 9 6 6 24 32 pos0_9.tsv q9_6_024.tsv >>  NEWlogBiggerest.txt
;$(OUTDIR)/mps_demonstration.exe 9 6 7 24 32 pos0_9.tsv q9_7_024.tsv >> NEWlogBIGgerest.txt
;$(OUTDIR)/mps_demonstration.exe 9 6 5 28 32 pos0_9.tsv q9_028.tsv >> NEWlogBiggest.txt


MoveBiggest.txt :$(OUTDIR)/mps_demonstration.exe
;$(OUTDIR)/mps_demonstration.exe 9 6 5 24 32 pos_k0_9.tsv q9k0.tsv > MoveLog.txt
;$(OUTDIR)/mps_demonstration.exe 9 6 5 24 32 pos_k1_9.tsv q9k1.tsv >> MoveLog.txt
;$(OUTDIR)/mps_demonstration.exe 9 6 5 24 32 pos_k2_9.tsv q9k2.tsv >> MoveLog.txt
;$(OUTDIR)/mps_demonstration.exe 9 6 5 24 32 pos_k3_9.tsv q9k3.tsv >> MoveLog.txt
;$(OUTDIR)/mps_demonstration.exe 9 6 5 24 32 pos_k5_9.tsv q9k5.tsv >> MoveLog.txt



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
;rm -f $(OUTDIR)/*.exe
