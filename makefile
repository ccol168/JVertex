LIBS:=`root-config --libs`
INCS:=`root-config --cflags`

do :
	make JUNO_PMTs.o
	make JVertex.o

JUNO_PMTs.o: JUNO_PMTs.cxx JUNO_PMTs.h 
	g++ -c JUNO_PMTs.cxx -o JUNO_PMTs.o ${INCS} ${LIBS}
JVertex.o: JVertex.cxx 
	g++ -c JVertex.cxx -o JVertex.o -lSpectrum ${INCS} ${LIBS}

