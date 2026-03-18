LIBS:=`root-config --libs`
INCS:=`root-config --cflags`

do :
	make Apply_JVertex.x

Apply_JVertex.x : JUNO_PMTs.o JVertex.o Apply_JVertex.o
	g++ -o Apply_JVertex.x JUNO_PMTs.o JVertex.o Apply_JVertex.o ${INCS} ${LIBS}

Apply_JVertex.o: Apply_JVertex.cxx 
	g++ -c Apply_JVertex.cxx -o Apply_JVertex.o ${INCS} ${LIBS}
JUNO_PMTs.o: JUNO_PMTs.cxx JUNO_PMTs.h 
	g++ -c JUNO_PMTs.cxx -o JUNO_PMTs.o ${INCS} ${LIBS}
JVertex.o: JVertex.cxx JVertex.h
	g++ -c JVertex.cxx -o JVertex.o -lSpectrum ${INCS} ${LIBS}

