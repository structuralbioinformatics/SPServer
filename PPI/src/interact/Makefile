C    = cc  
OBJSI =  interact.o readpdb.o allocate.o
OBJSF =  fold.o readpdb.o allocate.o
LIBS =  -lm 
.SUFFIXES: .c .o
.c.o:
	$(C) -c $*.c 
interact: $(OBJSI) 
	$(C) -o interact $(OBJSI) $(LIBS)
fold: $(OBJSF) 
	$(C) -o fold $(OBJSF) $(LIBS)
clean: 
	\rm *.o 
