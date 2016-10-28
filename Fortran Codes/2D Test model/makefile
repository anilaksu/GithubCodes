FCOMP = gfortran
FCFLAGS = -g -O2
LDFLAGS = -lblas
PROGRAM =  2DSEM_FISModel
SRCS = 2DFEMTest.f90 MeshTransform.f90
OBJECTS = $(SRCS:.f90=.o)

all: $(PROGRAM)

$(PROGRAM): $(OBJECTS)
	$(FCOMP) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.f90
	$(FCOMP) $(FCFLAGS) -c $<

.PHONY: clean veryclean

clean:
	rm -f *.o *.mod *.MOD

veryclean: clean
	rm -f *~ $(PROGRAM)


