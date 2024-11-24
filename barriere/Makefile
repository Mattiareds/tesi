PROG = md_prova

SRCS = module_parameters.f90 module_cluster.f90 module_structures.f90\
       md_main.f90 leggi_input.f90 \
       init.f90 bigvoi.f90 voisin.f90 \
	force_rgl.f90 therma.f90 vel.f90 \
	gauss.f90 write_avg.f90 write_photo.f90 ranmar_bh.f \
	growth.f90

OBJS = module_parameters.o module_cluster.o module_structures.o\
       md_main.o leggi_input.o \
       init.o bigvoi.o voisin.o \
	force_rgl.o therma.o vel.o \
	gauss.o write_avg.o write_photo.o ranmar_bh.o \
	growth.o

LIBS =  

FC= gfortran -O3
FFLAGS= -std=legacy
F90 = gfortran -O3
F90FLAGS = -std=legacy  #serve per trattare tutte le variabili 
#		 reali in doppia precisione 
LDFLAGS =  

all: $(PROG)

$(PROG): $(OBJS)
	$(F90) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	rm -f $(PROG) $(OBJS) *.mod

.SUFFIXES: $(SUFFIXES) .f90 

.f90.o:
	$(F90) $(F90FLAGS) -c  $<
