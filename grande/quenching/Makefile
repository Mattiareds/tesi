PROG = md_prova

SRCS = module_parameters.f90 module_cluster.f90\
       md_main.f90 leggi_input.f90 leggi_pot_lj.f90 leggi_pot_morse.f90 leggi_pot_gupta.f90\
       init.f90 bigvoi.f90 voisin.f90 \
	force_lj.f90 force_morse.f90 force_gupta.f90 therma.f90 vel.f90 \
	gauss.f90 write_avg.f90 write_photo.f90 ranmar_bh.f \
	growth.f90
	

OBJS = module_parameters.o module_cluster.o\
       md_main.o leggi_input.o leggi_pot_lj.o leggi_pot_morse.o leggi_pot_gupta.o\
       init.o bigvoi.o voisin.o \
	force_lj.o force_morse.o force_gupta.o therma.o vel.o \
	gauss.o write_avg.o write_photo.o ranmar_bh.o \
	growth.o
	

LIBS =  

#FC= gfortran  -g -O0 -fbacktrace
FC= gfortran -O3 
FFLAGS= 
#F90 = gfortran  -g -O0 -fbacktrace
F90 = gfortran -O3
F90FLAGS = 
LDFLAGS =  

all: $(PROG)

$(PROG): $(OBJS)
	$(F90) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	rm -f $(PROG) $(OBJS) *.mod

.SUFFIXES: $(SUFFIXES) .f90 

.f90.o:
	$(F90) $(F90FLAGS) -c  $<
