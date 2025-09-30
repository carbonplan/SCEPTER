# Start of the makefile
# Defining variables

FC            = gfortran
# FC            = ifort

CPFLAGS       = 
CPFLAGS       += -Dno_intr_findloc # need to use in cluster
# CPFLAGS       += -Dshow_PSDiter # showing iteration process during PSD calculation
# CPFLAGS       += -Dparallel_ON # testing parallelization
# CPFLAGS       += -Dnpar_in=8 # number of threads for parallelization 
# CPFLAGS       += -Dnpar_in=1 # number of threads for parallelization 
CPFLAGS       += -Dparpsd_chk # checking parallelization results
# CPFLAGS       += -Dksld_chk # checking rate consts for sld species
CPFLAGS       += -DolddustPSD # using old PSD for dust (not user input but prescribed one)
# CPFLAGS       += -Derrmtx_printout # 
CPFLAGS       += -Dmod_basalt_cmp # using basalt composition defined in <basalt_define.h>
# CPFLAGS       += -Ddef_flx_save_alltime # flux reported each integration (costs lots of bites)
# CPFLAGS       += -Dfull_flux_report # output all cumulative flux
# CPFLAGS       += -Ddisp_lim # limiting the display of results
# CPFLAGS       += -Ddiss_only # not allowing precipitation of minerals
# CPFLAGS       += -Dlim_minsld # limiting mineral lowest conc. 
# CPFLAGS       += -Dporoiter # do iteration for porosity  
# CPFLAGS       += -Dcalcw_full # fully coupled w calcuation  
# CPFLAGS       += -Ddispiter # showing PSD flux in each iteration   
# CPFLAGS       += -DdispPSDiter # showing PSD flux in each iteration   
# CPFLAGS       += -Dcalcporo_full # fully coupled poro calcuation  
# CPFLAGS       += -Diwtypein=0 # uplift type 0--cnst w, 1-- cnst poro*w, 2-- cnst (1-poro)*w, 3--- w-flexible, if not defined 0 is taken

ifeq ($(FC),gfortran)
  # CFLAGS        = -fcheck=all -g -O3  
  # CFLAGS        = -Wall -O3 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace
  # CFLAGS        = -Wall -O3 -g -fcheck=all -fbacktrace
  CFLAGS        = -fimplicit-none  -Wall  -Wline-truncation  -Wcharacter-truncation  -Wsurprising  \
	  -Waliasing  -Wimplicit-interface  -Wunused-parameter  -fwhole-file  -fcheck=all  -std=gnu  -pedantic  -fbacktrace -O3

endif

ifeq ($(FC),ifort)
  CFLAGS        = -O3 -heap-arrays -g -traceback -check bounds -fp-stack-check -gen-interfaces -warn interfaces -check arg_temp_created 
endif 

# LDFLAGS       = -L/usr/local/lib
LDFLAGS       = 

LIBS          = -lopenblas

ifneq (,$(findstring -Dmod_basalt_cmp,$(CPFLAGS)))
  # Found -Dmod_basalt_cmp
  INC          = -I/home/tykukla/SCEPTER/data 
else
  # Not found
  INC          = 
endif

ifneq (,$(findstring -Dparallel_ON,$(CPFLAGS)))
  # Found -Dparallel_ON
  CFLAGS        += -fopenmp
else
  # Not found
endif

OBJS          = scepter.o 
SRC           = scepter.f90  
                            
PROGRAM       = scepter

# - TK added to handle alternative scepter.f90 files -----------------------
OBJS_A      = scepter_rateA.o
SRC_A       = scepter_rateA.f90
PROGRAM_A   = scepter_rateA

rateA:           $(PROGRAM_A)

$(PROGRAM_A): $(OBJS_A)
	$(FC) $(OBJS_A) -o $(PROGRAM_A) -cpp $(CPFLAGS) $(CFLAGS) $(LIBS) $(LDFLAGS) $(INC)

$(OBJS_A):   $(SRC_A)
	$(FC) $(SRC_A) -c -cpp $(CPFLAGS) $(CFLAGS) $(LIBS) $(LDFLAGS) $(INC)

# ---
OBJS_B      = scepter_rateB.o
SRC_B       = scepter_rateB.f90
PROGRAM_B   = scepter_rateB

rateB:           $(PROGRAM_B)

$(PROGRAM_B): $(OBJS_B)
	$(FC) $(OBJS_B) -o $(PROGRAM_B) -cpp $(CPFLAGS) $(CFLAGS) $(LIBS) $(LDFLAGS) $(INC)

$(OBJS_B):   $(SRC_B)
	$(FC) $(SRC_B) -c -cpp $(CPFLAGS) $(CFLAGS) $(LIBS) $(LDFLAGS) $(INC)

# ---
OBJS_C      = scepter_rateC.o
SRC_C       = scepter_rateC.f90
PROGRAM_C   = scepter_rateC

rateC:           $(PROGRAM_C)

$(PROGRAM_C): $(OBJS_C)
	$(FC) $(OBJS_C) -o $(PROGRAM_C) -cpp $(CPFLAGS) $(CFLAGS) $(LIBS) $(LDFLAGS) $(INC)

$(OBJS_C):   $(SRC_C)
	$(FC) $(SRC_C) -c -cpp $(CPFLAGS) $(CFLAGS) $(LIBS) $(LDFLAGS) $(INC)

# ---
OBJS_D      = scepter_rateD.o
SRC_D       = scepter_rateD.f90
PROGRAM_D   = scepter_rateD

rateD:           $(PROGRAM_D)

$(PROGRAM_D): $(OBJS_D)
	$(FC) $(OBJS_D) -o $(PROGRAM_D) -cpp $(CPFLAGS) $(CFLAGS) $(LIBS) $(LDFLAGS) $(INC)

$(OBJS_D):   $(SRC_D)
	$(FC) $(SRC_D) -c -cpp $(CPFLAGS) $(CFLAGS) $(LIBS) $(LDFLAGS) $(INC)
# --------------------------------------------------------------------------

all:            $(PROGRAM)

$(PROGRAM):     $(OBJS)
	$(FC) $(OBJS) -o $(PROGRAM) -cpp $(CPFLAGS) $(CFLAGS) $(LIBS) $(LDFLAGS) $(INC)

$(OBJS):        $(SRC) 
	$(FC) $(SRC) -c -cpp $(CPFLAGS) $(CFLAGS) $(LIBS) $(LDFLAGS) $(INC)

clean:;         rm -f *.o  *~ $(PROGRAM) $(PROGRAM_A) $(PROGRAM_B) $(PROGRAM_C) $(PROGRAM_D)
blank:;         truncate -s 0 *.out
cleanall:;         rm -f *.o *.out *~ $(PROGRAM) $(PROGRAM_A) $(PROGRAM_B) $(PROGRAM_C) $(PROGRAM_D)
