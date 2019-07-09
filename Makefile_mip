GUROBI_HOME = /opt/gurobi/8.1.0/linux64

# ---------------------------------------------------------------------
# Compiler selection 
# ---------------------------------------------------------------------

#GPP = /act/gcc-4.7.2/bin/g++
GPP = /usr/bin/g++

# ---------------------------------------------------------------------
# Compiler options 
# ---------------------------------------------------------------------

OPT = -m64 -O3 -fopenmp

# ---------------------------------------------------------------------
# Link options and libraries
# ---------------------------------------------------------------------

INC_GUROBI     = -I$(GUROBI_HOME)/include	

LIB_COMMON     = -lm -lpthread
LIBDIR_GUROBI  = $(GUROBI_HOME)/lib
LIB_GUROBI     = -L$(LIBDIR_GUROBI) -lgurobi_c++ -lgurobi81 

# ------------------------------------------------------------

clean :
	/bin/rm -rf *.o *_c *_c++ *.class *.log *.rlp *.lp *.bas *.ilp
	/bin/rm -rf *~ *.mps *.ord *.sos *.sav *.net *.msg *.clp

	$(GPP) $(OPT) $(INC_GUROBI) -o prog prog.cpp $(LIB_COMMON) $(LIB_GUROBI)
# ------------------------------------------------------------
