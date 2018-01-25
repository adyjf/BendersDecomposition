#--
#-- Nom exe ... 
EXE = Benders

#-
#-- Liste des sources a compiler .. 
LIST_OF_OBJ_FILES = main.o

# ---------------------------------------------------------------------
# Compiler selection 
# ---------------------------------------------------------------------
CCC = g++

# ---------------------------------------------------------------------
# Compiler options 
# ---------------------------------------------------------------------

CCOPT :=
CCOPT += -g -O3 -fomit-frame-pointer -pipe -DNDEBUG -pedantic-errors -Wparentheses -Wreturn-type -Wcast-qual -Wall -Wpointer-arith -Wwrite-strings -Wconversion
CCOPT +=  -O -fPIC -fexceptions -DIL_STD -std=c++11

# ---------------------------------------------------------------------
# COIN-OR SETTING: Link options and libraries: C++
# ---------------------------------------------------------------------

COININCDIR = # .. A COMPLETER
COINLIBDIR = # .. A COMPLETER

COIN_LIBS = -L$(COINLIBDIR) -lCbc -lCgl -lOsiClp -lOsi -lClp -lCoinUtils \
	        -lm -llapack -lblas\
	 
COIN_INCS = -I$(COININCDIR) 

#------------------------------------------------------------
# AutoCompiling
#------------------------------------------------------------

TO_CLEAN = *~ $(LIST_OF_OBJ_FILES)

LIBS = $(COIN_LIBS)
INCL = $(COIN_INCS) $(PATH_TO_MY_HPP)

#------------------------------------------------------------
# AutoCompiling
#------------------------------------------------------------

all: $(EXE)

.SUFFIXES: .cpp .c .o .obj

$(EXE): $(LIST_OF_OBJ_FILES)
	bla=;\
	for file in $(LIST_OF_OBJ_FILES); do bla="$$bla `echo $$file`"; done; \
	$(CCC) $(CCOPT) -o $@ $$bla $(LIBS)

clean:
	rm -rf $(TO_CLEAN)

mrproper: clean
	rm -rf $(EXE)
	
doc:
	doxygen ../MakefileConf/doxygen.conf
	
.cpp.o:
	$(CCC) $(CCOPT) $(INCL) -c -o $@ `test -f '$<' || echo '$(SRCDIR)/'`$<

.cpp.obj:
	$(CCC) $(CCOPT) $(INCL) -c -o $@ `if test -f '$<'; then echo '$<'; else echo '$(SRCDIR)/$<'; fi`

.c.o:
	$(CC) $(COPT) $(INCL) -c -o $@ `test -f '$<' || echo '$(SRCDIR)/'`$<

.c.obj:
	$(CC) $(COPT) $(INCL) -c -o $@ `if test -f '$<'; then echo '$<'; else echo '$(SRCDIR)/$<'; fi`
