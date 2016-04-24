# set executable name
EXECNAME   	= EMC

###############################################
#	compiler
#
#CXX                     = g++-5
CXX                     = clang++-3.8
#OPT                     = -O3 -fopenmp -march=native -mno-avx -DEIGEN_NO_DEBUG -DNDEBUG
OPT                     = -fopenmp=libiomp5 -std=c++11
LDFLAGS                 = -Wextra
#LDFLAGS                 = -Wall -Werror -std=c++11 -Wno-deprecated-declarations
CXXFLAGS    =   $(LDFLAGS) $(OPT)
###############################################

# Directory structure
SOURCEDIR= source
OBJDIR= obj

SRC		:= $(wildcard $(SOURCEDIR)/*.cpp)
OBJS	:= $(addprefix $(OBJDIR)/,$(notdir $(SRC:.cpp=.o)))


# Compilation
all: $(EXECNAME)

$(EXECNAME): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

$(OBJDIR)/%.o: $(SOURCEDIR)/%.cpp
	@mkdir -p $(OBJDIR)
	$(CXX) $(CXXFLAGS) -c -o $@ $<



# debugging compilation
db:
	$(CXX) $(SRC) $(LDFLAGS) $(DEBUG) -o $(EXECNAME).db

# profiling compilation
prof:
	$(CXX) $(SRC) $(LDFLAGS) $(PROF) -o $(EXECNAME).prof

#depend: .depend
#.depend: $(SRC)
#	rm -f ./.depend
#	$(CXX) -MM $^  >> ./.depend;

#.PHONY: clean
.PHONY: clean
clean:
	$(RM) -r $(OBJDIR)
#clean:
#	$(RM) $(OBJS)

#dist-clean: clean
#	$(RM) *~ ./.depend $(EXECNAME).prof $(EXECNAME).db $(EXECNAME).db.dSYM $(EXECNAME)

#include .depend
