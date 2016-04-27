# check the hostname
HN := $(shell /bin/hostname -s)

# default compiler settings
CXX         =   g++
DEBUG       =   -g  -Wall -Werror # -ansi
PROF        =   -pg	-fno-omit-frame-pointer -O2 -DNDEBUG -fno-inline-functions 
LDFLAGS     =  -fopenmp -std=c++11 -Wno-deprecated-declarations
OPT			=  -O3 -march=native -mno-avx -DEIGEN_NO_DEBUG -DNDEBUG
CXXFLAGS	=  $(LDFLAGS) $(OPT)


###############################################
#	Change compiler depending on the host machine
#
ifneq (,$(findstring ThinkPad,$(HN)))
CXX         = g++
OPT         = -O3 -march=native -mno-avx -DEIGEN_NO_DEBUG -DNDEBUG
LDFLAGS     = -Wall -Werror -fopenmp -std=c++11 -Wno-deprecated-declarations
CXXFLAGS    =   $(LDFLAGS) $(OPT)
endif


# set executable name
EXECNAME   	= EMC

###############################################
#	clang compiler
#
#CXX                     = g++
#CXX                     = clang++-3.8
#OPT                     = -fopenmp=libiomp5 -std=c++11
#LDFLAGS                 = -Wextra
#CXXFLAGS    =   $(LDFLAGS) $(OPT)
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

.PHONY: clean

clean:
	$(RM) $(OBJS)

