
# Command-line parameter to switch gdb on/off
# Use with -B to force a rebuild
GDB =

# Subdirectories with header files
INCLUDES = graph util test

# Default compiler and linker flags
override CXXFLAGS += -Wall $(MODE) -std=c++0x -frounding-math -fopenmp
override LDFLAGS += -lCGAL -lgmp -lboost_graph -lboost_thread -lboost_unit_test_framework -fopenmp -lrt

# Automatically find all sources and use implicit make rules
SRCS = $(shell find * -name \*.cpp)
OBJS = $(SRCS:.cpp=.o)
DEPS = $(OBJS:.o=.d)
EXEC = $(OBJS:.o=.exe)

RM = rm -f

# Using target-specific variables would be nicer,
#   but they force you to build the whole source tree
ifeq ($(GDB), off)
	MODE = -O4
else
	ifeq ($(GDB), on)
		MODE = -O0 -g -ggdb
	else
		MODE = -O2 -g
	endif
endif

.PHONY: all clean cleanall cleanobjs rules

all: $(EXEC)

# Use .exe suffix for easier integration with Git
%.exe : %.o
	$(CXX) -o $@ $< $(LDFLAGS)

cleanall: clean
	$(RM) **/*~ *.dot

clean: cleanobjs
	$(RM) $(EXEC)

cleanobjs:
	$(RM) $(OBJS) $(DEPS)

# Files, flags and rules for auto dependencies
.PRECIOUS: %.d %.o
override CXXFLAGS += -MP -MMD $(INCLUDES:%=-I%)
-include $(DEPS)

# Show variables for Makefile debugging
rules:
	@echo SRCS = $(SRCS)
	@echo EXEC = $(EXEC)
	@echo OBJS = $(OBJS)
	@echo DEPS = $(DEPS)
	@echo PWD = $(PWD)
	@echo CC = `$(CC) --version | head -n 1`
	@echo CXX = `$(CXX) --version | head -n 1`
	@echo CXXFLAGS = $(CXXFLAGS)
	@echo LDFLAGS = $(LDFLAGS)
