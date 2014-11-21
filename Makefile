
# Files to build
EXEC = eval.exe gen.exe prototype/speed.exe prototype/complexity.exe
SRCS = $(EXEC:.exe=.cpp)

# Subdirectories with header files
INCLUDES = graph util test

# Changes below this line shouldn't be needed

OBJS = $(SRCS:.cpp=.o)
DEPS = $(SRCS:.cpp=.d)

RM = rm -f

# Target-specific build mode
MODE = -O2 -g
release : MODE = -O2
debug : MODE = -O0 -g -ggdb

CXXFLAGS = -Wall $(MODE) -std=c++0x

# Automatic dependencies
.PRECIOUS: %.d %.o
CXXFLAGS += -MP -MMD $(INCLUDES:%=-I%)

# Targets

.PHONY : all release debug clean clean-intermediate rules

all : $(EXEC)
release : $(EXEC)
debug : $(EXEC)

%.exe : %.o
	$(CXX) -o $@ $< $(LDFLAGS)

# Automatic dependencies
-include $(DEPS)

clean: clean-intermediate
	$(RM) $(EXEC)

clean-intermediate:
	$(RM) $(OBJS) $(DEPS)

rules:
	@echo EXEC = $(EXEC)
	@echo OBJS = $(OBJS)
	@echo SRCS = $(SRCS)
	@echo PWD = $(PWD)
	@echo CC = $(CC)
	@echo CXX = $(CXX)
	@$(CXX) --version
	@echo CXXFLAGS = $(CXXFLAGS)
	@echo LDFLAGS = $(LDFLAGS)
