#
#  Update 2020-01-29:
#  o Uses required argument of serial, omp, or mpi.
#  o Use VPATH for finding cpp file in different directories -- this simplifies rules
#  o Abort if gsl-config isn't available
#  o Fixed (INTEL) compiler search 0=found | 1=notfound ; make conditional simple (ifeq 0|1)
#  o Also use conditional for GCC
#  o For objects, use basename to get base file name
#  o Clean up directory prefix (shorten variable names and group)
#  o Simplified obj and bin rule logic and readability.
#  o put rules in canonical order
#  o Now has PROF for profiling. (This is by default overrided with empty PROF.)
#  o Now uses INCS. CXXFLAGS is used for C++ specific options.
#  o Make executables with suffixes ( nerdss_serial | nerdss_mpi | nerdss_omp).
#  --  a bit cleaner                                            Kent milfeld@tacc.utexas.edu
#
# TODO: use function to create VPATH
# TODO: Fix MPI after learning purpose
# TODO: Make rules for *.hpp's
#
# Set terminal width to 220 to avoid viewing wrapped lines in output. A width of 200 avoids most wrapping.
#
# Update 2025-08-25:
# o dded new PHONY: "debug" and "profile". 
# o Use `make serial debug` to debug with gdb
# o use `make serial profile` to profile
#

BDIR   = bin
SDIR   = src
EDIR   = EXEs

PROF   =

.PHONY: serial mpi debug profile asan ubsan sanitizers static-analysis clean syntax

# ---------------- REQUIREMENTS: gsl and directories
hasGSL = $(shell type gsl-config >/dev/null 2>&1; echo $$?)
ifeq ($(hasGSL),1)
$(error " GSL must be installed, and gsl-config must be in path.")
else
$(shell mkdir -p bin)
$(shell mkdir -p obj)
endif

# ---------------- EXECUTABLE SETUP
INCLUDE_FOLDERS = boundary_conditions classes error math parser reactions system_setup trajectory_functions io

BUILD_GOALS = $(filter serial mpi debug profile asan ubsan sanitizers,$(MAKECMDGOALS))
SERIAL_GOALS = $(filter serial debug profile asan ubsan sanitizers,$(MAKECMDGOALS))
SANITIZER_GOALS = $(filter asan ubsan sanitizers,$(MAKECMDGOALS))

ifneq (,$(SERIAL_GOALS))
	_EXEC = nerdss
	ENTRYPOINT = nerdss
endif

ifneq (,$(filter mpi,$(MAKECMDGOALS)))
	_EXEC = nerdss_mpi
	ENTRYPOINT = nerdss_mpi
	DEFS = -Dmpi_
	INCLUDE_FOLDERS += debug io_mpi mpi
endif

ifneq (,$(filter mpi,$(MAKECMDGOALS)))
ifneq (,$(SANITIZER_GOALS))
$(error Sanitizer targets are currently supported for serial builds only.)
endif
endif

ifneq (,$(filter debug,$(MAKECMDGOALS)))
	ENABLE_DEBUG = true
endif

ifneq (,$(filter profile,$(MAKECMDGOALS)))
	ENABLE_PROFILING = true
endif

ifneq (,$(SANITIZER_GOALS))
	ENABLE_DEBUG = true
endif

SRCS = $(foreach dir,$(INCLUDE_FOLDERS),$(wildcard $(SDIR)/$(dir)/*.cpp))

OS    := $(shell uname)
INTEL = $(shell type icpc  >/dev/null 2>&1; echo $$?)
GCC   = $(shell type g++   >/dev/null 2>&1; echo $$?)

INCS    = $(shell gsl-config --cflags) -Iinclude
CXXFLAGS = -std=c++0x
LIBS     = $(shell gsl-config --libs)
SANITIZER_FLAGS =
SANITIZER_NAME =
EXEC_SUFFIX =
BUILD_VARIANT = release

# ---------------- COMPILER SETUP
override PROF   =

ifeq ($(GCC),0)
	CC      = g++
	ifneq (,$(filter mpi,$(MAKECMDGOALS)))
		CC = mpicxx
	endif
	CFLAGS  = -O3 # use -O2 if profiling is confused by optimization
endif

ifeq ($(INTEL),0)
	CC      = icpc
	ifneq (,$(filter mpi,$(MAKECMDGOALS)))
		CC = mpicxx
	endif
	CFLAGS  = -O3 # use -O2 if profiling is confused by optimization
endif

# ---------------- Feature toggles
# Set debug flags if ENABLE_DEBUG is true
ifdef ENABLE_DEBUG
	BUILD_VARIANT = debug
	CFLAGS = -g -O0 -fno-omit-frame-pointer
	CXXFLAGS += -DDEBUG
endif

ifneq (,$(filter asan,$(MAKECMDGOALS)))
	SANITIZER_FLAGS += -fsanitize=address
	SANITIZER_NAME = asan
endif

ifneq (,$(filter ubsan,$(MAKECMDGOALS)))
	SANITIZER_FLAGS += -fsanitize=undefined
	SANITIZER_NAME = ubsan
endif

ifneq (,$(filter sanitizers,$(MAKECMDGOALS)))
	SANITIZER_FLAGS += -fsanitize=address -fsanitize=undefined
	SANITIZER_NAME = sanitizers
endif

ifneq (,$(filter asan,$(MAKECMDGOALS)))
ifneq (,$(filter ubsan,$(MAKECMDGOALS)))
	SANITIZER_NAME = sanitizers
endif
endif

ifneq (,$(SANITIZER_FLAGS))
	BUILD_VARIANT = $(SANITIZER_NAME)
	EXEC_SUFFIX = _$(SANITIZER_NAME)
	CFLAGS = -g -O0 -fno-omit-frame-pointer $(SANITIZER_FLAGS)
	LIBS += $(SANITIZER_FLAGS)
endif

# Set profiling flags if ENABLE_PROFILING is true
ifdef ENABLE_PROFILING
	BUILD_VARIANT := $(BUILD_VARIANT)-profile
	PROF += -pg
	CFLAGS += -DENABLE_PROFILING 
	LIBS += $(shell pkg-config --libs libprofiler)
endif

ODIR = obj/$(BUILD_VARIANT)
EXEC = $(patsubst %,$(BDIR)/%$(EXEC_SUFFIX),$(_EXEC))

# ---------------- OBJECT FILES
OBJS = $(patsubst $(SDIR)/%.cpp,$(ODIR)/%.o,$(SRCS))

# ---------------- RULES
syntax:
	@echo "------------------------------------"
	@printf '\033[31m%s\033[0m\n' " USAGE: make serial|mpi [debug] [profile]"
	@echo "------------------------------------"
	exit 0

ifneq ($(strip $(BUILD_GOALS)),)
$(BUILD_GOALS): $(EXEC)
	@echo "Finished making (re-)building $(MAKECMDGOALS) version, $(EXEC)."
endif

$(EXEC): $(OBJS)
	@echo "Compiling $(EDIR)/$(ENTRYPOINT).cpp"
	$(CC) $(CFLAGS) $(CXXFLAGS) $(INCS) $(PROF) -o $@ $(EDIR)/$(ENTRYPOINT).cpp $(OBJS) $(LIBS) $(PLANG)
	@echo "------------"

$(ODIR)/%.o: $(SDIR)/%.cpp
	@echo "Compiling $< to $@"
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) $(CXXFLAGS) $(INCS) $(PROF) -c $< -o $@ $(PLANG) $(DEFS)
	@echo "------------"

clean:
	rm -rf obj bin

static-analysis:
	tools/run_static_analysis.sh


# Reference: https://www.gnu.org/software/make/manual/html_node/Quick-Reference.html
#            https://www.gnu.org/software/make/
#            https://www.cmcrossroads.com/article/basics-vpath-and-vpath
#            https://www.gnu.org/software/make/manual/html_node/Implicit-Variables.html
