# --- Compiler ---
CXX = __CXX__
# Note: MAKE is implicitly defined, no need to export unless overriding

# --- Paths to Dependencies ---
EIGEN3_PATH    = __EIGEN3__
# Eigen3: The path where you can find "Eigen/", "signature_of_eigen3_matrix_library" and "unsupported/".

# --- Project Structure ---
SRCDIR      = src
INCDIR      = include
OBJDIR      = obj
LIBDIR      = lib
SOURCES     = $(shell find $(SRCDIR) -name '*.cpp')
# Generate corresponding object file paths in OBJDIR, preserving subdirectory structure
OBJECTS     = $(patsubst $(SRCDIR)/%.cpp,$(OBJDIR)/%.o,$(SOURCES))
# Generate dependency files (optional but good practice for header changes)
DEPS        = $(OBJECTS:.o=.d)

# --- Build Flags ---
# General Compiler Flags (apply to compilation steps)
# Using -isystem for external libraries suppresses warnings from their headers
# Added -MMD -MP to generate dependency files (.d)
CPPFLAGS    = -isystem $(EIGEN3_PATH) \
              -DEIGEN_INITIALIZE_MATRICES_BY_ZERO \
              -MMD -MP # Generate dependency files

CXXFLAGS    = -Wall -Wextra -Wpedantic -fPIC -march=native -O3 -std=c++2a

# --- Main Rules ---

# Default target: compile all sources and expose the header
.PHONY: all
all: $(OBJECTS) | $(INCDIR) $(LIBDIR)
	@echo "Compilation done!"
	ar -rv $(LIBDIR)/libmwfn.a $^
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -shared -o $(LIBDIR)/libmwfn.so $^
	@echo "The lib files libmwfn.a and libmwfn.so are put in $(LIBDIR)/."
	cat src/NecessaryHeaders.h > $(INCDIR)/libmwfn.h
	cat src/MwfnShell.h >> $(INCDIR)/libmwfn.h
	cat src/MwfnCenter.h >> $(INCDIR)/libmwfn.h
	cat src/MwfnOrbital.h >> $(INCDIR)/libmwfn.h
	cat src/Mwfn.h >> $(INCDIR)/libmwfn.h
	sed -i "s/EigenMatrix/Eigen::MatrixXd/g" $(INCDIR)/libmwfn.h
	sed -i "s/EigenVector/Eigen::VectorXd/g" $(INCDIR)/libmwfn.h
	@echo "The header libmwfn.h is put in $(INCDIR)/."

# --- Compilation Rule ---

# Pattern rule to compile .cpp files from SRCDIR into .o files in OBJDIR
# This handles source files in subdirectories of SRCDIR as well.
$(OBJDIR)/%.o: $(SRCDIR)/%.cpp | $(OBJDIR) # Use order-only prerequisite for OBJDIR
	@echo "Compiling $< -> $@ ..."
	@mkdir -p $(@D) # Create subdirectory in obj/ if it doesn't exist
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@
# $< expands to the first prerequisite (the .cpp file)
# $@ expands to the target (the .o file)
# $(@D) expands to the directory part of the target

# --- Utility Rules ---

# Rule to create the include directory
$(INCDIR):
	@echo "Creating directory $@"
	@mkdir -p $@

# Rule to create the lib directory
$(LIBDIR):
	@echo "Creating directory $@"
	@mkdir -p $@

# Rule to create the object directory
# This is triggered by the order-only prerequisite in the compilation rule
$(OBJDIR):
	@echo "Creating directory $@"
	@mkdir -p $@

# Rule to clean up generated files
.PHONY: clean
clean:
	@echo "Cleaning..."
	rm -rf $(INCDIR)  # Remove header
	rm -rf $(LIBDIR)  # Remove libraries
	rm -rf $(OBJDIR)  # Remove object directory and all its contents (.o, .d files)

# Include dependency files, if they exist
# This makes Make automatically recompile files if included headers change
-include $(DEPS)
