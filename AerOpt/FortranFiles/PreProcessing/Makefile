# Start of the makefile

# Defining variables
ifeq ($(OS),Windows_NT)
EXT = .exe
else
EXT =
endif
TARGET = PreProcessing$(EXT)
OBJECTS = Aggl2d_def.o Aggl2d_Aggl.o Aggl2d_main.o

# Find all sources
#SOURCES:=$(wildcard *.f90)
# Make targets out of them
#OBJECTS:=$(SOURCES:%.f90=%.o) #NOTE Order is important sometimes

# Set default output
default :
	make -s comment
comment :
	echo " "
	echo "   Makefile for Aggl2d        "
	echo "   -------------------------- "
	echo "   clean  to remove object files "              
	echo "   ifort  to compile and link with intel compiler"
	echo "   gfort  to compile and link with gfortran compiler"
	echo " "

ifort : COMP = ifort
ifort : FLAGS = -r8 -O3 -static -static-intel
ifort : $(OBJECTS)
	$(COMP) $(FLAGS) $(OBJECTS) -o $(TARGET)

gfort : COMP = gfortran
gfort : FLAGS = -fdefault-real-8 -O3 -static -static-libgfortran
gfort : $(OBJECTS)
	$(COMP) $(FLAGS) $(OBJECTS) -o $(TARGET)

# List of object files
%.o : %.f90
	$(COMP) $(FLAGS) -c $<

# Cleaning everything
clean : 
	find .. -name '*.o' -delete
	find .. -name '*.mod' -delete
# End of the makefile
