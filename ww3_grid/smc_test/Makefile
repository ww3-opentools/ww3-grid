
FC=ifort
FFLAGS=-r8 -mcmodel=medium
TARGET_DIR=./

.SUFFIXES: .F90 .f90 .o

OBJS=smcProps.o

$(TARGET_DIR)/genSides: ${OBJS}
	$(FC) $(FFLAGS) ${OBJS} -o $(TARGET_DIR)/smcProps

## General rule for building objects from f90 source files:
.F90.o:
	$(FC) $(FFLAGS) -c $< -o $@

.f90.o:
	$(FC) $(FFLAGS) -c $< -o $@

clean:
	rm -f *.o *.mod $(TARGET_DIR)/smcProps
