CC = gcc
LFLAGS = -lnetcdf `pkg-config --cflags --libs glib-2.0 gmodule-2.0` -lm `pkg-config --cflags --libs python3`
CFLAGS = -g $(LFLAGS)
COMPILE = $(CC) $(CFLAGS) -c
OBJFILES = bathy_interp.o
HEADERS =  
SRCS_C = bathy_interp.c
SRCS = $(SRCS_C) $(HEADERS)


# defaults
resolution=1
bathy_src=etopo1

all: bathy_interp.exe

bathy_interp.exe: $(OBJFILES)
	$(CC) $(OBJFILES) $(LFLAGS) -o bathy_interp.exe

bathy_interp.o: bathy_interp.c
	$(COMPILE) bathy_interp.c

clean:
	rm -f *.o

run: bathy_interp.exe
	./bathy_interp.exe 0 360 2 -70 70 2 in.nc out.nc

glob: bathy_interp.exe
	# ./bathy_interp.exe 0 360 $(resolution) -80 80 $(resolution) /source/gridgen/noaa/reference_data/$(bathy_src).nc glob-$(bathy_src)-$(resolution).nc  
	# ./bathy_interp.exe 0 360 $(resolution) -80 80 $(resolution) ./$(bathy_src).nc glob-$(bathy_src)-$(resolution).nc  
	./bathy_interp.exe 

debug: bathy_interp.exe
	# ./bathy_interp.exe 0 360 $(resolution) -80 80 $(resolution) /source/gridgen/noaa/reference_data/$(bathy_src).nc glob-$(bathy_src)-$(resolution).nc  
	valgrind ./bathy_interp.exe 0 360 $(resolution) -80 80 $(resolution) ./$(bathy_src).nc glob-$(bathy_src)-$(resolution).nc  

nz: bathy_interp.exe
	./bathy_interp.exe 165 180 $(resolution) -50 -30 $(resolution) ../../example_data/$(bathy_src).nc nz-$(bathy_src)-$(resolution).nc  
	ncview nz-$(bathy_src)-$(resolution).nc

view: 
	ncview /tmp/glob.nc

TAGS: $(SRCS)
	etags $(SRCS)
tags: $(SRCS)
	etags $(SRCS)
