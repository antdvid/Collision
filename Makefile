CXX=mpicxx -g -fopenmp -Wno-unused-result
F77=mpif77 
F77_LIBS=

libext = 
incs =   -I/usr/local/pkg/mpich2/include -D__MPI__  -I/usr/local/pkg/hdf/include -DUSE_HDF   -D__HYPRE__  
libincs =  -L/usr/lib -L/usr/lib -L/usr/local/pkg/hdf/lib  -L/usr/local/pkg/mpich2/lib 
libs = -lgd -lmfhdf -ldf -ljpeg -lz  -lmpich -lpthread

collid3d: collid3d.o dcollid.o cdinit.o
	$(CXX) $^ -I../include -L../lib/x86_64 $(libincs) -lFronTier -lm $(libs) -lCGAL_Core -lCGAL_ImageIO -lCGAL -lgmp -frounding-math -o collid3d 

%.o: %.cpp
	${CXX} $< -c -I../include $(incs) -frounding-math
clean:
	rm -rf *.o collid3d 
tagsfile:
	ctags *.cpp ../src/*/*.[chf]

-include ../devel-deps.inc
