#CXX=mpicxx -g -fopenmp -Wno-unused-result
CXX=mpicxx -g -pedantic -Wall -Wextra -Wno-undef -Wno-comment -Wno-unused-parameter -Wno-long-long -Wctor-dtor-privacy -Wdisabled-optimization -Wformat=2 -Winit-self -Wlogical-op -Wmissing-declarations -Wmissing-include-dirs -Wnoexcept -Wold-style-cast -Woverloaded-virtual -Wredundant-decls -Wshadow -Wsign-conversion -Wsign-promo -Wstrict-null-sentinel -Wstrict-overflow=5 -Wswitch-default -Werror -Wno-old-style-cast -Wno-redundant-decls
F77=mpif77 
F77_LIBS=

libext = 
incs =   -I/usr/local/pkg/mpich2/include -D__MPI__  -I/usr/local/pkg/hdf/include -DUSE_HDF   -D__HYPRE__  
libincs =  -L/usr/lib -L/usr/lib -L/usr/local/pkg/hdf/lib  -L/usr/local/pkg/mpich2/lib 
libs = -lgd -lmfhdf -ldf -ljpeg -lz  -lmpich -lpthread

test_collid: test_collid.o dcollid.o cdinit.o dcollid3d.o
	$(CXX) $^ -I../include -L../lib/x86_64 $(libincs) -lFronTier -lm $(libs) -lCGAL_Core -lCGAL_ImageIO -lCGAL -lgmp -frounding-math -o test_collid 

%.o: %.cpp
	${CXX} $< -c -I../include $(incs) -frounding-math
clean:
	rm -rf *.o test_collid 
tagsfile:
	ctags *.cpp ../src/*/*.[chf]

-include ../devel-deps.inc
