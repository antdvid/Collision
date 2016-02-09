CXX=mpicxx -O3 -Wno-unused-result
#comment the following if debugging needed
#CXX=mpicxx -g -pedantic -Wall -Wextra -Wno-undef -Wno-comment -Wno-unused-parameter -Wno-long-long -Wctor-dtor-privacy -Wdisabled-optimization -Wformat=2 -Winit-self -Wlogical-op -Wmissing-declarations -Wmissing-include-dirs -Wnoexcept -Wold-style-cast -Woverloaded-virtual -Wredundant-decls -Wshadow -Wsign-conversion -Wsign-promo -Wstrict-null-sentinel -Wstrict-overflow=5 -Wswitch-default -Werror -Wno-old-style-cast -Wno-redundant-decls
F77=mpif77 
F77_LIBS=

USE_VTK=
VTK_LIB=
VTK_INCL=
#uncomment the following if vtk lib is used
#USE_VTK=-D__VTK__
#VTK_LIB=/usr/lib/x86_64-linux-gnu/libvtk*.so
#VTK_INCL=/usr/include/vtk-6.0

incs =   -I/usr/local/pkg/mpich2/include -D__MPI__  -I/usr/local/pkg/hdf/include $(USE_VTK) -DUSE_HDF   -D__HYPRE__ -I../include 
libincs =  -L/usr/lib -L/usr/lib -L/usr/local/pkg/hdf/lib  -L/usr/local/pkg/mpich2/lib -L../lib/x86_64
libs = -lgd -lmfhdf -ldf -ljpeg -lz  -lmpich -lpthread $(VTK_LIB)

test: test.o dcollid.o cdinit.o dcollid3d.o
	$(CXX) $^ $(libincs) -lFronTier -lm $(libs) -lCGAL_Core -lCGAL_ImageIO -lCGAL -lgmp -frounding-math -o test 

vtk.o: vtk.cpp
	$(CXX) vtk.cpp -c $(incs) -I$(VTK_INCL)

%.o: %.cpp
	${CXX} $< -c -I../include $(incs) -frounding-math

-include ../devel-deps.inc

clean:
	rm -rf *.o test 
tagsfile:
	ctags *.cpp ../src/*/*.[chf]

