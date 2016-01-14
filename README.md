# C++ library for robust collision handling
##Introduction
This is a library implenmented to handle fabric-fabric or fabric-rigid body collision.
The fabric surface is modeled with triangle mesh and its interior dynmaics is simulated with spring model (not included in this library).
For convenience, the data structure SURFACE, TRI, POINT are from FronTier library (a computational fluid dynamic library with interface tracking capability from Stony Brook University). However, people can define their own 
data structure as long as the SURFACE has a link list for triangles, TRI has pointers to three POINTs and POINT has coordinates(x,y,z).
##Numerical Test
The basic idea of the implementation is originated from the paper "Robust Treatment of Collisions, Contact and Fricition for Cloth Animation", but we modify the method to couple it with FronTier easily. A few numerical test are shown below:
<img style="float: left;" src="http://guest.ams.sunysb.edu/~zgao/work/collision/fall-sphere.gif" width="230">
<img style="float: left;" src="http://guest.ams.sunysb.edu/~zgao/work/collision/fall-body.gif" width="230">
<img style="float: left;" src="http://guest.ams.sunysb.edu/~zgao/work/collision/fall-box.gif" width="230">

