# Contour Parametrization via Anisotropic Mean Curvature Flow (AMCF)

This module is an implementation of the Mean Curvature Flow with constraints. The available code was used to generate the examples in [Link]().

The mean curvature flow (MCF) deforms a hypersurface in the normal direction with a speed proportional to its mean curvature. An example of a constrained MCF is the evolution of a closed planar curve enclosing an unique object. This code evolves a curve drawn around a single object picture, the evolution stops at the points which reach the object.

## Python Dependencies 

- pycuda
- PIL
- pylab
- numpy
- matplotlib

## Run 

- Edit line 18 in "ImageRec_v4.py" to select the image.
- Run "ImageRec_v4.py".
