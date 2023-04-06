# Contour Parametrization via Anisotropic Mean Curvature Flow (AMCF)

This module is an implementation of the Mean Curvature Flow with constraints. The available code was used to generate the examples  [here](http://arxiv.org/abs/1803.03724).

The mean curvature flow (MCF) deforms a hypersurface in the normal direction with a speed proportional to its mean curvature. An example of a constrained MCF is the evolution of a closed planar curve enclosing an unique object. This code evolves a curve drawn around a single object picture, the evolution stops at the points which reach the object.

#### Examples

![convex2](https://user-images.githubusercontent.com/36924228/37216935-c83e5b6e-2381-11e8-864b-15df43a98197.gif)![bean2](https://user-images.githubusercontent.com/36924228/37216917-b2243f7e-2381-11e8-8d95-18383fa96c08.gif)


## Python Dependencies 

- pycuda
- PIL
- pylab
- numpy
- matplotlib

## Run 

- Edit line 18 in "ImageRec_v4.py" to select the image.
- Run "ImageRec_v4.py".

##### Note:
- <s>For Python 3.X, 'print' syntax must be modified.</s>
- A new version for Python 3 is now included: **ImageRec_v4_PY3.py**.
