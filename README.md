# Hermite interpolation symbolic coefficient generator
 
## Implemented orders:

* 3D interpolation at 2nd-order accuracy (in `3d_2o`)
* 3D interpolation at 3rd-order accuracy (in `3d_3o`)

## Validation notes:

Interpolation coefficients validated to agree perfectly with the trusted, `Maple`-generated codes in [`AEILocalInterp`](https://bitbucket.org/cactuscode/numerical/src/master/AEILocalInterp/)

Validation dependencies: NRPy, installed via `pip install nrpy` , [`GNU sed`](https://www.gnu.org/software/sed/).

Instructions:
 1. Go into appropriate directory
 1. Run ./validate.sh
