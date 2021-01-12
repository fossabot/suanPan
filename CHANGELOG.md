# Changelog

## version 1.x

### User

1. add phase field enabled elements: `DCP3`, `DCP4`, `DC3D4`, `DC3D8` elements
2. add support to record nodal damping/inertial force `DF` and `IF`
3. add regularized `Yeoh` model for compressible rubbers
4. improve stability of `RambergOsgood` model
5. add `LeeNewmarkFull` damping model, improve performance of `LeeNewmark` damping model
6. add shared memory `SuperLU` solver
7. add `Spike` solver for banded matrices
8. add displacement based beam element with end moment release: `B21EL` and `B21EH` elements
9. correct name of `Kelvin` model

### Developer

1. Validation of material/section is moved to element initialization. For all derived elements from `MaterialElement` and `SectionElement`, there is no need to validate material/section. As long as the element is active, material/section must have the correct type.

## version 1.0

1. initial release
