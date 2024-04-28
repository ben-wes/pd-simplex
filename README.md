# simplex~
Pure Data (Pd) external for signal rate simplex noise with 1d, 2d, 3d and 4d multichannel input

The initial code was an experiment generated via chatgpt. after some research, it was pretty obviously a quite precise copy of https://github.com/SRombauts/SimplexNoise/blob/master/src/SimplexNoise.cpp though.

The current version is based on https://github.com/stegu/perlin-noise/blob/master/src/simplexnoise1234.c

## Installation
Compilation requires 
~~~
git clone --recurse-submodules https://github.com/ben-wes/simplex.git
cd simplex
make install
~~~

## Usage
* the left inlet expects a multichannel signal with up to 4 dimensions for the noise sampling position.
* the right inlet can be used to control the persistence value with a signal

### Creation args
* optional `-n` flag at the beginning activates normalization of the octaves' sum to keep values in the `[-1..1]` range
* first numerical argument sets the number of octaves (defaults to 1)
* second numerical argument sets the persistence (while right inlet is not connected)

### Messages
* `[seed <int>(` generates different (deterministic) permutation tables for the simplex function
* `[normalize 0/1(` (de)activates normalization
* `[octaves <int>(` dynamically changes the number of octaves
