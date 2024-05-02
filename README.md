# simplex~
Pure Data (Pd) external for signal rate simplex noise with 1d, 2d, 3d and 4d multichannel input

The initial code was an experiment generated via ChatGPT4. The fascinating part of this was that with a few feedback loops, it managed to create a working Pd external. After some research, the noise part of the external turned out to be a pretty exact copy of [this SimplexNoise.cpp](https://github.com/SRombauts/SimplexNoise/blob/master/src/SimplexNoise.cpp) though.

The current version is based on https://github.com/stegu/perlin-noise/blob/master/src/simplexnoise1234.c with minor changes.

## Installation
Compilation requires [Makefile.pdlibbuilder](https://github.com/pure-data/pd-lib-builder/), which is included as a submodule here.
~~~
git clone --recurse-submodules https://github.com/ben-wes/simplex.git
cd simplex
make install
~~~

## Usage
* left inlet expects a multichannel signal with up to 4 dimensions for the noise sampling position. The noise algorithm (and therefore its performance) is selected according to the number of channels.
* right inlet can be used to control the persistence value with a signal (or float inputs).
* output is the sampled noise in signal rate.

#### Creation arguments
* optional `-n` flag at the beginning activates normalization of the octaves' sum to keep values in the [-1..1] range
* first numerical argument sets the number of octaves (defaults to 1)
* second numerical argument sets initial persistence (defaults to 0.5) 

#### Messages
* `[seed <int>(` generates different (deterministic) permutation tables for the simplex function
* `[normalize 0/1(` (de)activates normalization
* `[octaves <int>(` dynamically changes the number of octaves

## Additional information
#### Octaves
In case of more than 1 octave, additional noise octaves get sampled, i.e. downscaled instances of the noise space (higher noise frequencies) are added. This downscaling is achieved by multiplying the coordinate with `2^(octave-1)` - for example:
* 2d coordinate `(1, -3)` samples from `2*(1, -3) = (2,  -6)` for 2nd octave
* 2d coordinate `(1, -3)` samples from `4*(1, -3) = (4, -12)` for 3rd octave

#### Persistence
The persistence value determines the influence of successive octaves on the final output. It is a multiplier applied to each octave’s amplitude. Lower persistence values cause the amplitudes to decrease rapidly with each octave, leading to a smoother noise pattern. Conversely, higher persistence values maintain stronger amplitudes in higher octaves, resulting in more detailed and rougher (high frequency) patterns - for example:
* persistence of `0.5` yields octaves' amplitudes `1, 0.5, 0.25, 0.125, 0.0625, ...`
* persistence of `0.9` yields octaves' amplitudes `1, 0.9, 0.81, 0.729, 0.6561, ...`

If normalization is activated (through creation argument or message), the sum of all octaves' amplitudes gets normalized to 1.
