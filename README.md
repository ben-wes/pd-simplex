# simplex~
Pure Data (Pd) external for signal rate simplex noise sampling with 1d, 2d, 3d and 4d multichannel input

The current version's noise implementation is based on https://github.com/stegu/perlin-noise/blob/master/src/sdnoise1234.c by Stefan Gustavson.

![simplex~-help.pd screenshot](simplex%7E-help.pd.png)

## Build and Installation Instructions
To compile and install, follow these steps. Note that [Makefile.pdlibbuilder](https://github.com/pure-data/pd-lib-builder/) is included as submodule, requiring the `--recursive` parameter.
~~~
git clone --recursive https://github.com/ben-wes/pd-simplex.git
cd pd-simplex
make install
~~~

## Usage
### Inlets and Outlets
* **1st inlet:** expects a multichannel signal with up to 4 dimensions for the noise sampling position (the noise algorithm, and hence its performance, is determined by the number of channels)
* **2nd inlet:** can be used to control the persistence value with a signal or float inputs
* **outlet:** outputs the sampled noise value at signal rate

### Creation Arguments
#### Optional flags
* `-n` flag activates normalization of the octaves' sum to keep values in the [-1..1] range
* `-s <int>` flag initializes the noise's permutation table with a given seed
* `-d` flag activates derivatives multichannel output on additional 2nd outlet

#### Optional numerical arguments (must be written after the flags)
* **first arg** sets the number of octaves (defaults to `1`, max. is `24`)
* **second arg** sets initial persistence (defaults to `0.5`) 

### Messages
* `[seed <int>(` generates different (deterministic) permutation tables for the simplex function
* `[normalize 0/1(` (de)activates normalization (off by default)
* `[octaves <int>(` dynamically changes the number of octaves
* `[persistence <float>(` sets persistence (while no signal is connected)
* `[coeffs <float> <float> ... <float>(` let's you customize the octave scaling (default is `1, 2, 4, ...`) 

---

## Additional Information
### Octaves
In case of more than 1 octave, additional noise octaves get sampled, i.e. downscaled instances of the noise space (higher noise frequencies) are added. This downscaling is achieved by multiplying the sampling coordinate with `2^(octave-1)` - for example:
* input `(1, -3)` samples from `2*(1, -3) = (2,  -6)` for 2nd octave
* input `(1, -3)` samples from `4*(1, -3) = (4, -12)` for 3rd octave

### Persistence
The persistence value determines the influence of successive octaves on the final output. It is a multiplier applied to each octave’s amplitude. Lower persistence values cause the amplitudes to decrease rapidly with each octave, leading to a smoother noise pattern. Conversely, higher persistence values maintain stronger amplitudes in higher octaves, resulting in more detailed and rougher (high frequency) patterns - for example:
* persistence of `0.5` yields octaves' amplitudes `1, 0.5, 0.25, 0.125, 0.0625, ...`
* persistence of `0.9` yields octaves' amplitudes `1, 0.9, 0.81, 0.729, 0.6561, ...`

### Normalization
If normalization is activated (via creation argument or message), the sum of all octaves' amplitudes is normalized to 1.

---

## History

The initial code of this external was an experiment generated with ChatGPT4. The fascinating part of this was, that with a few feedback loops (feeding errors back to it), it managed to create a working Pd external. After some more research around simplex noise, the noise algorithm of the external turned out to be a pretty exact copy of [this SimplexNoise.cpp](https://github.com/SRombauts/SimplexNoise/blob/master/src/SimplexNoise.cpp) though (including variable names and comments).
