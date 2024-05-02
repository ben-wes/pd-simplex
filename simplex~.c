#include "m_pd.h"
#include <math.h>
#include <time.h>
#include <stdlib.h>

#define max(a,b) ( ((a) > (b)) ? (a) : (b) )
#define fastfloor(x) ( ((int)(x)<=(x)) ? ((int)x) : (((int)x)-1) )

#define PERSISTENCE 0.5f

#define F2   0.36602540378f // 0.5*(sqrt(3.0)-1.0);

#define G2   0.2113248654f  // (3.0-sqrt(3.0))/6.0;
#define G2_2 0.42264973081f

#define F3   0.33333333333f // 1.0/3.0;
#define G3   0.16666666666f // 1.0/6.0
#define G3_2 0.33333333333f
#define G3_3 0.5f

#define F4   0.30901699437f // (sqrt(5.0)-1.0)/4.0;
#define G4   0.13819660112f // (5.0-sqrt(5.0))/20.0;
#define G4_2 0.27639320225f
#define G4_3 0.41458980337f
#define G4_4 0.5527864045f

static t_class *simplex_tilde_class;

typedef struct _simplex_tilde {
    t_object x_obj;
    t_inlet *x_inlet_persistence;
    int x_normalize;
    int x_octaves;
    int x_perm[512];
} t_simplex_tilde;

// simplex noises code from https://github.com/stegu/perlin-noise/blob/master/src/simplexnoise1234.c

static const signed char grad[] = {1,2,3,4,5,6,7,8,-1,-2,-3,-4,-5,-6,-7,-8};

static inline t_float grad1( int hash, t_float x ) {
    return grad[hash & 15] * x;
}

static inline t_float grad2( int hash, t_float x, t_float y ) {
    int h = hash & 7;      // Convert low 3 bits of hash code
    t_float u = h<4 ? x : y;  // into 8 simple gradient directions,
    t_float v = h<4 ? y : x;  // and compute the dot product with (x,y).
    return ((h&1)? -u : u) + ((h&2)? -2.0f*v : 2.0f*v);
}

static inline t_float grad3( int hash, t_float x, t_float y , t_float z ) {
    int h = hash & 15;     // Convert low 4 bits of hash code into 12 simple
    t_float u = h<8 ? x : y; // gradient directions, and compute dot product.
    t_float v = h<4 ? y : h==12||h==14 ? x : z; // Fix repeats at h = 12 to 15
    return ((h&1)? -u : u) + ((h&2)? -v : v);
}

static inline t_float grad4( int hash, t_float x, t_float y, t_float z, t_float t ) {
    int h = hash & 31;      // Convert low 5 bits of hash code into 32 simple
    t_float u = h<24 ? x : y; // gradient directions, and compute dot product.
    t_float v = h<16 ? y : z;
    t_float w = h<8 ? z : t;
    return ((h&1)? -u : u) + ((h&2)? -v : v) + ((h&4)? -w : w);
}

// A lookup table to traverse the simplex around a given point in 4D.
static const unsigned char simplex4[64][4] = {
    {0,1,2,3}, {0,1,3,2}, {0,0,0,0}, {0,2,3,1}, {0,0,0,0}, {0,0,0,0}, {0,0,0,0}, {1,2,3,0},
    {0,2,1,3}, {0,0,0,0}, {0,3,1,2}, {0,3,2,1}, {0,0,0,0}, {0,0,0,0}, {0,0,0,0}, {1,3,2,0},
    {0,0,0,0}, {0,0,0,0}, {0,0,0,0}, {0,0,0,0}, {0,0,0,0}, {0,0,0,0}, {0,0,0,0}, {0,0,0,0},
    {1,2,0,3}, {0,0,0,0}, {1,3,0,2}, {0,0,0,0}, {0,0,0,0}, {0,0,0,0}, {2,3,0,1}, {2,3,1,0},
    {1,0,2,3}, {1,0,3,2}, {0,0,0,0}, {0,0,0,0}, {0,0,0,0}, {2,0,3,1}, {0,0,0,0}, {2,1,3,0},
    {0,0,0,0}, {0,0,0,0}, {0,0,0,0}, {0,0,0,0}, {0,0,0,0}, {0,0,0,0}, {0,0,0,0}, {0,0,0,0},
    {2,0,1,3}, {0,0,0,0}, {0,0,0,0}, {0,0,0,0}, {3,0,1,2}, {3,0,2,1}, {0,0,0,0}, {3,1,2,0},
    {2,1,0,3}, {0,0,0,0}, {0,0,0,0}, {0,0,0,0}, {3,1,0,2}, {0,0,0,0}, {3,2,0,1}, {3,2,1,0}
};

static const unsigned char simplex3[8][6] = {
    // The entries correspond to (i1, j1, k1, i2, j2, k2) for each condition
    {0,0,1,0,1,1}, // ZYX 0
    {0,0,0,0,0,0},
    {0,1,0,0,1,1}, // YZX 2:   2
    {0,1,0,1,1,0}, // YXZ 3:   2 1
    {0,0,1,1,0,1}, // ZXY 4: 4
    {1,0,0,1,0,1}, // XZY 5: 4   1
    {0,0,0,0,0,0},
    {1,0,0,1,1,0}  // XYZ 7: 4 2 1
};

void shuffle(int *array, int n, unsigned int *seed) {
    srand(*seed);
    for (int i = n - 1; i > 0; i--) {
        int j = rand() % (i + 1);
        int temp = array[i];
        array[i] = array[j];
        array[j] = temp;
    }
}

void initPermutation(t_simplex_tilde *x, unsigned int seed) {
    int basePermutation[256];
    for (int i = 0; i < 256; i++) {
        basePermutation[i] = i;
    }
    shuffle(basePermutation, 256, &seed);
    for (int i = 0; i < 256; i++) {
        x->x_perm[i] = basePermutation[i];
        x->x_perm[i+256] = basePermutation[i];
    }
}

// 1D simplex noise
t_float snoise1(t_simplex_tilde *x, t_float x_in) {
    int i0 = fastfloor(x_in);
    int i1 = i0 + 1;
    t_float x0 = x_in - i0;
    t_float x1 = x0 - 1.0f;

    t_float n0, n1;

    t_float t0 = 1.0f - x0*x0;
    t0 *= t0;
    n0 = t0 * t0 * grad1(x->x_perm[i0 & 0xff], x0);

    t_float t1 = 1.0f - x1*x1;
    t1 *= t1;
    n1 = t1 * t1 * grad1(x->x_perm[i1 & 0xff], x1);
    return 0.395 * (n0 + n1);
}

// 2D simplex noise
t_float snoise2(t_simplex_tilde *x, t_float x_in, t_float y_in) {
    t_float n0, n1, n2; // Noise contributions from the three corners

    // Skew the input space to determine which simplex cell we're in
    t_float s = (x_in+y_in)*F2; // Hairy factor for 2D
    t_float xs = x_in + s;
    t_float ys = y_in + s;
    int i = fastfloor(xs);
    int j = fastfloor(ys);

    t_float t = (float)(i+j)*G2;
    t_float X0 = i-t; // Unskew the cell origin back to (x,y) space
    t_float Y0 = j-t;
    t_float x0 = x_in-X0; // The x,y distances from the cell origin
    t_float y0 = y_in-Y0;

    // Offsets for second (middle) corner of simplex in (i,j) coords
    int i1 = x0>y0;
    int j1 = !i1;
    // lower triangle, XY order: (0,0)->(1,0)->(1,1)
    // upper triangle, YX order: (0,0)->(0,1)->(1,1)

    t_float x1 = x0 - i1 + G2; // Offsets for middle corner in (x,y) unskewed coords
    t_float y1 = y0 - j1 + G2;
    t_float x2 = x0 - 1.0f + G2_2; // Offsets for last corner in (x,y) unskewed coords
    t_float y2 = y0 - 1.0f + G2_2;

    // Wrap the integer indices at 256, to avoid indexing x->x_perm[] out of bounds
    int ii = i & 0xff;
    int jj = j & 0xff;

    // Calculate the contribution from the three corners
    t_float t0 = 0.5f - x0*x0-y0*y0;
    t0 = max(t0, 0.0f);
    t0 *= t0;
    n0 = t0 * t0 * grad2(x->x_perm[ii+x->x_perm[jj]], x0, y0); 

    t_float t1 = 0.5f - x1*x1-y1*y1;
    t1 = max(t1, 0.0f);
    t1 *= t1;
    n1 = t1 * t1 * grad2(x->x_perm[ii+i1+x->x_perm[jj+j1]], x1, y1);

    t_float t2 = 0.5f - x2*x2-y2*y2;
    t2 = max(t2, 0.0f);
    t2 *= t2;
    n2 = t2 * t2 * grad2(x->x_perm[ii+1+x->x_perm[jj+1]], x2, y2);

    // Add contributions from each corner to get the final noise value.
    // The result is scaled to return values in the interval [-1,1].
    return 40.0f * (n0 + n1 + n2); // TODO: The scale factor is preliminary!
}

// 3D simplex noise
t_float snoise3(t_simplex_tilde *x, t_float x_in, t_float y_in, t_float z_in) {
    t_float n0, n1, n2, n3; // Noise contributions from the four corners

    // Skew the input space to determine which simplex cell we're in
    t_float s = (x_in+y_in+z_in)*F3; // Very nice and simple skew factor for 3D
    t_float xs = x_in+s;
    t_float ys = y_in+s;
    t_float zs = z_in+s;
    int i = fastfloor(xs);
    int j = fastfloor(ys);
    int k = fastfloor(zs);

    t_float t = (float)(i+j+k)*G3; 
    t_float X0 = i-t; // Unskew the cell origin back to (x,y,z) space
    t_float Y0 = j-t;
    t_float Z0 = k-t;
    t_float x0 = x_in-X0; // The x,y,z distances from the cell origin
    t_float y0 = y_in-Y0;
    t_float z0 = z_in-Z0;

    // For the 3D case, the simplex shape is a slightly irregular tetrahedron.
    // Determine which simplex we are in.
    int c = (x0 >= y0) * 4 + (y0 >= z0) * 2 + (x0 >= z0);

    int i1 = simplex3[c][0];
    int j1 = simplex3[c][1];
    int k1 = simplex3[c][2];
    int i2 = simplex3[c][3];
    int j2 = simplex3[c][4];
    int k2 = simplex3[c][5];

    t_float x1 = x0 - i1 + G3; // Offsets for second corner in (x,y,z) coords
    t_float y1 = y0 - j1 + G3;
    t_float z1 = z0 - k1 + G3;
    t_float x2 = x0 - i2 + G3_2; // Offsets for third corner in (x,y,z) coords
    t_float y2 = y0 - j2 + G3_2;
    t_float z2 = z0 - k2 + G3_2;
    t_float x3 = x0 - 1.0f + G3_3; // Offsets for last corner in (x,y,z) coords
    t_float y3 = y0 - 1.0f + G3_3;
    t_float z3 = z0 - 1.0f + G3_3;

    // Wrap the integer indices at 256, to avoid indexing x->x_perm[] out of bounds
    int ii = i & 0xff;
    int jj = j & 0xff;
    int kk = k & 0xff;

    // Calculate the contribution from the four corners
    t_float t0 = 0.5f - x0*x0 - y0*y0 - z0*z0;
    t0 = max(t0, 0.0f);
    t0 *= t0;
    n0 = t0 * t0 * grad3(x->x_perm[ii+x->x_perm[jj+x->x_perm[kk]]], x0, y0, z0);

    t_float t1 = 0.5f - x1*x1 - y1*y1 - z1*z1;
    t1 = max(t1, 0.0f);
    t1 *= t1;
    n1 = t1 * t1 * grad3(x->x_perm[ii+i1+x->x_perm[jj+j1+x->x_perm[kk+k1]]], x1, y1, z1);

    t_float t2 = 0.5f - x2*x2 - y2*y2 - z2*z2;
    t2 = max(t2, 0.0f);
    t2 *= t2;
    n2 = t2 * t2 * grad3(x->x_perm[ii+i2+x->x_perm[jj+j2+x->x_perm[kk+k2]]], x2, y2, z2);

    t_float t3 = 0.5f - x3*x3 - y3*y3 - z3*z3;
    t3 = max(t3, 0.0f);
    t3 *= t3;
    n3 = t3 * t3 * grad3(x->x_perm[ii+1+x->x_perm[jj+1+x->x_perm[kk+1]]], x3, y3, z3);

    // Add contributions from each corner to get the final noise value.
    // The result is scaled to stay just inside [-1,1]
    return 72.0f * (n0 + n1 + n2 + n3);
}


// 4D simplex noise
t_float snoise4(t_simplex_tilde *x, t_float x_in, t_float y_in, t_float z_in, t_float w_in) {
    t_float n0, n1, n2, n3, n4; // Noise contributions from the five corners

    // Skew the (x,y,z,w) space to determine which cell of 24 simplices we're in
    t_float s = (x_in + y_in + z_in + w_in) * F4; // Factor for 4D skewing
    t_float xs = x_in + s;
    t_float ys = y_in + s;
    t_float zs = z_in + s;
    t_float ws = w_in + s;
    int i = fastfloor(xs);
    int j = fastfloor(ys);
    int k = fastfloor(zs);
    int l = fastfloor(ws);

    t_float t = (i + j + k + l) * G4; // Factor for 4D unskewing
    t_float X0 = i - t; // Unskew the cell origin back to (x,y,z,w) space
    t_float Y0 = j - t;
    t_float Z0 = k - t;
    t_float W0 = l - t;

    t_float x0 = x_in - X0;  // The x,y,z,w distances from the cell origin
    t_float y0 = y_in - Y0;
    t_float z0 = z_in - Z0;
    t_float w0 = w_in - W0;

    // For the 4D case, the simplex is a 4D shape I won't even try to describe.
    // To find out which of the 24 possible simplices we're in, we need to
    // determine the magnitude ordering of x0, y0, z0 and w0.
    // The method below is a good way of finding the ordering of x,y,z,w and
    // then find the correct traversal order for the simplex we're in.
    // First, six pair-wise comparisons are performed between each possible pair
    // of the four coordinates, and the results are used to add up binary bits
    // for an integer index.
    int c1 = (x0 > y0) << 5;
    int c2 = (x0 > z0) << 4;
    int c3 = (y0 > z0) << 3;
    int c4 = (x0 > w0) << 2;
    int c5 = (y0 > w0) << 1;
    int c6 = (z0 > w0);
    int c = c1 + c2 + c3 + c4 + c5 + c6;

    int i1, j1, k1, l1; // The integer offsets for the second simplex corner
    int i2, j2, k2, l2; // The integer offsets for the third simplex corner
    int i3, j3, k3, l3; // The integer offsets for the fourth simplex corner

    // simplex4[c] is a 4-vector with the numbers 0, 1, 2 and 3 in some order.
    // Many values of c will never occur, since e.g. x>y>z>w makes x<z, y<w and x<w
    // impossible. Only the 24 indices which have non-zero entries make any sense.
    // We use a thresholding to set the coordinates in turn from the largest magnitude.
    // The number 3 in the "simplex" array is at the position of the largest coordinate.
    i1 = simplex4[c][0] > 2;
    j1 = simplex4[c][1] > 2;
    k1 = simplex4[c][2] > 2;
    l1 = simplex4[c][3] > 2;
    // The number 2 in the "simplex" array is at the second largest coordinate.
    i2 = simplex4[c][0] > 1;
    j2 = simplex4[c][1] > 1;
    k2 = simplex4[c][2] > 1;
    l2 = simplex4[c][3] > 1;
    // The number 1 in the "simplex" array is at the second smallest coordinate.
    i3 = simplex4[c][0] > 0;
    j3 = simplex4[c][1] > 0;
    k3 = simplex4[c][2] > 0;
    l3 = simplex4[c][3] > 0;
    // The fifth corner has all coordinate offsets = 1, so no need to look that up.

    t_float x1 = x0 - i1 + G4; // Offsets for second corner in (x,y,z,w) coords
    t_float y1 = y0 - j1 + G4;
    t_float z1 = z0 - k1 + G4;
    t_float w1 = w0 - l1 + G4;
    t_float x2 = x0 - i2 + G4_2; // Offsets for third corner in (x,y,z,w) coords
    t_float y2 = y0 - j2 + G4_2;
    t_float z2 = z0 - k2 + G4_2;
    t_float w2 = w0 - l2 + G4_2;
    t_float x3 = x0 - i3 + G4_3; // Offsets for fourth corner in (x,y,z,w) coords
    t_float y3 = y0 - j3 + G4_3;
    t_float z3 = z0 - k3 + G4_3;
    t_float w3 = w0 - l3 + G4_3;
    t_float x4 = x0 - 1.0f + G4_4; // Offsets for last corner in (x,y,z,w) coords
    t_float y4 = y0 - 1.0f + G4_4;
    t_float z4 = z0 - 1.0f + G4_4;
    t_float w4 = w0 - 1.0f + G4_4;

    // Wrap the integer indices at 256, to avoid indexing x->x_perm[] out of bounds
    int ii = i & 0xff;
    int jj = j & 0xff;
    int kk = k & 0xff;
    int ll = l & 0xff;

    // Calculate the contribution from the five corners
    t_float t0 = 0.5f - x0*x0 - y0*y0 - z0*z0 - w0*w0;
    t0 = max(t0, 0.0f);
    t0 *= t0;
    n0 = t0 * t0 * grad4(x->x_perm[ii+x->x_perm[jj+x->x_perm[kk+x->x_perm[ll]]]], x0, y0, z0, w0);

    t_float t1 = 0.5f - x1*x1 - y1*y1 - z1*z1 - w1*w1;
    t1 = max(t1, 0.0f);
    t1 *= t1;
    n1 = t1 * t1 * grad4(x->x_perm[ii+i1+x->x_perm[jj+j1+x->x_perm[kk+k1+x->x_perm[ll+l1]]]], x1, y1, z1, w1);

    t_float t2 = 0.5f - x2*x2 - y2*y2 - z2*z2 - w2*w2;
    t2 = max(t2, 0.0f);
    t2 *= t2;
    n2 = t2 * t2 * grad4(x->x_perm[ii+i2+x->x_perm[jj+j2+x->x_perm[kk+k2+x->x_perm[ll+l2]]]], x2, y2, z2, w2);

    t_float t3 = 0.5f - x3*x3 - y3*y3 - z3*z3 - w3*w3;
    t3 = max(t3, 0.0f);
    t3 *= t3;
    n3 = t3 * t3 * grad4(x->x_perm[ii+i3+x->x_perm[jj+j3+x->x_perm[kk+k3+x->x_perm[ll+l3]]]], x3, y3, z3, w3);

    t_float t4 = 0.5f - x4*x4 - y4*y4 - z4*z4 - w4*w4;
    t4 = max(t4, 0.0f);
    t4 *= t4;
    n4 = t4 * t4 * grad4(x->x_perm[ii+1+x->x_perm[jj+1+x->x_perm[kk+1+x->x_perm[ll+1]]]], x4, y4, z4, w4);

    // Sum up and scale the result to cover the range [-1,1]
    return 62.0f * (n0 + n1 + n2 + n3 + n4);
}

static t_float generate_noise(t_simplex_tilde *x, t_float *pos, t_float persistence, int channel_count) {
    t_float result = 0, persistence_sum = 0, scale, coeff;

    for (int octave = 0; octave < x->x_octaves; octave++) {
        scale = pow(2, octave);
        coeff = pow(persistence, octave);
        persistence_sum += fabs(coeff);

        switch (channel_count) {
            case 1:
                result += snoise1(x, pos[0] * scale) * coeff;
                break;
            case 2:
                result += snoise2(x, pos[0] * scale, pos[1] * scale) * coeff;
                break;
            case 3:
                result += snoise3(x, pos[0] * scale, pos[1] * scale, pos[2] * scale) * coeff;
                break;
            case 4:
                result += snoise4(x, pos[0] * scale, pos[1] * scale, pos[2] * scale, pos[3] * scale) * coeff;
                break;
        }
    }
    return x->x_normalize ? result / persistence_sum : result;
}

static t_int *simplex_tilde_perform(t_int *w) {
    t_simplex_tilde *x = (t_simplex_tilde *)(w[1]);
    t_int nsamples = (t_int)(w[2]);
    t_int nchans = (t_int)(w[3]);
    t_sample *in_coord = (t_sample *)(w[4]);
    t_sample *in_persistence = (t_sample *)(w[5]);
    t_sample *out = (t_sample *)(w[6]);

    for (int i = 0; i < nsamples; i++) {
        t_float pos[4] = {0};
        for (int channel = 0; channel < nchans; channel++)
            pos[channel] = in_coord[nsamples * channel + i];
        *out++ = generate_noise(x, pos, in_persistence[i], nchans);
    }
    return w+7;
}

void simplex_tilde_dsp(t_simplex_tilde *x, t_signal **sp) {
    signal_setmultiout(&sp[2], 1);
    dsp_add(simplex_tilde_perform, 6, x, (t_int)sp[0]->s_n, (t_int)sp[0]->s_nchans, sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec);
}

static void simplex_tilde_octaves(t_simplex_tilde *x, t_floatarg f){
    x->x_octaves = fastfloor(f);
}

static void simplex_tilde_normalize(t_simplex_tilde *x, t_floatarg f){
    x->x_normalize = fastfloor(f);
}

static void simplex_tilde_seed(t_simplex_tilde *x, t_floatarg f){
    initPermutation(x, fastfloor(f));
}

void *simplex_tilde_new(t_symbol *s, int ac, t_atom *av) {
    t_float persistence;

    t_simplex_tilde *x = (t_simplex_tilde *)pd_new(simplex_tilde_class);
    x->x_normalize = 0;
    if (ac && av->a_type == A_SYMBOL) {
        if (atom_getsymbol(av) == gensym("-n"))
            x->x_normalize = 1;
        else
            pd_error(x, "[simplex~]: invalid argument");
        ac--, av++;
    }
    x->x_octaves = ac-- ? atom_getint(av++) : 1;
    persistence = ac ? atom_getfloat(av) : PERSISTENCE;
    x->x_inlet_persistence = inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal);
        pd_float((t_pd *)x->x_inlet_persistence, persistence);
    outlet_new(&x->x_obj, &s_signal);

    initPermutation(x, time(NULL));
    return x;
}

void *simplex_tilde_free(t_simplex_tilde *x){
    inlet_free(x->x_inlet_persistence);
    return (void *)x;
}

void simplex_tilde_setup(void) {
    simplex_tilde_class = class_new(gensym("simplex~"), (t_newmethod)simplex_tilde_new,
        (t_method)simplex_tilde_free, sizeof(t_simplex_tilde), CLASS_MULTICHANNEL, A_GIMME, 0);
    class_addmethod(simplex_tilde_class, nullfn, gensym("signal"), 0);
    class_addmethod(simplex_tilde_class, (t_method)simplex_tilde_dsp, gensym("dsp"), 0);
    class_addmethod(simplex_tilde_class, (t_method)simplex_tilde_octaves, gensym("octaves"), A_FLOAT, 0);
    class_addmethod(simplex_tilde_class, (t_method)simplex_tilde_seed, gensym("seed"), A_FLOAT, 0);
    class_addmethod(simplex_tilde_class, (t_method)simplex_tilde_normalize, gensym("normalize"), A_FLOAT, 0);
}
