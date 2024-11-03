// simplex algorithms by Stefan Gustavson (slightly adapted) from https://github.com/stegu/perlin-noise/blob/master/src/sdnoise1234.c

#ifndef SIMPLEX_COMMON_H
#define SIMPLEX_COMMON_H

#include "m_pd.h"
#include <math.h>
#include <time.h>
#include <stdlib.h>

// Constants
#define DEFAULT_PERSISTENCE 0.5
#define MAX_DIMENSIONS 4
#define MAX_OCTAVES 1024

// Skew factors
#define F2   0.36602540378 // (sqrt(3) - 1) / 2
#define G2   0.2113248654  // (3 - sqrt(3)) / 6
#define G2_2 0.42264973081

#define F3   0.33333333333 // 1 / 2
#define G3   0.16666666666 // 1 / 6
#define G3_2 0.33333333333
#define G3_3 0.5

#define F4   0.30901699437 // (sqrt(5) - 1) / 4
#define G4   0.13819660112 // (5 - sqrt(5)) / 20
#define G4_2 0.27639320225
#define G4_3 0.41458980337
#define G4_4 0.5527864045

// Utility macros
#define max(a,b) ( ((a) > (b)) ? (a) : (b) )
#define min(a,b) ( ((a) < (b)) ? (a) : (b) )
#define clamp(a,b,c) (min(max((a), (b)), (c)))
#define fastfloor(x) ( ((t_int)(x)<=(x)) ? ((t_int)x) : (((t_int)x)-1) )

// Generic config struct used by both externals
typedef struct _simplex_config {
    t_float *octave_factors;
    int octaves;
    int normalize;
    unsigned char *perm;
} t_simplex_config;

// Lookup tables
static t_float grad2lut[8][2] = {
    {-1,-1}, { 1, 0}, {-1, 0}, { 1, 1},
    {-1, 1}, { 0,-1}, { 0, 1}, { 1,-1}
};

/*
 * Gradient directions for 3D.
 * These vectors are based on the midpoints of the 12 edges of a cube.
 * A larger array of random unit length vectors would also do the job,
 * but these 12 (including 4 repeats to make the array length a power
 * of two) work better. They are not random, they are carefully chosen
 * to represent a small, isotropic set of directions.
 */

static t_float grad3lut[16][3] = {
    { 1, 0, 1}, { 0, 1, 1}, // 12 cube edges
    {-1, 0, 1}, { 0,-1, 1},
    { 1, 0,-1}, { 0, 1,-1},
    {-1, 0,-1}, { 0,-1,-1},
    { 1,-1, 0}, { 1, 1, 0},
    {-1, 1, 0}, {-1,-1, 0},
    { 1, 0, 1}, {-1, 0, 1}, // 4 repeats to make 16
    { 0, 1,-1}, { 0,-1,-1}
};

static t_float grad4lut[32][4] = {
    { 0, 1, 1, 1}, { 0, 1, 1,-1}, { 0, 1,-1, 1}, { 0, 1,-1,-1},
    { 0,-1, 1, 1}, { 0,-1, 1,-1}, { 0,-1,-1, 1}, { 0,-1,-1,-1},
    { 1, 0, 1, 1}, { 1, 0, 1,-1}, { 1, 0,-1, 1}, { 1, 0,-1,-1},
    {-1, 0, 1, 1}, {-1, 0, 1,-1}, {-1, 0,-1, 1}, {-1, 0,-1,-1},
    { 1, 1, 0, 1}, { 1, 1, 0,-1}, { 1,-1, 0, 1}, { 1,-1, 0,-1},
    {-1, 1, 0, 1}, {-1, 1, 0,-1}, {-1,-1, 0, 1}, {-1,-1, 0,-1},
    { 1, 1, 1, 0}, { 1, 1,-1, 0}, { 1,-1, 1, 0}, { 1,-1,-1, 0},
    {-1, 1, 1, 0}, {-1, 1,-1, 0}, {-1,-1, 1, 0}, {-1,-1,-1, 0}
};

static unsigned char simplex[64][4] = {
    {0,1,2,3}, {0,1,3,2}, {0,0,0,0}, {0,2,3,1}, {0,0,0,0}, {0,0,0,0}, {0,0,0,0}, {1,2,3,0},
    {0,2,1,3}, {0,0,0,0}, {0,3,1,2}, {0,3,2,1}, {0,0,0,0}, {0,0,0,0}, {0,0,0,0}, {1,3,2,0},
    {0,0,0,0}, {0,0,0,0}, {0,0,0,0}, {0,0,0,0}, {0,0,0,0}, {0,0,0,0}, {0,0,0,0}, {0,0,0,0},
    {1,2,0,3}, {0,0,0,0}, {1,3,0,2}, {0,0,0,0}, {0,0,0,0}, {0,0,0,0}, {2,3,0,1}, {2,3,1,0},
    {1,0,2,3}, {1,0,3,2}, {0,0,0,0}, {0,0,0,0}, {0,0,0,0}, {2,0,3,1}, {0,0,0,0}, {2,1,3,0},
    {0,0,0,0}, {0,0,0,0}, {0,0,0,0}, {0,0,0,0}, {0,0,0,0}, {0,0,0,0}, {0,0,0,0}, {0,0,0,0},
    {2,0,1,3}, {0,0,0,0}, {0,0,0,0}, {0,0,0,0}, {3,0,1,2}, {3,0,2,1}, {0,0,0,0}, {3,1,2,0},
    {2,1,0,3}, {0,0,0,0}, {0,0,0,0}, {0,0,0,0}, {3,1,0,2}, {0,0,0,0}, {3,2,0,1}, {3,2,1,0}
};

// Common initialization functions
static void init_permutation_with_seed(unsigned char *perm, unsigned int seed) {
    int i;
    unsigned char basePermutation[256];
    for (i = 0; i < 256; i++) {
        basePermutation[i] = i;
    }
    srand(seed);
    for (i = 255; i > 0; i--) {
        int j = rand() % (i + 1);
        unsigned char temp = basePermutation[i];
        basePermutation[i] = basePermutation[j];
        basePermutation[j] = temp;
    }
    for (i = 0; i < 256; i++) {
        perm[i] = basePermutation[i];
        perm[i + 256] = basePermutation[i];
    }
}

static void init_permutation(unsigned char *perm) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    init_permutation_with_seed(perm, (unsigned int)ts.tv_nsec);
}

static inline void init_octave_factors(t_float *octave_factors, int octaves) {
    for (int octave = 0; octave < octaves; octave++)
        octave_factors[octave] = (t_float)(1 << octave);
}

static void grad1(t_int hash, t_float *gx) {
    t_int h = hash & 15;
    *gx = 1 + (h & 7); 
    if (h & 8) *gx = -(*gx); 
}

static void grad2(t_int hash, t_float *gx, t_float *gy) {
    t_int h = hash & 7;
    *gx = grad2lut[h][0];
    *gy = grad2lut[h][1];
}

static void grad3(t_int hash, t_float *gx, t_float *gy, t_float *gz) {
    t_int h = hash & 15;
    *gx = grad3lut[h][0];
    *gy = grad3lut[h][1];
    *gz = grad3lut[h][2];
}

static void grad4(t_int hash, t_float *gx, t_float *gy, t_float *gz, t_float *gw) {
    t_int h = hash & 31;
    *gx = grad4lut[h][0];
    *gy = grad4lut[h][1];
    *gz = grad4lut[h][2];
    *gw = grad4lut[h][3];
}

// 1D simplex noise
static t_float snoise1(t_float *pos, t_float sc, t_float coeff, unsigned char *perm, t_float *derivatives) {
    t_float x = sc * pos[0];

    t_int i0 = fastfloor(x);
    t_int i1 = i0 + 1;
    t_float x0 = x - i0;
    t_float x1 = x0 - 1;
    t_float gx0, gx1;
    t_float n0, n1;
    t_float t1, t20, t40, t21, t41, x21;
    t_float x20 = x0 * x0;
    t_float t0 = 1 - x20;
    t20 = t0 * t0;
    t40 = t20 * t20;
    grad1(perm[i0 & 0xff], &gx0);
    n0 = t40 * gx0 * x0;
    x21 = x1 * x1;
    t1 = 1 - x21;
    t21 = t1 * t1;
    t41 = t21 * t21;
    grad1(perm[i1 & 0xff], &gx1);
    n1 = t41 * gx1 * x1;
    if (derivatives) {
        t_float d;
        d = t20 * t0 * gx0 * x20;
        d += t21 * t1 * gx1 * x21;
        d *= -8;
        d += t40 * gx0 + t41 * gx1;
        d *= 0.25;
        derivatives[0] += coeff * d;
    }
    return coeff * 0.25 * (n0 + n1);
}

// 2D simplex noise
static t_float snoise2(t_float *pos, t_float sc, t_float coeff, unsigned char *perm, t_float *derivatives) {
    t_float x = sc * pos[0], y = sc * pos[1];

    t_float n0, n1, n2;
    t_float gx0, gy0, gx1, gy1, gx2, gy2;
    t_float t0, t1, t2, x1, x2, y1, y2;
    t_float t20, t40, t21, t41, t22, t42;
    t_float temp0, temp1, temp2, noise;
    t_float s = (x + y) * F2;
    t_float xs = x + s;
    t_float ys = y + s;
    t_int ii, i = fastfloor(xs);
    t_int jj, j = fastfloor(ys);
    t_float t = (t_float)(i + j) * G2;
    t_float X0 = i - t;
    t_float Y0 = j - t;
    t_float x0 = x - X0;
    t_float y0 = y - Y0;
    t_int i1, j1;
    if (x0 > y0) { i1 = 1; j1 = 0; }
    else { i1 = 0; j1 = 1; }
    x1 = x0 - i1 + G2;
    y1 = y0 - j1 + G2;
    x2 = x0 - 1 + G2_2;
    y2 = y0 - 1 + G2_2;
    ii = i & 0xff;
    jj = j & 0xff;
    t0 = 0.5 - x0 * x0 - y0 * y0;
    if (t0 < 0) t40 = t20 = t0 = n0 = gx0 = gy0 = 0;
    else {
        grad2(perm[ii + perm[jj]], &gx0, &gy0);
        t20 = t0 * t0;
        t40 = t20 * t20;
        n0 = t40 * (gx0 * x0 + gy0 * y0);
    }
    t1 = 0.5 - x1 * x1 - y1 * y1;
    if (t1 < 0) t21 = t41 = t1 = n1 = gx1 = gy1 = 0;
    else {
        grad2(perm[ii + i1 + perm[jj + j1]], &gx1, &gy1);
        t21 = t1 * t1;
        t41 = t21 * t21;
        n1 = t41 * (gx1 * x1 + gy1 * y1);
    }
    t2 = 0.5 - x2 * x2 - y2 * y2;
    if (t2 < 0) t42 = t22 = t2 = n2 = gx2 = gy2 = 0;
    else {
        grad2(perm[ii + 1 + perm[jj + 1]], &gx2, &gy2);
        t22 = t2 * t2;
        t42 = t22 * t22;
        n2 = t42 * (gx2 * x2 + gy2 * y2);
    }
    noise = 40 * (n0 + n1 + n2);
    if (derivatives) {
        t_float d[2];
        temp0 = t20 * t0 * (gx0 * x0 + gy0 * y0);
        d[0] = temp0 * x0;
        d[1] = temp0 * y0;
        temp1 = t21 * t1 * (gx1 * x1 + gy1 * y1);
        d[0] += temp1 * x1;
        d[1] += temp1 * y1;
        temp2 = t22 * t2 * (gx2 * x2 + gy2 * y2);
        d[0] += temp2 * x2;
        d[1] += temp2 * y2;
        d[0] *= -8;
        d[1] *= -8;
        d[0] += t40 * gx0 + t41 * gx1 + t42 * gx2;
        d[1] += t40 * gy0 + t41 * gy1 + t42 * gy2;
        d[0] *= 40;
        d[1] *= 40;
        derivatives[0] += coeff * d[0];
        derivatives[1] += coeff * d[1];
    }
    return coeff * noise;
}

// 3D simplex noise
static t_float snoise3(t_float *pos, t_float sc, t_float coeff, unsigned char *perm, t_float *derivatives) {
    t_float x = sc * pos[0], y = sc * pos[1], z = sc * pos[2];

    t_float n0, n1, n2, n3;
    t_float noise;
    t_float gx0, gy0, gz0, gx1, gy1, gz1;
    t_float gx2, gy2, gz2, gx3, gy3, gz3;
    t_float x1, y1, z1, x2, y2, z2, x3, y3, z3;
    t_float t0, t1, t2, t3, t20, t40, t21, t41, t22, t42, t23, t43;
    t_float temp0, temp1, temp2, temp3;
    t_float s = (x + y + z) * F3;
    t_float xs = x + s;
    t_float ys = y + s;
    t_float zs = z + s;
    t_int ii, i = fastfloor(xs);
    t_int jj, j = fastfloor(ys);
    t_int kk, k = fastfloor(zs);
    t_float t = (t_float)(i + j + k) * G3; 
    t_float X0 = i - t;
    t_float Y0 = j - t;
    t_float Z0 = k - t;
    t_float x0 = x - X0;
    t_float y0 = y - Y0;
    t_float z0 = z - Z0;
    t_int i1, j1, k1;
    t_int i2, j2, k2;
    if (x0 >= y0) {
        if (y0 >= z0) {
            i1 = 1; j1 = 0; k1 = 0; i2 = 1; j2 = 1; k2 = 0;
        } else if (x0 >= z0) {
            i1 = 1; j1 = 0; k1 = 0; i2 = 1; j2 = 0; k2 = 1;
        } else {
            i1 = 0; j1 = 0; k1 = 1; i2 = 1; j2 = 0; k2 = 1;
        }
    } else {
        if (y0 < z0) {
            i1 = 0; j1 = 0; k1 = 1; i2 = 0; j2 = 1; k2 = 1;
        } else if (x0 < z0) {
            i1 = 0; j1 = 1; k1 = 0; i2 = 0; j2 = 1; k2 = 1;
        } else {
            i1 = 0; j1 = 1; k1 = 0; i2 = 1; j2 = 1; k2 = 0;
        }
    }
    x1 = x0 - i1 + G3;
    y1 = y0 - j1 + G3;
    z1 = z0 - k1 + G3;
    x2 = x0 - i2 + G3_2;
    y2 = y0 - j2 + G3_2;
    z2 = z0 - k2 + G3_2;
    x3 = x0 - 1 + G3_3;
    y3 = y0 - 1 + G3_3;
    z3 = z0 - 1 + G3_3;
    ii = i & 0xff;
    jj = j & 0xff;
    kk = k & 0xff;
    t0 = 0.5 - x0 * x0 - y0 * y0 - z0 * z0;
    if (t0 < 0) n0 = t0 = t20 = t40 = gx0 = gy0 = gz0 = 0;
    else {
        grad3(perm[ii + perm[jj + perm[kk]]], &gx0, &gy0, &gz0);
        t20 = t0 * t0;
        t40 = t20 * t20;
        n0 = t40 * (gx0 * x0 + gy0 * y0 + gz0 * z0);
    }
    t1 = 0.5 - x1 * x1 - y1 * y1 - z1 * z1;
    if (t1 < 0) n1 = t1 = t21 = t41 = gx1 = gy1 = gz1 = 0;
    else {
        grad3(perm[ii + i1 + perm[jj + j1 + perm[kk + k1]]], &gx1, &gy1, &gz1);
        t21 = t1 * t1;
        t41 = t21 * t21;
        n1 = t41 * (gx1 * x1 + gy1 * y1 + gz1 * z1);
    }
    t2 = 0.5 - x2 * x2 - y2 * y2 - z2 * z2;
    if (t2 < 0) n2 = t2 = t22 = t42 = gx2 = gy2 = gz2 = 0;
    else {
        grad3(perm[ii + i2 + perm[jj + j2 + perm[kk + k2]]], &gx2, &gy2, &gz2);
        t22 = t2 * t2;
        t42 = t22 * t22;
        n2 = t42 * (gx2 * x2 + gy2 * y2 + gz2 * z2);
    }
    t3 = 0.5 - x3 * x3 - y3 * y3 - z3 * z3;
    if (t3 < 0) n3 = t3 = t23 = t43 = gx3 = gy3 = gz3 = 0;
    else {
        grad3(perm[ii + 1 + perm[jj + 1 + perm[kk + 1]]], &gx3, &gy3, &gz3);
        t23 = t3 * t3;
        t43 = t23 * t23;
        n3 = t43 * (gx3 * x3 + gy3 * y3 + gz3 * z3);
    }
    noise = 72 * (n0 + n1 + n2 + n3);
    if (derivatives) {
        t_float d[3];
        temp0 = t20 * t0 * (gx0 * x0 + gy0 * y0 + gz0 * z0);
        d[0] = temp0 * x0;
        d[1] = temp0 * y0;
        d[2] = temp0 * z0;
        temp1 = t21 * t1 * (gx1 * x1 + gy1 * y1 + gz1 * z1);
        d[0] += temp1 * x1;
        d[1] += temp1 * y1;
        d[2] += temp1 * z1;
        temp2 = t22 * t2 * (gx2 * x2 + gy2 * y2 + gz2 * z2);
        d[0] += temp2 * x2;
        d[1] += temp2 * y2;
        d[2] += temp2 * z2;
        temp3 = t23 * t3 * (gx3 * x3 + gy3 * y3 + gz3 * z3);
        d[0] += temp3 * x3;
        d[1] += temp3 * y3;
        d[2] += temp3 * z3;
        d[0] *= -8;
        d[1] *= -8;
        d[2] *= -8;
        d[0] += t40 * gx0 + t41 * gx1 + t42 * gx2 + t43 * gx3;
        d[1] += t40 * gy0 + t41 * gy1 + t42 * gy2 + t43 * gy3;
        d[2] += t40 * gz0 + t41 * gz1 + t42 * gz2 + t43 * gz3;
        d[0] *= 72;
        d[1] *= 72;
        d[2] *= 72;
        derivatives[0] += coeff * d[0];
        derivatives[1] += coeff * d[1];
        derivatives[2] += coeff * d[2];
    }
    return coeff * noise;
}

// 4D simplex noise
static t_float snoise4(t_float *pos, t_float sc, t_float coeff, unsigned char *perm, t_float *derivatives) {
    t_float x = sc * pos[0], y = sc * pos[1], z = sc * pos[2], w = sc * pos[3];

    t_float n0, n1, n2, n3, n4;
    t_float noise;
    t_float gx0, gy0, gz0, gw0, gx1, gy1, gz1, gw1;
    t_float gx2, gy2, gz2, gw2, gx3, gy3, gz3, gw3, gx4, gy4, gz4, gw4;
    t_float t20, t21, t22, t23, t24;
    t_float t40, t41, t42, t43, t44;
    t_float x1, y1, z1, w1, x2, y2, z2, w2, x3, y3, z3, w3, x4, y4, z4, w4;
    t_float t0, t1, t2, t3, t4;
    t_float temp0, temp1, temp2, temp3, temp4;
    t_float s = (x + y + z + w) * F4;
    t_float xs = x + s;
    t_float ys = y + s;
    t_float zs = z + s;
    t_float ws = w + s;
    t_int ii, i = fastfloor(xs);
    t_int jj, j = fastfloor(ys);
    t_int kk, k = fastfloor(zs);
    t_int ll, l = fastfloor(ws);
    t_float t = (i + j + k + l) * G4;
    t_float X0 = i - t;
    t_float Y0 = j - t;
    t_float Z0 = k - t;
    t_float W0 = l - t;
    t_float x0 = x - X0;
    t_float y0 = y - Y0;
    t_float z0 = z - Z0;
    t_float w0 = w - W0;
    t_int c1 = (x0 > y0) << 5;
    t_int c2 = (x0 > z0) << 4;
    t_int c3 = (y0 > z0) << 3;
    t_int c4 = (x0 > w0) << 2;
    t_int c5 = (y0 > w0) << 1;
    t_int c6 = (z0 > w0);
    t_int c = c1 | c2 | c3 | c4 | c5 | c6;
    t_int i1, j1, k1, l1;
    t_int i2, j2, k2, l2;
    t_int i3, j3, k3, l3;
    i1 = simplex[c][0] > 2;
    j1 = simplex[c][1] > 2;
    k1 = simplex[c][2] > 2;
    l1 = simplex[c][3] > 2;
    i2 = simplex[c][0] > 1;
    j2 = simplex[c][1] > 1;
    k2 = simplex[c][2] > 1;
    l2 = simplex[c][3] > 1;
    i3 = simplex[c][0] > 0;
    j3 = simplex[c][1] > 0;
    k3 = simplex[c][2] > 0;
    l3 = simplex[c][3] > 0;
    x1 = x0 - i1 + G4;
    y1 = y0 - j1 + G4;
    z1 = z0 - k1 + G4;
    w1 = w0 - l1 + G4;
    x2 = x0 - i2 + G4_2;
    y2 = y0 - j2 + G4_2;
    z2 = z0 - k2 + G4_2;
    w2 = w0 - l2 + G4_2;
    x3 = x0 - i3 + G4_3;
    y3 = y0 - j3 + G4_3;
    z3 = z0 - k3 + G4_3;
    w3 = w0 - l3 + G4_3;
    x4 = x0 - 1 + G4_4;
    y4 = y0 - 1 + G4_4;
    z4 = z0 - 1 + G4_4;
    w4 = w0 - 1 + G4_4;
    ii = i & 0xff;
    jj = j & 0xff;
    kk = k & 0xff;
    ll = l & 0xff;
    t0 = 0.5 - x0 * x0 - y0 * y0 - z0 * z0 - w0 * w0;
    if (t0 < 0) n0 = t0 = t20 = t40 = gx0 = gy0 = gz0 = gw0 = 0;
    else {
        t20 = t0 * t0;
        t40 = t20 * t20;
        grad4(perm[ii + perm[jj + perm[kk + perm[ll]]]], &gx0, &gy0, &gz0, &gw0);
        n0 = t40 * (gx0 * x0 + gy0 * y0 + gz0 * z0 + gw0 * w0);
    }
    t1 = 0.5 - x1 * x1 - y1 * y1 - z1 * z1 - w1 * w1;
    if (t1 < 0) n1 = t1 = t21 = t41 = gx1 = gy1 = gz1 = gw1 = 0;
    else {
        t21 = t1 * t1;
        t41 = t21 * t21;
        grad4(perm[ii + i1 + perm[jj + j1 + perm[kk + k1 + perm[ll + l1]]]], &gx1, &gy1, &gz1, &gw1);
        n1 = t41 * (gx1 * x1 + gy1 * y1 + gz1 * z1 + gw1 * w1);
    }
    t2 = 0.5 - x2 * x2 - y2 * y2 - z2 * z2 - w2 * w2;
    if (t2 < 0) n2 = t2 = t22 = t42 = gx2 = gy2 = gz2 = gw2 = 0;
    else {
        t22 = t2 * t2;
        t42 = t22 * t22;
        grad4(perm[ii + i2 + perm[jj + j2 + perm[kk + k2 + perm[ll + l2]]]], &gx2, &gy2, &gz2, &gw2);
        n2 = t42 * (gx2 * x2 + gy2 * y2 + gz2 * z2 + gw2 * w2);
    }
    t3 = 0.5 - x3 * x3 - y3 * y3 - z3 * z3 - w3 * w3;
    if (t3 < 0) n3 = t3 = t23 = t43 = gx3 = gy3 = gz3 = gw3 = 0;
    else {
        t23 = t3 * t3;
        t43 = t23 * t23;
        grad4(perm[ii + i3 + perm[jj + j3 + perm[kk + k3 + perm[ll + l3]]]], &gx3, &gy3, &gz3, &gw3);
        n3 = t43 * (gx3 * x3 + gy3 * y3 + gz3 * z3 + gw3 * w3);
    }
    t4 = 0.5 - x4 * x4 - y4 * y4 - z4 * z4 - w4 * w4;
    if (t4 < 0) n4 = t4 = t24 = t44 = gx4 = gy4 = gz4 = gw4 = 0;
    else {
        t24 = t4 * t4;
        t44 = t24 * t24;
        grad4(perm[ii + 1 + perm[jj + 1 + perm[kk + 1 + perm[ll + 1]]]], &gx4, &gy4, &gz4, &gw4);
        n4 = t44 * (gx4 * x4 + gy4 * y4 + gz4 * z4 + gw4 * w4);
    }
    noise = 62 * (n0 + n1 + n2 + n3 + n4);
    if (derivatives) {
        t_float d[4];
        temp0 = t20 * t0 * (gx0 * x0 + gy0 * y0 + gz0 * z0 + gw0 * w0);
        d[0] = temp0 * x0;
        d[1] = temp0 * y0;
        d[2] = temp0 * z0;
        d[3] = temp0 * w0;
        temp1 = t21 * t1 * (gx1 * x1 + gy1 * y1 + gz1 * z1 + gw1 * w1);
        d[0] += temp1 * x1;
        d[1] += temp1 * y1;
        d[2] += temp1 * z1;
        d[3] += temp1 * w1;
        temp2 = t22 * t2 * (gx2 * x2 + gy2 * y2 + gz2 * z2 + gw2 * w2);
        d[0] += temp2 * x2;
        d[1] += temp2 * y2;
        d[2] += temp2 * z2;
        d[3] += temp2 * w2;
        temp3 = t23 * t3 * (gx3 * x3 + gy3 * y3 + gz3 * z3 + gw3 * w3);
        d[0] += temp3 * x3;
        d[1] += temp3 * y3;
        d[2] += temp3 * z3;
        d[3] += temp3 * w3;
        temp4 = t24 * t4 * (gx4 * x4 + gy4 * y4 + gz4 * z4 + gw4 * w4);
        d[0] += temp4 * x4;
        d[1] += temp4 * y4;
        d[2] += temp4 * z4;
        d[3] += temp4 * w4;
        d[0] *= -8;
        d[1] *= -8;
        d[2] *= -8;
        d[3] *= -8;
        d[0] += t40 * gx0 + t41 * gx1 + t42 * gx2 + t43 * gx3 + t44 * gx4;
        d[1] += t40 * gy0 + t41 * gy1 + t42 * gy2 + t43 * gy3 + t44 * gy4;
        d[2] += t40 * gz0 + t41 * gz1 + t42 * gz2 + t43 * gz3 + t44 * gz4;
        d[3] += t40 * gw0 + t41 * gw1 + t42 * gw2 + t43 * gw3 + t44 * gw4;
        d[0] *= coeff * 62;
        d[1] *= coeff * 62;
        d[2] *= coeff * 62;
        d[3] *= coeff * 62;
        derivatives[0] += coeff * d[0];
        derivatives[1] += coeff * d[1];
        derivatives[2] += coeff * d[2];
        derivatives[3] += coeff * d[3];
    }
    return coeff * noise;
}

static inline t_float generate_noise(t_simplex_config *cfg,
									 t_float *pos, 
                                     t_float persistence,
                                     int func_index, 
                                     t_float *derivatives) {
    t_float result = 0;
    t_float coeff = 1;
    t_float normalize_factor = 1;
    t_float scale;
    t_float abs_persistence = fabs(persistence);
    int i;

    if (cfg->normalize) {
        if (abs_persistence == 1) {
            normalize_factor = 1.0f / cfg->octaves;
        } else {
            normalize_factor = abs_persistence - 1;
            normalize_factor /= pow(abs_persistence, cfg->octaves) - 1;
        }
    }
    
    static t_float (*noise_func[])(t_float *, t_float, t_float, unsigned char *, t_float *) = {
        snoise1, snoise2, snoise3, snoise4
    };
    
    if (derivatives) {
        for (i = 0; i < MAX_DIMENSIONS; i++)
            derivatives[i] = 0;
    }
    
    for (i = 0; i < cfg->octaves; i++) {
        if (i) coeff *= persistence;
        scale = cfg->octave_factors[i];
        result += noise_func[func_index](pos, scale, coeff, cfg->perm, derivatives);
    }
    
    return result * normalize_factor;
}

#endif // SIMPLEX_COMMON_H