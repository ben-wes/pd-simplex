#include "m_pd.h"
#include <math.h>

#define max(a,b) ( ((a) > (b)) ? (a) : (b) )

#define PERSISTENCE 0.5

#define F2 0.36602540378 // 0.5*(sqrt(3.0)-1.0);
#define G2 0.2113248654  // (3.0-sqrt(3.0))/6.0;
#define F3 0.33333333333 // 1.0/3.0;
#define G3 0.16666666666 // 1.0/6.0
#define F4 0.30901699437 // (sqrt(5.0)-1.0)/4.0;
#define G4 0.13819660112 // (5.0-sqrt(5.0))/20.0;

static t_class *simplex_tilde_class;

typedef struct _simplex_tilde {
    t_object x_obj;
    t_inlet *x_inlet_persistence;
    int x_octaves;
    int x_perm[512];
} t_simplex_tilde;

// static inline uint8_t hash(t_simplex_tilde *x, int i) {
//     return x->x_p[(uint8_t)i]; // Assuming x->p is your permutation array within t_simplex_tilde
// }

// static inline t_float grad(int hash, t_float x, t_float y, t_float z) {
//     int h = hash & 15;
//     t_float u = h < 8 ? x : y;
//     t_float v = h < 4 ? y : h == 12 || h == 14 ? x : z;
//     return ((h & 1) ? -u : u) + ((h & 2) == 0 ? v : -v);
// }

// // 3D simplex noise function
// t_float simplex(t_simplex_tilde *x, t_float xin, t_float yin, t_float zin) {
//     t_float n0, n1, n2, n3;  // Noise contributions from the four simplex corners

//     // Skewing and unskewing factors for 3 dimensions
//     static t_float F3 = 1.0 / 3.0;
//     static t_float G3 = 1.0 / 6.0;

//     // Skew the input space to determine which simplex cell we're in
//     t_float s = (xin + yin + zin) * F3;  // Very nice and simple skew factor for 3D
//     int i = fastfloor(xin + s);
//     int j = fastfloor(yin + s);
//     int k = fastfloor(zin + s);
//     t_float t = (i + j + k) * G3;
//     t_float X0 = i - t;  // Unskew the cell origin back to (x, y, z) space
//     t_float Y0 = j - t;
//     t_float Z0 = k - t;
//     t_float x0 = xin - X0;  // The x, y, z distances from the cell origin
//     t_float y0 = yin - Y0;
//     t_float z0 = zin - Z0;
//     // Determine which simplex we are in.
//     int i1, j1, k1;  // Offsets for second corner of simplex in (i, j, k) coords
//     int i2, j2, k2;  // Offsets for third corner of simplex in (i, j, k) coords
//     if (x0 >= y0) {
//         if (y0 >= z0) {
//             i1 = 1; j1 = 0; k1 = 0; i2 = 1; j2 = 1; k2 = 0;  // X Y Z order
//         } else if (x0 >= z0) {
//             i1 = 1; j1 = 0; k1 = 0; i2 = 1; j2 = 0; k2 = 1;  // X Z Y order
//         } else {
//             i1 = 0; j1 = 0; k1 = 1; i2 = 1; j2 = 0; k2 = 1;  // Z X Y order
//         }
//     } else {  // x0 < y0
//         if (y0 < z0) {
//             i1 = 0; j1 = 0; k1 = 1; i2 = 0; j2 = 1; k2 = 1;  // Z Y X order
//         } else if (x0 < z0) {
//             i1 = 0; j1 = 1; k1 = 0; i2 = 0; j2 = 1; k2 = 1;  // Y Z X order
//         } else {
//             i1 = 0; j1 = 1; k1 = 0; i2 = 1; j2 = 1; k2 = 0;  // Y X Z order
//         }
//     }
//     // A step of (1,0,0) in (i,j,k) means a step of (1-c,-c,-c) in (x,y,z),
//     // a step of (0,1,0) in (i,j,k) means a step of (-c,1-c,-c) in (x,y,z), and
//     // a step of (0,0,1) in (i,j,k) means a step of (-c,-c,1-c) in (x,y,z), where
//     // c = 1/6.

//     t_float x1 = x0 - i1 + G3; // Offsets for second corner in (x,y,z) coords
//     t_float y1 = y0 - j1 + G3;
//     t_float z1 = z0 - k1 + G3;
//     t_float x2 = x0 - i2 + 2.0 * G3; // Offsets for third corner in (x,y,z) coords
//     t_float y2 = y0 - j2 + 2.0 * G3;
//     t_float z2 = z0 - k2 + 2.0 * G3;
//     t_float x3 = x0 - 1.0 + 3.0 * G3; // Offsets for last corner in (x,y,z) coords
//     t_float y3 = y0 - 1.0 + 3.0 * G3;
//     t_float z3 = z0 - 1.0 + 3.0 * G3;

//     int gi0 = hash(x, i + hash(x, j + hash(x, k)));
//     int gi1 = hash(x, i + i1 + hash(x, j + j1 + hash(x, k + k1)));
//     int gi2 = hash(x, i + i2 + hash(x, j + j2 + hash(x, k + k2)));
//     int gi3 = hash(x, i + 1 + hash(x, j + 1 + hash(x, k + 1)));

//     // Calculate the contribution from the four corners
//     t_float t0 = 0.6 - x0*x0 - y0*y0 - z0*z0;
//     if (t0 < 0) {
//         n0 = 0.0;
//     } else {
//         t0 *= t0;
//         n0 = t0 * t0 * grad(gi0, x0, y0, z0);
//     }
//     t_float t1 = 0.6 - x1*x1 - y1*y1 - z1*z1;
//     if (t1 < 0) {
//         n1 = 0.0;
//     } else {
//         t1 *= t1;
//         n1 = t1 * t1 * grad(gi1, x1, y1, z1);
//     }
//     t_float t2 = 0.6 - x2*x2 - y2*y2 - z2*z2;
//     if (t2 < 0) {
//         n2 = 0.0;
//     } else {
//         t2 *= t2;
//         n2 = t2 * t2 * grad(gi2, x2, y2, z2);
//     }
//     t_float t3 = 0.6 - x3*x3 - y3*y3 - z3*z3;
//     if (t3 < 0) {
//         n3 = 0.0;
//     } else {
//         t3 *= t3;
//         n3 = t3 * t3 * grad(gi3, x3, y3, z3);
//     }
//     // Add contributions from each corner to get the final noise value.
//     // The result is scaled to stay just inside [-1,1]
//     return 32.0 * (n0 + n1 + n2 + n3);
// }

// noise functions for 2d/3d/4d copied from https://github.com/weswigham/simplex/blob/master/c/src/simplex.c

int gradients3d[12][3] = {{1,1,0},{-1,1,0},{1,-1,0},{-1,-1,0},
{1,0,1},{-1,0,1},{1,0,-1},{-1,0,-1},
{0,1,1},{0,-1,1},{0,1,-1},{0,-1,-1}};
int gradients4d[32][4] = {{0,1,1,1}, {0,1,1,-1}, {0,1,-1,1}, {0,1,-1,-1},
{0,-1,1,1}, {0,-1,1,-1}, {0,-1,-1,1}, {0,-1,-1,-1},
{1,0,1,1}, {1,0,1,-1}, {1,0,-1,1}, {1,0,-1,-1},
{-1,0,1,1}, {-1,0,1,-1}, {-1,0,-1,1}, {-1,0,-1,-1},
{1,1,0,1}, {1,1,0,-1}, {1,-1,0,1}, {1,-1,0,-1},
{-1,1,0,1}, {-1,1,0,-1}, {-1,-1,0,1}, {-1,-1,0,-1},
{1,1,1,0}, {1,1,-1,0}, {1,-1,1,0}, {1,-1,-1,0},
{-1,1,1,0}, {-1,1,-1,0}, {-1,-1,1,0}, {-1,-1,-1,0}};
int simplex[64][4] = {
{0,1,2,3},{0,1,3,2},{0,0,0,0},{0,2,3,1},{0,0,0,0},{0,0,0,0},{0,0,0,0},{1,2,3,0},
{0,2,1,3},{0,0,0,0},{0,3,1,2},{0,3,2,1},{0,0,0,0},{0,0,0,0},{0,0,0,0},{1,3,2,0},
{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
{1,2,0,3},{0,0,0,0},{1,3,0,2},{0,0,0,0},{0,0,0,0},{0,0,0,0},{2,3,0,1},{2,3,1,0},
{1,0,2,3},{1,0,3,2},{0,0,0,0},{0,0,0,0},{0,0,0,0},{2,0,3,1},{0,0,0,0},{2,1,3,0},
{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
{2,0,1,3},{0,0,0,0},{0,0,0,0},{0,0,0,0},{3,0,1,2},{3,0,2,1},{0,0,0,0},{3,1,2,0},
{2,1,0,3},{0,0,0,0},{0,0,0,0},{0,0,0,0},{3,1,0,2},{0,0,0,0},{3,2,0,1},{3,2,1,0}};

unsigned int lcg_next(unsigned int *seed) {
    const unsigned int a = 1664525;
    const unsigned int c = 1013904223;
    *seed = (a * (*seed) + c); // Modulo 2^32 is implicit due to unsigned int overflow
    return *seed;
}

void shuffle(int *array, int n, unsigned int seed) {
    for (int i = n - 1; i > 0; i--) {
        unsigned int rand = lcg_next(&seed);
        int j = rand % (i + 1);
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
    shuffle(basePermutation, 256, seed);
    for (int i = 0; i < 256; i++) {
        x->x_perm[i] = basePermutation[i];
        x->x_perm[i+256] = basePermutation[i];
    }
}

static inline int fastfloor(float x) {
    return x > 0 ? (int)x : (int)x - 1;
}

t_float Dot2D(int tbl[],t_float x,t_float y)
{
    return tbl[0]*x + tbl[1]*y; 
}

t_float Dot3D(int tbl[],t_float x,t_float y,t_float z)
{
    return tbl[0]*x + tbl[1]*y + tbl[2]*z;
}

t_float Dot4D(int tbl[],t_float x,t_float y,t_float z,t_float w) 
{
    return tbl[0]*x + tbl[1]*y + tbl[2]*z + tbl[3]*w;
}

t_float Noise2D(t_simplex_tilde *x, t_float x_in, t_float y_in)
{
    t_float n0, n1, n2; // Noise contributions from the three corners
    // Skew the input space to determine which simplex cell we're in
    t_float s = (x_in+y_in)*F2; // Hairy factor for 2D
    int i = fastfloor(x_in+s);
    int j = fastfloor(y_in+s);
    t_float t = (i+j)*G2;
    t_float X0 = i-t; // Unskew the cell origin back to (x,y) space
    t_float Y0 = j-t;
    t_float x0 = x_in-X0; // The x,y distances from the cell origin
    t_float y0 = y_in-Y0;
    
    // For the 2D case, the simplex shape is an equilateral triangle.
    // Determine which simplex we are in.
    int i1, j1; // Offsets for second (middle) corner of simplex in (i,j) coords
    if(x0>y0){
        i1=1; 
        j1=0;  // lower triangle, XY order: (0,0)->(1,0)->(1,1)
    }
    else {
        i1=0;
        j1=1; // upper triangle, YX order: (0,0)->(0,1)->(1,1)
    }
    
    // A step of (1,0) in (i,j) means a step of (1-c,-c) in (x,y), and
    // a step of (0,1) in (i,j) means a step of (-c,1-c) in (x,y), where
    // c = (3-sqrt(3))/6

    t_float x1 = x0 - i1 + G2; // Offsets for middle corner in (x,y) unskewed coords
    t_float y1 = y0 - j1 + G2;
    t_float x2 = x0 - 1.0 + 2.0 * G2; // Offsets for last corner in (x,y) unskewed coords
    t_float y2 = y0 - 1.0 + 2.0 * G2;

    // Work out the hashed gradient indices of the three simplex corners
    int ii = i & 255;
    int jj = j & 255;
    int gi0 = x->x_perm[ii+x->x_perm[jj]] % 12;
    int gi1 = x->x_perm[ii+i1+x->x_perm[jj+j1]] % 12;
    int gi2 = x->x_perm[ii+1+x->x_perm[jj+1]] % 12;

    // Calculate the contribution from the three corners
    t_float t0 = 0.5 - x0*x0-y0*y0;
    if (t0<0){
        n0 = 0.0;
    }
    else{
        t0 = t0 * t0;
        n0 = t0 * t0 * Dot2D(gradients3d[gi0], x0, y0); // (x,y) of Gradients3D used for 2D gradient
    }
    
    t_float t1 = 0.5 - x1*x1-y1*y1;
    if (t1<0){
        n1 = 0.0;
    }
    else{
        t1 = t1*t1;
        n1 = t1 * t1 * Dot2D(gradients3d[gi1], x1, y1);
    }
    
    t_float t2 = 0.5 - x2*x2-y2*y2;
    if (t2<0){
        n2 = 0.0;
    }
    else{
        t2 = t2*t2;
        n2 = t2 * t2 * Dot2D(gradients3d[gi2], x2, y2);
    }
    
    // Add contributions from each corner to get the final noise value.
    // The result is scaled to return values in the localerval [-1,1].
    t_float ret = (70.0 * (n0 + n1 + n2));
    return ret;
}

t_float Noise3D(t_simplex_tilde *x, t_float x_in, t_float y_in, t_float z_in)
{
    t_float n0, n1, n2, n3; // Noise contributions from the four corners
    
    // Skew the input space to determine which simplex cell we're in

    t_float s = (x_in+y_in+z_in)*F3; // Very nice and simple skew factor for 3D
    int i = fastfloor(x_in+s);
    int j = fastfloor(y_in+s);
    int k = fastfloor(z_in+s);
    
    t_float t = (i+j+k)*G3;
    
    t_float X0 = i-t; // Unskew the cell origin back to (x,y,z) space
    t_float Y0 = j-t;
    t_float Z0 = k-t;
    
    t_float x0 = x_in-X0; // The x,y,z distances from the cell origin
    t_float y0 = y_in-Y0;
    t_float z0 = z_in-Z0;
    
    // For the 3D case, the simplex shape is a slightly irregular tetrahedron.
    // Determine which simplex we are in.
    int i1, j1, k1; // Offsets for second corner of simplex in (i,j,k) coords
    int i2, j2, k2; // Offsets for third corner of simplex in (i,j,k) coords
    
    if (x0>=y0){
        if (y0>=z0){
            i1=1; j1=0; k1=0; i2=1; j2=1; k2=0; // X Y Z order
        }
        else if (x0>=z0){
            i1=1; j1=0; k1=0; i2=1; j2=0; k2=1; // X Z Y order
        }
        else{
            i1=0; j1=0; k1=1; i2=1; j2=0; k2=1;  // Z X Y order
        }
    }
    else{ // x0<y0
        if (y0<z0){
            i1=0; j1=0; k1=1; i2=0; j2=1; k2=1; // Z Y X order
        }
        else if (x0<z0){
            i1=0; j1=1; k1=0; i2=0; j2=1; k2=1; // Y Z X order
        }
        else{ 
            i1=0; j1=1; k1=0; i2=1; j2=1; k2=0; // Y X Z order
        }
    }
    
    // A step of (1,0,0) in (i,j,k) means a step of (1-c,-c,-c) in (x,y,z),
    // a step of (0,1,0) in (i,j,k) means a step of (-c,1-c,-c) in (x,y,z), and
    // a step of (0,0,1) in (i,j,k) means a step of (-c,-c,1-c) in (x,y,z), where
    // c = 1/6.
    
    t_float x1 = x0 - i1 + G3; // Offsets for second corner in (x,y,z) coords
    t_float y1 = y0 - j1 + G3;
    t_float z1 = z0 - k1 + G3;
    
    t_float x2 = x0 - i2 + 2.0*G3; // Offsets for third corner in (x,y,z) coords
    t_float y2 = y0 - j2 + 2.0*G3;
    t_float z2 = z0 - k2 + 2.0*G3;
    
    t_float x3 = x0 - 1.0 + 3.0*G3; // Offsets for last corner in (x,y,z) coords
    t_float y3 = y0 - 1.0 + 3.0*G3;
    t_float z3 = z0 - 1.0 + 3.0*G3;
    
    // Work out the hashed gradient indices of the four simplex corners
    int ii = i & 255;
    int jj = j & 255;
    int kk = k & 255;
    
    int gi0 = x->x_perm[ii+x->x_perm[jj+x->x_perm[kk]]] % 12;
    int gi1 = x->x_perm[ii+i1+x->x_perm[jj+j1+x->x_perm[kk+k1]]] % 12;
    int gi2 = x->x_perm[ii+i2+x->x_perm[jj+j2+x->x_perm[kk+k2]]] % 12;
    int gi3 = x->x_perm[ii+1+x->x_perm[jj+1+x->x_perm[kk+1]]] % 12;
    
    // Calculate the contribution from the four corners
    t_float t0 = 0.5 - x0*x0 - y0*y0 - z0*z0;
    
    if (t0<0){
        n0 = 0.0;
    }
    else {
        t0 = t0*t0;
        n0 = t0 * t0 * Dot3D(gradients3d[gi0], x0, y0, z0);
    }
    
    t_float t1 = 0.5 - x1*x1 - y1*y1 - z1*z1;
    
    if (t1<0){ 
        n1 = 0.0;
    }
    else{
        t1 = t1*t1;
        n1 = t1 * t1 * Dot3D(gradients3d[gi1], x1, y1, z1);
    }
    
    t_float t2 = 0.5 - x2*x2 - y2*y2 - z2*z2;
    
    if (t2<0){
        n2 = 0.0;
    }
    else{
        t2 = t2*t2;
        n2 = t2 * t2 * Dot3D(gradients3d[gi2], x2, y2, z2);
    }
    
    t_float t3 = 0.5 - x3*x3 - y3*y3 - z3*z3;
    
    if (t3<0){
        n3 = 0.0;
    }
    else{
        t3 = t3*t3;
        n3 = t3 * t3 * Dot3D(gradients3d[gi3], x3, y3, z3);
    }
    
    // Add contributions from each corner to get the final noise value.
    // The result is scaled to stay just inside [-1,1]
    t_float retval = 32.0*(n0 + n1 + n2 + n3);
    return retval;
}

t_float Noise4D(t_simplex_tilde *x, t_float x_in, t_float y_in, t_float z_in, t_float w_in)
{
    // The skewing and unskewing factors are hairy again for the 4D case
    t_float n0, n1, n2, n3, n4; // Noise contributions from the five corners
    // Skew the (x,y,z,w) space to determine which cell of 24 simplices we're in
    t_float s = (x_in + y_in + z_in + w_in) * F4; // Factor for 4D skewing
    int i = fastfloor(x_in + s);
    int j = fastfloor(y_in + s);
    int k = fastfloor(z_in + s);
    int l = fastfloor(w_in + s);
    t_float t = (i + j + k + l) * G4; // Factor for 4D unskewing
    t_float X0 = i - t; // Unskew the cell origin back to (x,y,z,w) space
    t_float Y0 = j - t;
    t_float Z0 = k - t;
    t_float W0 = l - t;
    t_float x0 = x_in - X0; // The x,y,z,w distances from the cell origin
    t_float y0 = y_in - Y0;
    t_float z0 = z_in - Z0;
    t_float w0 = w_in - W0;
    // For the 4D case, the simplex is a 4D shape I won't even try to describe.
    // To find out which of the 24 possible simplices we're in, we need to
    // determine the magnitude ordering of x0, y0, z0 and w0.
    // The method below is a good way of finding the ordering of x,y,z,w and
    // then find the correct traversal order for the simplex we’re in.
    // First, six pair-wise comparisons are performed between each possible pair
    // of the four coordinates, and the results are used to add up binary bits
    // for an localeger index.
    int c1 = (x0 > y0) ? 32 : 1;
    int c2 = (x0 > z0) ? 16 : 1;
    int c3 = (y0 > z0) ? 8 : 1;
    int c4 = (x0 > w0) ? 4 : 1;
    int c5 = (y0 > w0) ? 2 : 1;
    int c6 = (z0 > w0) ? 1 : 1;
    int c = c1 + c2 + c3 + c4 + c5 + c6;
    int i1, j1, k1, l1; // The localeger offsets for the second simplex corner
    int i2, j2, k2, l2; // The localeger offsets for the third simplex corner
    int i3, j3, k3, l3; // The localeger offsets for the fourth simplex corner
    
    // simplex[c] is a 4-vector with the numbers 0, 1, 2 and 3 in some order.
    // Many values of c will never occur, since e.g. x>y>z>w makes x<z, y<w and x<w
    // impossible. Only the 24 indices which have non-zero entries make any sense.
    // We use a thresholding to set the coordinates in turn from the largest magnitude.
    // The number 3 in the "simplex" array is at the position of the largest coordinate.
    
    i1 = simplex[c][0]>=3 ? 1 : 0;
    j1 = simplex[c][1]>=3 ? 1 : 0;
    k1 = simplex[c][2]>=3 ? 1 : 0;
    l1 = simplex[c][3]>=3 ? 1 : 0;
    // The number 2 in the "simplex" array is at the second largest co:dinate.
    i2 = simplex[c][0]>=2 ? 1 : 0;
    j2 = simplex[c][1]>=2 ? 1 : 0;
    k2 = simplex[c][2]>=2 ? 1 : 0;
    l2 = simplex[c][3]>=2 ? 1 : 0;
    // The number 1 in the "simplex" array is at the second smallest co:dinate.
    i3 = simplex[c][0]>=1 ? 1 : 0;
    j3 = simplex[c][1]>=1 ? 1 : 0;
    k3 = simplex[c][2]>=1 ? 1 : 0;
    l3 = simplex[c][3]>=1 ? 1 : 0;
    // The fifth corner has all coordinate offsets = 1, so no need to look that up.
    t_float x1 = x0 - i1 + G4; // Offsets for second corner in (x,y,z,w) coords
    t_float y1 = y0 - j1 + G4;
    t_float z1 = z0 - k1 + G4;
    t_float w1 = w0 - l1 + G4;
    t_float x2 = x0 - i2 + 2.0*G4; // Offsets for third corner in (x,y,z,w) coords
    t_float y2 = y0 - j2 + 2.0*G4;
    t_float z2 = z0 - k2 + 2.0*G4;
    t_float w2 = w0 - l2 + 2.0*G4;
    t_float x3 = x0 - i3 + 3.0*G4; // Offsets for fourth corner in (x,y,z,w) coords
    t_float y3 = y0 - j3 + 3.0*G4;
    t_float z3 = z0 - k3 + 3.0*G4;
    t_float w3 = w0 - l3 + 3.0*G4;
    t_float x4 = x0 - 1.0 + 4.0*G4; // Offsets for last corner in (x,y,z,w) coords
    t_float y4 = y0 - 1.0 + 4.0*G4;
    t_float z4 = z0 - 1.0 + 4.0*G4;
    t_float w4 = w0 - 1.0 + 4.0*G4;
    
    // Work out the hashed gradient indices of the five simplex corners
    int ii = i & 255;
    int jj = j & 255;
    int kk = k & 255;
    int ll = l & 255;
    int gi0 = x->x_perm[ii+x->x_perm[jj+x->x_perm[kk+x->x_perm[ll]]]] % 32;
    int gi1 = x->x_perm[ii+i1+x->x_perm[jj+j1+x->x_perm[kk+k1+x->x_perm[ll+l1]]]] % 32;
    int gi2 = x->x_perm[ii+i2+x->x_perm[jj+j2+x->x_perm[kk+k2+x->x_perm[ll+l2]]]] % 32;
    int gi3 = x->x_perm[ii+i3+x->x_perm[jj+j3+x->x_perm[kk+k3+x->x_perm[ll+l3]]]] % 32;
    int gi4 = x->x_perm[ii+1+x->x_perm[jj+1+x->x_perm[kk+1+x->x_perm[ll+1]]]] % 32;
    
    
    // Calculate the contribution from the five corners
    t_float t0 = 0.5 - x0*x0 - y0*y0 - z0*z0 - w0*w0;
    if (t0<0){
        n0 = 0.0;
    }
    else{
        t0 = t0*t0;
        n0 = t0 * t0 * Dot4D(gradients4d[gi0], x0, y0, z0, w0);
    }
    
    t_float t1 = 0.5 - x1*x1 - y1*y1 - z1*z1 - w1*w1;
    if (t1<0){
        n1 = 0.0;
    }
    else{
        t1 = t1*t1;
        n1 = t1 * t1 * Dot4D(gradients4d[gi1], x1, y1, z1, w1);
    }
    
    t_float t2 = 0.5 - x2*x2 - y2*y2 - z2*z2 - w2*w2;
    if (t2<0){
        n2 = 0.0;
    }
    else{
        t2 = t2*t2;
        n2 = t2 * t2 * Dot4D(gradients4d[gi2], x2, y2, z2, w2);
    }
    
    t_float t3 = 0.5 - x3*x3 - y3*y3 - z3*z3 - w3*w3;
    if (t3<0){
        n3 = 0.0;
    }
    else {
        t3 = t3*t3;
        n3 = t3 * t3 * Dot4D(gradients4d[gi3], x3, y3, z3, w3);
    }
    
    t_float t4 = 0.5 - x4*x4 - y4*y4 - z4*z4 - w4*w4;
    if (t4<0){
        n4 = 0.0;
    }
    else{
        t4 = t4*t4;
        n4 = t4 * t4 * Dot4D(gradients4d[gi4], x4, y4, z4, w4);
    }
    
    // Sum up and scale the result to cover the range [-1,1]
    t_float retval = 27.0 * (n0 + n1 + n2 + n3 + n4);
    return retval;
}



static t_int *simplex_tilde_perform(t_int *w) {
    int i;
    t_simplex_tilde *x = (t_simplex_tilde *)(w[1]);
    t_int n = (t_int)(w[2]);
    t_int nchans = (t_int)(w[3]);
    t_sample *in_coord = (t_sample *)(w[4]);
    t_sample *in_persistence = (t_sample *)(w[5]);
    t_sample *out = (t_sample *)(w[6]);
    if (nchans == 2) {
        for(i = 0; i < n; i++) {
            t_float x_coord = in_coord[i];
            t_float y_coord = in_coord[n + i];
            t_float persistence = in_persistence[i];
            t_float result = Noise2D(x, x_coord, y_coord) * persistence;
            *out++ = result;
        }
    } else if (nchans == 3) {
        for(i = 0; i < n; i++) {
            t_float x_coord = in_coord[i];
            t_float y_coord = in_coord[n + i];
            t_float z_coord = in_coord[2*n + i];
            t_float persistence = in_persistence[i];
            t_float result = Noise3D(x, x_coord, y_coord, z_coord) * persistence;
            *out++ = result;
        }
    } else if (nchans == 4) {
        for(i = 0; i < n; i++) {
            t_float x_coord = in_coord[i];
            t_float y_coord = in_coord[n + i];
            t_float z_coord = in_coord[2*n + i];
            t_float w_coord = in_coord[3*n + i];
            t_float persistence = in_persistence[i];
            t_float result = Noise4D(x, x_coord, y_coord, z_coord, w_coord) * persistence;
            *out++ = result;
        }
    } else {
        for(i = 0; i < n; i++)
            *out++ = 0;
    }
    return(w+7);
}

void simplex_tilde_dsp(t_simplex_tilde *x, t_signal **sp) {
    signal_setmultiout(&sp[2], 1);
    dsp_add(simplex_tilde_perform, 6, x, (t_int)sp[0]->s_n, (t_int)sp[0]->s_nchans, sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec);
}

static void simplex_tilde_octaves(t_simplex_tilde *x, t_floatarg f){
    x->x_octaves = fastfloor(f);
}

static void simplex_tilde_seed(t_simplex_tilde *x, t_floatarg f){
    initPermutation(x, fastfloor(f));
}

void *simplex_tilde_new(t_symbol *s, int ac, t_atom *av) {
    float persistence;

    t_simplex_tilde *x = (t_simplex_tilde *)pd_new(simplex_tilde_class);
    x->x_octaves = (ac >= 1) ? max(1, atom_getintarg(0, ac, av)) : 1;
    persistence = (ac >= 2) ? atom_getfloatarg(1, ac, av) : PERSISTENCE;
    x->x_inlet_persistence = inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal); // persistence
        pd_float((t_pd *)x->x_inlet_persistence, persistence);
    outlet_new(&x->x_obj, &s_signal);

    initPermutation(x, (int)0);
    return(x);
}

void *simplex_tilde_free(t_simplex_tilde *x){
    inlet_free(x->x_inlet_persistence);
    return(void *)x;
}

void simplex_tilde_setup(void) {
    simplex_tilde_class = class_new(gensym("simplex~"), (t_newmethod)simplex_tilde_new,
        (t_method)simplex_tilde_free, sizeof(t_simplex_tilde), CLASS_MULTICHANNEL, A_GIMME, 0);
    class_addmethod(simplex_tilde_class, nullfn, gensym("signal"), 0);
    class_addmethod(simplex_tilde_class, (t_method)simplex_tilde_dsp, gensym("dsp"), 0);
    class_addmethod(simplex_tilde_class, (t_method)simplex_tilde_octaves, gensym("octaves"), A_FLOAT, 0);
    class_addmethod(simplex_tilde_class, (t_method)simplex_tilde_seed, gensym("seed"), A_FLOAT, 0);
}
