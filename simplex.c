// original external created by Ben Wesch, 2024
// control rate version adapted by Cyrille Henry, 2024

#include "m_pd.h"
#include "simplex_common.h"

static t_class *simplex_class;

typedef struct _simplex {
    t_object x_obj;
    t_float *derivatives;
    t_outlet *derivatives_outlet;
    t_float persistence;
    t_float octave_factors[MAX_OCTAVES];
    t_atom derivative_list[MAX_DIMENSIONS]; 
    t_float last_value;
    t_int last_dim;
    int normalize;
    int octaves;
    int dimensions;
    unsigned char perm[512];
} t_simplex;

void simplex_list(t_simplex *x, t_symbol *s, int ac, t_atom *av) {
    t_simplex_config cfg = {
        .octave_factors = x->octave_factors,
        .octaves = x->octaves,
        .normalize = x->normalize,
        .perm = x->perm
    };
    t_float out = 0.;
    t_float pos[MAX_DIMENSIONS] = {0};
    int dim = ac;
    if(dim==0) {
        if (x->derivatives) {
            outlet_list(x->derivatives_outlet, &s_list, x->last_dim, x->derivative_list);
        }
        outlet_float(x->x_obj.ob_outlet, x->last_value);
    	return;
    }

    if(dim > MAX_DIMENSIONS) {
        pd_error(x, "too much data, limiting to %d dimensions", MAX_DIMENSIONS);
        dim = MAX_DIMENSIONS;
    }

    x->last_dim = dim;

    for(int i=0; i<dim; i++) {
        pos[i] = atom_getfloat(av+i);
    }

    out = generate_noise(&cfg, pos, x->persistence, dim-1, x->derivatives);
    x->last_value = out;

    if (x->derivatives) {
        for (int i = 0; i < dim; i++) {
            SETFLOAT(&(x->derivative_list[i]), x->derivatives[i]);
        }
        outlet_list(x->derivatives_outlet, &s_list, dim, x->derivative_list);
    }

    outlet_float(x->x_obj.ob_outlet, out);
    (void)s;
}

static void simplex_octaves(t_simplex *x, t_floatarg f) {
    x->octaves = (int)clamp(f, 1, MAX_OCTAVES);
    init_octave_factors(x->octave_factors, x->octaves);
}

static void simplex_persistence(t_simplex *x, t_floatarg f) {
    x->persistence = f;
}

static void simplex_normalize(t_simplex *x, t_symbol *s, int ac, t_atom *av) {
    // activate normalization if no argument is given or anything that is not false
    x->normalize = !ac || atom_getfloat(av);
    (void)s;
}

static void simplex_seed(t_simplex *x, t_symbol *s, int ac, t_atom *av) {
    if (ac)
        init_permutation_with_seed(x->perm, atom_getfloat(av));
    else
        init_permutation(x->perm);
    (void)s;
}

static void simplex_coeffs(t_simplex *x, t_symbol *s, int ac, t_atom *av) {
    int i;
    x->octaves = clamp(ac, 1, MAX_OCTAVES);
    for (i = 0; i < x->octaves; i++) {
        x->octave_factors[i] = atom_getfloat(av);
        av++;
    }
    (void)s;
}

static void simplex_free(t_simplex *x) {
    if (x->derivatives) {
        freebytes(x->derivatives, MAX_DIMENSIONS * sizeof(t_float));
        outlet_free(x->derivatives_outlet);
    }
}

static void *simplex_new(t_symbol *s, int ac, t_atom *av) {
    t_simplex *x = (t_simplex *)pd_new(simplex_class);
    outlet_new(&x->x_obj, &s_float);
    x->normalize = 0;
    x->octaves = 1;
    x->dimensions = 0;
    x->last_value = 0.;
    x->last_dim = 0;
    int maj = 0, min = 0, bug = 0;
    sys_getversion(&maj, &min, &bug);
    init_permutation(x->perm);
    while (ac && av->a_type == A_SYMBOL) {
        if (atom_getsymbol(av) == gensym("-n"))
            x->normalize = 1;
        else if (atom_getsymbol(av) == gensym("-d"))
            x->derivatives = (t_float *)getbytes(MAX_DIMENSIONS * sizeof(t_float));
        else if (atom_getsymbol(av) == gensym("-dim")) {
            ac--, av++;
            x->dimensions = clamp((int)atom_getint(av), 1, MAX_DIMENSIONS);
        } else if (atom_getsymbol(av) == gensym("-s")) {
            ac--, av++;
            init_permutation_with_seed(x->perm, (unsigned int)atom_getint(av));
        } else
            pd_error(x, "[simplex]: invalid argument");
        ac--, av++;
    }
    if (ac) {
        x->octaves = clamp(atom_getint(av), 1, MAX_OCTAVES);
        ac--, av++;
    }
    x->persistence = ac ? atom_getfloat(av) : DEFAULT_PERSISTENCE;
    init_octave_factors(x->octave_factors, x->octaves);

    if (x->derivatives) {
      x->derivatives_outlet = outlet_new(&x->x_obj, 0);
    }
    (void)s;
    return x;
}

void simplex_setup(void) {
    simplex_class = class_new(
        gensym("simplex"),
        (t_newmethod)simplex_new,
        (t_method)simplex_free,
        sizeof(t_simplex), CLASS_DEFAULT, A_GIMME, 0);
    class_addmethod(simplex_class, (t_method)simplex_seed, gensym("seed"), A_GIMME, 0);
    class_addmethod(simplex_class, (t_method)simplex_normalize, gensym("normalize"), A_GIMME, 0);
    class_addmethod(simplex_class, (t_method)simplex_coeffs, gensym("coeffs"), A_GIMME, 0);
    class_addmethod(simplex_class, (t_method)simplex_octaves, gensym("octaves"), A_FLOAT, 0);
    class_addmethod(simplex_class, (t_method)simplex_persistence, gensym("persistence"), A_FLOAT, 0);
    class_addmethod(simplex_class, (t_method)simplex_list, &s_list, A_GIMME, 0);
}
