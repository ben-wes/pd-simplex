// external created by Ben Wesch, 2024

#include "m_pd.h"
#include "simplex_common.h"
#ifdef _WIN32
#include <windows.h>
#else
#include <dlfcn.h>
#endif

typedef void (*t_signal_setmultiout)(t_signal **, int); 
t_signal_setmultiout g_signal_setmultiout;
static t_class *simplex_tilde_class;

typedef struct _simplex_tilde {
    t_object x_obj;
    t_inlet *inlet_persistence;
    t_sample **coordinate_vector;
    t_sample *in_persistence;
    t_sample *out_value;
    t_sample **derivatives_vector;
    t_float *derivatives;
    t_float octave_factors[MAX_OCTAVES];
    t_sample f;
    int multichannel;
    int normalize;
    int octaves;
    int dimensions;
    unsigned char perm[512];
} t_simplex_tilde;

static t_int *simplex_tilde_perform(t_int *w) {
    t_simplex_tilde *x = (t_simplex_tilde *)(w[1]);
    t_simplex_config cfg = {
        .octave_factors = x->octave_factors,
        .octaves = x->octaves,
        .normalize = x->normalize,
        .perm = x->perm
    };
    t_int n_samples = w[2];
    int sample, dim;
    for (sample = 0; sample < n_samples; sample++) {
        t_float pos[MAX_DIMENSIONS] = {0};
        for (dim = 0; dim < x->dimensions; dim++)
            pos[dim] = x->coordinate_vector[dim][sample];
        x->out_value[sample] = generate_noise(&cfg, pos, x->in_persistence[sample], 
                                              x->dimensions - 1, x->derivatives);
        if (x->derivatives) {
            for (dim = 0; dim < x->dimensions; dim++)
                x->derivatives_vector[dim][sample] = x->derivatives[dim];
        }
    }
    return w+3;
}

void simplex_tilde_dsp(t_simplex_tilde *x, t_signal **sp) {
    t_int n_samples = (t_int)sp[0]->s_n;
    int dim;

    if (x->multichannel) {
        x->dimensions = (int)min(MAX_DIMENSIONS, sp[0]->s_nchans);
        for (dim = 0; dim < x->dimensions; dim++)
            x->coordinate_vector[dim] = sp[0]->s_vec + n_samples * dim;
        x->in_persistence = sp[1]->s_vec;
        g_signal_setmultiout(&sp[2], 1);
        x->out_value = sp[2]->s_vec;
        if (x->derivatives) {
            g_signal_setmultiout(&sp[3], x->dimensions);
            for (dim = 0; dim < x->dimensions; dim++)
                x->derivatives_vector[dim] = sp[3]->s_vec + n_samples * dim;
        }
    } else {
        for (dim = 0; dim < x->dimensions; dim++)
            x->coordinate_vector[dim] = sp[dim]->s_vec;
        x->in_persistence = sp[x->dimensions]->s_vec;
        if (g_signal_setmultiout) g_signal_setmultiout(&sp[x->dimensions + 1], 1);
        x->out_value = sp[x->dimensions + 1]->s_vec;
        if (x->derivatives) {
            for (dim = 0; dim < x->dimensions; dim++) {
                if (g_signal_setmultiout) g_signal_setmultiout(&sp[x->dimensions + 2 + dim], 1);
                x->derivatives_vector[dim] = sp[x->dimensions + 2 + dim]->s_vec;
            }
        }
    }
    dsp_add(simplex_tilde_perform, 2, x, n_samples);
}

static void simplex_tilde_octaves(t_simplex_tilde *x, t_floatarg f) {
    x->octaves = (int)clamp(f, 1, MAX_OCTAVES);
    init_octave_factors(x->octave_factors, x->octaves);
}

static void simplex_tilde_persistence(t_simplex_tilde *x, t_floatarg f) {
    pd_float((t_pd *)x->inlet_persistence, f);
}

static void simplex_tilde_normalize(t_simplex_tilde *x, t_symbol *s, int ac, t_atom *av) {
    // activate normalization if no argument is given or anything that is not false
    x->normalize = !ac || atom_getfloat(av);
    (void)s;
}

static void simplex_tilde_seed(t_simplex_tilde *x, t_symbol *s, int ac, t_atom *av) {
    if (ac)
        init_permutation_with_seed(x->perm, atom_getfloat(av));
    else
        init_permutation(x->perm);
    (void)s;
}

static void simplex_tilde_coeffs(t_simplex_tilde *x, t_symbol *s, int ac, t_atom *av) {
    int i;
    x->octaves = clamp(ac, 1, MAX_OCTAVES);
    for (i = 0; i < x->octaves; i++) {
        x->octave_factors[i] = atom_getfloat(av);
        av++;
    }
    (void)s;
}

static void simplex_tilde_free(t_simplex_tilde *x) {
    inlet_free(x->inlet_persistence);
    freebytes(x->coordinate_vector, MAX_DIMENSIONS * sizeof(t_sample *));
    if (x->derivatives) {
        freebytes(x->derivatives, MAX_DIMENSIONS * sizeof(t_float));
        freebytes(x->derivatives_vector, MAX_DIMENSIONS * sizeof(t_sample *));
    }
}

static void *simplex_tilde_new(t_symbol *s, int ac, t_atom *av) {
    t_simplex_tilde *x = (t_simplex_tilde *)pd_new(simplex_tilde_class);
    t_float persistence;
    x->normalize = 0;
    x->octaves = 1;
    x->dimensions = 0;
    int i;
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
            pd_error(x, "[simplex~]: invalid argument");
        ac--, av++;
    }
    if (ac) {
        x->octaves = clamp(atom_getint(av), 1, MAX_OCTAVES);
        ac--, av++;
    }
    persistence = ac ? atom_getfloat(av) : DEFAULT_PERSISTENCE;
    init_octave_factors(x->octave_factors, x->octaves);
    x->multichannel = !x->dimensions && g_signal_setmultiout;
    x->coordinate_vector = (t_sample **)getbytes(MAX_DIMENSIONS * sizeof(t_sample *));
    if (!g_signal_setmultiout && !x->dimensions) {
        x->dimensions = 1;
        pd_error(x, "[simplex~]: no multichannel support in Pd %i.%i-%i", maj, min, bug);
        post("[simplex~]: use '-dim <count>' flag to set dimension count (default: 1)");
    }
    if (!x->multichannel) { // add inlets for non-multichannel mode
        for (i=0; i < x->dimensions-1; i++)
            inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal);
    }
    x->inlet_persistence = inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal);
    pd_float((t_pd *)x->inlet_persistence, persistence);

    outlet_new(&x->x_obj, &s_signal); // sampled value outlet
    if (x->derivatives) {
        x->derivatives_vector = (t_sample **)getbytes(MAX_DIMENSIONS * sizeof(t_sample *));
        outlet_new(&x->x_obj, &s_signal); // derivatives outlet
        if (!x->multichannel) { // add outlets for non-multichannel mode
            for (i=0; i < x->dimensions-1; i++)
                outlet_new(&x->x_obj, &s_signal);
        }
    }
    (void)s;
    return x;
}

void simplex_tilde_setup(void) {
// multichannel handling copied from https://github.com/Spacechild1/vstplugin/blob/3f0ed8a800ea238bf204a2ead940b2d1324ac909/pd/src/vstplugin~.cpp#L4122-L4136
#ifdef _WIN32
    // get a handle to the module containing the Pd API functions.
    // NB: GetModuleHandle("pd.dll") does not cover all cases.
    HMODULE module;
    if (GetModuleHandleEx(
            GET_MODULE_HANDLE_EX_FLAG_FROM_ADDRESS | GET_MODULE_HANDLE_EX_FLAG_UNCHANGED_REFCOUNT,
            (LPCSTR)&pd_typedmess, &module)) {
        g_signal_setmultiout = (t_signal_setmultiout)(void *)GetProcAddress(
            module, "signal_setmultiout");
    }
#else
    // search recursively, starting from the main program
    g_signal_setmultiout = (t_signal_setmultiout)dlsym(
        dlopen(NULL, RTLD_NOW), "signal_setmultiout");
#endif

    simplex_tilde_class = class_new(
        gensym("simplex~"),
        (t_newmethod)simplex_tilde_new,
        (t_method)simplex_tilde_free,
        sizeof(t_simplex_tilde), CLASS_MULTICHANNEL, A_GIMME, 0);

    CLASS_MAINSIGNALIN(simplex_tilde_class, t_simplex_tilde, f);
    class_addmethod(simplex_tilde_class, (t_method)simplex_tilde_dsp, gensym("dsp"), 0);
    class_addmethod(simplex_tilde_class, (t_method)simplex_tilde_seed, gensym("seed"), A_GIMME, 0);
    class_addmethod(simplex_tilde_class, (t_method)simplex_tilde_normalize, gensym("normalize"), A_GIMME, 0);
    class_addmethod(simplex_tilde_class, (t_method)simplex_tilde_coeffs, gensym("coeffs"), A_GIMME, 0);
    class_addmethod(simplex_tilde_class, (t_method)simplex_tilde_octaves, gensym("octaves"), A_FLOAT, 0);
    class_addmethod(simplex_tilde_class, (t_method)simplex_tilde_persistence, gensym("persistence"), A_FLOAT, 0);
}
