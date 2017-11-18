#include <stdio.h>

typedef struct {
    double x, y, z, w;
} xyzw;

typedef struct {
    unsigned char r, g, b, a;
} rgba;

typedef struct {
    int len;
    double *vals;
} vec;

typedef struct {
    xyzw pos;
    vec val;
} point;

typedef struct {
    double depth;
    vec varying;
} fragment;

typedef point (*vshader)(const vec uniform, const vec attrs);
typedef rgba  (*fshader)(const vec uniform, const vec varying);

typedef struct {
    vshader vertex_shader;
    fshader fragment_shader;

    int num_vertices;
    vec *vertices;
    int num_triangles;
    int *triangles;

    vec vertex_uniform;
    vec fragment_uniform;
} scene;

typedef struct {
    int w;
    int h;
    uint *data;
} image;

void render(const scene scene, image image);

void render_main(
        FILE *out, const uint frames,
        scene scene, const int w, const int h);

void render_main_default(scene scene);

#define MIN(a,b) \
  ({ __typeof__ (a) _a = (a); \
      __typeof__ (b) _b = (b); \
    _a < _b ? _a : _b; })
#define MAX(a,b) \
  ({ __typeof__ (a) _a = (a); \
      __typeof__ (b) _b = (b); \
    _a > _b ? _a : _b; })

