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

void render(
        const vshader vertex_shader, const fshader fragment_shader,
        const vec vertex_uniform, const vec fragment_uniform,
        const int num_vertices, const vec *vertices,
        const int num_triangles, const int *triangles,
        const int w, const int h, uint *img);

void render_main(
        FILE *out, const uint frames,
        const vshader vertex_shader, const fshader fragment_shader,
        const int num_vertices, const vec *vertices,
        const int num_triangles, const int *triangles,
        const int w, const int h);

void render_main_default(
        const vshader vertex_shader, const fshader fragment_shader,
        const int num_vertices, const vec *vertices,
        const int num_triangles, const int *triangles);

#define MIN(a,b) \
  ({ __typeof__ (a) _a = (a); \
      __typeof__ (b) _b = (b); \
    _a < _b ? _a : _b; })
#define MAX(a,b) \
  ({ __typeof__ (a) _a = (a); \
      __typeof__ (b) _b = (b); \
    _a > _b ? _a : _b; })

