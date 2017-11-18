#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "renderer.h"

#define HYPO_LEN 0.5773502691896257
#define LEG_LEN 0.28867513459481287
#define HEIGHT 0.816496580927726
#define A(a, b, c, d, e, f, g) {7, (double[]){a, b, c, d, e, f, g}}
int tetra_num_vertices = 4;
vec tetra_vertices[4] = {
    A(   0,      0, -HYPO_LEN, 1, 0, 0, 0),
    A(-0.5,      0,   LEG_LEN, 0, 1, 0, 0),
    A( 0.5,      0,   LEG_LEN, 0, 0, 1, 0),
    A(   0, HEIGHT,         0, 0, 0, 0, 1),
};
int tetra_num_triangles = 4;
int tetra_triangles[4 * 3] = {
    0, 1, 2,
    0, 1, 3,
    1, 2, 3,
    2, 0, 3,
};

point tetra_vshader(const vec uniform, const vec attrs) {
    double b = uniform.vals[0];
    double *v = attrs.vals;
    point p;
    p.pos.x = v[0] * cos(b) + v[2] * sin(b);
    double ty = p.pos.y = v[1],
           tz = v[2] * cos(b) - v[0] * sin(b);
    p.pos.y = ty * cos(b / 2) - tz * sin(b / 2);
    p.pos.z = ty * sin(b / 2) + tz * cos(b / 2);
    p.pos.w = (p.pos.z + 3) / 2;
    p.val.len = 4;
    p.val.vals = calloc(4, sizeof(double));
    p.val.vals[0] = v[3];
    p.val.vals[1] = v[4];
    p.val.vals[2] = v[5];
    p.val.vals[3] = v[6];
    return p;
}

rgba tetra_fshader(const vec uniform, const vec varying) {
    int count = 0;
    for (int i = 0; i < 4; i++)
        if (varying.vals[i] > 0.03) count++;
    return count > 2 ?
        (rgba){0xFF, 0xFF, 0xFF, 0xFF} : (rgba){0xFF, 0, 0, 0xFF};
}

int main(void) {
    scene scene = {
        tetra_vshader, tetra_fshader,
        tetra_num_vertices, tetra_vertices,
        tetra_num_vertices, tetra_triangles};
    render_main_default(scene);
    return 0;
}

