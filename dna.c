#include <stdlib.h>
#include <math.h>
#include "renderer.h"
#include "phong.h"

void dna_transform(double time, double v[3]) {
    double x = v[0] / 2, y = v[2] / -2, z = v[1] / -2;
    v[0] = x * cos(time) + z * sin(time);
    v[1] = y;
    v[2] = z * cos(time) - x * sin(time);
}

point dna_vshader(const vec uniform, const vec attrs) {
    double time = uniform.vals[0];
    double *v = attrs.vals;
    point p;
    p.val.len = 6;
    p.val.vals = calloc(p.val.len, sizeof(double));
    p.val.vals[0] = v[0];
    p.val.vals[1] = v[1];
    p.val.vals[2] = v[2];
    p.val.vals[3] = v[5];
    p.val.vals[4] = v[6];
    p.val.vals[5] = v[7];
    dna_transform(time, p.val.vals);
    dna_transform(time, p.val.vals + 3);
    normalize(p.val.vals + 3);
    p.pos.x = p.val.vals[0];
    p.pos.y = p.val.vals[1];
    p.pos.z = p.val.vals[2];
    p.pos.w = (p.pos.z + 3) / 2;
    return p;
}

rgba dna_fshader(const vec uniform, const vec varying) {
    double *v = varying.vals,
           *normal = v + 3;
    double intensity = phong(v, normal, m) + phong(v, normal, n);
    intensity = MAX(0, MIN(1, intensity));
    return (rgba){0xFF * intensity, 0, 0, 0xFF};
}

#define VERTEX_DAT_ENTRIES 155000
int main(void) {
    int dna_num_vertices = 19375;
    int dna_num_triangles = 31799;

    double *dna_vertex_dat = calloc(VERTEX_DAT_ENTRIES, sizeof(double));
    vec *dna_vertices = calloc(19375, sizeof(vec));
    int *dna_triangles = calloc(31799 * 3, sizeof(int));

    FILE *f = fopen("dna_vertices", "r");
    fread(dna_vertex_dat, VERTEX_DAT_ENTRIES, sizeof(double), f);
    fclose(f);
    f = fopen("dna_triangles", "r");
    fread(dna_triangles, dna_num_triangles * 3, sizeof(int), f);
    fclose(f);
    for (int i = 0; i < 19375; i++) {
        dna_vertices[i] = (vec){8, dna_vertex_dat + i * 8};
    }

    render_main_default(
            dna_vshader, dna_fshader,
            dna_num_vertices, dna_vertices,
            dna_num_triangles, dna_triangles);
    return 0;
}

