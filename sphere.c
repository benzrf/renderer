#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "renderer.h"
#include "phong.h"

void sphere_coords(double theta, double phi, double *vert) {
    double theta_ = 2 * M_PI * theta, phi_ = 2 * M_PI * phi;
    vert[0] = sin(theta_) * cos(phi_) / 2;
    vert[1] = sin(theta_) * sin(phi_) / 2;
    vert[2] = cos(theta_) / 2;
}

void sphere_vertices_triangles(
    int res_t, int res_p,
    vec **vertices,
    int **triangles
) {
    vec *curV = *vertices = calloc((res_t - 1) * res_p, sizeof(vec));
    for (int theta = 1; theta < res_t; theta++) {
        for (int phi = 0; phi < res_p; phi++) {
            float th = (float)theta / (2 * (float)res_t),
                  ph = (float)phi / ((float)res_p);
            curV->len = 5;
            curV->vals = calloc(curV->len, sizeof(double));
            sphere_coords(th, ph, curV->vals);
            curV->vals[3] = theta;
            curV->vals[4] = phi;
            curV += 1;
        }
    }
    int *curT = *triangles = calloc((res_t - 2) * res_p * 2 * 3, sizeof(int));
    for (int theta = 0; theta < res_t - 2; theta++) {
        for (int phi = 0; phi < res_p; phi++) {
            int me = phi + res_p * theta,
                right = (phi - 1 + res_p) % res_p + res_p * theta,
                down = me + res_t,
                down_left = (phi + 1) % res_p + res_p * (theta + 1);
            curT[0] = me;
            curT[1] = right;
            curT[2] = down;
            curT[3] = me;
            curT[4] = down;
            curT[5] = down_left;
            curT += 6;
        }
    }
}

point sphere_vshader(const vec uniform, const vec attrs) {
    double b = uniform.vals[0];
    double *v = attrs.vals;
    point p;
    p.pos.x = v[0] * cos(b) + v[2] * sin(b);
    p.pos.y = v[1];
    p.pos.z = v[2] * cos(b) - v[0] * sin(b);
    p.pos.w = (p.pos.z + 3) / 2;
    p.val.len = 5;
    p.val.vals = calloc(p.val.len, sizeof(double));
    p.val.vals[0] = p.pos.x;
    p.val.vals[1] = p.pos.y;
    p.val.vals[2] = p.pos.z;
    p.val.vals[3] = fmod(v[3], 3) / 2;
    p.val.vals[4] = fmod(v[4], 3) / 2;
    return p;
}

uint *sphere_tex;
int stw, sth;
rgba sphere_fshader(const vec uniform, const vec varying) {
    double *v = varying.vals,
           normal[3] = T(v[_]);
    normalize(normal);
    double intensity = phong(v, normal, m) + phong(v, normal, n);
    intensity = MAX(0, MIN(1, intensity));
    int tx = MAX(0, MIN(stw - 1, v[3] * stw)),
        ty = MAX(0, MIN(sth - 1, (1 - v[4]) * sth));
    uint texel = sphere_tex[ty * stw + tx];
    return (rgba){
        (texel         & 0xFF) * intensity,
        ((texel >> 8)  & 0xFF) * intensity,
        ((texel >> 16) & 0xFF) * intensity,
        ((texel >> 24) & 0xFF)};
}

uint *load_pam(FILE *f, int *w, int *h) {
    fseek(f, 9, SEEK_CUR);
    fscanf(f, "%d", w);
    fseek(f, 8, SEEK_CUR);
    fscanf(f, "%d", h);
    fseek(f, 46, SEEK_CUR);
    uint *dat = calloc(*w * *h, sizeof(uint));
    fread(dat, sizeof(uint), *w * *h, f);
    return dat;
}

#define RES_T 15
#define RES_P 15
int main(void) {
    int sphere_num_vertices = (RES_T - 1) * RES_P,
        sphere_num_triangles = (RES_T - 2) * RES_P * 2;
    vec *sphere_vertices;
    int *sphere_triangles;
    sphere_vertices_triangles(
            RES_T, RES_P, &sphere_vertices, &sphere_triangles);
    FILE *f = fopen("gear.pam", "r");
    sphere_tex = load_pam(f, &stw, &sth);
    fclose(f);

    render_main_default(
            sphere_vshader, sphere_fshader,
            sphere_num_vertices, sphere_vertices,
            sphere_num_triangles, sphere_triangles);
    return 0;
}

