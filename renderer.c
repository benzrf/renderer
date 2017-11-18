#include <stdlib.h>
#include <stdio.h>
#include "renderer.h"
#include "constants.h"

void barycentric(
        double ax, double ay,
        double bx, double by,
        double cx, double cy,
        double x, double y,
        double *u, double *v, double *w) {
    double v0x = bx - ax, v0y = by - ay,
           v1x = cx - ax, v1y = cy - ay,
           v2x = x - ax, v2y = y - ay,
           d00 = v0x * v0x + v0y * v0y,
           d01 = v0x * v1x + v0y * v1y,
           d11 = v1x * v1x + v1y * v1y,
           d20 = v2x * v0x + v2y * v0y,
           d21 = v2x * v1x + v2y * v1y,
           denom = d00 * d11 - d01 * d01;
    *v = (d11 * d20 - d01 * d21) / denom;
    *w = (d00 * d21 - d01 * d20) / denom;
    *u = 1 - *v - *w;
}

#define SX(x) (((x + 1) / 2) * w)
#define SY(y) (((y + 1) / 2) * h)
void rasterize(
        const int w, const int h, fragment *raster,
        const point a, const point b, const point c) {
    double pax = a.pos.x / a.pos.w, pay = a.pos.y / a.pos.w,
           pbx = b.pos.x / b.pos.w, pby = b.pos.y / b.pos.w,
           pcx = c.pos.x / c.pos.w, pcy = c.pos.y / c.pos.w;
    int minX = MAX(0, MIN(SX(pax), MIN(SX(pbx), SX(pcx)))),
        minY = MAX(0, MIN(SY(pay), MIN(SY(pby), SY(pcy)))),
        maxX = MIN(w - 1, MAX(SX(pax), MAX(SX(pbx), SX(pcx)))),
        maxY = MIN(h - 1, MAX(SY(pay), MAX(SY(pby), SY(pcy))));
    for (int ix = minX; ix <= maxX; ix++) {
        for (int iy = minY; iy <= maxY; iy++) {
            double x = (double)ix * 2 / w - 1, y = (double)iy * 2 / h - 1;

            double pu, pv, pw;
            barycentric(pax, pay, pbx, pby, pcx, pcy,
                    x, y, &pu, &pv, &pw);
            if (!(pu >= 0 && pu <= 1 && pv >= 0 && pv <= 1 &&
                        pw >= 0 && pw <= 1)) continue;

            double contraction = 1/(pu/a.pos.w + pv/b.pos.w + pw/c.pos.w);
            double ru, rv, rw;
            barycentric(a.pos.x, a.pos.y, b.pos.x, b.pos.y, c.pos.x, c.pos.y,
                    x * contraction, y * contraction, &ru, &rv, &rw);
            double depth = ru * a.pos.z + rv * b.pos.z + rw * c.pos.z;
            if (depth >= raster[ix + iy * w].depth) continue;

            vec varying;
            varying.len = a.val.len;
            varying.vals = calloc(varying.len, sizeof(double));
            for (int i = 0; i < varying.len; i++)
                varying.vals[i] = ru * a.val.vals[i] +
                    rv * b.val.vals[i] + rw * c.val.vals[i];

            double *old = raster[ix + iy * w].varying.vals;
            raster[ix + iy * w] = (fragment){depth, varying};
            free(old);
        }
    }
}

void render(const scene scene, image image) {
    int w = image.w, h = image.h;
    fragment *raster = calloc(w * h, sizeof(fragment));
    for (int y = 0; y < h; y++)
        for (int x = 0; x < w; x++)
            raster[x + y * w] = (fragment){1.0/0.0, {0, NULL}};

    point vpoints[scene.num_vertices];
    for (int v = 0; v < scene.num_vertices; v++)
        vpoints[v] = scene.vertex_shader(
                scene.vertex_uniform, scene.vertices[v]);

    point a, b, c;
    for (int t = 0; t < scene.num_triangles; t++) {
        a = vpoints[scene.triangles[t * 3]];
        b = vpoints[scene.triangles[t * 3 + 1]];
        c = vpoints[scene.triangles[t * 3 + 2]];
        rasterize(w, h, raster, a, b, c);
    }

    for (int v = 0; v < scene.num_vertices; v++)
        free(vpoints[v].val.vals);

    fragment f;
    rgba col;
    for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++) {
            f = raster[x + y * w];
            if (f.varying.vals == NULL)
                col = (rgba){0, 0, 0, 0xFF};
            else {
                col = scene.fragment_shader(
                        scene.fragment_uniform, f.varying);
                free(f.varying.vals);
            }
            image.data[(h - y - 1) * w + x] =
                col.a << 24 | col.b << 16 |
                col.g << 8  | col.r;
        }
    }
    free(raster);
}

void render_main(
        FILE *out, const uint frames,
        scene scene, const int w, const int h) {
    image image = {w, h};
    image.data = calloc(w * h, sizeof(uint));
    for (int t = 0; !frames || t < frames; t++) {
        scene.vertex_uniform =
            scene.fragment_uniform = (vec){1, (double[]){t / 30.}};
        render(scene, image);
        fwrite(image.data, w * h, sizeof(uint), out);
    }
    free(image.data);
}

void render_main_default(scene scene) {
    render_main(stdout, FRAMES, scene, W, H);
}

