#include <stdlib.h>
#include <stdio.h>
#include <math.h>

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

#define MIN(a,b) \
  ({ __typeof__ (a) _a = (a); \
      __typeof__ (b) _b = (b); \
    _a < _b ? _a : _b; })
#define MAX(a,b) \
  ({ __typeof__ (a) _a = (a); \
      __typeof__ (b) _b = (b); \
    _a > _b ? _a : _b; })
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

void render(
        const vshader vertex_shader, const fshader fragment_shader,
        const vec vertex_uniform, const vec fragment_uniform,
        const int num_vertices, const vec *vertices,
        const int num_triangles, const int *triangles,
        const int w, const int h, uint *img) {
    fragment *raster = calloc(sizeof(fragment), w * h);
    for (int y = 0; y < h; y++)
        for (int x = 0; x < w; x++)
            raster[x + y * w] = (fragment){1.0/0.0, {0, NULL}};

    point vpoints[num_vertices];
    for (int v = 0; v < num_vertices; v++)
        vpoints[v] = vertex_shader(vertex_uniform, vertices[v]);

    point a, b, c;
    for (int t = 0; t < num_triangles; t++) {
        a = vpoints[triangles[t * 3]];
        b = vpoints[triangles[t * 3 + 1]];
        c = vpoints[triangles[t * 3 + 2]];
        rasterize(w, h, raster, a, b, c);
    }

    for (int v = 0; v < num_vertices; v++)
        free(vpoints[v].val.vals);

    fragment f;
    rgba col;
    for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++) {
            f = raster[x + y * w];
            if (f.varying.vals == NULL)
                col = (rgba){0, 0, 0, 0xFF};
            else {
                col = fragment_shader(fragment_uniform, f.varying);
                free(f.varying.vals);
            }
            img[(h - y - 1) * w + x] =
                col.a << 24 | col.b << 16 |
                col.g << 8  | col.r;
        }
    }
    free(raster);
}


#define HYPO_LEN 0.5773502691896257
#define LEG_LEN 0.28867513459481287
#define HEIGHT 0.816496580927726
#define A(a, b, c, d, e, f, g) {7, (double[]){a, b, c, d, e, f, g}}
const int tetra_num_vertices = 4;
const vec tetra_vertices[4] = {
    A(   0,      0, -HYPO_LEN, 1, 0, 0, 0),
    A(-0.5,      0,   LEG_LEN, 0, 1, 0, 0),
    A( 0.5,      0,   LEG_LEN, 0, 0, 1, 0),
    A(   0, HEIGHT,         0, 0, 0, 0, 1),
};
const int tetra_num_triangles = 4;
const int tetra_triangles[4 * 3] = {
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


static void sphere_coords(double theta, double phi, double *vert) {
    double theta_ = 2 * M_PI * theta, phi_ = 2 * M_PI * phi;
    vert[0] = sin(theta_) * cos(phi_) / 2;
    vert[1] = sin(theta_) * sin(phi_) / 2;
    vert[2] = cos(theta_) / 2;
}

static void sphere_vertices_triangles(
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

// note to self: CAREFUL WITH THE NON-PARENTHESIZED ARGUMENTS
// useful with "-normal", but...
#define DOT(a, b) (a[0] * b[0] + a[1] * b[1] + a[2] * b[2])
#define T(e) {({int _ = 0; e;}), ({int _ = 1; e;}), ({int _ = 2; e;})}

void normalize(double v[3]) {
    double mag = sqrt(DOT(v, v));
    v[0] /= mag;
    v[1] /= mag;
    v[2] /= mag;
}

const double k_s = 0.4, k_d = 0.4, k_a = 0.2;
const int a = 16;
double phong(const double v[3], const double normal[3], const double s[3]) {
    double l_s[3] = T(s[_] - v[_]);
    normalize(l_s);
    double f = DOT(l_s, normal),
           r_s[3] = T(2 * f * normal[_] - l_s[_]);
    normalize(r_s);
    double intensity = k_a + k_d * DOT(l_s, normal) +
        k_s * pow(DOT(r_s, -normal), a);
    return intensity;
}

static uint *sphere_tex;
static int stw, sth;
const double m[3] = {-0.8, -0.1, -1},
      n[3] = {1, 1, -0.2};
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


const int dna_num_vertices = 19375;
vec dna_vertices[19375];
double dna_vertex_dat[155000];
const int dna_num_triangles = 31799;
int dna_triangles[31799 * 3];

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


const int w = 720, h = 720;

const int res_t = 15, res_p = 15;

static int sphere_init = 0;
static int sphere_num_vertices;
static vec *sphere_vertices;
static int sphere_num_triangles;
static int *sphere_triangles;

static int dna_init = 0;

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

void render_main(vec vertex_uniform, vec fragment_uniform, uint *img) {
    if (!dna_init) {
        FILE *f = fopen("dna_vertices", "r");
        fread(dna_vertex_dat, 1, sizeof(dna_vertex_dat), f);
        fclose(f);
        f = fopen("dna_triangles", "r");
        fread(dna_triangles, 1, sizeof(dna_triangles), f);
        fclose(f);
        for (int i = 0; i < 19375; i++) {
            dna_vertices[i] = (vec){8, dna_vertex_dat + i * 8};
        }
        dna_init = 1;
    }
    render(dna_vshader, dna_fshader,
            vertex_uniform, fragment_uniform,
            dna_num_vertices, dna_vertices,
            dna_num_triangles, dna_triangles,
            w, h, img);
    /*
    if (!sphere_init) {
        sphere_num_vertices = (res_t - 1) * res_p;
        sphere_num_triangles = (res_t - 2) * res_p * 2;
        sphere_vertices_triangles(
                res_t, res_p, &sphere_vertices, &sphere_triangles);
        FILE *f = fopen("gear.pam", "r");
        sphere_tex = load_pam(f, &stw, &sth);
        fclose(f);
        sphere_init = 1;
    }
    render(sphere_vshader, sphere_fshader,
            vertex_uniform, fragment_uniform,
            sphere_num_vertices, sphere_vertices,
            sphere_num_triangles, sphere_triangles,
            w, h, img);
            */
    /*
    render(tetra_vshader, tetra_fshader,
            vertex_uniform, fragment_uniform,
            tetra_num_vertices, tetra_vertices,
            tetra_num_triangles, tetra_triangles,
            w, h, img);
            */
}

const uint frames = 0;
int main(int argc, const char *argv[]) {
    uint *img = calloc(w * h, sizeof(uint));
    FILE *f = stdout; //fopen("/dev/null", "w");
    for (int t = 0; !frames || t < frames; t++) {
        vec uni = {1, (double[]){t / 30.}};
        render_main(uni, uni, img);
        fwrite(img, w * h, sizeof(uint), f);
    }
    fclose(f);
    free(img);
    return 0;
}

