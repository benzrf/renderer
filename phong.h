#include <math.h>

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

const double m[3] = {-0.8, -0.1, -1}, n[3] = {1, 1, -0.2};

