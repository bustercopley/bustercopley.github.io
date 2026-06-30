#version 300 es
precision mediump float;
in vec2 U; // Texture coordinate
in vec2 q; // Flattened cone coordinate
in vec4 N; // World space normal
in mat2 g; // Pullback of screen metric to texture space
in mat2 h; // Pullback of screen metric to flattened cone
in float scale;
in float angle;
in float alpha;
out vec4 color;
uniform vec4 lightDirection;
uniform float phi;

#define PI 3.14159265358979323844
#define RATIO (5.0 / 3.0)

float deth;
float detg;

// Colors
vec4 Z = vec4(1.0, 1.0, 1.0, 0.5); // Background
vec4 A = vec4(0.0, 0.7, 0.0, 0.8); // Horizontal lines
vec4 B = vec4(0.0, 0.0, 1.0, 0.8); // Vertical lines
vec4 C = vec4(0.8, 0.0, 0.0, 0.8); // Radial lines and tangent circle
vec4 D = vec4(0.1, 0.1, 0.1, 0.9); // Background arrow
vec4 F = vec4(0.0, 0.0, 0.0, 1.0); // Foreground arrow

// Apply metric tensor to tangent vectors
float hdot(highp vec2 u, highp vec2 v) {
  return dot(u, h * v);
}

// Apply volume form (squared)
highp float hvolsq(vec2 u, vec2 v) {
  highp float a = u.x * v.y - u.y * v.x;
  return deth * (a * a);
}

// Line segment helper
float segment(vec2 u, vec2 v, vec2 w) {
  return max(hvolsq(u, v) / hdot(u, u),
             scale * max(hdot(u, -v), hdot(u, w)));
}

// Screen distance from arrow figure in flattened-cone space
float arrow(vec2 A, vec2 dir) {
  vec2 u = (0.27 * scale) * vec2(dir.x, dir.y);
  vec2 u1 = (0.06 * scale) * vec2(dir.x - dir.y, dir.y + dir.x);
  vec2 u2 = (0.06 * scale) * vec2(dir.x + dir.y, dir.y - dir.x);

  vec2 v = scale * (q - A);
  vec2 w = v - u;
  vec2 v1 = w + u1;
  vec2 v2 = w + u2;

  return min(min(min(hdot(v, v),            // Disk center A
                     hdot(w, w)),           // Disk center B = A + u
                 min(hdot(v1, v1),          // Disk center C = B + u1
                     hdot(v2, v2))),        // Disk center D = B + u2
             min(segment(u, v, w),          // Line segment AB
                 min(segment(u1, v1, w),    // Line segment BC
                     segment(u2, v2, w)))); // Line segment CD
}

// Smoothed characteristic function of a figure
float chi(float dsq, float width) {
  return clamp(width - 0.06 * dsq, 0.0, 1.0);
}

void main() {
  detg = g[0][0] * g[1][1] - g[0][1] * g[1][0];
  deth = h[0][0] * h[1][1] - h[0][1] * h[1][0];

  // Horizontal distance from nearest meridian, in texture space
  float ax = scale * (U.x - round(12.0 * U.x) / 12.0);
  // Vertical distance from tangent circle, in texture space
  float ay = scale * (U.y - (1.0 / RATIO));
  // Ruled lines in flattened cone space
  float bx = scale * (q.x - round(6.0 * q.x) / 6.0);
  float by = scale * (q.y - round(6.0 * q.y) / 6.0);

  // Perpendicular distances from the lines in screen space (to first order)
  float aX = (ax * ax) * (detg / g[1][1]);
  float aY = (ay * ay) * (detg / g[0][0]);
  float bX = (bx * bx) * (deth / h[1][1]);
  float bY = (by * by) * (deth / h[0][0]);

  // Two nearest background arrows, plus the one that peeps over the cut
  float S = floor(12.0 * U.x) * (angle / 12.0);
  float T = ceil(12.0 * U.x) * (angle / 12.0);
  float R = min(min(arrow(vec2(cos(S), sin(S)), vec2(0.0, 1.0)),
                    arrow(vec2(cos(T), sin(T)), vec2(0.0, 1.0))),
                arrow(vec2(1.0, 0.0), vec2(sin(angle), cos(angle))));

  // Foreground arrow
  float Q = arrow(vec2(cos(phi), sin(phi)), vec2(0.0, 1.0));

  // Overlays
  float a = chi(bX, 0.85);           // Horizontal lines
  float b = chi(bY, 0.85);           // Vertical lines
  float c = chi(min(aX, aY), 1.25);  // Radial lines and tangent circle
  float d = chi(R, 3.5);             // Background arrow stroke
  float e = chi(R, 1.25);            // Background arrow fill
  float f = chi(Q, 1.25);            // Foreground arrow

  // Lighting
  float l = max(0.0, min(1.0, -dot(N, lightDirection)));

  // Compositing
  vec4 X = mix(mix(mix(Z, A, a), B, b), C, c);
  vec4 Y = mix(mix(X, D, d), X, e);
  color = mix(vec4(l * Y.xyz, max(alpha, Y.w)), F, f);
}
