#version 300 es
in vec2 u;  // Texture coordinate
out vec2 U; // Texture coordinate
out vec2 q; // Flattened cone coordinate
out vec4 N; // World space normal
out mat2 g; // Pullback of screen metric to texture space
out mat2 h; // Pullback of screen metric to flattened cone
out float scale;
out float angle;
out float alpha;
uniform mat4 projection;
uniform mat4 modelview;
uniform float width, height;
uniform float theta;
const bool flattened = false;

#define PI 3.14159265358979323844
#define RATIO (5.0 / 3.0)

void main() {
  U = u;
  alpha = 0.0;
  scale = min(width, height);
  float x0 = width / scale;
  float y0 = height / scale;

  vec2 v = vec2(2.0 * PI * u.x, RATIO * u.y); // Angle, radial distance
  float cth = cos(theta);
  float sth = sin(theta);
  float r = cth / sth;
  float s = cth * r;
  float sx = sin(v.x);
  float cx = cos(v.x);

  // Model space
  vec4 w = vec4(v.y * cth * sx, sth + s * (1.0 - v.y), v.y * cth * cx, 1.0);
  vec4 x = modelview * w;                         // World space
  vec4 y = projection * x;                        // Clip space
  q = v.y * vec2(cos(sth * v.x), sin(sth * v.x)); // Flattened space

  // Normal in world space
  N = normalize(modelview * vec4(sx * cth, sth, cx * cth, 0.0));

  angle = (2.0 * PI) * sth;
  float qinv = 1.0 / v.y;
  float qinvsq = qinv * qinv;
  mat2 dudq = mat2(-q.y * qinvsq / angle, q.x * qinv / RATIO,
                   q.x * qinvsq / angle, q.y * qinv / RATIO);

  mat2x4 dwdu = (RATIO * cth) *
    mat2x4((2.0 * PI) * (u.y * cx), 0.0, (-2.0 * PI) * (u.y * sx), 0.0,
           sx, -r, cx, 0.0);

  mat2x4 dydu = projection * modelview * dwdu;

  // Pushforward of texture coordinate frame to screen space
  float pdiv = 1.0 / y.w;
  mat2 dpdu = pdiv * mat2(x0 * (dydu[0].x - pdiv * y.x * dydu[0].w),
                          y0 * (dydu[0].y - pdiv * y.y * dydu[0].w),
                          x0 * (dydu[1].x - pdiv * y.x * dydu[1].w),
                          y0 * (dydu[1].y - pdiv * y.y * dydu[1].w));

  // Pushforward of flattened cone frame to screen space
  mat2 dpdq = dpdu * dudq;

  // Pullback of screen space metric to texture space
  g = mat2(dot(dpdu[0], dpdu[0]), dot(dpdu[0], dpdu[1]),
           dot(dpdu[1], dpdu[0]), dot(dpdu[1], dpdu[1]));

  // Pullback of screen space metric to flattened cone
  h = mat2(dot(dpdq[0], dpdq[0]), dot(dpdq[0], dpdq[1]),
           dot(dpdq[1], dpdq[0]), dot(dpdq[1], dpdq[1]));

  gl_Position = y;
  //gl_Position = vec4(2.0 * u.x - 1.0, 1.0 - 2.0 * u.y, 0.0, 1.0);
}
