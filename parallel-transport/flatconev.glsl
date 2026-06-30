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

#define PI 3.14159265358979323844
#define RATIO (5.0 / 3.0)
#define SCALE 0.25

void main() {
  U = u;
  alpha = 1.0;
  scale = SCALE * min(width, height);
  float x0 = scale / width;
  float y0 = scale / height;
  angle = (2.0 * PI) * sin(theta);

  vec2 v = vec2(angle * u.x, RATIO * u.y); // Angle, radial distance
  q = v.y * vec2(cos(v.x), sin(v.x));      // Vertex in flattened cone space
  N = vec4(0.0, 0.0, 10.0, 0.0);           // Normal in world space

  // Pullback of screen space metric to texture space
  g = mat2((angle * angle) * (v.y * v.y), 0.0, 0.0, RATIO * RATIO);

  // Pullback of screen space metric to flattened cone space
  h = mat2(1.0, 0.0, 0.0, 1.0);

  gl_Position = vec4((q.x - 1.5) * x0, (q.y + 1.5) * y0, 0.0, 1.0);
}
