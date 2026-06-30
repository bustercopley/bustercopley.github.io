#version 300 es
precision mediump float;
in vec2 U;              // Texture coordinate
in vec2 p;              // Screen coordinate
in vec4 N;              // World space normal
in float g00, g01, g11; // Pullback of screen space metric to texture space
in float scale;
out vec4 color;
uniform sampler2D sampler;
uniform vec4 lightDirection;

void main() {
  // Overlays
  // Horizontal distance from nearest meridian, in texture space
  float x = scale * (U.x - round(12.0 * U.x) / 12.0);
  // Vertical distance from nearest latitude line, in texture space
  float y = scale * (U.y - round(max(1.0, min(7.0, 8.0 * U.y))) / 8.0);
  // Perpendicular distances from the lines in screen space (to first order)
  float sx = (x * x) * (g00 - g01 * g01 / g11);
  float sy = (y * y) * (g11 - g01 * g01 / g00);

  // Smoothed characteristic function of the union of the overlays
  float e = clamp(1.0 - 0.25 * sqrt(min(sx, sy)), 0.0, 1.0);

  // Lighting
  float l = max(0.0, min(1.0, -dot(N, lightDirection)));
  color = vec4(mix(texture(sampler, U).xyz, vec3(0.3, 0.8, 1.0), e) * l, 1.0);
}
