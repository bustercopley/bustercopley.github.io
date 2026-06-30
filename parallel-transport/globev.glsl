#version 300 es
in vec2 u;  // Texture coordinate
out vec2 U; // Texture coordinate
out vec2 p; // Screen coordinate
out vec4 N; // World space normal
out float g00, g01, g11; // Pullback of screen space metric to texture space
out float scale;
uniform mat4 projection;
uniform mat4 modelview;
uniform float width, height;

#define PI 3.14159265358979323844

void main() {
  scale = min(width, height);
  float x0 = width / scale;
  float y0 = height / scale;

  U = u;                                   // Texture coordinate
  vec2 v = vec2(2.0 * PI * u.x, PI * u.y); // Longitude, colatitude (radians)
  vec2 sinv = sin(v);
  vec2 cosv = cos(v);

  vec4 w = vec4(sinv.y * sinv.x, cosv.y, sinv.y * cosv.x, 1.0); // Model space
  vec4 x = modelview * w;                                       // World space
  vec4 y = projection * x;                                      // Clip space
  float pdiv = 1.0 / y.w;
  p = 0.5 * pdiv * vec2(width * y.x, height * y.y);             // Screen space

  // Normal in world space
  N = normalize(x - modelview[3]);

  mat2x4 dwdu = mat2x4(
    (2.0 * PI) * sinv.y * cosv.x, 0.0, (-2.0 * PI) * sinv.y * sinv.x, 0.0,
    PI * cosv.y * sinv.x, -PI * sinv.y, PI * cosv.y * cosv.x, 0.0
  );

  mat2x4 dydu = projection * modelview * dwdu;

  // Pushforward of texture coordinate frame to screen space
  mat2 dpdu = pdiv * mat2(x0 * (dydu[0].x - pdiv * y.x * dydu[0].w),
                          y0 * (dydu[0].y - pdiv * y.y * dydu[0].w),
                          x0 * (dydu[1].x - pdiv * y.x * dydu[1].w),
                          y0 * (dydu[1].y - pdiv * y.y * dydu[1].w));

  // Pullback of screen space metric to texture space
  g00 = dot(dpdu[0], dpdu[0]);
  g01 = dot(dpdu[0], dpdu[1]);
  g11 = dot(dpdu[1], dpdu[1]);

  gl_Position = y;
  //gl_Position = vec4(2.0 * u.x - 1.0, 1.0 - 2.0 * u.y, 0.0, 1.0);
}
