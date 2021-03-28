attribute vec2 x;
attribute vec2 u;
varying vec2 U;

void main() {
  U = u;
  gl_Position = vec4(x.x, x.y, 0.0, 1.0);
}
