// Copyright 2014-2020 Richard Copley

// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at

//     http://www.apache.org/licenses/LICENSE-2.0

// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

precision mediump float;
varying vec3 X, N, H;

// Lighting.
uniform vec3 l[4], d[4];     // Position, diffuse reflectance
const float s = .6, e = 32.; // Specular reflectance and exponent

// Line-width parameters.
const float a0 = 0.002; // Inner distance (100% black)
const float a1 = 0.004; // Outer distance (0% black)
const float m0 = 1. / (a1 - a0);
const float c0 = -a0 * m0;

void main() {
  vec3 c, n = normalize(N);
  float xi = inversesqrt(dot(X, X)), alpha = 0.9;
  if (!gl_FrontFacing) {
    alpha = 1.;
    n = -n;
  }
  // There are four lights.
  for (int i = 0; i != 4; ++i) {
    vec3 L = normalize(l[i] - X);
    float k = dot(L, n);
    float j = max(0., dot(L - 2. * k * n, X) * xi);
    c += max(0., k) * d[i] + pow(j, e) * s;
  }
  // Line shading.
  c *= clamp(c0 + m0 * abs(H[0]), 0., 1.);
  float s = clamp(c0 + m0 * abs(H[1]), 0., 1.);
  c[1] = min(c[1], s);
  c[2] = min(c[2], s);
  // Workaround for bad implementations.
  c = min(c, vec3(.9, .9, .9));
  gl_FragColor = vec4(c, alpha);
}
