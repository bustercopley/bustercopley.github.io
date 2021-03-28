// Copyright 2014-2019 Richard Copley

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

// Lights.
uniform vec3 a, l [4], d [4], s [4];
uniform float e [4];

float z;
vec3 v, n;

float amplify (float d)
{
  d = clamp (1.25 - 128. * d, 0., 1.);
  return (4./3.) * (exp2 (-2. * d * d) - .25);
}

vec3 light (vec3 L, vec3 d, vec3 s, float e) {
  vec3 l = normalize (L - X);
  float k = dot (n, l);
  vec3 r = 2. * k * n - l;
  float j = max (0., dot (r, v));
  return max (0., k) * d + pow (j, e) * s;
}

void main ()
{
  z = amplify (min (H [0], min (H [1], H [2])));
  vec4 c;
  if (gl_FrontFacing) {
    v = normalize (-X);
    n = normalize (N);

    // There are four lights.
    c = vec4 (z * (a
                   + light (l [0], d [0], s [0], e [0])
                   + light (l [1], d [1], s [1], e [1])
                   + light (l [2], d [2], s [2], e [2])
                   + light (l [3], d [3], s [3], e [3])), .95);
  }

  if (! gl_FrontFacing) {
    v = normalize (-X);
    n = normalize (-N);
    vec3 d = vec3 (.3);
    vec3 s = vec3 (.7);
    c = vec4 (z * (a
                   + light (l [0], d, s, e [0])
                   + light (l [1], d, s, e [1])
                   + light (l [2], d, s, e [2])
                   + light (l [3], d, s, e [3])), 1.);
  }
  gl_FragColor = c;
}
