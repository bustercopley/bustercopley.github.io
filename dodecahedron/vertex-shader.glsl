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

attribute vec3 x, n, h;
varying vec3 X, N, H;
uniform mat4 m, p;

void main() {
  H = h;
  N = (m * vec4(n, 0)).xyz;
  vec4 P = m * vec4(x, 1);
  X = P.xyz;
  gl_Position = p * P;
}
