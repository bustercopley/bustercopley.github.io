// -*- coding: utf-8; -*-

// Copyright 2014-2019, Richard Copley <buster at buster dot me dot uk>

// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at

//     http://www.apache.org/licenses/LICENSE-2.0

// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

function Body(stdlib, foreign, heap) {
  "use asm";

  var f = stdlib.Math.fround;
  var sqrt = stdlib.Math.sqrt;

  var data = new stdlib.Float32Array(heap);

  // Data layout:
  // Byte 0:   4 vectors of 4 floats (64 bytes): scratch space
  // Byte 128: 7 vectors of 4 floats (112 bytes): primary body
  // Byte 256: 7 vectors of 4 floats (112 bytes): secondary body

  // Field byte offsets from start of body:
  var PX = 0;  // Linear position vector (4 floats, xyz0).
  var PU = 16; // Angular position vector (4 floats, xyz0).
  var PW = 32; // Angular velocity vector (4 floats, xyz0).
  var PM = 64; // Modelview matrix (16 floats, OpenGL).

  // fpoly, gpoly, hpoly are minimax cubic polynomials on [0, (pi/2)^2],
  // computed by the Remes algorithm, for the functions:
  //   f(x) = sin(s)/s,               equioscillating error +-5.82769815e-7;
  //   g(x) = (1-cos(s))/s^2,         equioscillating error +-5.88282774e-8;
  //   h(x) = (1-(s/2)*cot(s/2))/s^2, equioscillating error +-5.48014034e-9;
  // where s = sqrt(x)).

  function fpoly(x) {
    x = f(x);
    var c0 = f(9.99999225e-1);
    var c1 = f(-1.66656822e-1);
    var c2 = f(8.31325911e-3);
    var c3 = f(-1.85243567e-4);
    return f(c0 + f(x * f(c1 + f(x * f(c2 + f(x * c3))))));
  }

  function gpoly(x) {
    x = f(x);
    var c0 = f(4.99999911e-1);
    var c1 = f(-4.16656733e-2);
    var c2 = f(1.38686667e-3);
    var c3 = f(-2.34776126e-5);
    return f(c0 + f(x * f(c1 + f(x * f(c2 + f(x * c3))))));
  }

  function hpoly(x) {
    x = f(x);
    var c0 = f(8.33333284e-2);
    var c1 = f(1.38897949e-3);
    var c2 = f(3.28886745e-5);
    var c3 = f(9.39777692e-7);
    return f(c0 + f(x * f(c1 + f(x * f(c2 + f(x * c3))))));
  }

  // Dot product of two three-dimensional vectors.
  function dot(pa, pb) {
    pa = pa | 0;
    pb = pb | 0;
    var a = f(0);
    a = f(data[pa + 0 >> 2] * data[pb + 0 >> 2]);
    a = f(a + f(data[pa + 4 >> 2] * data[pb + 4 >> 2]));
    a = f(a + f(data[pa + 8 >> 2] * data[pb + 8 >> 2]));
    return a;
  }

  // Compute the OpenGL modelview matrix from the position x and orientation u.
  function matrix(byteOffset) {
    byteOffset = byteOffset | 0;
    var xsq = f(0);
    var x = f(0);
    var x1 = f(0);
    var x1sq = f(0);
    var s = f(0);
    var c = f(0);
    var u0 = f(0);
    var u1 = f(0);
    var u2 = f(0);
    var u0sq = f(0);
    var u1sq = f(0);
    var u2sq = f(0);
    var symm0 = f(0);
    var symm1 = f(0);
    var symm2 = f(0);
    var skew0 = f(0);
    var skew1 = f(0);
    var skew2 = f(0);
    var px = 0;
    var pu = 0;
    var pm = 0;

    px = (byteOffset + PX) | 0;
    pu = (byteOffset + PU) | 0;
    pm = (byteOffset + PM) | 0;

    u0 = f(data[pu + 0 >> 2]);
    u1 = f(data[pu + 4 >> 2]);
    u2 = f(data[pu + 8 >> 2]);

    u0sq = f(u0 * u0);
    u1sq = f(u1 * u1);
    u2sq = f(u2 * u2);

    // Calculate s = sin(x) / x, c = (1 - cos(x)) / xsq.
    xsq = f(f(u0sq + u1sq) + u2sq);
    if (xsq < f(2.46740110)) {
      // Quadrant 1 (0 <= x < pi/2).
      s = f(fpoly(xsq)); // sin(x) / x
      c = f(gpoly(xsq)); // (1 - cos(x)) / xsq
    }
    else {
      // Quadrants 2, 3 and 4 (pi/2 <= x < 2pi).
      x = f(sqrt(xsq));
      x1 = f(f(x - f(3.14159274)) + f(8.74227766e-8)); // x - pi
      x1sq = f(x1 * x1);
      // fpoly(x1sq) = sin(x1) / x1 = -sin(x) / x1
      // gpoly(x1sq) = (1 - cos(x1)) / x1sq = (1 + cos(x)) / x1sq
      s = f(f(f(-x1) * f(fpoly(x1sq))) / x);           // sin(x) / x
      c = f(f(f(2) - f(x1sq * f(gpoly(x1sq)))) / xsq); // (1 - cos(x)) / xsq
    }

    skew0 = f(s * u0);
    skew1 = f(s * u1);
    skew2 = f(s * u2);

    symm0 = f(c * f(u1 * u2));
    symm1 = f(c * f(u2 * u0));
    symm2 = f(c * f(u0 * u1));

    data[pm + 0 >> 2] = f(f(1) - f(c * f(u1sq + u2sq)));
    data[pm + 4 >> 2] = f(symm2 + skew2);
    data[pm + 8 >> 2] = f(symm1 - skew1);
    data[pm + 12 >> 2] = f(0);

    data[pm + 16 >> 2] = f(symm2 - skew2);
    data[pm + 20 >> 2] = f(f(1) - f(c * f(u2sq + u0sq)));
    data[pm + 24 >> 2] = f(symm0 + skew0);
    data[pm + 28 >> 2] = f(0);

    data[pm + 32 >> 2] = f(symm1 + skew1);
    data[pm + 36 >> 2] = f(symm0 - skew0);
    data[pm + 40 >> 2] = f(f(1) - f(c * f(u0sq + u1sq)));
    data[pm + 44 >> 2] = f(0);

    data[pm + 48 >> 2] = data[px + 0 >> 2];
    data[pm + 52 >> 2] = data[px + 4 >> 2];
    data[pm + 56 >> 2] = data[px + 8 >> 2];
    data[pm + 60 >> 2] = f(1);
  }

  // Given vectors u and w in R³, calculate the total derivative of u(t)
  // w.r.t. t at t=0, where u(t) satisfies f(u(t)) = f(wt) f(u),
  // where f: R³ → SO(R³) is the Rodrigues function f(v) = exp(φ(v)),
  // where φ: R³ → so(R³) is the intertwinor from the identity
  // representation of the Lie algebra so(R³) to the adjoint representation.
  function tangent(p, u0, u1, u2, w0, w1, w2) {
    p = p | 0;
    u0 = f(u0);
    u1 = f(u1);
    u2 = f(u2);
    w0 = f(w0);
    w1 = f(w1);
    w2 = f(w2);
    var xsq = f(0);
    var x = f(0);
    var x1 = f(0);
    var x1sq = f(0);
    var g = f(0);
    var h = f(0);
    var k = f(0);
    var uwh = f(0);
    var c = f(0);
    var s = f(0);

    // Compute g = (x/2)cot(x/2) and h = (1-g)/x^2.
    xsq = f(f(f(u0 * u0) + f(u1 * u1)) + f(u2 * u2));
    if (xsq < f(2.46740110)) {
      // Quadrant 1 (0 <= x < pi/2).
      h = f(hpoly(xsq));
      g = f(f(1) - f(xsq * h));
    }
    else {
      // Quadrants 2, 3 and 4 (pi/2 <= x < 2pi).
      x = f(sqrt(xsq));
      // (pi/2) - (x/2)
      x1 = f(f(-0.5) * f(f(x - f(3.14159274)) + f(8.74227766e-8)));
      x1sq = f(x1 * x1);
      // fpoly(x1sq) = sin(x1) / x1 = cos(x/2) / x1
      // gpoly(x1sq) = (1 - cos(x1)) / x1sq = (1 - sin(x/2)) / x1sq
      c = f(x1 * f(fpoly(x1sq)));             // cos(x/2)
      s = f(f(1) - f(x1sq * f(gpoly(x1sq)))); // sin(x/2)
      g = f(f(0.5) * f(f(x * c) / s));        // (x/2)cot(x/2)
      h = f(f(f(1) - g) / xsq);
    }

    // <u, w>h
    uwh = f(h * f(f(f(u0 * w0) + f(u1 * w1)) + f(u2 * w2)));

    // u'(t) = <u,w>hu + gw - 0.5[u,w]
    data[p + 0 >> 2] =
      f(f(f(uwh * u0) + f(g * w0)) - f(f(0.5) * f(f(u1 * w2) - f(u2 * w1))));
    data[p + 4 >> 2] =
      f(f(f(uwh * u1) + f(g * w1)) - f(f(0.5) * f(f(u2 * w0) - f(u0 * w2))));
    data[p + 8 >> 2] =
      f(f(f(uwh * u2) + f(g * w2)) - f(f(0.5) * f(f(u0 * w1) - f(u1 * w0))));
  }

  // Return false if |w*dt| < 1.0e-6.
  // Otherwise, set u <- unhat(log(exp(hat(w*dt))exp(hat(u)))).
  function bch4(byteOffset, dt) {
    // One step of the classical fourth-order Runge-Kutta method.
    byteOffset = byteOffset | 0;
    dt = f(dt);
    var pa = 0;
    var pb = 0;
    var pc = 0;
    var pd = 0;
    var pu = 0;
    var pw = 0;
    var xsq = f(0);
    var wdt_sq = f(0);
    var k = f(0);
    var u0 = f(0);
    var u1 = f(0);
    var u2 = f(0);
    var v0 = f(0);
    var v1 = f(0);
    var v2 = f(0);
    var wdt0 = f(0);
    var wdt1 = f(0);
    var wdt2 = f(0);

    pa = 0;
    pb = 16;
    pc = 32;
    pd = 48;
    pu = (byteOffset + PU) | 0;
    pw = (byteOffset + PW) | 0;

    wdt0 = f(data[pw + 0 >> 2] * dt);
    wdt1 = f(data[pw + 4 >> 2] * dt);
    wdt2 = f(data[pw + 8 >> 2] * dt);

    wdt_sq = f(f(f(wdt0 * wdt0) + f(wdt1 * wdt1)) + f(wdt2 * wdt2));
    if (wdt_sq < f(1.0e-12)) { return 0; }

    u0 = f(data[pu + 0 >> 2]);
    u1 = f(data[pu + 4 >> 2]);
    u2 = f(data[pu + 8 >> 2]);
    tangent(pa, u0, u1, u2, wdt0, wdt1, wdt2);

    v0 = f(u0 + f(f(0.5) * data[pa + 0 >> 2]));
    v1 = f(u1 + f(f(0.5) * data[pa + 4 >> 2]));
    v2 = f(u2 + f(f(0.5) * data[pa + 8 >> 2]));
    tangent(pb, v0, v1, v2, wdt0, wdt1, wdt2);

    v0 = f(u0 + f(f(0.5) * data[pb + 0 >> 2]));
    v1 = f(u1 + f(f(0.5) * data[pb + 4 >> 2]));
    v2 = f(u2 + f(f(0.5) * data[pb + 8 >> 2]));
    tangent(pc, v0, v1, v2, wdt0, wdt1, wdt2);

    v0 = f(u0 + data[pc + 0 >> 2]);
    v1 = f(u1 + data[pc + 4 >> 2]);
    v2 = f(u2 + data[pc + 8 >> 2]);
    tangent(pd, v0, v1, v2, wdt0, wdt1, wdt2);

    // u <- u + (1/6)a + (2/6)b + (2/6)c + (1/6)d.
    v0 = f(data[pb + 0 >> 2] + data[pc + 0 >> 2]);
    v1 = f(data[pb + 4 >> 2] + data[pc + 4 >> 2]);
    v2 = f(data[pb + 8 >> 2] + data[pc + 8 >> 2]);

    v0 = f(f(v0 + v0) + f(data[pa + 0 >> 2] + data[pd + 0 >> 2]));
    v1 = f(f(v1 + v1) + f(data[pa + 4 >> 2] + data[pd + 4 >> 2]));
    v2 = f(f(v2 + v2) + f(data[pa + 8 >> 2] + data[pd + 8 >> 2]));

    u0 = f(u0 + f(f(0.16666667) * v0));
    u1 = f(u1 + f(f(0.16666667) * v1));
    u2 = f(u2 + f(f(0.16666667) * v2));

    // If |u| is above a threshold strictly between pi and 2*pi (here sqrt(10)),
    // replace u with a smaller element that represents the same rotation.
    xsq = f(f(f(u0 * u0) + f(u1 * u1)) + f(u2 * u2));
    if (xsq > f(10)) {
      k = f(f(1) - f(f(6.28318531) / f(sqrt(xsq))));
      u0 = f(k * u0);
      u1 = f(k * u1);
      u2 = f(k * u2);
    }

    data[pu + 0 >> 2] = u0;
    data[pu + 4 >> 2] = u1;
    data[pu + 8 >> 2] = u2;

    return 1;
  }

  return { advance : bch4, matrix : matrix };
}
