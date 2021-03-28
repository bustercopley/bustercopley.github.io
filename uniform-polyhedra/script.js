// -*- coding: utf-8; -*-

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

function Body(stdlib, foreign, heap) {
  "use asm";

  // asm.js obliges us to alias these stdlib functions.
  var f = stdlib.Math.fround;
  var sqrt = stdlib.Math.sqrt;

  var data = new stdlib.Float32Array(heap);

  // Data layout:
  // Byte 0:   4 vectors of 4 floats (64 bytes): scratch space
  // Byte 128: 7 vectors of 4 floats (112 bytes): primary body
  // Byte 256: 7 vectors of 4 floats (112 bytes): secondary body

  // Field byte offsets from start of body:
  var PX = 0;   // Linear position vector (4 floats, xyz0).
  var PU = 16;  // Angular position vector (4 floats, xyz0).
  var PW = 32;  // Angular velocity vector (4 floats, xyz0).
  var PM = 64;  // Modelview matrix (16 floats, OpenGL).

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
    pa = pa|0;
    pb = pb|0;
    var a = f (0);
    a = f(data[pa + 0 >> 2] * data[pb + 0 >> 2]);
    a = f(a + f(data[pa + 4 >> 2] * data[pb + 4 >> 2]));
    a = f(a + f(data[pa + 8 >> 2] * data[pb + 8 >> 2]));
    return a;
  }

  // Compute the OpenGL modelview matrix from the position x and orientation u.
  function matrix(byteOffset) {
    byteOffset = byteOffset|0;
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

    px = (byteOffset + PX)|0;
    pu = (byteOffset + PU)|0;
    pm = (byteOffset + PM)|0;

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
      s = f(fpoly(xsq));  // sin(x) / x
      c = f(gpoly(xsq));  // (1 - cos(x)) / xsq
    }
    else {
      // Quadrants 2, 3 and 4 (pi/2 <= x < 2pi).
      x = f(sqrt(xsq));
      x1 = f(f(x - f(3.14159274)) + f(8.74227766e-8));  // x - pi
      x1sq = f(x1 * x1);
      // fpoly(x1sq) = sin(x1) / x1 = -sin(x) / x1
      // gpoly(x1sq) = (1 - cos(x1)) / x1sq = (1 + cos(x)) / x1sq
      s = f(f(f(-x1) * f(fpoly(x1sq))) / x);            // sin(x) / x
      c = f(f(f(2) - f(x1sq * f(gpoly(x1sq)))) / xsq);  // (1 - cos(x)) / xsq
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
    p = p|0;
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
      c = f(x1 * f(fpoly(x1sq)));                 // cos(x/2)
      s = f(f(1) - f(x1sq * f(gpoly(x1sq))));     // sin(x/2)
      g = f(f(0.5) * f(f(x * c) / s));            // (x/2)cot(x/2)
      h = f(f(f(1) - g) / xsq);
    }

    // <u, w>h
    uwh = f(h * f(f(f(u0 * w0) + f(u1 * w1)) + f(u2 * w2)));

    // u'(t) = <u,w>hu + gw - 0.5[u,w]
    data[p + 0 >> 2] = f(f(f(uwh * u0) + f(g * w0)) - f(f(0.5) * f(f(u1 * w2) - f(u2 * w1))));
    data[p + 4 >> 2] = f(f(f(uwh * u1) + f(g * w1)) - f(f(0.5) * f(f(u2 * w0) - f(u0 * w2))));
    data[p + 8 >> 2] = f(f(f(uwh * u2) + f(g * w2)) - f(f(0.5) * f(f(u0 * w1) - f(u1 * w0))));
  }

  // Return false if |w*dt| < 1.0e-6.
  // Otherwise, set u <- unhat(log(exp(hat(w*dt))exp(hat(u)))).
  function bch4(byteOffset, dt) {
    // One step of the classical fourth-order Runge-Kutta method.
    byteOffset = byteOffset|0;
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
    pu = (byteOffset + PU)|0;
    pw = (byteOffset + PW)|0;

    wdt0 = f(data[pw + 0 >> 2] * dt);
    wdt1 = f(data[pw + 4 >> 2] * dt);
    wdt2 = f(data[pw + 8 >> 2] * dt);

    wdt_sq = f(f(f(wdt0 * wdt0) + f(wdt1 * wdt1)) + f(wdt2 * wdt2));
    if (wdt_sq < f(1.0e-12)) {
      return 0;
    }

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

  return { advance: bch4, matrix: matrix };
}

(function () {
  "use strict";

  // Create asm.js Body object.
  var heap = new ArrayBuffer(0x10000);  // 64 kilobytes (the minimum).
  var body = Body(window, null, heap);

  // Bodies (byte offsets from start of heap):
  var MAIN = 128; // Main body.

  // Fields (byte offsets from start of body):
  var PX = 0;   // Linear position.
  var PU = 16;  // Angular position.
  var PW = 32;  // Angular velocity.
  var PM = 64;  // Modelview matrix.

  // Views into the main body.
  var x = new Float32Array(heap, MAIN + PX, 4);  // Centre of sphere
  var u = new Float32Array(heap, MAIN + PU, 4);  // Angular position
  var w = new Float32Array(heap, MAIN + PW, 4);  // Angular momentum
  var m = new Float32Array(heap, MAIN + PM, 16); // Modelview matrix

  // DOM manipulation helper functions.
  function removeChildren(element) {
    while (element.firstChild) {
      element.removeChild(element.firstChild);
    }
  }

  function appendText(node, text) {
    var textNode = document.createTextNode(text);
    node.appendChild(textNode);
  }

  function appendPreText(node, pretext) {
    if (pretext) {
      var pre = document.createElement("pre");
      appendText(pre, pretext);
      node.appendChild(pre);
    }
  }

  function showError(text, pretext) {
    removeChildren(document.body);
    var div = document.createElement("div");
    appendText(div, text);
    appendPreText(div, pretext);
    document.body.appendChild(div);
  }

  // Viewing frustum dimensions. (The "t" stands for "tank".)
  var tz = 4;  // distance from the eye to the screen (i.e., to the front of the tank)
  var td = 3;  // distance from the screen to the rear clipping plane (the tank depth)
  var ts = 2;  // half of the shorter side of the front wall of the tank

  // Material colour and lighting parameters.
  var ambientReflection = [0.0, 0.0, 0.0];
  var lightPosition = [ 0, 0, -tz + 3,
                        -4, 3, -tz,
                        +4, 3, -tz,
                        0, -5, -tz ];
  var diffuseReflectance = [ 0.6, 0.6, 0.6,
                             0.4, 0.4, 0.4,
                             0.4, 0.4, 0.4,
                             0.4, 0.4, 0.4 ];
  var specularReflection = [ 0.6, 0.6, 0.6,
                             0.6, 0.6, 0.6,
                             0.6, 0.6, 0.6,
                             0.6, 0.6, 0.6 ];
  var specularExponent = [ 32, 32, 32, 32 ];

  // Operations on a three dimensional inner product space.
  function dot(u, v) {
    return u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
  }

  function cross(u, v) {
    var r = Array(3);
    r[0] = u[1] * v[2] - u[2] * v[1];
    r[1] = u[2] * v[0] - u[0] * v[2];
    r[2] = u[0] * v[1] - u[1] * v[0];
    return r;
  }

  function rotate(x) {
    var result = new Float32Array(3);
    result[0] = m[0] * x[0] + m[4] * x[1] + m[8] * x[2];
    result[1] = m[1] * x[0] + m[5] * x[1] + m[9] * x[2];
    result[2] = m[2] * x[0] + m[6] * x[1] + m[10] * x[2];
    return result;
  }

  function unrotate(x) {
    var result = new Float32Array(3);
    result[0] = m[0] * x[0] + m[1] * x[1] + m[2] * x[2];
    result[1] = m[4] * x[0] + m[5] * x[1] + m[6] * x[2];
    result[2] = m[8] * x[0] + m[9] * x[1] + m[10] * x[2];
    return result;
  }

  // Set the angular velocity of the body so that the direction vector anchor
  // (in model space) will point in direction target (in world space)
  // after time dt.
  function spinTo(anchor, target, dt) {
    var current = rotate(anchor);
    var crs = cross(current, target);
    // |crs| = sin(θ), where θ is the angle between current and target.
    var sinTheta = Math.sqrt(dot(crs, crs));
    // The required value is in the direction of crs and has magnitude θ/dt.
    var xsq = sinTheta * sinTheta;
    var thetaOverSinTheta = xsq < 1.0e-8 ? 1 - xsq / 6 : Math.asin(sinTheta) / sinTheta;
    var lambda = thetaOverSinTheta / dt;
    w[0] = lambda * crs[0];
    w[1] = lambda * crs[1];
    w[2] = lambda * crs[2];
    w[3] = 0.0;
  }

  function pick(p, rect) {
    // p is a pair of coordinates in client space (the mouse position).
    // rect is the bounding rectangle of the canvas, in client space.
    // X is a triple of coordinates in world space (the centre of the object).
    // Let L be the line in world space whose points all project to p in client space.
    // Let T be the point in world space where L intersects the unit sphere centre X,
    // or the nearer of the two points if the line intersects the sphere twice,
    // or null if they do not intersect.
    // Return the unit vector T - X (or null if T is null).

    var w = rect.right - rect.left;
    var h = rect.bottom - rect.top;
    var s0 = rect.left + 0.5 * w;
    var s1 = rect.top + 0.5 * h;
    var ss = ts / Math.min(w, h);

    // Find the point in world space where L intersects the near clipping plane.
    var q = [ss * (p[0] - s0), -ss * (p[1] - s1), -tz];

    // Since 0, Q and T are on the line L, assuming Q ≠ 0, there is
    // some real t such that T = 0 + t * Q.
    // Since T is on the unit sphere centre X, we have |T - X| = 1.
    // Combining these two facts we get a quadratic equation in t.
    var qsq = dot(q, q);
    var xsq = dot(x, x);
    var qx = dot(q, x);

    // The quadratic equation is qsq*t² - 2*qx*t + xsq-1 = 0.
    // Its solutions are t = (qx ± √(qx² - qsq*(xsq-1))) / qsq.
    // If the discriminant is negative then L doesn't intersect the sphere.
    var scrim = qx * qx - qsq * (xsq - 1);
    if (scrim < 0) {
      return null;
    }

    // Take the smaller of the two solutions to get the intersection closer to Q.
    var t = (qx - Math.sqrt(scrim)) / qsq;

    // Return T - X.
    return [t * q[0] - x[0], t * q[1] - x[1], t * q[2] - x[2]];
  }

  function randomVectorInBall(x, r) {
    var rsq = r * r;
    var xsq;
    do {
      x[0] = 2 * r * Math.random() - r;
      x[1] = 2 * r * Math.random() - r;
      x[2] = 2 * r * Math.random() - r;
      xsq = x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
    }
    while (xsq >= rsq);
  }

  function randomVectorOnSphere(x, r) {
    // Marsaglia (1972); see http://mathworld.wolfram.com/SpherePointPicking.html.
    // Get a vector in the unit disc.
    var xsq;
    do {
      x[0] = 2 * Math.random() - 1;
      x[1] = 2 * Math.random() - 1;
      xsq = x[0] * x[0] + x[1] * x[1];
    }
    while (xsq >= 1);
    // Apply Marsaglia's transformation.
    var t = 2 * r * Math.sqrt(1 - xsq);
    x[0] = t * x[0];
    x[1] = t * x[1];
    x[2] = r - 2 * r * xsq;
  }

  var DISABLE_REQUEST_CACHING = false;
  function download(url, type, handler, data) {
    var r = new XMLHttpRequest();
    r.onload = function () {
      if (this.status !== 200) {
        showError(this.status.toString() + " " +
                  this.statusText + ", requested \"" +
                  url + "\".");
      }
      else {
        handler(this.response, data);
      }
    };
    r.open("GET", url + (DISABLE_REQUEST_CACHING ? "?v=2" : ""), true);
    r.responseType = type;
    r.send();
  }

  var previousHue = -1;

  function randomColor(sat) {
    var color = new Float32Array(3);
    var h;
    var i;
    var j;
    var d;
    do {
      h = 6 * Math.random();
      d = Math.abs(h - previousHue);
      d = Math.min(d, 6 - d);
    }
    while (d < 1);
    previousHue = h;
    i = h > 4 ? h - 4 : h + 2;
    j = h > 2 ? h - 2 : h + 4;
    color[0] = (1 - sat) + sat * (h < 2 ? 1 : h < 3 ? 3 - h : h < 5 ? 0 : h - 5);
    color[1] = (1 - sat) + sat * (i < 2 ? 1 : i < 3 ? 3 - i : i < 5 ? 0 : i - 5);
    color[2] = (1 - sat) + sat * (j < 2 ? 1 : j < 3 ? 3 - j : j < 5 ? 0 : j - 5);
    return color;
  }

  // Graphics.
  var gl;
  var program;
  var vertexCount;

  function initializeGraphics(vshader, fshader) {
    if (!gl) {
      showError("WebGL context creation failed.");
      return false;
    }

    program = initShaderProgram(gl, vshader, fshader);
    if (!program) {
      return false;
    }

    gl.clearColor(1, 1, 1, 1);
    gl.enable(gl.CULL_FACE);
    gl.enable(gl.BLEND);
    gl.blendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA);

    gl.uniform3fv(gl.getUniformLocation(program, "l"), lightPosition);
    gl.uniform3fv(gl.getUniformLocation(program, "a"), ambientReflection);
    gl.uniform3fv(gl.getUniformLocation(program, "s"), specularReflection);
    gl.uniform1fv(gl.getUniformLocation(program, "e"), specularExponent);

    resize(true);
    return true;
  }

  function resize(force) {
    try {
      var ratio = devicePixelRatio;
      var rect = gl.canvas.getBoundingClientRect();
      var w = Math.floor(ratio * rect.right) - Math.floor(ratio * rect.left);
      var h = Math.floor(ratio * rect.bottom) - Math.floor(ratio * rect.top);

      if (force || gl.canvas.width !== w || gl.canvas.height !== w) {
        gl.canvas.width = w;
        gl.canvas.height = h;
        w = gl.drawingBufferWidth;
        h = gl.drawingBufferHeight;

        var scaleFactor = 0.5 / Math.min(w, h);
        var tw = ts * w * scaleFactor;
        var th = ts * h * scaleFactor;
        var tt = 2 * tz / td;
        var projectionMatrix = [
          tz / tw, 0, 0, 0,
          0, tz / th, 0, 0,
          0, 0, tt + 1, -1,
          0, 0, tz * (tt + 2), 0
        ];

        gl.viewport(0, 0, w, h);
        setMatrixUniform(gl, program, "p", projectionMatrix);
      }
    }
    catch (e) {
      showError(e.toString(), e.stack);
    }
  }

  function loadVertexBuffer(vbuffer) {
    var vertexBuffer = gl.createBuffer();
    var attributes = new Float32Array(vbuffer);
    gl.bindBuffer(gl.ARRAY_BUFFER, vertexBuffer);
    gl.bufferData(gl.ARRAY_BUFFER, attributes, gl.STATIC_DRAW);

    function attrib(name, offset) {
      var location = gl.getAttribLocation(program, name);
      gl.vertexAttribPointer(location, 3, gl.FLOAT, false, 0, offset);
      gl.enableVertexAttribArray(location);
    }

    var chunkSize = attributes.byteLength / 3;
    attrib("x", 0 * chunkSize);
    attrib("n", 1 * chunkSize);
    attrib("h", 2 * chunkSize);

    vertexCount = chunkSize / (3 * 4);
  }

  function setColor(color) {
    var diffuseReflection = new Float32Array(12);
    var k;
    var i;
    for (k = 0; k !== 4; k += 1) {
      for (i = 0; i !== 3; i += 1) {
        diffuseReflection[3 * k + i] = color[i] * diffuseReflectance[3 * k + i];
      }
    }
    gl.uniform3fv(gl.getUniformLocation(program, "d"), diffuseReflection);
  }

  function paint() {
    resize();
    setMatrixUniform(gl, program, "m", m);
    gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);

    gl.cullFace(gl.FRONT);
    gl.drawArrays(gl.TRIANGLES, 0, vertexCount);

    gl.cullFace(gl.BACK);
    gl.drawArrays(gl.TRIANGLES, 0, vertexCount);
  }

  function initShaderProgram(gl, vshader, fshader) {
    function addShader(source, type) {
      function checkShader() {
        if (gl.getShaderParameter(shader, gl.COMPILE_STATUS)) {
          return true;
        }
        var s = "Shader compilation failed:";
        switch (type) {
        case gl.VERTEX_SHADER:
          s = "Vertex " + s;
          break;
        case gl.FRAGMENT_SHADER:
          s = "Fragment " + s;
          break;
        }
        showError(s, gl.getShaderInfoLog(shader));
        return false;
      }
      var shader = gl.createShader(type);
      gl.shaderSource(shader, source);
      gl.compileShader(shader);
      if (checkShader()) {
        gl.attachShader(program, shader);
        return true;
      }
    }

    function checkProgram() {
      if (gl.getProgramParameter(program, gl.LINK_STATUS)) {
        return true;
      }
      showError("Shader program linking failed:", gl.getProgramInfoLog(program));
      return false;
    }

    var program = gl.createProgram();
    if (addShader(vshader, gl.VERTEX_SHADER) && addShader(fshader, gl.FRAGMENT_SHADER)) {
      gl.linkProgram(program);
      if (checkProgram()) {
        gl.useProgram(program);
        return program;
      }
    }
    return null;
  }

  function setMatrixUniform(gl, program, name, matrix) {
    var location = gl.getUniformLocation(program, name);
    gl.uniformMatrix4fv(location, false, matrix);
  }

  function getEventPoint(event) {
    return[ event.clientX, event.clientY ];
  }

  function makeSwipable(element, handleStart, handleEnd, handleMove) {
    function getEventPoint(event) {
      return [event.clientX, event.clientY];
    }

    function resetMouse() {
      document.removeEventListener("mouseup", onMouseup, false);
      element.removeEventListener("mousemove", onMousemove, false);
    }

    function resetTouch() {
      document.removeEventListener("touchend", onTouchend, false);
      element.removeEventListener("touchmove", onTouchmove, false);
    }

    function onMousedown(event) {
      if (event.button === 0) {
        handleStart(getEventPoint(event));
        document.addEventListener("mouseup", onMouseup, false);
        element.addEventListener("mousemove", onMousemove, false);
        if (event.currentTarget === element) {
          // Inhibit selecting nearby text on double click.
          event.preventDefault();
        }
      }
      else {
        resetMouse();
      }
    }

    function onMouseup(event) {
      resetMouse();
      handleEnd();
    }

    function onMousemove(event) {
      handleMove(getEventPoint(event));
    }

    var touchId = null;

    function onTouchstart(event) {
      var e;
      touchId = null;
      if (event.touches.length === 1) {
        e = event.touches[0];
        touchId = e.identifier;
        handleStart(getEventPoint(e));
        document.addEventListener("touchend", onTouchend, false);
        element.addEventListener("touchmove", onTouchmove, false);
        if (event.currentTarget === element) {
          // Inhibit selecting nearby text on long click.
          event.preventDefault();
        }
      }
    }

    function onTouchend(event) {
      touchId = null;
      resetTouch();
      handleEnd();
    }

    function onTouchmove(event) {
      var i;
      var e;
      if (touchId !== null && event.touches.length === 1) {
        e = event.touches[0];
        if (touchId === e.identifier) {
          handleMove(getEventPoint(e));
        }
        else {
          onTouchend(null);
        }
        // Inhibit scroll and zoom during swipe.
        event.preventDefault();
      }
      else {
        onTouchend(null);
      }
    }

    element.addEventListener("mousedown", onMousedown, false);
    element.addEventListener("touchstart", onTouchstart, false);
  }

  function setup() {
    var canvas = null;
    var updatePending = false;
    var lastUpdateTime = null;
    var pickAnchor = null;
    var eventPoint = null;
    var loadState = 0;

    // Network resources.
    var vshader = null;
    var fshader = null;
    var modelIndex = 0;
    var models= [
      { id: "2331", caption: " 1: Tetrahedron" },
      //{ id: "2332", caption: " 1: Tetrahedron" },
      { id: "2330", caption: " 2: Octahedron" },
      //{ id: "2431", caption: " 2: Octahedron" },
      { id: "2432", caption: " 3: Cube" },
      { id: "2337", caption: " 4: Icosahedron" },
      //{ id: "2531", caption: " 4: Icosahedron" },
      { id: "2532", caption: " 5: Dodecahedron" },
      { id: "2334", caption: " 6: Truncated tetrahedron" },
      //{ id: "2335", caption: " 6: Truncated tetrahedron" },
      { id: "2336", caption: " 7: Truncated octahedron" },
      //{ id: "2435", caption: " 7: Truncated octahedron" },
      { id: "2434", caption: " 8: Truncated cube" },
      { id: "2535", caption: " 9: Truncated icosahedron" },
      { id: "2534", caption: "10: Truncated dodecahedron" },
      { id: "2333", caption: "11: Cuboctahedron" },
      //{ id: "2430", caption: "11: Cuboctahedron" },
      { id: "2530", caption: "12: Icosidodecahedron" },
      { id: "2433", caption: "13: Rhombicuboctahedron" },
      { id: "2533", caption: "14: Rhombicosidodecahedron" },
      { id: "2436", caption: "15: Rhombitruncated cuboctahedron" },
      { id: "2536", caption: "16: Rhombitruncated icosidodecahedron" },
      { id: "2437", caption: "17: Snub cube" },
      { id: "2537", caption: "18: Snub dodecahedron" },
    ];

    function sendModelRequest(index) {
      models[index].vertices = null;
      var url = "data/" + models[index].id + ".dat";
      download(url, "arraybuffer", modelResponse, index);
    }

    function modelResponse(response, index) {
      models[index].vertices = response;
      updateLoadState();
      // If we're not already waiting for a response and we haven't
      // already got all the models, start loading another one now.
      var i = 0, j = 0;
      while (i != models.length && models[i].vertices !== undefined) ++i;
      while (j != models.length && models[j].vertices !== null) ++j;
      if (i != models.length && j == models.length) {
        sendModelRequest(i);
      }
    }

    function nextModel(direction) {
      // Change modelIndex.
      var N = models.length;
      modelIndex = (modelIndex + N + direction) % N;
      // Show caption.
      var captionDiv = document.getElementById("caption");
      removeChildren(captionDiv);
      appendText(captionDiv, models[modelIndex].caption);
      // Start loading unless already started loading.
      if (models[modelIndex].vertices === undefined) {
        sendModelRequest(modelIndex);
      }
      // Display when loaded.
      loadState = Math.min(loadState, 4);
      updateLoadState();
    }

    function pick1() {
      return eventPoint && pick(eventPoint, canvas.getBoundingClientRect(), x);
    }

    function swipeStart(point) {
      eventPoint = point;
      var pickPoint = pick1();
      pickAnchor = pickPoint && unrotate(pickPoint);
    }

    function swipeMove(point) {
      eventPoint = point;
      updateLoadState();
    }

    function swipeEnd() {
      eventPoint = null;
    }

    // Ensure exactly one update is pending.
    function invalidate() {
      if (!updatePending) {
        updatePending = true;
        requestAnimationFrame(update);
      }
    }

    // Update the view.
    function update(updateTime) {
      updatePending = false;
      var updateInterval = Math.max(1, Math.min(30, updateTime - lastUpdateTime));
      var pickPoint = pick1();
      lastUpdateTime = updateTime;
      if (!pickPoint) {
        pickAnchor = null;
      }
      else if (!pickAnchor) {
        pickAnchor = unrotate(pickPoint);
        w[0] = 0.0;
        w[1] = 0.0;
        w[2] = 0.0;
      }
      else {
        spinTo(pickAnchor, pickPoint, updateInterval);
      }

      var inMotion = body.advance(MAIN, updateInterval) !== 0;
      body.matrix(MAIN);

      try {
        paint();
        // Ensure an update is pending if the ball is in motion.
        if (inMotion) {
          invalidate();
        }
      }
      catch (e) {
        showError(e.toString(), e.stack);
      }
    }

    function updateLoadState() {
      try {
        if (loadState === 0) {
          loadState = 1;
          addEventListener("resize", function (ignore) {
            eventPoint = null;
            updateLoadState();
          }, false);

          canvas = document.getElementById("canvas");
          if (canvas) {
            makeSwipable(canvas, swipeStart, swipeEnd, swipeMove);
            loadState = 2;
          }
        }

        if (loadState === 2 && vshader && fshader) {
          loadState = 3;
          gl = canvas.getContext("experimental-webgl");
          if (initializeGraphics(vshader, fshader)) {
            x[2] = -tz - 0.5 * td;
            randomVectorInBall(u, 3.14159265);
            randomVectorOnSphere(w, 0.001);
            loadState = 4;
          }
        }

        if (loadState === 4 && models[modelIndex].vertices) {
          loadState = 5;
          loadVertexBuffer(models[modelIndex].vertices);
          setColor(randomColor(0.5));
          loadState = 6;
        }

        if (loadState === 6) {
          invalidate();
        }
      } catch (e) {
        showError(e.toString(), e.stack);
      }
    }

    function start() {
      function setupButton(id, handler) {
        var button = document.getElementById(id);
        if (button) {
          button.addEventListener("click", handler);
        }
      }
      download("vertex-shader.glsl", "text", function (response) { vshader = response; updateLoadState(); });
      download("fragment-shader.glsl", "text", function (response) { fshader = response; updateLoadState(); });
      nextModel(0);
      setupButton("previous", function (ignore) { nextModel(-1); });
      setupButton("next", function (ignore) { nextModel(+1); });
    };

    start();
  }

  document.addEventListener("DOMContentLoaded", function () {
    setup();
  }, true);
}());
