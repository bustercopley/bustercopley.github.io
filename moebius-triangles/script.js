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

(function () {
  "use strict";

  // Create asm.js Body object.
  var heap = new ArrayBuffer(0x10000);  // 64 kilobytes (the minimum).
  var body = Body(window, null, heap);

  // Bodies (byte offsets from start of heap):
  var MAIN = 128; // Main body.
  var AUX = 256;  // Auxiliary body.

  // Fields (byte offsets from start of body):
  var PX = 0;   // Linear position.
  var PU = 16;  // Angular position.
  var PW = 32;  // Angular velocity.
  var PM = 64;  // Modelview matrix.

  // Views into the main body.
  var x = new Float32Array(heap, MAIN + PX, 4);  // Centre of sphere.
  var u = new Float32Array(heap, MAIN + PU, 4);  // Angular position.
  var w = new Float32Array(heap, MAIN + PW, 4);  // Angular momentum.
  var m = new Float32Array(heap, MAIN + PM, 16); // Modelview matrix.

  // Operations on three dimensional vectors.
  function dot(u, v) {
    return u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
  }

  function cross(u, v) {
    return [
      u[1] * v[2] - u[2] * v[1],
      u[2] * v[0] - u[0] * v[2],
      u[0] * v[1] - u[1] * v[0]];
  }

  function rotate(x) {
    return [
      m[0] * x[0] + m[4] * x[1] + m[8] * x[2],
      m[1] * x[0] + m[5] * x[1] + m[9] * x[2],
      m[2] * x[0] + m[6] * x[1] + m[10] * x[2]];
  }

  function unrotate(x) {
    return [
      m[0] * x[0] + m[1] * x[1] + m[2] * x[2],
      m[4] * x[0] + m[5] * x[1] + m[6] * x[2],
      m[8] * x[0] + m[9] * x[1] + m[10] * x[2]];
  }

  function normalize(u) {
    var s = Math.hypot(u[0], u[1], u[2]);
    return [u[0] / s, u[1] / s,  u[2] / s];
  }

  function scale(u, scale) {
    return [scale * u[0], scale * u[1], scale * u[2]];
  }

  function sq(x) {
    return x * x;
  }

  // The system represents the tiling in model space, with no boundary circle.
  var systemDarts = [], systemNodes = [], systemCircles = [];
  var systemSigma = [], systemAlpha = [];

  var systemParameters = [
    [2, 3, 3], [4, 3, 2], [2, 3, 5],
    [2, 2, 2], [2, 2, 3], [2, 2, 4], [2, 2, 5], [2, 2, 6], [2, 2, 7],
  ];

  function setSystemIndex(index) {
    makeSystem(systemParameters[index]);
  }

  function moebiusOrder(p, q, r) {
    return (2 * p * q * r) / (q * r + r * p + p * q - p * q * r);
  }

  // Find coordinates in model space of three vertices of a spherical
  // triangle whose angles are angles π/p, π/q, π/r. (Assume p = 2.)
  function getTriangle(p, q, r) {
    // The angles of the spherical triangle.
    var A = Math.PI / p;
    var B = Math.PI / q;
    var C = Math.PI / r;
    var sinA = Math.sin(A), cosA = Math.cos(A);
    var sinB = Math.sin(B), cosB = Math.cos(B);
    var sinC = Math.sin(C), cosC = Math.cos(C);
    // Cosines of the edge lengths, by the spherical law of cosines.
    var cosb = (cosB + cosC * cosA) / (sinC * sinA);
    var cosc = (cosC + cosA * cosB) / (sinA * sinB);
    // Sines of the edge lengths, by the trigonometrical identity.
    var sinb = Math.sqrt(1.0 - cosb * cosb);
    var sinc = Math.sqrt(1.0 - cosc * cosc);
    // Vertices of the spherical triangle.
    return [
      [1.0, 0.0, 0.0],
      [cosc, sinc, 0.0],
      [cosb, cosA * sinb, -sinA * sinb],
    ];
  }

  // Return a function that rotates a vector 'v' through 'angle' about 'axis'.
  function Rotator(axis, order) {
    var u = new Float32Array(heap, AUX + PU, 4);  // View of angular position.
    var m = new Float32Array(heap, AUX + PM, 16); // View of modelview matrix.
    var angle = 6.283185307179586 / order;
    u[0] = angle * axis[0];
    u[1] = angle * axis[1];
    u[2] = angle * axis[2];
    u[3] = 0;
    body.matrix(AUX);
    return function (v) {
      var result = new Array(3);
      result[0] = m[0] * v[0] + m[4] * v[1] + m[8] * v[2];
      result[1] = m[1] * v[0] + m[5] * v[1] + m[9] * v[2];
      result[2] = m[2] * v[0] + m[6] * v[1] + m[10] * v[2];
      return result;
    }
  }

  function sigmaInverse(dart, order) {
    if (dart % order < 2) {
      dart = (dart ^ 1) + order;
    }
    return (dart - 2) | 0;
  }

  function makeSystem(p) {
    var n = moebiusOrder(p[0], p[1], p[2]);
    var i, j, m, t, d;
    var a = Array(3), b = Array(3), next = Array(3);
    var mu = [0, 1, 2], even = true;

    systemAlpha.length = 6 * n;
    systemSigma.length = 6 * n;
    systemDarts.length = 6 * n;
    systemNodes.length = n + 2;
    systemCircles.length = 0;

    for (i = 0; i !== 6 * n; ++i) {
      systemAlpha[i] = undefined;
      systemDarts[i] = {
        node: undefined,
        shaded: undefined,
        pole: undefined,
        reverse: undefined,
        argument: undefined,
      };
    }

    var triangle = getTriangle(p[0], p[1], p[2]);
    var dart = 0, node = 0, predecessor;
    for (i = 0; i !== 3; ++i) {
      a[i] = dart;
      next[i] = dart + 2 * p[i];
      systemNodes[node] = triangle[i];
      for (m = 0; m !== n / p[i]; ++m) {
        predecessor = dart + 2 * p[i] - 1;
        for (j = 0; j !== p[i]; ++j) {
          systemSigma[predecessor] = dart;
          systemSigma[1 ^ predecessor] = 1 ^ dart;
          predecessor = dart;
          systemDarts[dart].node = node;
          systemDarts[dart].shaded = !(j & 1);
          systemDarts[1 ^ dart].node = node;
          systemDarts[1 ^ dart].shaded = !((j & 1) ^ (p[i] & 1));
          dart += 2;
        }
        ++node;
      }
      a[i] = systemSigma[a[i]];
    }
    var rotate = Rotator(systemNodes[0], p[0]);

    while (true) {
      systemAlpha[systemSigma[a[0]]] = a[1];
      systemAlpha[systemSigma[a[1]]] = a[2];
      systemAlpha[systemSigma[a[2]]] = a[0];

      if (systemAlpha[a[2]] === undefined) {
        b[1] = sigmaInverse(a[2], 2 * p[mu[2]]);
        if (systemAlpha[b[1]] !== undefined) {
          b[0] = sigmaInverse(systemAlpha[b[1]], 2 * p[mu[0]]);
          b[2] = systemSigma[a[1]];
          systemAlpha[systemSigma[b[0]]] = b[1];
          systemAlpha[systemSigma[b[1]]] = b[2];
          systemAlpha[systemSigma[b[2]]] = b[0];
        }
      }

      if (systemAlpha[systemSigma[systemSigma[a[0]]]] === undefined) {
        a[0] = systemSigma[a[0]];
        t = mu[1];
        mu[1] = mu[2];
        mu[2] = t;
      }
      else if (systemAlpha[systemSigma[systemSigma[a[1]]]] === undefined) {
        a[0] = systemSigma[a[1]];
        t = mu[0];
        mu[0] = mu[1];
        mu[1] = t;
        rotate = Rotator(systemNodes[systemDarts[a[0]].node], p[mu[0]]);
      }
      else break;
      even = !even;

      d = systemAlpha[systemSigma[systemSigma[a[0]]]];
      if (d !== undefined) {
        a[1] = systemSigma[systemAlpha[systemSigma[d]]];
      }
      else {
        a[1] = next[mu[1]];
        next[mu[1]] += 2 * p[mu[1]];
        d = systemAlpha[systemSigma[systemAlpha[a[0]]]];
        systemNodes[systemDarts[a[1]].node] = rotate(
          systemNodes[systemDarts[d].node]);
        if (even) a[1] = systemSigma[a[1]];
      }
      a[2] = sigmaInverse(systemAlpha[a[0]], 2 * p[mu[2]]);
    }

    // Find the γ-orbits and calculate poles and arguments.
    for (i = 0; i !== 6 * n; ++i) {
      if (systemDarts[i].shaded && systemDarts[i].pole === undefined) {
        var x = systemNodes[systemDarts[i].node];
        var y = systemNodes[systemDarts[systemAlpha[i]].node];
        var e2 = normalize(cross(x, y));
        var e0 = x; // Already a unit vector.
        var e1 = cross(e2, e0);
        var argument = 0.0;
        var pole = systemCircles.push({dart: i, e2: e2, e0: e0}) - 1;
        d = i;
        do {
          systemDarts[d].pole = pole;
          systemDarts[d].reverse = false;
          systemDarts[d].argument = argument;
          systemDarts[1 ^ d].pole = pole;
          systemDarts[1 ^ d].reverse = true;
          systemDarts[1 ^ d].argument = argument;
          d = 1 ^ systemAlpha[d];
          x = systemNodes[systemDarts[d].node];
          argument = Math.atan2(dot(x, e1), dot(x, e0));
        } while (d !== i);
      }
    }
  }

  // The zoom parameter is from 0 to 1 inclusive. It controls the distance x₂
  // from the eye (at the origin) to the centre of the object sphere and the
  // distance zs from the eye to the surface of projection.

  var zs;              // zs: z coordinate of the surface of projection.
  var stereographic;   // True if eye is on sphere.
  var zoomPower = 4.0;
  var targetRadius = 0.67;
  var zoomThreshold = 0.875;
  var near = 0.0005;
  var far = 20.0;
  var thresholdCentreDistance = centreDistance(zoomThreshold);

  // Calculate z-coordinate of centre of object sphere from zoom parameter.
  function centreDistance(zoom) {
    return -(1 + near) - (far - near) * Math.pow(1 - zoom, zoomPower);
  }

  function setZoom(zoom) {
    stereographic = zoom === 1;
    x[2] = centreDistance(zoom);
    var z = zoom < zoomThreshold ? x[2] : thresholdCentreDistance;
    // We have r = -zs / √(x₂² - 1) (see renderModel),
    // whence zs = -r √(x₂² - 1).
    zs = -targetRadius * Math.sqrt(z * z - 1);
  }

  // Circles in world space.

  // For r ∊ R and c, n ∊ R³, where n is a unit vector, define figures
  //   Σ(r², c)     the sphere radius r, centre c,
  //   Π(c, n)      the plane normal to n passing through c,
  //   Λ(c)         the set of points X such that OXC is a right angle,
  //   Γ(r², c, n)  the great circle Σ(r², c) ∩ Π(n, c),
  //   Δ(r², c)     the boundary circle Σ(r², c) ∩ Λ(c),
  //   Ω(r², c, n)  the figure Γ(r², c, n) ∩ Δ(r², c).
  // It turns out Λ(c) is the sphere centre ½c, radius ‖½c‖.

  // Σ(r², c)    = {x ∊ R³: ‖x - c‖² = r²}.
  // Π(c, n)     = {x ∊ R³: <n, x - c> = 0}.
  // Λ(c)        = {x ∊ R³: <x, x - c> = 0}.
  // Γ(r², c, n) = Σ(r², c) ∩ Π(c, n).
  // Δ(r², c)    = Σ(r², c) ∩ Λ(c).
  // Ω(r², c, n) = Σ(r², c) ∩ Π(n, c) ∩ Λ(c).

  // Return r′², c′, n′ such that Γ(r′², c′, n′) = Δ(r², c).
  function getBoundary(rsq, c) {
    // Since x is in Σ ∩ Λ we have ‖x‖² = <c, x> = ‖c‖² - r² (*).
    // Draw the triangle OXC, with a right angle at X and |XC| = r,
    // and drop a perpendicular from X to meet OC in a right angle at B.
    // Let h = |BX|. Resolving x along c gives b = (<c, x>/‖c‖²)c,
    // and by similarity of triangles OBX and OXC, h²=(‖x‖²/‖c‖²)r².
    // Then by (*), b = (1 - r²/‖c‖²)c and h² = (1 - r²/‖c‖²)r².
    var n = normalize(c);
    var s = 1.0 - rsq / dot(c, c);
    return {rsq: s * rsq, centre: scale(c, s), normal: n};
  }

  // Return the intersection Ω of the great circle Γ = Γ(r², c, n)
  // and the boundary circle Δ = Σ ∩ Λ as an array of two points in
  // world space, or undefined if Γ and Δ do not intersect.
  function boundaryIntersects(rsq, c, n) {
    // Choose an orthogonal basis {e₀, e₁, e₂} such that n = e₀
    // and c = c₀e₀ + c₁e₁, and suppose that x = x₀e₀ + x₁e₁ + x₂e₂
    // is an element of Ω(r², c, n) = Π(n, c) ∩ Σ(r², c) ∩ Λ(c). Then
    // <n, x - c> = 0 and ‖x - c‖² = r² and <x, x - c> = 0.
    // Deduce that <n, x> = <n, c> and <c, x> = ‖c‖² - r², or x₀ = c₀
    // and c₀x₀ + c₁x₁ = c₀² + c₁² - r², whence x₁ = (c₁² - r²) / c₁.
    // To make all this explicit, we have
    //   c₀ = <c, n>,
    //   c₁e₁ = c - c₀e₀ = c - <c, n>n,
    //   c₁² = ‖c₁e₁‖² = ‖c‖² - <c, n>²,
    //   e₁ = (c - c₀e₀) / c₁ = (c - <c, n>n) / c₁,
    //   e₂ = e₀ × e₁ = n × c / c₁.
    // Therefore
    //   x = x₀e₀ + x₁e₁ + x₂e₂
    //     = c₀n + [(c₁² - r²) / c₁](c - <c, n>n) / c₁ + x₂e₂
    //     = <c, n>n + [1 - r² / c₁²] (c - <c, n>n) + x₂e₂
    //     = (r² / c₁²)<c, n>n + [1 - r² / c₁²]c + x₂e₂.
    // We have ‖x - c‖² = r², or
    //   ‖(r² / c₁²)<c, n>n + [1 - r² / c₁²]c + x₂e₂ - c‖² = r²,
    // but e₂ is orthogonal to n and c, so by the Pythagorean identity,
    //   x₂² + ‖(r² / c₁²)<c, n>n - (r² / c₁²)c‖² = r², or
    //   x₂² + (r⁴ / c₁⁴)‖c - <c, n>n‖² = r², or
    //   x₂² + r⁴ / c₁² = r², or
    //   x₂² = r²(1 - r² / c₁²), whence
    //   x₂ = ± √[r²(1 - r² / c₁²)], if the square root exists.
    // Putting it all together, we have x = un + vc ± wn × c, where
    //   u = (r² / c₁²)<c, n>,
    //   v = 1 - r² / c₁²,
    //   w = √[(r² / c₁²)(1 - r² / c₁²)].
    var c0 = dot(c, n);
    var csq = dot(c, c);
    var c1sq = csq - c0 * c0;
    if (c1sq > rsq) {
      var qsq = rsq / c1sq;
      var u = qsq * c0;
      var v = 1 - qsq;
      var w = Math.sqrt(qsq * v);
      var d0 = u * n[0] + v * c[0];
      var d1 = u * n[1] + v * c[1];
      var d2 = u * n[2] + v * c[2];
      var e0 = w * (n[1] * c[2] - n[2] * c[1]);
      var e1 = w * (n[2] * c[0] - n[0] * c[2]);
      var e2 = w * (n[0] * c[1] - n[1] * c[0]);
      return [[d0 - e0, d1 - e1, d2 - e2], [d0 + e0, d1 + e1, d2 + e2]];
    }
  }

  // Ellipses in user space.

  // Project a circle Γ = Γ(r², c, n) to an ellipse in user space.
  // Return the centre [u₀, u₁], the semi-axes r₀ and r₁ and the rotation φ.
  function projectCircle(rsq, c, n) {
    // Compute the coefficients of the quadratic equation
    // Au² + Bv² + Cw² + 2Dvw + 2Ewu + 2Fuv = 0
    // describing the ellipse in device coordinates. See
    // <https://www.math.utah.edu/~treiberg/Perspect/Perspect.htm#circle>.
    var nc0 = n[0] * c[0];
    var nc1 = n[1] * c[1];
    var nc2 = n[2] * c[2];
    var nc = nc0 + nc1 + nc2;
    var ssq = c[0] * c[0] + c[1] * c[1] + c[2] * c[2] - rsq;
    var A = n[0] * n[0] * ssq + nc * (nc1 + nc2 - nc0);
    var B = n[1] * n[1] * ssq + nc * (nc2 + nc0 - nc1);
    var C = n[2] * n[2] * ssq + nc * (nc0 + nc1 - nc2);
    var D = n[1] * n[2] * ssq - nc * (n[1] * c[2] + n[2] * c[1]);
    var E = n[2] * n[0] * ssq - nc * (n[2] * c[0] + n[0] * c[2]);
    var F = n[0] * n[1] * ssq - nc * (n[0] * c[1] + n[1] * c[0]);

    // Compute the centre [u₀, u₁], the semi-axes r₀ and r₁
    // and the rotation φ from the quadratic equation. See
    // <http://mathworld.wolfram.com/Ellipse.html>.
    var discriminant = A * B - F * F;
    if (discriminant) {
      var fdbe = (F * D - B * E) / discriminant;
      var fead = (F * E - A * D) / discriminant;
      var u0 = zs * fdbe;
      var u1 = zs * fead;
      var alpha = -2 * (fdbe * E + fead * D + C);
      var y = -2 * F;
      var x = B - A;
      var delta = Math.hypot(x, y);
      var phi = 0.5 * Math.atan2(y, x);
      var q0 = A + B + delta;
      var r0 = -0.5 * zs * Math.sqrt(alpha * q0 / discriminant);
      var r1 = -zs * Math.sqrt(alpha / q0);
      if (r1 >= 0.001) {
        return {c: [u0, u1], r: [r0, r1], phi: phi};
      }
    }
  }

  // Find the argument of the point p along the ellipse e.
  function argumentInEllipse(p, e) {
    var c = Math.cos(e.phi);
    var s = Math.sin(e.phi);
    var x = c * (p[0] - e.c[0]) + s * (p[1] - e.c[1]);
    var y = -s * (p[0] - e.c[0]) + c * (p[1] - e.c[1]);
    return Math.atan2(y / e.r[1], x / e.r[0]);
  }

  // The scheme represents the tiling in world space, with boundary circle.
  var schemeCircles = [], schemeNodes = [], schemeDarts = [], schemeFaces = [];
  var schemeAlpha = [], schemeSigma = [];

  // Let Γ be the circle specified by (rsq, centre, normal).
  // Find a right-handed orthogonal basis {e₀, e₁, e₂} with e₂ == ±normal.
  // Then the path γ given by γ(θ) = (√rsq)((cosθ)e₀ + (sinθ)e₁)
  // is a parameterization of Γ. Choose the sense of e₂ so that γ traverses Γ
  // in the positive (anticlockwise) sense as seen from the origin.
  function addSchemeCircle(rsq, centre, normal) {
    var reverse = dot(centre, normal) >= 0.0;
    var n = scale(normal, reverse ? -1.0 : 1.0);
    var ellipse = projectCircle(rsq, centre, n);
    return schemeCircles.push({
      rsq: rsq,
      centre: centre,
      e0: undefined,
      e1: undefined,
      e2: n,
      reverse: reverse,
      ellipse: ellipse,
      meetVertices: undefined,
      dart: undefined,
    }) - 1;
  }

  function addSchemeDart(node, pole, shaded, reverse, front) {
    return schemeDarts.push({
      node: node,
      pole: pole,
      argument: undefined,
      shaded: shaded,
      reverse: reverse,
      front: front,
      ellipseArgument: undefined,
    }) - 1;
  }

  function split(d0, d1, d2, d3) {
    schemeAlpha[d3] = d0;
    schemeAlpha[d0] = d3;
    schemeAlpha[d2] = d1;
    schemeAlpha[d1] = d2;
  }

  function perspective(v) {
    return [zs * v[0] / v[2], zs * v[1] / v[2]];
  }

  function between(t, begin, end) {
    return t - end < ((end < begin) - (t < begin)) * 6.283185307179586;
  }

  // Find a dart crossing front to back.
  function findCrossingDart(circle, reverse) {
    var e = circle.dart;
    // First, search for a front dart whose twin is a rear dart.
    if (schemeDarts[e].reverse !== reverse) {
      e = schemeAlpha[e];
    }
    var d = e;
    do {
      if (schemeDarts[d].front && !schemeDarts[schemeAlpha[d]].front) {
        return d;
      }
      d = 1 ^ schemeAlpha[d];
    } while (d !== e);

    // If that failed, the boundary intersects the great circle twice between
    // two adjacent rear nodes, creating a front-facing digonal face. Search
    // for the dart that crosses their midpoint (on the great circle). This
    // selects the correct dart even if one intersection is very near a node.
    var midpoint = [
      circle.meetVertices[0][0] + circle.meetVertices[1][0] - 2.0 * circle.centre[0],
      circle.meetVertices[0][1] + circle.meetVertices[1][1] - 2.0 * circle.centre[1],
      circle.meetVertices[0][2] + circle.meetVertices[1][2] - 2.0 * circle.centre[2],
    ];
    var x = dot(midpoint, circle.e0);
    var y = dot(midpoint, circle.e1);
    var arg = Math.atan2(y, x);

    // Search darts in ascending order of argument.
    e = circle.dart;
    do {
      d = e;
      e = 1 ^ schemeAlpha[d];
    } while (!between(arg, schemeDarts[d].argument, schemeDarts[e].argument));
    if (schemeDarts[d].reverse !== reverse) {
      d = schemeAlpha[d];
    }
    return d;
  }

  function makeScheme(rsq) {
    var i, n;
    var d, e;
    var v, p;
    var circle, boundary;
    var r = Math.sqrt(rsq);

    var maxSchemeDartCount = systemDarts.length + 8 * systemCircles.length;
    schemeAlpha.length = maxSchemeDartCount;
    schemeSigma.length = maxSchemeDartCount;

    // Add the system circles.
    schemeCircles.length = 0;
    for (i = 0; i !== systemCircles.length; ++i) {
      addSchemeCircle(rsq, x, rotate(systemCircles[i].e2));
      circle = schemeCircles[i];
      circle.e0 = rotate(systemCircles[i].e0);
      circle.e1 = cross(circle.e2, circle.e0);
      circle.meetVertices =
        boundaryIntersects(circle.rsq, circle.centre, circle.e2);
      circle.dart = systemCircles[i].dart;
      if (circle.reverse) {
        circle.dart = 1 ^ circle.dart;
      }
    }

    // Add the boundary circle.
    var b = getBoundary(rsq, x);
    var boundary = addSchemeCircle(b.rsq, b.centre, b.normal);

    // Add system nodes and mark front nodes.
    circle = schemeCircles[boundary];
    schemeNodes.length = systemNodes.length;
    var bcn = dot(circle.centre, circle.e2);
    for (i = 0; i !== systemNodes.length; ++i) {
      v = rotate(systemNodes[i]);
      v[0] = x[0] + r * v[0];
      v[1] = x[1] + r * v[1];
      v[2] = x[2] + r * v[2];
      schemeNodes[i] = {
        point: perspective(v),
        front: dot(v, circle.e2) > bcn,
      };
    }

    // Add system darts.
    schemeDarts.length = systemDarts.length;
    for (d = 0; d !== systemDarts.length; ++d) {
      var dart = systemDarts[d];
      var circle = schemeCircles[dart.pole];
      var reverse = circle.reverse;
      var front = schemeNodes[dart.node].front;
      var argument = (reverse ? -1 : 1) * dart.argument;
      schemeDarts[d] = {
        node: dart.node,
        pole: dart.pole,
        argument: argument,
        shaded: dart.shaded,
        reverse: dart.reverse !== reverse,
        front: front,
        ellipseArgument: undefined,
      };
      schemeAlpha[d] = systemAlpha[d];
      schemeSigma[d] = systemSigma[d];
    }

    var b0, b1, b2, b3, c0, c1, c2, c3;

    // Break great circles at boundary intersections.
    for (i = 0; i !== boundary; ++i) {
      circle = schemeCircles[i];
      if (circle.meetVertices) {
        for (n = 0; n !== 2; ++n) {
          var reverse = n === 1;
          var node = schemeNodes.push({
            point: perspective(circle.meetVertices[n]),
            front: false,
          }) - 1;
          c0 = findCrossingDart(circle, reverse);
          c1 = schemeAlpha[c0];
          c2 = addSchemeDart(node, i, schemeDarts[c0].shaded, reverse, false);
          c3 = addSchemeDart(node, i, !schemeDarts[c0].shaded, !reverse, true);
          split(c0, c1, c2, c3);
        }
      }
    }

    // Break the boundary circle and stitch up.
    if (c2 !== undefined) {
      c0 = c2;
      b0 = schemeCircles[boundary].dart = schemeDarts.length;
      b1 = 1 ^ b0;
      do {
        dart = schemeDarts[c2];
        c3 = 1 ^ c2;
        b2 = addSchemeDart(dart.node, boundary, !dart.shaded, false, false);
        b3 = addSchemeDart(dart.node, boundary, dart.shaded, true, true);
        split(b0, b1, b2, b3);
        schemeSigma[c2] = b2;
        schemeSigma[b2] = c3;
        schemeSigma[c3] = b3;
        schemeSigma[b3] = c2;
        b1 = b3;
        // Advance c2 to the next break.
        c2 = schemeAlpha[c2];
        c2 = schemeAlpha[schemeSigma[c2]];
        if (c2 < systemDarts.length) {
          c2 = schemeAlpha[schemeSigma[c2]];
          if (c2 < systemDarts.length) {
            c2 = schemeAlpha[schemeSigma[c2]];
          }
        }
      } while (c2 !== c0);
    }

    // Fill in the ellipseArgument field.
    for (i = 0; i !== schemeCircles.length; ++i) {
      var ellipse = schemeCircles[i].ellipse;
      n = schemeCircles[i].dart;
      if (ellipse !== undefined && n !== undefined) {
        e = 1 ^ n;
        d = n;
        do {
          p = schemeNodes[schemeDarts[d].node].point;
          var theta = argumentInEllipse(p, ellipse);
          schemeDarts[d].ellipseArgument = theta;
          schemeDarts[e].ellipseArgument = theta;
          e = schemeAlpha[d];
          d = 1 ^ e;
        } while (d !== n);
      }
    }

    // Build faces (orbits of φ = ασ).
    schemeFaces.length = 0;
    var dartDone = Array(schemeDarts.length);
    for (n = 0; n !== schemeDarts.length; ++n) {
      if (!dartDone[n]) {
        var outer = !schemeDarts[n].front;
        d = n;
        do {
          dartDone[d] = true;
          outer = outer && schemeDarts[d].reverse;
          d = schemeSigma[schemeAlpha[d]];
        } while (d !== n);
        schemeFaces.push({dart: n, outer: outer});
      }
    }
  }

  // [rear, front]. Edge doesn't support "#rgba".
  var faceStyle = ["#ff4", "rgba(0,255,255,0.67)"];
  var edgeStyle = ["#788", "#000"];

  function renderModel(canvas) {
    var width = canvas.width;
    var height = canvas.height;
    var scale = Math.min(width, height) / 2;
    var context = canvas.getContext("2d");
    var boundary = systemCircles.length;
    var hasOuterFace = false;

    function resetTransform() {
      context.setTransform(scale, 0, 0, scale, 0.5 * width, 0.5 * height);
    }

    // Add an elliptical arc.
    function addEllipse(ellipse, t0, t1, reverse) {
      var r0 = scale * ellipse.r[0];
      var r1 = scale * ellipse.r[1];
      var x = scale * ellipse.c[0];
      var y = scale * ellipse.c[1];
      var c = Math.cos(ellipse.phi);
      var s = Math.sin(ellipse.phi);
      context.setTransform(r0 * c, r0 * s, r1 * -s, r1 * c,
                           0.5 * width + x, 0.5 * height + y);
      context.arc(0, 0, 1, t0, t1, reverse);
    }

    makeScheme(1.0);

    context.setTransform(1, 0, 0, 1, 0, 0);
    context.clearRect(0, 0, width, height);
    context.lineWidth = 0.5;
    context.lineJoin = "bevel";
    context.lineCap = "butt";
    // Fill and then stroke the faces (rear faces first).
    for (var facing = 0; facing !== 2; ++facing) {
      var front = facing === 1;
      context.fillStyle = faceStyle[facing];
      context.strokeStyle = edgeStyle[facing | stereographic];
      // Hide front faces in stereographic mode.
      if (stereographic && front) {
        continue;
      }
      for (var i = 0; i !== schemeFaces.length; ++i) {
        var d = schemeFaces[i].dart;
        if (schemeDarts[d].shaded && schemeDarts[d].front === front) {
          context.beginPath();
          // If this is an outer face, the area to be filled is the
          // inside of the boundary ellipse minus the inside of the
          // face, so add the boundary ellipse to the path.
          if (schemeFaces[i].outer) {
            ellipse = schemeCircles[boundary].ellipse;
            hasOuterFace = true;
            addEllipse(ellipse, 0, 6.2832, false);
          }
          resetTransform();
          var p0 = schemeNodes[schemeDarts[d].node].point;
          context.moveTo(p0[0], p0[1]);
          do {
            // Add an elliptical arc or straight line to the subpath.
            var e = schemeAlpha[d];
            var p1 = schemeNodes[schemeDarts[e].node].point;
            var ellipse = schemeCircles[schemeDarts[d].pole].ellipse;
            if (ellipse && Math.hypot(p1[0] - p0[0], p1[1] - p0[1]) > 0.025) {
              addEllipse(ellipse,
                         schemeDarts[d].ellipseArgument,
                         schemeDarts[e].ellipseArgument,
                         schemeDarts[d].reverse);
            }
            else {
              resetTransform();
              context.lineTo(p1[0], p1[1]);
            }
            p0 = p1;
            d = schemeSigma[e];
          } while (d !== schemeFaces[i].dart);
          // Stroke the face.
          context.closePath();
          context.setTransform(1, 0, 0, 1, 0, 0);
          context.fill("evenodd");
          context.stroke();
        }
      }
    }

    if (!stereographic) {
      // Start a new path and add the boundary ellipse.
      ellipse = schemeCircles[boundary].ellipse;
      context.beginPath();
      addEllipse(ellipse, 0, 6.2832, false);
      if (hasOuterFace) {
        // If there is a shaded rear outer face, the whole interior of
        // the boundary ellipse is a shaded front-facing face.
        context.fill();
      }
      // Stroke the boundary ellipse.
      context.setTransform(1, 0, 0, 1, 0, 0);
      context.lineWidth = 1.0;
      context.stroke();
    }

    // Garbage (as long as no one draws the same scheme twice).
    schemeCircles.length = 0;
    schemeNodes.length = 0;
    schemeDarts.length = 0;
    schemeFaces.length = 0;
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
      touchId = null;
      if (event.touches.length === 1) {
        var e = event.touches[0];
        touchId = e.identifier;
        document.addEventListener("touchend", onTouchend, false);
        element.addEventListener("touchmove", onTouchmove, false);
        if (handleStart(getEventPoint(e))) {
          // Rule out the start of a gesture.
          if (event.currentTarget === element) {
            event.preventDefault();
          }
        }
      }
    }

    function onTouchend(event) {
      touchId = null;
      resetTouch();
      handleEnd();
    }

    function onTouchmove(event) {
      if (touchId !== null && event.touches.length === 1) {
        var e = event.touches[0];
        if (touchId === e.identifier) {
          handleMove(getEventPoint(e));
          return;
        }
      }
      onTouchend(null);
    }

    element.addEventListener("mousedown", onMousedown, false);
    element.addEventListener("touchstart", onTouchstart, false);
  }

  // Sphere point picking.

  function pick(p, rect) {
    // p is a pair of coordinates in client space (the mouse position).
    // rect is the bounding rectangle of the canvas, in client space.
    // X is a triple of coordinates in world space (the centre of the object).
    // Let L be the line in world space whose points all project to p
    // in client space.
    // Let T be point in world space where L intersects the unit
    // sphere centre X, or the nearer of the two points if the line
    // intersects the sphere twice, or null if they do not intersect.
    // Return the unit vector T - X (or null if T is null).

    var w = rect.right - rect.left;
    var h = rect.bottom - rect.top;
    var s0 = rect.left + 0.5 * w;
    var s1 = rect.top + 0.5 * h;
    var ss = 0.5 * Math.min(w, h);

    // Find the point Q in world space where L intersects the near
    // clipping plane.
    var q = [(p[0] - s0) / ss, (p[1] - s1) / ss, zs];

    // Since 0, Q and T are on the line L, and Q ≠ 0 (since zs ≠ 0), there is
    // some real t such that T = tQ.
    // And since T is on the unit sphere centre X, we have ‖T - X‖² = 1.
    // Combining these two facts we get
    //   ‖tQ - X‖² = 1,
    // which is a quadratic equation in t, namely
    //   ‖Q‖²t² - 2<Q,X>t + ‖X‖² - 1 = 0,
    // whose solutions are t = (<Q,X> ± √(<Q,X>² - ‖Q‖²(‖X‖²-1))) / ‖Q‖².
    // If the discriminant is negative then L doesn't intersect the sphere.
    var qsq = dot(q, q);
    var xsq = dot(x, x);
    var qx = dot(q, x);
    var scrim = qx * qx + qsq - qsq * xsq;
    if (scrim >= 0) {
      // Take the smaller of the two solutions to get the intersection
      // closer to Q, unless we're in stereographic projection.
      var t = (qx - (stereographic ? -1 : 1) * Math.sqrt(scrim)) / qsq;
      // Return T - X.
      return [t * q[0] - x[0], t * q[1] - x[1], t * q[2] - x[2]];
    }
  }

  // Set the angular velocity of the body so that the direction vector anchor
  // (in model space) will point in direction target (in world space)
  // after time dt.
  function spinTo(anchor, target, dt) {
    var current = rotate(anchor);
    var crs = cross(current, target);
    // ‖crs‖ = sin(θ), where θ is the angle between current and target.
    var sinSquaredTheta = dot(crs, crs);
    var thetaOverSinTheta;
    if (sinSquaredTheta < 1.0e-8) {
      thetaOverSinTheta = 1 + sinSquaredTheta / 6;
    }
    else {
      var sinTheta = Math.sqrt(sinSquaredTheta);
      thetaOverSinTheta = Math.asin(sinTheta) / sinTheta;
    }
    // The required value is in the direction of crs and has magnitude θ/dt.
    var lambda = thetaOverSinTheta / dt;
    w[0] = lambda * crs[0];
    w[1] = lambda * crs[1];
    w[2] = lambda * crs[2];
    w[3] = 0.0;
  }

  // Random vectors.

  function randomVectorInBall(x, r) {
    var rsq = r * r;
    var xsq;
    do {
      x[0] = 2 * r * Math.random() - r;
      x[1] = 2 * r * Math.random() - r;
      x[2] = 2 * r * Math.random() - r;
      xsq = x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
    } while (xsq >= rsq);
  }

  function randomVectorOnSphere(x, r) {
    // Get a vector in the unit disc. Marsaglia (1972); see
    // <http://mathworld.wolfram.com/SpherePointPicking.html>.
    var xsq;
    do {
      x[0] = 2 * Math.random() - 1;
      x[1] = 2 * Math.random() - 1;
      xsq = x[0] * x[0] + x[1] * x[1];
    } while (xsq >= 1.0);
    // Apply Marsaglia's transformation.
    var t = 2 * r * Math.sqrt(1 - xsq);
    x[0] = t * x[0];
    x[1] = t * x[1];
    x[2] = r - 2 * r * xsq;
  }

  // Initialization, input and animation.

  function setup() {
    var updatePending = false;
    var lastUpdateTime = null;
    var pickAnchor = null;
    var eventPoint = null;

    var canvas = document.getElementById("diagram");
    var zoomInput = document.getElementById("zoom");
    var buttons = document.querySelectorAll("#controls>button");

    function pick1() {
      return eventPoint && pick(eventPoint, canvas.getBoundingClientRect(), x);
    }

    function swipeStart(point) {
      eventPoint = point;
      var pickPoint = pick1();
      pickAnchor = pickPoint && unrotate(pickPoint);
      // If point hits the object, return true to inhibit touch actions.
      return pickPoint;
    }

    function swipeMove(point) {
      eventPoint = point;
      invalidate();
    }

    function swipeEnd() {
      eventPoint = null;
    }

    function zoomRangeInput(e) {
      setZoom(+zoomInput.value);
      invalidate();
    }

    function buttonClick(e) {
      if (e.target.hasAttribute("data-index")) {
        setSystemIndex(+e.target.getAttribute("data-index"));
        invalidate();
      }
    }

    function resize() {
      var context = canvas.getContext("2d");
      var rect = canvas.getBoundingClientRect();
      var ratio = window.devicePixelRatio;
      var width =
          Math.round(ratio * rect.right) - Math.round(ratio * rect.left);
      var height =
          Math.round(ratio * rect.bottom) - Math.round(ratio * rect.top);
      if (canvas.width !== width) {
        canvas.width = width;
      }
      if (canvas.height !== height) {
        canvas.height = height;
      }
      invalidate();
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
      var updateInterval =
          Math.max(1, Math.min(30, updateTime - lastUpdateTime));
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
      renderModel(canvas);
      // Ensure an update is pending if the ball is in motion.
      if (inMotion) {
        invalidate();
      }
    }

    // Initialize the model.
    setSystemIndex(3);
    randomVectorInBall(u, 3.14159265);
    randomVectorOnSphere(w, 0.001);

    // Set event handlers.
    makeSwipable(canvas, swipeStart, swipeEnd, swipeMove);
    window.addEventListener("resize", resize);
    zoomInput.addEventListener("input", zoomRangeInput);
    for (var i = 0; i !== buttons.length; ++i) {
      buttons[i].setAttribute("data-index", i);
      buttons[i].addEventListener("click", buttonClick);
    }

    // Start animating.
    zoomRangeInput();
    resize();
  }

  document.addEventListener("DOMContentLoaded", setup, true);
}());
