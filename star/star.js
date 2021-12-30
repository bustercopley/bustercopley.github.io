// -*- coding: utf-8; -*-

// Copyright 2021, Richard Copley <buster at buster dot me dot uk>

// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at

//     http://www.apache.org/licenses/LICENSE-2.0

// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

(function() {
"use strict";

let context;
let width;
let height;
let m = 7;
let n = 3;
let input;
const trace = new Array(0);

document.addEventListener("DOMContentLoaded", function() {
  context = document.getElementById("canvas").getContext("2d");
  input = document.getElementById("input");
  input.addEventListener("input", inputEvent);
  inputEvent();
  requestAnimationFrame(render);
});

function gcd(a, b) { return (b | 0) ? gcd(b | 0, (a | 0) % (b | 0)) : (a | 0); }

// Validate and apply change of input
function inputEvent() {
  let state = 0, i = 0, j = 0;
  for (let c of input.value) {
    if (state === 0 && c === '{') {
      ++state;
    } else if (state === 1 && ('0' <= c && c <= '9')) {
      if (i >= 100) {
        state == 5;
      } else {
        i = 10 * i + (c - '0');
      }
    } else if (state === 1 && c === '/') {
      ++state;
    } else if (state === 2 && ('0' <= c && c <= '9')) {
      if (j >= 100) {
        state == 5;
      } else {
        j = 10 * j + (c - '0');
      }
    } else if (state === 2 && c === '}') {
      ++state;
    } else if (c != ' ') {
      state = undefined;
    }
  }
  let msg;
  if (state === 3) {
    if (1 <= j && j < i && gcd(i, j) === 1) {
      m = i;
      n = j;
      trace.length = 0;
    } else if (!i || !j) {
      msg = "specify two nonzero numbers"
    } else if (j > i) {
      msg = "𝑛 ≥ 𝑚";
    } else if (gcd(i, j) != 1) {
      msg = "𝑚 and 𝑛 are not coprime";
    } else {
      msg = "𝑚 and 𝑛 are not valid";
    }
  } else if (state < 3) {
    msg = "incomplete";
  } else if (state == 5) {
    msg = "too big!";
  } else {
    msg = "syntax error";
  }

  if (msg) {
    input.style["background"] = "pink";
    document.getElementById("error").firstChild.nodeValue =
        "Using {" + m + "/" + n + "} (" + msg + ")";
  } else {
    input.style["background"] = "";
    document.getElementById("error").firstChild.nodeValue = " ";
  }
}

function render(time) {
  const omega = 0.0005; // kiloradians per second
  const tau = 0.0005;   // kilostages per second
  const ratio = 0.8;

  const rect = context.canvas.getBoundingClientRect();
  const dpr = window.devicePixelRatio;
  if (rect.width !== width * dpr || rect.height != height * dpr) {
    width = rect.width * dpr;
    height = rect.height * dpr;
    context.canvas.width = width;
    context.canvas.height = height;
    trace.length = 0;
  }

  context.lineWidth = 2.0;
  context.resetTransform();
  context.clearRect(0.0, 0.0, width, height);

  const side = Math.min(width, height);
  const theta = omega * time;
  const phi = theta * (1.0 - m / n);

  const c = new DOMPoint(0.5 * width, 0.5 * height);
  const r0 = 0.4 * side;
  const r1 = r0 * n / m;

  function fade(t, t0, t1) {
    return Math.max(0.0, Math.min(1.0, t - t0, t1 - t));
  }

  function getRGBA(rgb, alpha) { return "rgba(" + rgb + "," + alpha + ")"; }

  const t = (tau * time) % 26.0;
  const A = fade(t, 2.0, 5.0) + fade(t, 12.0, 17.0);
  const B = fade(t, 6.0, 9.0) + fade(t, 14.0, 17.0);
  const C = fade(t, 10.0, 17.0);
  const D = fade(t, 18.0, 23.0);
  const E = fade(t, 20.0, 23.0);

  // Comments assume symbol is "{7/3}" for concreteness
  // 7-pointed star
  if (C) {
    if (!trace.length) {
      const steps = 300;
      for (let k = 0; k != steps; ++k) {
        const theta = k * n * 2.0 * Math.PI / steps;
        const phi = theta * (1.0 - m / n);
        const p = new DOMPoint(c.x + (r0 - r1) * Math.cos(theta),
                               c.y + (r0 - r1) * Math.sin(theta));
        const q = new DOMPoint(p.x + ratio * r1 * Math.cos(phi),
                               p.y + ratio * r1 * Math.sin(phi));
        trace.push(q);
      }
    }

    context.beginPath();
    for (let q of trace) {
      context.lineTo(q.x, q.y);
    }
    context.closePath();
    context.strokeStyle = getRGBA("80,80,80", C);
    context.stroke();
  }

  // Outer large circle
  if (D) {
    context.beginPath();
    context.arc(c.x, c.y, r0, 0.0, 2.0 * Math.PI);
    context.strokeStyle = getRGBA("0,0,255", D);
    context.stroke();
  }

  for (let i = 0; i != m - n; ++i) {
    const theta = omega * time + 2.0 * Math.PI * i / (m - n);
    const p = new DOMPoint(c.x + (r0 - r1) * Math.cos(theta),
                           c.y + (r0 - r1) * Math.sin(theta));

    const phi = (theta - 2.0 * Math.PI * i / (m - n)) * (1.0 - m / n);
    // One small circle (or all 4 small circles if E)
    if ((!i && D) || E) {
      context.beginPath();
      context.arc(p.x, p.y, r1, 0.0, 2.0 * Math.PI);
      context.fillStyle = getRGBA("0,0,0", 0.5 * ((!i && D) || E));
      context.fill();
      // One small radius
      const q =
          new DOMPoint(p.x + r1 * Math.cos(phi), p.y + r1 * Math.sin(phi));
      context.beginPath();
      context.moveTo(p.x, p.y);
      context.lineTo(q.x, q.y);
      context.strokeStyle = getRGBA("0,0,0", (!i && D) || E);
      context.stroke();
    }

    // Four triangles
    if (A) {
      context.beginPath();
      for (let j = 0; j != n; ++j) {
        const psi = phi + 2.0 * Math.PI * j / n;
        const q = new DOMPoint(p.x + ratio * r1 * Math.cos(psi),
                               p.y + ratio * r1 * Math.sin(psi));
        context.lineTo(q.x, q.y);
      }
      context.closePath();
      context.strokeStyle = getRGBA("0,180,0", A);
      context.stroke();
    }
  }

  // Three squares
  if (B) {
    context.strokeStyle = getRGBA("0,0,240", B);
    for (let j = 0; j != n; ++j) {
      context.beginPath();
      for (let i = 0; i != m - n; ++i) {
        const theta = omega * time + 2.0 * Math.PI * i / (m - n);
        const phi = (theta - 2.0 * Math.PI * i / (m - n)) * (1.0 - m / n);
        const psi = phi + 2.0 * Math.PI * j / n;
        const p = new DOMPoint(c.x + (r0 - r1) * Math.cos(theta),
                               c.y + (r0 - r1) * Math.sin(theta));
        const q = new DOMPoint(p.x + ratio * r1 * Math.cos(psi),
                               p.y + ratio * r1 * Math.sin(psi));
        context.lineTo(q.x, q.y);
      }
      context.closePath();
      context.stroke();
    }
  }

  // Spots (always drawn)
  context.fillStyle = "orange";
  for (let i = 0; i != m - n; ++i) {
    const theta = omega * time + 2.0 * Math.PI * i / (m - n);
    const p = new DOMPoint(c.x + (r0 - r1) * Math.cos(theta),
                           c.y + (r0 - r1) * Math.sin(theta));
    for (let j = 0; j != n; ++j) {
      const phi = (theta - 2.0 * Math.PI * i / (m - n)) * (1.0 - m / n);
      const psi = phi + 2.0 * Math.PI * j / n;
      const q = new DOMPoint(p.x + ratio * r1 * Math.cos(psi),
                             p.y + ratio * r1 * Math.sin(psi));
      context.beginPath();
      context.arc(q.x, q.y, 15.0, 0.0, 2.0 * Math.PI);
      context.fill();
    }
  }

  requestAnimationFrame(render);
}
}());
