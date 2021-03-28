// -*- coding: utf-8; -*-

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

(function () {
  "use strict";

  // Create asm.js Body object.
  var heap = new ArrayBuffer(0x10000); // 64 kilobytes (the minimum).
  var body = Body(window, null, heap);

  // Bodies (byte offsets from start of heap):
  var MAIN = 128; // Main body.

  // Fields (byte offsets from start of body):
  var PX = 0;  // Linear position.
  var PU = 16; // Angular position.
  var PW = 32; // Angular velocity.
  var PM = 64; // Modelview matrix.

  // Views into the main body.
  var x = new Float32Array(heap, MAIN + PX, 4);  // Centre of sphere
  var u = new Float32Array(heap, MAIN + PU, 4);  // Angular position
  var w = new Float32Array(heap, MAIN + PW, 4);  // Angular momentum
  var m = new Float32Array(heap, MAIN + PM, 16); // Modelview matrix

  // Viewing frustum dimensions. (The "t" stands for "tank".)
  var tz = 6;   // Distance from eye to screen
  var ts = 2.5; // Half of the shorter side of the front wall of the screen
  var z = -tz;

  // Material colour and lighting parameters.
  var lightPosition = [
    -1, -2, z + 2.25,
    +1, +2, z + 2.25,
    -3, +4, z + 1.50,
    +3, -4, z + 1.50
  ];
  var diffuseReflectance = [
    0.4, 0.4, 0.4,
    0.7, 0.7, 0.7,
    0.4, 0.4, 0.4,
    0.5, 0.5, 0.5
  ];

  // Operations on a three dimensional inner product space.
  function dot(u, v) { return u[0] * v[0] + u[1] * v[1] + u[2] * v[2]; }

  function cross(u, v) {
    var r = Array(3);
    r[0] = u[1] * v[2] - u[2] * v[1];
    r[1] = u[2] * v[0] - u[0] * v[2];
    r[2] = u[0] * v[1] - u[1] * v[0];
    return r;
  }

  // Multiply a three dimensional vector by the rotation part of the modelview
  // matrix.
  function rotate(x) {
    var result = new Float32Array(3);
    result[0] = m[0] * x[0] + m[4] * x[1] + m[8] * x[2];
    result[1] = m[1] * x[0] + m[5] * x[1] + m[9] * x[2];
    result[2] = m[2] * x[0] + m[6] * x[1] + m[10] * x[2];
    return result;
  }

  // Multiply a three dimensional vector by the inverse of the rotation part of
  // the modelview matrix.
  function unrotate(x) {
    var result = new Float32Array(3);
    result[0] = m[0] * x[0] + m[1] * x[1] + m[2] * x[2];
    result[1] = m[4] * x[0] + m[5] * x[1] + m[6] * x[2];
    result[2] = m[8] * x[0] + m[9] * x[1] + m[10] * x[2];
    return result;
  }

  // Set the angular velocity of the body so that the direction vector anchor
  // (in model space) will point in direction target (in world space) after
  // time dt.
  function spinTo(anchor, target, dt) {
    var current = rotate(anchor);
    var crs = cross(current, target);
    // sin²(θ) = ⏸crs⏸², where θ is the angle between current and target.
    var ssq = dot(crs, crs);
    var thetaOverSinTheta;
    if (ssq < 1.0e-8) { thetaOverSinTheta = 1 - ssq / 6; }
    else {
      var s = Math.sqrt(ssq);
      thetaOverSinTheta = Math.asin(s) / s;
    }
    // The required value is in the direction of crs and has magnitude θ/dt.
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
    // L: the line in world space whose points all project to p in client space.
    // T: the point in world space where L intersects the unit sphere centre X,
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
    if (scrim < 0) { return null; }

    // Take the smaller of the solutions to get the intersection closer to Q.
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
    } while (xsq >= rsq);
  }

  function randomVectorOnSphere(x, r) {
    // Get a vector on the unit sphere in R³ (Marsaglia (1972)), see [1].
    // [1] <http://mathworld.wolfram.com/SpherePointPicking.html>

    // Get a vector in the unit disc.
    var xsq;
    do {
      x[0] = 2 * Math.random() - 1;
      x[1] = 2 * Math.random() - 1;
      xsq = x[0] * x[0] + x[1] * x[1];
    } while (xsq >= 1);

    // Apply Marsaglia's transformation.
    var t = 2 * r * Math.sqrt(1 - xsq);
    x[0] = t * x[0];
    x[1] = t * x[1];
    x[2] = r - 2 * r * xsq;
  }

  function showError(text) {
    document.getElementById("status").firstChild.nodeValue = text;
  }

  var DISABLE_REQUEST_CACHING = false;
  function download(url, type, handler, data) {
    var r = new XMLHttpRequest();
    r.onload = function () {
      if (this.status !== 200) {
        showError(url + ": " + this.status.toString() + " " + this.statusText);
      }
      else {
        handler(this.response, data);
      }
    };
    r.open("GET", url + (DISABLE_REQUEST_CACHING ? "?v=2" : ""), true);
    r.responseType = type;
    r.send();
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
    if (!program) { return false; }

    gl.clearColor(1, 1, 1, 1);
    gl.enable(gl.CULL_FACE);
    gl.enable(gl.BLEND);
    gl.blendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA);

    gl.uniform3fv(gl.getUniformLocation(program, "l"), lightPosition);

    return resize(true);
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

        var tw = ts * w;
        var th = ts * h;
        var tt = 2 * tz * Math.min(w, h);
        // This projection matrix to give z = 0 in clip coordinates, in order to
        // avoid clipping in the z dimension. We have no use for the z coordinate
        // as long as depth testing is not enabled.
        var projectionMatrix = [
          tt / tw, 0, 0,  0,
          0, tt / th, 0,  0,
          0,    0,    0, -1,
          0,    0,    0,  0
        ];

        gl.viewport(0, 0, w, h);
        setMatrixUniform(gl, program, "p", projectionMatrix);
      }
      return true;
    }
    catch (e) {
      showError(e.toString());
      return false;
    }
  }

  function loadVertexBuffer(vbuffer) {
    var vertexBuffer = gl.createBuffer();
    var attributes = new Float32Array(vbuffer);
    gl.bindBuffer(gl.ARRAY_BUFFER, vertexBuffer);
    gl.bufferData(gl.ARRAY_BUFFER, attributes, gl.STATIC_DRAW);

    function attrib(name, offset, size) {
      var location = gl.getAttribLocation(program, name);
      gl.vertexAttribPointer(location, size, gl.FLOAT, false, 0, offset);
      gl.enableVertexAttribArray(location);
    }

    var chunkSize = attributes.byteLength / 3;
    attrib("x", 0 * chunkSize, 3);
    attrib("n", 1 * chunkSize, 3);
    attrib("h", 2 * chunkSize, 3);

    vertexCount = chunkSize / (3 * 4);

    // mesh.dat:
    //   12 faces
    //   5 triangles per face
    //   3 triangle vertices per triangle
    //   3 attributes per vertex chunks
    //   3 floats per attribute
    //   4 bytes per float
    // Total bytes 6480
  }

  function setColor(color) {
    var color = [0.85, 0.7, 0.2];
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
        if (gl.getShaderParameter(shader, gl.COMPILE_STATUS)) { return true; }
        var s = "Shader compilation failed";
        switch (type) {
        case gl.VERTEX_SHADER: s = "Vertex " + s; break;
        case gl.FRAGMENT_SHADER: s = "Fragment " + s; break;
        }
        console.log(gl.getShaderInfoLog(shader));
        showError(s);
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
      if (gl.getProgramParameter(program, gl.LINK_STATUS)) { return true; }
      showError("Shader program linking failed");
    }

    var program = gl.createProgram();
    if (!addShader(vshader, gl.VERTEX_SHADER)) return;
    if (!addShader(fshader, gl.FRAGMENT_SHADER)) return;
    gl.linkProgram(program);
    if (checkProgram()) {
      gl.useProgram(program);
      return program;
    }
  }

  function setMatrixUniform(gl, program, name, matrix) {
    var location = gl.getUniformLocation(program, name);
    gl.uniformMatrix4fv(location, false, matrix);
  }

  // HTML event handling.
  function makeSwipable(element, handleStart, handleEnd, handleMove) {
    function getEventPoint(event) { return [event.clientX, event.clientY]; }

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

    function onMousemove(event) { handleMove(getEventPoint(event)); }

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
        if (touchId === e.identifier) { handleMove(getEventPoint(e)); }
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
    var model;

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
      invalidate();
    }

    function swipeEnd() { eventPoint = null; }

    // Ensure exactly one update is pending.
    function invalidate() {
      if (!updatePending) {
        updatePending = true;
        requestAnimationFrame(update);
      }
    }

    // Update the view.
    function update(time) {
      updatePending = false;
      var updateInterval = Math.max(1, Math.min(30, time - lastUpdateTime));
      var pickPoint = pick1();
      lastUpdateTime = time;
      if (!pickPoint) { pickAnchor = null; }
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
        if (inMotion) { invalidate(); }
      }
      catch (e) {
        showError(e.toString());
      }
    }

    function updateLoadState() {
      try {
        if (vshader && fshader && model) {
          canvas = document.getElementById("canvas");
          if (canvas) {
            makeSwipable(canvas, swipeStart, swipeEnd, swipeMove);
            gl = canvas.getContext("experimental-webgl");
            if (initializeGraphics(vshader, fshader)) {
              x[2] = z;
              randomVectorInBall(u, 3.14159265);
              randomVectorOnSphere(w, 0.0005);
              loadVertexBuffer(model);
              setColor();
              invalidate()
            }
          }
        }
      } catch (e) {
        showError(e.toString());
      }
    }

    function start() {
      function setupButton(id, handler) {
        var button = document.getElementById(id);
        if (button) { button.addEventListener("click", handler); }
      }
      download("vertex-shader.glsl", "text", function(response) {
        vshader = response;
        updateLoadState();
      });
      download("fragment-shader.glsl", "text", function(response) {
        fshader = response;
        updateLoadState();
      });
      download("attributes.dat", "arraybuffer", function (response) {
        model = response;
        updateLoadState();
      });
      updateLoadState();
    };

    start();
  }

  document.addEventListener("DOMContentLoaded", function() { setup(); }, true);
}());
