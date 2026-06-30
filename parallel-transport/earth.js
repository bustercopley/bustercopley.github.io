(function() {
"use strict";

let gl;
let rectangles;
let angle = 0.0;
let previousTime;
let width;
let height;

const r = 1.0;               // scale factor
const theta = Math.PI / 4.0; // latitude where cone is tangent to globe
const omega = 0.00025;       // angular speed
const z0 = -1.00;            // z coordinate of screen and of near clip plane
const z1 = -40.00;           // z coordinate of far clip plane
const position = new DOMPoint(0.0, -1.0, -30.0, 1.0); // center of globe

document.addEventListener("DOMContentLoaded", function() {
  gl = document.getElementById("canvas").getContext("webgl2");
  gl.clearColor(0.0, 0.0, 0.0, 1.0);
  gl.enable(gl.CULL_FACE);
  gl.cullFace(gl.BACK);
  gl.enable(gl.BLEND);
  gl.blendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA);
  loadAll();
});

async function loadAll() {
  // Image credit: Natural Earth III <http://www.shadedrelief.com/natural3/>
  rectangles = await Promise.all([
    loadRectangle(gl, 48, 96, "globev.glsl", "globef.glsl", "earth-1k.jpg"),
    loadRectangle(gl, 48, 96, "conev.glsl", "conef.glsl"),
    loadRectangle(gl, 48, 96, "flatconev.glsl", "conef.glsl"),
  ]);
  // Now everything is loaded, light the scene
  const theta = 23.4999 * (Math.PI / 180.0); // Northern Summer ...
  const phi = 0.25 * Math.PI;                // ... from the right
  rectangles.map(rectangle => rectangle.setLight([
    -Math.cos(theta) * Math.sin(phi),
    -Math.sin(theta),
    -Math.cos(theta) * Math.cos(phi),
    0.0,
  ]));
  // Start animating
  requestAnimationFrame(render);
}

async function loadRectangle(gl, m, n, vurl, furl, image) {
  const [vtext, ftext, imageData] = await Promise.all([
    fetch(vurl).then(response => response.text()),
    fetch(furl).then(response => response.text()),
    loadImage(image),
  ]);
  return createRectangle(gl, m, n, vurl, vtext, furl, ftext, imageData);
}

async function loadImage(url) {
  if (url) {
    const image = new Image;
    image.src = url;
    await image.decode();
    const canvas = document.createElement("canvas");
    const context = canvas.getContext("2d");
    const width = image.width;
    const height = image.height;
    canvas.width = width;
    canvas.height = height;
    context.drawImage(image, 0, 0, width, height);
    return context.getImageData(0, 0, width, height);
  }
}

function createRectangle(gl, m, n, vurl, vtext, furl, ftext, imageData) {
  const program = gl.createProgram();
  const vertexArray = gl.createVertexArray();
  let texture;

  const attributeLocations = {
    u: undefined,
  };

  const uniformLocations = {
    sampler: undefined,
    projection: undefined,
    modelview: undefined,
    lightDirection: undefined,
    width: undefined,
    height: undefined,
    theta: undefined,
    phi: undefined,
  };

  addShader(gl, program, gl.VERTEX_SHADER, vurl, vtext);
  addShader(gl, program, gl.FRAGMENT_SHADER, furl, ftext);
  gl.linkProgram(program);
  if (!gl.getProgramParameter(program, gl.LINK_STATUS)) {
    throw new Error("Shader link failed: " + gl.getProgramInfoLog(program));
  }
  gl.useProgram(program);

  for (const name in attributeLocations) {
    const location = gl.getAttribLocation(program, name);
    attributeLocations[name] = location;
  }

  for (const name in uniformLocations) {
    const location = gl.getUniformLocation(program, name);
    uniformLocations[name] = location;
  }

  if (imageData) {
    texture = gl.createTexture();
    gl.activeTexture(gl.TEXTURE0);
    gl.uniform1i(uniformLocations.sampler, 0);
    gl.bindTexture(gl.TEXTURE_2D, texture);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.REPEAT);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER,
                     gl.LINEAR_MIPMAP_LINEAR);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.LINEAR);
    gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, imageData.width, imageData.height,
                  0, gl.RGBA, gl.UNSIGNED_BYTE, imageData.data);
    gl.generateMipmap(gl.TEXTURE_2D);
  }

  gl.bindVertexArray(vertexArray);
  gl.bindBuffer(gl.ARRAY_BUFFER, gl.createBuffer());
  gl.bufferData(gl.ARRAY_BUFFER, makeRectangleMesh(m, n), gl.STATIC_DRAW);
  gl.vertexAttribPointer(attributeLocations.u, 2, gl.FLOAT, false, 2 * 4, 0);
  gl.enableVertexAttribArray(attributeLocations.u);

  function setModelview(p, q, r) {
    // Given
    //   p, a vector (the position),
    //   q, a quaternion (the orientation),
    //   r, a scalar (the radius),
    // compute the matrix taking x to rqxq⁻¹ + p.
    const modelview = new Float32Array([
      r * (1.0 - 2.0 * (q.y * q.y + q.z * q.z)),
      r * (2.0 * (q.x * q.y + q.w * q.z)),
      r * (2.0 * (q.x * q.z - q.w * q.y)),
      0.0,
      r * (2.0 * (q.x * q.y - q.w * q.z)),
      r * (1.0 - 2.0 * (q.x * q.x + q.z * q.z)),
      r * (2.0 * (q.y * q.z + q.w * q.x)),
      0.0,
      r * (2.0 * (q.x * q.z + q.w * q.y)),
      r * (2.0 * (q.y * q.z - q.w * q.x)),
      r * (1.0 - 2.0 * (q.x * q.x + q.y * q.y)),
      0.0,
      p.x,
      p.y,
      p.z,
      1.0,
    ]);

    gl.useProgram(program);
    gl.uniformMatrix4fv(uniformLocations.modelview, false, modelview);
    gl.uniform1f(uniformLocations.phi, angle * Math.sin(theta));
    // theta is used in conev.glsl (latitude where cone is tangent to sphere)
    gl.uniform1f(uniformLocations.theta, theta);
  }

  function setViewport(width, height) {
    const scale = 10.0 * Math.min(width, height);
    const x0 = width / scale;
    const y0 = height / scale;
    // No surprises: this is glFrustum(-x0, x0, -y0, y0, -z0, -z1).
    // This is overkill. The z coordinate could just as well be zero, in clip
    // space (since no clipping in fact occurs) and in NDC (since depth testing
    // is not enabled).
    const projection = new Float32Array([
      -z0 / x0, 0.0, 0.0, 0.0,                //
      0.0, -z0 / y0, 0.0, 0.0,                //
      0.0, 0.0, -(z0 + z1) / (z1 - z0), -1.0, //
      0.0, 0.0, 2 * z0 * z1 / (z1 - z0), 0.0
    ]);
    gl.useProgram(program);
    gl.uniformMatrix4fv(uniformLocations.projection, false, projection);
    gl.uniform1f(uniformLocations.width, width);
    gl.uniform1f(uniformLocations.height, height);
  }

  function setLight(direction) {
    gl.useProgram(program);
    gl.uniform4fv(uniformLocations.lightDirection, direction);
  }

  function render() {
    gl.bindVertexArray(vertexArray);
    gl.useProgram(program);
    if (texture !== undefined) {
      gl.activeTexture(gl.TEXTURE0);
      gl.bindTexture(gl.TEXTURE_2D, texture);
    }
    for (let i = 0; i != m; ++i) {
      gl.drawArrays(gl.TRIANGLE_STRIP, i * 2 * (n + 1), 2 * (n + 1));
    }
  }

  return Object.freeze({
    setModelview: setModelview,
    setViewport: setViewport,
    setLight: setLight,
    render: render,
  });
}

function makeRectangleMesh(m, n) {
  const vertexAttribBuffer = new Float32Array(m * 2 * (n + 1) * 2);
  let offset = 0;
  for (let i = 0; i != m; ++i) {
    for (let j = 0; j != n + 1; ++j) {
      const u = j / n;
      const v0 = i / m;
      const v1 = (i + 1) / m;
      vertexAttribBuffer[offset++] = u;
      vertexAttribBuffer[offset++] = v0;
      vertexAttribBuffer[offset++] = u;
      vertexAttribBuffer[offset++] = v1;
    }
  }
  return vertexAttribBuffer;
}

function addShader(gl, program, type, url, text) {
  const shader = gl.createShader(type);
  gl.shaderSource(shader, text);
  gl.compileShader(shader);
  if (!gl.getShaderParameter(shader, gl.COMPILE_STATUS)) {
    throw new Error((type == gl.VERTEX_SHADER ? "Vertex" : "Fragment") +
                    " shader " + url + ": " + gl.getShaderInfoLog(shader));
  }
  gl.attachShader(program, shader);
  return true;
}

function render(time) {
  const rect = gl.canvas.getBoundingClientRect();
  const dpr = window.devicePixelRatio;
  if (rect.width !== width * dpr || rect.height != height * dpr) {
    width = rect.width * dpr;
    height = rect.height * dpr;
    gl.canvas.width = width;
    gl.canvas.height = height;
    gl.viewport(0.0, 0.0, width, height);
    rectangles.map(rectangle => rectangle.setViewport(width, height));
  }

  if (previousTime) {
    angle = (angle + omega * (time - previousTime)) % (2.0 * Math.PI);
  }
  previousTime = time;

  const c = Math.cos(0.5 * angle);
  const s = Math.sin(0.5 * angle);
  const q = new DOMPoint(0.0, -s, 0.0, c);

  rectangles.map(rectangle => rectangle.setModelview(position, q, r));
  gl.clear(gl.COLOR_BUFFER_BIT);
  rectangles.map(rectangle => rectangle.render());
  requestAnimationFrame(render);
}
}());
