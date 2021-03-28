// Copyright 2020-2021 Richard Copley

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

  const MAX_SIZE = 1024;
  const MAX_TILES = 2 * MAX_SIZE * MAX_SIZE;
  const MAX_VERTICES = 4 * MAX_TILES;

  let gl;
  let context;
  let trim;
  let width;
  let height;
  let resizePending = true;
  let updatePending = false;
  let paused = false;
  let startTime;
  let pauseTime;
  let rendered;
  let size;
  let rate = 1; // milliseconds per step
  let newRate;
  let markCount;
  let maxTextureSize;

  const shaders = new Array(2);
  const FLOATS = Float32Array.BYTES_PER_ELEMENT;
  const vertexAttribBuffer = new Float32Array(4 * MAX_VERTICES);
  const indexBuffer = new Uint32Array(6 * MAX_TILES); // 2 triangles per tile
  const grid = new Uint8Array(MAX_SIZE * MAX_SIZE);
  const colors = [
    new Float32Array([0.412, 0.612, 1.000, 1]), // 0: blue
    new Float32Array([0.839, 0.698, 0.317, 1]), // 1: yellowy beige
    new Float32Array([0.847, 0.455, 0.376, 1]), // 2: red
    new Float32Array([0.455, 0.843, 0.373, 1]), // 3: green
    new Float32Array([0.976, 0.651, 0.216, 1]), // 4: orange yellow
    new Float32Array([0.376, 0.376, 0.376, 1]), // 5: dark grey
    new Float32Array([0.000, 0.000, 0.000, 1]), // 6: black
    new Float32Array([1.000, 1.000, 1.000, 1]), // 7: white
    new Float32Array([0.000, 0.000, 0.000, 0]), // 6: transparent
    new Float32Array([0.976, 0.918, 0.173, 1]), // 9: lemon yellow
  ];

  const tileColors = new Array(8);
  for (let n = 0; n != tileColors.length; ++n) {
    tileColors[n] = new Float32Array(4);
  }
  const edgeColor = new Float32Array(4);
  const outerEdgeColor = new Float32Array(4);
  const innerEdgeColor = new Float32Array(4);
  const arrowColors = [new Float32Array(4), new Float32Array(4)];
  const gridColor = new Float32Array(4);

  // fixed colors
  copyColor(tileColors[0], colors[0]); // tile 0 fill color
  copyColor(tileColors[1], colors[1]); // tile 1 fill color
  copyColor(tileColors[2], colors[2]); // tile 2 fill color
  copyColor(tileColors[3], colors[3]); // tile 3 fill color
  const gridStroke = colorToCSS(colors[6]);
  const edgeStroke = colorToCSS(colors[5]);

  for (let n = 0; n != MAX_TILES; ++n) {
    // element indices (one quad -> two triangles -> six vertices)
    indexBuffer[6 * n + 0] = 4 * n + 0;
    indexBuffer[6 * n + 1] = 4 * n + 1;
    indexBuffer[6 * n + 2] = 4 * n + 3;
    indexBuffer[6 * n + 3] = 4 * n + 3;
    indexBuffer[6 * n + 4] = 4 * n + 1;
    indexBuffer[6 * n + 5] = 4 * n + 2;
  }

  const transforms = [
    new DOMMatrix([1, 0, 0, 1, 0, 0]),
    new DOMMatrix([0, -1, 1, 0, 0, 0]),
    new DOMMatrix([0, 1, -1, 0, 0, 0]),
    new DOMMatrix([-1, 0, 0, -1, 0, 0]),
  ];

  // animate steps
  const SHOW_MOVE = 0;           // slide tiles along their proper direction
  const SHOW_MOVE_PAUSE = 1;     //
  const SHOW_EMPTY = 2;          // highlight unoccupied 2x2 blocks
  const SHOW_EMPTY_PAUSE = 3;    //
  const SHOW_INSERT = 4;         // add new pairs of tiles
  const SHOW_INSERT_PAUSE = 5;   //
  const SHOW_EXPAND = 6;         // expand the grid
  const SHOW_EXPAND_PAUSE = 7;   //
  const SHOW_ZAP = 8;            // highlight pairs of tiles that can't move
  const SHOW_ZAP_PAUSE = 9;      //
  const SHOW_DELETE = 10;        // delete the marked tiles
  const SHOW_DELETE_PAUSE = 11;  //
  const SHOW_STEPS = 12;         // total number of steps

  function fadeAlpha(data, color, alpha0, alpha1, theta) {
    data[0] = color[0];
    data[1] = color[1];
    data[2] = color[2];
    data[3] = alpha0 + theta * (alpha1 - alpha0);
  }

  function fade(data, color1, color2, theta) {
    data[0] = color1[0] + theta * (color2[0] - color1[0]);
    data[1] = color1[1] + theta * (color2[1] - color1[1]);
    data[2] = color1[2] + theta * (color2[2] - color1[2]);
    data[3] = color1[3] + theta * (color2[3] - color1[3]);
  }

  function copyColor(data, other) {
    data[0] = other[0];
    data[1] = other[1];
    data[2] = other[2];
    data[3] = other[3];
  }

  function colorToCSS(data) {
    return "rgba(" + ((data[0] * 255) | 0) + "," + ((data[1] * 255) | 0) + "," +
      ((data[2] * 255) | 0) + "," + data[3] + ")";
  }

  function showSize(text) {
    document.getElementById("size").firstChild.nodeValue = text;
  }

  let bufferIndices = new Array(2);
  function vertex(x, y, u, v, c) {
    // x, y: vertex position
    // u, v: texture coordinate within texture section
    // c: texture section selector (0 - 7)
    // marked (zapped/deleted/empty) tiles have c >> 2 == 1
    // draw marked tiles first (at index bufferIndices[1] < bufferIndices[0])
    vertexAttribBuffer[bufferIndices[c >> 2]++] = x;
    vertexAttribBuffer[bufferIndices[c >> 2]++] = y;
    vertexAttribBuffer[bufferIndices[c >> 2]++] = (2 * (c >> 2) + u) / 3;
    vertexAttribBuffer[bufferIndices[c >> 2]++] = (8 * (c & 3) + v) / 32
  }

  function render(time) {
    // resize canvas?
    const resized = resizePending;
    if (resizePending) {
      const rc = context.canvas.getBoundingClientRect();
      const pr = window.devicePixelRatio;
      width = Math.round(pr * rc.right) - Math.round(pr * rc.left);
      height = Math.round(pr * rc.bottom) - Math.round(pr * rc.top);
      context.canvas.width = width;
      context.canvas.height = height;
      gl.canvas.width = width;
      gl.canvas.height = height;
      gl.viewport(0, 0, width, height);
      resizePending = false;
    }

    // restart animation?
    if (startTime === undefined) {
      startTime = time;
      pauseTime = undefined;
      rendered = 0;
      grid.fill(0);
      size = 1;
      showSize(size + 1);
    }

    // re-render at paused time (after resize)?
    if (pauseTime !== undefined) {
      startTime += time - pauseTime;
      pauseTime = undefined;
    }

    // animation parameter
    let target = (time - startTime) / rate;

    // drop frames?
    if (target > rendered + SHOW_STEPS) {
      // skip ahead by a whole number of complete cycles
      const targetCycle = (target / SHOW_STEPS) | 0;
      const renderedCycle = (rendered / SHOW_STEPS) | 0;
      const cycles = Math.max(1, Math.min(8, targetCycle - renderedCycle));
      target = (renderedCycle + cycles) * SHOW_STEPS + SHOW_INSERT_PAUSE;
    } else if (target > rendered + 0.5) {
      // two complete steps per frame when we're skipping frames
      target = rendered + 2;
    }

    // stop when maximum size reached, after inserting new tiles
    const end = (MAX_SIZE - 3) * SHOW_STEPS + SHOW_INSERT_PAUSE;
    if (target >= end) {
      target = end;
      pause(true);
    }

    // change speed?
    if (newRate !== undefined) {
      rate = newRate;
    }
    newRate = undefined;

    startTime = time - target * rate;

    if (!resized && rendered > target && !(rendered & 1)) {
      // no changes to draw since last frame (odd steps are *_PAUSE)
      return;
    }

    // do all computations needed before rendering the target step
    let sizeChanged;
    let skip;
    for (; rendered <= target; ++rendered) {
      switch (rendered % SHOW_STEPS) {
      case SHOW_MOVE:
        markCount = 0;
        // expand grid and move tiles
        for (let i = 0; i != size; ++i) {
          for (let j = 0; j != size; ++j) {
            const tiles = grid[i * MAX_SIZE + j];
            for (let flip = 0; flip != 2; ++flip) {
              // is this a tile that has not been moved yet?
              if (tiles >> flip & 1) {
                // add tile in new position in expanded grid (in top 3 bits)
                switch (tiles >> 1 & 2 | flip) {
                case 0: // up (blue)
                  grid[i * MAX_SIZE + j] |= (1 << 5);
                  break;
                case 1: // down (green)
                  grid[(i + 1) * MAX_SIZE + j + 1] |= (2 << 5);
                  break;
                case 2: // left (yellow)
                  grid[(i + 1) * MAX_SIZE + j] |= (5 << 5);
                  break;
                case 3: // right (red)
                  grid[i * MAX_SIZE + j + 1] |= (6 << 5);
                  break;
                }
              }
            }
          }
        }
        ++size;
        // out with the old tiles and in with the new tiles
        for (let i = 0; i != size; ++i) {
          for (let j = 0; j != size; ++j) {
            grid[i * MAX_SIZE + j] >>= 5;
          }
        }
        break;
      case SHOW_EMPTY:
        // insert new tiles, highlighted
        {
          skip = rendered + SHOW_INSERT_PAUSE - SHOW_EMPTY <= target;
          const highlight1 = skip ? 0 : 8;
          const highlight2 = skip ? 0 : 16;
          for (let i = 1; i != size; ++i) {
            for (let j = 1; j != size; ++j) {
              const b = i * MAX_SIZE + j;
              const r = (i - 1) * MAX_SIZE + j;
              const l = i * MAX_SIZE + j - 1;
              const t = (i - 1) * MAX_SIZE + j - 1;
              if ((!(grid[b] & 3) || (grid[b] & 7) === 1) &&
                  (!(grid[r] & 3) || (grid[r] & 7) === 5) &&
                  (!(grid[l] & 3) || (grid[l] & 7) === 6) &&
                  (!(grid[t] & 3) || (grid[t] & 7) === 2)) {
                ++markCount;
                if (Math.random() < 0.5) {
                  grid[t] |= 1 | highlight1;
                  grid[b] |= 2 | highlight2;
                } else {
                  grid[l] |= 5 | highlight1;
                  grid[r] |= 6 | highlight2;
                }
              }
            }
          }
          if (skip) {
            markCount = 0;
          }
        }
        break;
      case SHOW_INSERT_PAUSE:
        // unhighlight tiles
        markCount = 0;
        if (!skip) {
          for (let i = 0; i != size; ++i) {
            for (let j = 0; j != size; ++j) {
              grid[i * MAX_SIZE + j] &= 7;
            }
          }
        }
        break;
      case SHOW_EXPAND:
        markCount = 0;
        sizeChanged = true;
        break;
      case SHOW_ZAP:
        {
          skip = rendered + SHOW_DELETE_PAUSE - SHOW_ZAP <= target;
          // highlight tiles to be zapped
          for (let i = 1; i != size; ++i) {
            for (let j = 1; j != size; ++j) {
              const n = i * MAX_SIZE + j;
              if ((grid[n] & 3) === 3) {
                grid[n] = skip ? 0 : (grid[n] | 24);
                ++markCount;
              }
            }
          }
          if (!markCount) {
            // nothing to zap, skip ahead
            target += 4;
            startTime -= 4 * rate;
          }
          if (skip) {
            markCount = 0;
          }
          break;
        }
      case SHOW_DELETE_PAUSE:
        // delete tiles of zapped points
        markCount = 0;
        if (!skip) {
          for (let i = 1; i != size; ++i) {
            for (let j = 1; j != size; ++j) {
              if (grid[i * MAX_SIZE + j] & 24) {
                grid[i * MAX_SIZE + j] = 0;
              }
            }
          }
        }
        break;
      }
    }

    if (sizeChanged) {
      showSize(size + 1);
    }

    const targetStep = (rendered - 1) % SHOW_STEPS;
    const theta = target - (rendered - 1);
    const expanding =
      targetStep < SHOW_EXPAND ? 0 :
      targetStep > SHOW_EXPAND ? 1 :
      theta;

    // 2d scale
    const side = Math.min(width, height);
    const scale = 0.5 * side / (size + expanding);
    const view = new DOMMatrix([scale, 0, 0, scale, 0.5 * width, 0.5 * height]);
    // gl scale
    const xscale = 2 * scale / width;
    const yscale = -2 * scale / height;

    const arrowAlpha = Math.max(0, Math.min(1, (scale - 8) / 16));
    const lineWidth = Math.max(0, Math.min(3, (scale - 4) / 16));

    // compute colors
    fadeAlpha(arrowColors[0], colors[5], 0, 1, arrowAlpha);
    switch(targetStep) {
    case SHOW_EMPTY:
      fadeAlpha(tileColors[4], colors[4], 0, 1, theta);
      copyColor(tileColors[5], tileColors[4]);
      copyColor(tileColors[6], tileColors[4]);
      copyColor(tileColors[7], tileColors[4]);
      fadeAlpha(outerEdgeColor, colors[5], 0, 1, theta);
      copyColor(innerEdgeColor, colors[8]);
      copyColor(arrowColors[1], colors[8]);
      break;
    case SHOW_EMPTY_PAUSE:
      copyColor(tileColors[4], colors[4]);
      copyColor(tileColors[5], tileColors[4]);
      copyColor(tileColors[6], tileColors[4]);
      copyColor(tileColors[7], tileColors[4]);
      copyColor(outerEdgeColor, colors[5]);
      copyColor(innerEdgeColor, colors[8]);
      copyColor(arrowColors[1], colors[8]);
      break;
    case SHOW_INSERT:
      fade(tileColors[4], colors[4], colors[0], theta);
      fade(tileColors[5], colors[4], colors[1], theta);
      fade(tileColors[6], colors[4], colors[2], theta);
      fade(tileColors[7], colors[4], colors[3], theta);
      copyColor(outerEdgeColor, colors[5]);
      fadeAlpha(innerEdgeColor, colors[5], 0, 1, theta);
      fadeAlpha(arrowColors[1], colors[5], 0, arrowAlpha, theta);
      break;
    case SHOW_ZAP:
      fade(tileColors[4], colors[0], colors[9], theta);
      fade(tileColors[5], colors[1], colors[9], theta);
      fade(tileColors[6], colors[2], colors[9], theta);
      fade(tileColors[7], colors[3], colors[9], theta);
      copyColor(outerEdgeColor, colors[5]);
      copyColor(innerEdgeColor, colors[5]);
      copyColor(arrowColors[1], arrowColors[0]);
      break;
    case SHOW_ZAP_PAUSE:
      copyColor(tileColors[4], colors[9]);
      copyColor(tileColors[5], colors[9]);
      copyColor(tileColors[6], colors[9]);
      copyColor(tileColors[7], colors[9]);
      copyColor(outerEdgeColor, colors[5]);
      copyColor(innerEdgeColor, colors[5]);
      copyColor(arrowColors[1], arrowColors[0]);
      break;
    case SHOW_DELETE:
      fadeAlpha(tileColors[4], colors[9], 1, 0, theta);
      copyColor(tileColors[5], tileColors[4]);
      copyColor(tileColors[6], tileColors[4]);
      copyColor(tileColors[7], tileColors[4]);
      fadeAlpha(outerEdgeColor, colors[5], 1, 0, theta);
      copyColor(innerEdgeColor, outerEdgeColor);
      fadeAlpha(arrowColors[1], colors[5], arrowAlpha, 0, theta);
      break;
    default:
      copyColor(tileColors[4], colors[0]);
      copyColor(tileColors[5], colors[1]);
      copyColor(tileColors[6], colors[2]);
      copyColor(tileColors[7], colors[3]);
      copyColor(outerEdgeColor, colors[5]);
      copyColor(innerEdgeColor, colors[5]);
      copyColor(arrowColors[1], arrowColors[0]);
      break;
    }

    // inflate each tile by half the line width
    const e = lineWidth / width;
    const f = -lineWidth / height;

    // store vertex attributes
    bufferIndices[0] = 32 * markCount; // 2 tiles * 4 vertices * 4 attributes
    bufferIndices[1] = 0;
    for (let i = 0; i != size; ++i) {
      for (let j = 0; j != size; ++j) {
        const px = -i + j;
        const py = i + j - size + 1;
        const tiles = grid[i * MAX_SIZE + j];
        if ((tiles & 5) === 1) {
          const c = ((tiles >> 1) & 4) | 0;
          const d = (targetStep === SHOW_MOVE) ? theta - 1 : 0;
          vertex(xscale * (px - 1) - e, yscale * (py - d) - f, 0, 3, c);
          vertex(xscale * (px + 1) + e, yscale * (py - d) - f, 1, 3, c);
          vertex(xscale * (px + 1) + e, yscale * (py + 1 - d) + f, 1, 5, c);
          vertex(xscale * (px - 1) - e, yscale * (py + 1 - d) + f, 0, 5, c);
        }
        if ((tiles & 5) === 5) {
          const c = ((tiles >> 1) & 4) | 1;
          const d = (targetStep === SHOW_MOVE) ? theta - 1 : 0;
          vertex(xscale * (px - d) - e, yscale * (py + 1) + f, 0, 3, c);
          vertex(xscale * (px - d) - e, yscale * (py - 1) - f, 1, 3, c);
          vertex(xscale * (px + 1 - d) + e, yscale * (py - 1) - f, 1, 5, c);
          vertex(xscale * (px + 1 - d) + e, yscale * (py + 1) + f, 0, 5, c);
        }
        if ((tiles & 6) === 6) {
          const c = ((tiles >> 2) & 4) | 2;
          const d = (targetStep === SHOW_MOVE) ? theta - 1 : 0;
          vertex(xscale * (px + d) + e, yscale * (py - 1) - f, 0, 3, c);
          vertex(xscale * (px + d) + e, yscale * (py + 1) + f, 1, 3, c);
          vertex(xscale * (px - 1 + d) - e, yscale * (py + 1) + f, 1, 5, c);
          vertex(xscale * (px - 1 + d) - e, yscale * (py - 1) - f, 0, 5, c);
        }
        if ((tiles & 6) === 2) {
          const c = ((tiles >> 2) & 4) | 3;
          const d = (targetStep === SHOW_MOVE) ? theta - 1 : 0;
          vertex(xscale * (px + 1) + e, yscale * (py + d) + f, 0, 3, c);
          vertex(xscale * (px - 1) - e, yscale * (py + d) + f, 1, 3, c);
          vertex(xscale * (px - 1) - e, yscale * (py - 1 + d) - f, 1, 5, c);
          vertex(xscale * (px + 1) + e, yscale * (py - 1 + d) - f, 0, 5, c);
        }
      }
    }

    context.resetTransform();
    context.clearRect(0, 0, width, height);

    // draw grid to 2d canvas
    if (lineWidth) {
      context.lineWidth = lineWidth;
      context.setTransform(view);
      context.beginPath();
      for (let i = 1; i != size; ++i) {
        context.rect(-i, -(size - i), 2 * i, 2 * (size - i));
      }
      context.moveTo(0, -size + 1);
      context.lineTo(0, size - 1);
      context.moveTo(-size + 1, 0);
      context.lineTo(size - 1, 0);
      context.resetTransform();
      context.strokeStyle = gridStroke;
      context.stroke();

      // draw outer part of expanded grid
      if (targetStep >= SHOW_EXPAND) {
        context.beginPath();
        for (let n = 0; n != 4; ++n) {
          context.setTransform(view.multiply(transforms[n]));
          context.moveTo(size - 1, 0);
          context.lineTo(size, 0);
          for (let i = 0; i != size; ++i) {
            context.lineTo(size - i, i + 1);
            context.lineTo(size - i - 1, i + 1);
          }
        }
        context.resetTransform();
        fadeAlpha(gridColor, colors[6], 0, 1, expanding);
        context.strokeStyle = colorToCSS(gridColor);
        context.stroke();
      }
    }

    // draw tile outlines and arrows to texture canvas
    const h = Math.ceil(Math.max(16, Math.min(maxTextureSize / 16, scale)));
    const w = 2 * h;
    const lw = lineWidth * (h / scale);
    trim.canvas.width = 3 * w;
    trim.canvas.height = 16 * h;
    trim.clearRect(0, 0, 3 * w, 16 * h);

    const outerEdgeStroke = colorToCSS(outerEdgeColor);
    const innerEdgeStroke = colorToCSS(innerEdgeColor);
    const arrowStrokes = [
      colorToCSS(arrowColors[0]), colorToCSS(arrowColors[1])];
    const transform = DOMMatrix.fromMatrix(view);

    for (let n = 0; n != 8; ++n) {
      const type = n & 3;
      const mark = n >> 2;
      trim.fillStyle = colorToCSS(tileColors[n]);
      trim.fillRect(1.5 * mark * w, 4 * type * h, 1.5 * w, 4 * h)
      if (lineWidth) {
        const left = 2 * mark * w;
        const top = (4 * type + 1.5) * h;
        // stroke long edge at head of arrow
        trim.fillStyle = mark ? outerEdgeStroke : edgeStroke;
        trim.fillRect(left, top - 2 * lw, w, 3 * lw);
        // stroke long edge at tail of arrow
        trim.fillStyle = mark ? innerEdgeStroke : edgeStroke;
        trim.fillRect(left, top + h - lw, w, 3 * lw);
        if (arrowAlpha > 0.4) {
          // stroke arrow
          trim.beginPath();
          trim.moveTo(left + 0.45 * w, top + 0.45 * h);
          trim.lineTo(left + 0.5 * w, top + 0.35 * h);
          trim.lineTo(left + 0.55 * w, top + 0.45 * h);
          trim.moveTo(left + 0.5 * w, top + 0.35 * h);
          trim.lineTo(left + 0.5 * w, top + 0.65 * h);
          trim.strokeStyle = arrowStrokes[mark];
          trim.lineWidth = lw;
          trim.stroke();
        }
      }
    }
    if (lw) {
      // stroke short edges (extend into space between tiles to avoid gaps)
      trim.fillStyle = edgeStroke;
      trim.fillRect(0, 0, lw, 16 * h);
      trim.fillRect(w - lw, 0, 3 * lw, 16 * h);
      trim.fillStyle = outerEdgeStroke;
      trim.fillRect(2 * w - 2 * lw, 0, 3 * lw, 16 * h);
      trim.fillRect(3 * w - lw, 0, lw, 16 * h);
    }

    gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, gl.RGBA, gl.UNSIGNED_BYTE,
      trim.canvas);

    // paint tile colors to gl context
    const vertexAttribs = vertexAttribBuffer.subarray(0, bufferIndices[0]);
    gl.clear(gl.COLOR_BUFFER_BIT);
    gl.bufferData(gl.ARRAY_BUFFER, vertexAttribs, gl.DYNAMIC_DRAW);
    // 4 floats per vertex, 4 vertices per quad, 1 quad per 2 triangles
    gl.drawElements(gl.TRIANGLES, bufferIndices[0] * 3 / 8, gl.UNSIGNED_INT, 0);

    // copy to canvas2d context
    context.drawImage(gl.canvas, 0, 0);
  }

  function addShader(gl, program, type, source) {
    const shader = gl.createShader(type);
    gl.shaderSource(shader, source);
    gl.compileShader(shader);
    if (!gl.getShaderParameter(shader, gl.COMPILE_STATUS)) {
      throw new Error((type == gl.VERTEX_SHADER ? "Vertex" : "Fragment") +
        " shader compilation failed: " + gl.getShaderInfoLog(shader));
    }
    gl.attachShader(program, shader);
    return true;
  }

  function makeProgram(gl, vshader, fshader) {
    const program = gl.createProgram();
    addShader(gl, program, gl.VERTEX_SHADER, vshader);
    addShader(gl, program, gl.FRAGMENT_SHADER, fshader);
    gl.linkProgram(program);
    if (!gl.getProgramParameter(program, gl.LINK_STATUS)) {
      throw new Error("Shader link failed: " + gl.getProgramInfoLog(program));
    }
    return program;
  }

  function initializeGraphics() {
    gl.getExtension("OES_element_index_uint");
    gl.clearColor(0.0, 0.0, 0.0, 0.0);
    maxTextureSize = gl.getParameter(gl.MAX_TEXTURE_SIZE);
    if (!maxTextureSize) {
      maxTextureSize = 2048;
    }

    const program = makeProgram(gl, shaders[0], shaders[1]);
    const x = gl.getAttribLocation(program, "x");
    const u = gl.getAttribLocation(program, "u");
    const sampler = gl.getUniformLocation(program, "sampler");
    gl.useProgram(program);

    gl.activeTexture(gl.TEXTURE0);
    gl.uniform1i(sampler, 0);
    const texture = gl.createTexture();
    gl.bindTexture(gl.TEXTURE_2D, texture);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.LINEAR);

    const size = vertexAttribBuffer.byteLength;
    gl.bindBuffer(gl.ARRAY_BUFFER, gl.createBuffer());
    gl.bufferData(gl.ARRAY_BUFFER, size, gl.DYNAMIC_DRAW);
    gl.vertexAttribPointer(x, 2, gl.FLOAT, false, 4 * FLOATS, 0 * FLOATS);
    gl.vertexAttribPointer(u, 2, gl.FLOAT, false, 4 * FLOATS, 2 * FLOATS);
    gl.enableVertexAttribArray(x);
    gl.enableVertexAttribArray(u);

    gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, gl.createBuffer());
    gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, indexBuffer, gl.STATIC_DRAW);
  }

  function update(time) {
    updatePending = false;
    render(time);
    if (paused) {
      pauseTime = time;
    } else {
      requestAnimationFrame(update);
      updatePending = true;
    }
  }

  function pause(state) {
    paused = state;
    if (paused) {
      document.getElementById("pause").classList.add("down");
    } else {
      document.getElementById("pause").classList.remove("down");
    }
    if (!paused && !updatePending) {
      requestAnimationFrame(update);
      updatePending = true;
    }
  }

  function start() {
    startTime = undefined;
    pause(false);
  }

  function resize() {
    resizePending = true;
    if (!updatePending) {
      requestAnimationFrame(update);
      updatePending = true;
    }
  }

  function restartClick() {
    start();
  }

  function pauseClick() {
    pause(!paused);
  }

  function rateInput() {
    const SLOWEST = 5000;
    const FASTEST = 1;
    const x = parseFloat(document.getElementById("rate").value);
    const y = SLOWEST * Math.pow(FASTEST / SLOWEST, x);
    // inverse: x = log(y / SLOWEST) / log(FASTEST / SLOWEST)
    const speedString = (y / 1000.0).toPrecision(3);
    document.getElementById("speed").firstChild.nodeValue = speedString;
    newRate = y;
  }

  function shaderLoaded() {
    if (shaders[0] && shaders[1]) {
      initializeGraphics();
      start();
    }
  }

  function loadShaders() {
    const urls = ["v.glsl", "f.glsl"];
    for (let n = 0; n != 2; ++n) {
      fetch(urls[n]).then(response => response.text()).then(data => {
        shaders[n] = data;
        shaderLoaded();
      });
    }
  }

  document.addEventListener("DOMContentLoaded", function () {
    addEventListener("resize", resize);
    const options = {premultipliedAlpha: false};
    context = document.getElementById("canvas").getContext("2d", options);
    trim = document.createElement("canvas").getContext("2d", options);
    gl = document.createElement("canvas").getContext("webgl", options);
    context.canvas.addEventListener("click", pauseClick);
    document.getElementById("restart").addEventListener("click", restartClick);
    document.getElementById("pause").addEventListener("click", pauseClick);
    document.getElementById("rate").addEventListener("input", rateInput);
    rateInput();
    loadShaders();
  });
}());
