// -*- coding: utf-8; -*-

// Copyright 2019, Richard Copley <buster at buster dot me dot uk>

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

  var GENERATION_COUNT = 200;

  document.addEventListener("DOMContentLoaded", function () {
    var canvas = document.getElementById("diagram");
    var statusText = document.getElementById("status").firstChild;
    var generationTime = 1000;
    var startTime, generation, data, status;
    var pauseClicked, pauseTime;
    var paintPending;

    function newLevel() {
      var level = new Array(5);
      for (var i = 0; i !== 5; ++i) {
        level[i] = new Array(5);
        for (var j = 0; j !== 5; ++j) {
          level[i][j] = 0;
        }
      }
      return level;
    }

    function getInput() {
      data = { 0: newLevel() };
      var input = document.getElementById("input").value.split('\n');
      var bugCount = 0;
      for (var i = 0; i !== 5; ++i) {
        for (var j = 0; j !== 5; ++j) {
          var value = (input[i][j] == '#' ? 1 : 0);
          data[0][i][j] = value;
          bugCount += value;
        }
      }
      return bugCount;
    }

    function addLevelStatus(level) {
      status = status + "\nLevel " + level.toString() + "\n";
      for (var i = 0; i !== 5; ++i) {
        for (var j = 0; j !== 5; ++j) {
          status = status + (data[level][i][j] === 1 ? "#" : ".");
        }
        status = status + "\n";
      }
    }

    function neighbours(level, i, j, maxLevel) {
      var n = 0, k;
      // above
      if (i === 0) {
        if (level > -maxLevel) {
          n += data[level - 1][1][2];
        }
      }
      else if (i === 3 && j === 2) {
        if (j === 2) {
          if (level < maxLevel) {
            for (k = 0; k !== 5; ++k) {
              n += data[level + 1][4][k];
            }
          }
        }
      }
      else {
        n += data[level][i - 1][j];
      }
      // below
      if (i === 4) {
        if (level > -maxLevel) {
          n += data[level - 1][3][2];
        }
      }
      else if (i === 1 && j === 2) {
        if (level < maxLevel) {
          for (k = 0; k !== 5; ++k) {
            n += data[level + 1][0][k];
          }
        }
      }
      else {
        n += data[level][i + 1][j];
      }
      // left
      if (j === 0) {
        if (level > -maxLevel) {
          n += data[level - 1][2][1];
        }
      }
      else if (j === 3 && i === 2) {
        if (level < maxLevel) {
          for (k = 0; k !== 5; ++k) {
            n += data[level + 1][k][4];
          }
        }
      }
      else {
        n += data[level][i][j - 1];
      }
      // right
      if (j === 4) {
        if (level > -maxLevel) {
          n += data[level - 1][2][3];
        }
      }
      else if (j === 1 && i === 2) {
        if (level < maxLevel) {
          for (k = 0; k !== 5; ++k) {
            n += data[level + 1][k][0];
          }
        }
      }
      else {
        n += data[level][i][j + 1];
      }

      return n;
    }

    function paint(updateTime) {
      paintPending = false;

      if (startTime === undefined) {
        startTime = updateTime;
        pauseClicked = undefined;
        pauseTime = undefined;
        document.getElementById("pause").removeAttribute("disabled");

        generation = 0;
        status = "Generation 0\nBug count " + getInput() + "\n";
        addLevelStatus(0);
        statusText.nodeValue = status;
      }

      if (pauseClicked) {
        pauseClicked = undefined;
        if (pauseTime) {
          startTime = startTime + updateTime - pauseTime;
          pauseTime = undefined;
        }
        else {
          pauseTime = updateTime;
        }
      }

      if (pauseTime) {
        updateTime = pauseTime;
      }

      var time = (updateTime - startTime) / generationTime;
      var generationNeeded = Math.min(time|0, GENERATION_COUNT);
      time = Math.min(time, generation + 1);
      var level, i, j, n, newData, value, bugCount;
      var maxLevel = (((generation + 1) / 2) | 0);

      while (generation < generationNeeded) {
        ++generation;
        bugCount = 0;
        maxLevel = (((generation + 1) / 2) | 0);
        if (generation % 2 === 1) {
          data[-maxLevel] = newLevel();
          data[maxLevel] = newLevel();
        }
        newData = {};
        for (level = -maxLevel; level != maxLevel + 1; ++level) {
          newData[level] = newLevel();
          for (i = 0; i !== 5; ++i) {
            for (j = 0; j !== 5; ++j) {
              if (i !== 2 || j !== 2) {
                n = neighbours(level, i, j, maxLevel);
                if (data[level][i][j] === 1) {
                  value = (n === 1 ? 1 : 0);
                }
                else {
                  value = (n === 1 || n === 2 ? 1 : 0);
                }
                bugCount = bugCount + value;
                newData[level][i][j] = value;
              }
            }
          }
        }
        data = newData;
        status = "Generation " + generation.toString() + "\n";
        status = status + "Bug count " + bugCount.toString() + "\n";
        if (maxLevel >= 1) {
          addLevelStatus(-maxLevel);
        }
        if (maxLevel >= 2) {
          addLevelStatus(-maxLevel + 1);
        }
        if (maxLevel >= 3) {
          addLevelStatus(-maxLevel + 2);
        }
        addLevelStatus(0);
        if (maxLevel >= 3) {
          addLevelStatus(maxLevel - 2);
        }
        if (maxLevel >= 2) {
          addLevelStatus(maxLevel - 1);
        }
        if (maxLevel >= 1) {
          addLevelStatus(maxLevel);
        }
        statusText.nodeValue = status;
      }

      var middleRadius = 0.5;
      var deltaRadius = middleRadius / (time + 3);
      var radii = new Array(3);
      var PI4 = Math.PI / 4;
      var w = canvas.width, h = canvas.height, s = Math.min(w, h);
      var context = canvas.getContext("2d");
      var vertices, vertex, u, v, outerRing, clockwise, k;
      // http://colorbrewer2.org/?type=qualitative&scheme=Paired&n=10
      var colors = ["#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99",
                    "#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a"];
      var colorIndex;
      var corners = new Array(4);

      context.setTransform(1, 0, 0, 1, 0, 0);
      context.clearRect(0, 0, w, h);
      context.setTransform(0.5 * s, 0, 0, 0.5 * s, 0.5 * w, 0.5 * h);

      for (level = -maxLevel; level != maxLevel + 1; ++ level) {
        colorIndex = 2 * (maxLevel - Math.abs(level)) % colors.length;
        for (k = 0; k != 3; ++ k) {
          radii[k] = middleRadius + (0.5 - 2 * level + k) * deltaRadius;
        }

        // Cosine rule a bunch of times
        var phi = 0.95 * Math.PI / 2;
        var cphi = Math.cos(phi), sphi = Math.sin(phi);
        var rr0 = radii[0] / radii[1];
        var rr1 = radii[1] / radii[2];
        var t1 = Math.acos(sphi * sphi * rr1 +
          cphi * Math.sqrt(1 - sphi * sphi * rr1 * rr1)) / PI4;
        var t3 = Math.acos(sphi * sphi * rr0 +
          cphi * Math.sqrt(1 - sphi * sphi * rr0 *rr0)) / PI4;
        var phi1 = phi - t3;
        var cphi1 = Math.cos(phi1), sphi1 = Math.sin(phi1);
        var t2 = t3 + Math.acos(sphi1 * sphi1 * rr1 +
          cphi1 * Math.sqrt(1 - sphi1 * sphi1 * rr1 * rr1)) / PI4;
        // Don't cross the streamers
        t2 = Math.min(t2, 0.95);

        var vertices = [ // [radius-index, angle / (PI / 4)]
          [[2, 5], [2, 5 + t1], [2, 5 + t2], [2, 7 - t2], [2, 7 - t1], [2, 7]],
          [[2, 5 - t1], [1, 5], [1, 5 + t3], [1, 7 - t3], [1, 7], [2, 7 + t1]],
          [[2, 5 - t2], [1, 5 - t3], [0, 5], [0, 7], [1, 7 + t3], [2, 7 + t2]],
          [[2, 3 + t2], [1, 3 + t3], [0, 3], [0, 1], [1, 1 - t3], [2, 1 - t2]],
          [[2, 3 + t1], [1, 3], [1, 3 - t3], [1, 1 + t3], [1, 1], [2, 1 - t1]],
          [[2, 3], [2, 3 - t1], [2, 3 - t2], [2, 1 + t2], [2, 1 + t1], [2, 1]]];

        for (i = 0; i !== 5; ++i) {
          for (j = 0; j !== 5; ++j) {
            if (i !== 2 || j !== 2) {
              if (data[level][i][j]) {
                outerRing = 0;
                corners[0] = vertices[i][j];
                corners[1] = vertices[i][j + 1];
                corners[2] = vertices[i + 1][j + 1];
                corners[3] = vertices[i + 1][j];
                u = corners[3];
                context.beginPath();
                context.moveTo(radii[u[0]] * Math.cos(PI4 * u[1]),
                               radii[u[0]] * Math.sin(PI4 * u[1]))
                for (k = 0; k != 4; ++k) {
                  v = corners[k];
                  if (v[0] == 2) {
                    outerRing = 1;
                  }
                  if (v[0] != u[0]) {
                    context.lineTo(radii[v[0]] * Math.cos(PI4 * v[1]),
                                   radii[v[0]] * Math.sin(PI4 * v[1]));
                  }
                  else {
                    clockwise = (8 + v[1] - u[1]) % 8 > 4;
                    context.arc(0, 0, radii[u[0]],
                                PI4 * u[1], PI4 * v[1], clockwise);
                  }
                  u = v;
                }
                context.closePath();
                context.fillStyle = colors[colorIndex + outerRing];
                context.fill();
              }
            }
          }
        }
      }
      if (generation >= GENERATION_COUNT) {
        document.getElementById("pause").setAttribute("disabled", "");
      }
      else if (!paintPending && !pauseTime) {
        paintPending = true;
        requestAnimationFrame(paint);
      }
    }

    function startClick(e) {
      startTime = undefined;
      generationTime = +e.target.getAttribute("data-generation-time");
      if (!paintPending) {
        paintPending = true;
        requestAnimationFrame(paint);
      }
    }

    function pauseClick() {
      pauseClicked = true;
      if (pauseTime && !paintPending) {
        paintPending = true;
        requestAnimationFrame(paint);
      }
    }

    function windowResize() {
      var rect = canvas.getBoundingClientRect();
      var dpr = window.devicePixelRatio;
      var width = Math.round(dpr * rect.right) - Math.round(dpr * rect.left);
      var height = Math.round(dpr * rect.bottom) - Math.round(dpr * rect.top);
      canvas.width = width;
      canvas.height = height;
      if (!paintPending) {
        paintPending = true;
        requestAnimationFrame(paint);
      }
    }

    document.getElementById("start0").addEventListener("click", startClick);
    document.getElementById("start1").addEventListener("click", startClick);
    document.getElementById("start2").addEventListener("click", startClick);
    document.getElementById("pause").addEventListener("click", pauseClick);
    window.addEventListener("resize", windowResize);
    window.dispatchEvent(new Event("resize"));
  }, true);
}());
