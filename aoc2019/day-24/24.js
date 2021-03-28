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

  var GENERATION_TIME = 500; // milliseconds
  var GENERATION_COUNT = 200;
  var CELLSPEC = [
    [[3, 3], [3, 2], [3, 1], [2, 1], [1, 1], [1, 2], [1, 3], [2, 3]],
    [[3, 4], [4, 4], [4, 3], [4, 2], [4, 1], [4, 0], [3, 0], [2, 0],
     [1, 0], [0, 0], [0, 1], [0, 2], [0, 3], [0, 4], [1, 4], [2, 4]]];


  document.addEventListener("DOMContentLoaded", function () {
    var canvas = document.getElementById("diagram");
    var statusText = document.getElementById("status").firstChild;
    var startTime, generation, data;

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
      var level = newLevel();
      var input = document.getElementById("input").value.split('\n');
      for (var i = 0; i !== 5; ++i) {
        for (var j = 0; j !== 5; ++j) {
          level[i][j] = (input[i][j] == '#' ? 1 : 0);
        }
      }
      return level;
    }

    function neighbours(level, i, j, minLevel, maxLevel) {
      var n = 0, k;
      // above
      if (i === 0) {
        if (level > minLevel) {
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
        if (level > minLevel) {
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
        if (level > minLevel) {
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
        if (level > minLevel) {
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
      if (startTime === undefined) {
        generation = 0;
        data = { 0: getInput() };
        startTime = updateTime;
      }
      var time = (updateTime - startTime) / GENERATION_TIME;

      var generationNeeded = Math.min(time|0, GENERATION_COUNT);
      var bugCount;
      time = Math.min(time, generation + 1);
      var level, i, j, n, newData, newValue;
      var minLevel = -(((generation + 1) / 2) | 0);
      var maxLevel = (((generation + 3) / 2) | 0);

      while (generation < generationNeeded) {
        ++generation;
        bugCount = 0;
        minLevel = -(((generation + 1) / 2) | 0);
        maxLevel = (((generation + 1) / 2) | 0);

        if (generation % 2 === 1) {
          data[minLevel] = newLevel();
          data[maxLevel] = newLevel();
        }

        newData = {};
        for (level = minLevel; level != maxLevel + 1; ++level) {
          newData[level] = newLevel();
          for (i = 0; i !== 5; ++i) {
            for (j = 0; j !== 5; ++j) {
              if (i !== 2 || j !== 2) {
                n = neighbours(level, i, j, minLevel, maxLevel);
                if (data[level][i][j] === 1) {
                  newValue = (n === 1 ? 1 : 0);
                }
                else {
                  newValue = (n === 1 || n === 2 ? 1 : 0);
                }
                bugCount = bugCount + newValue;
                newData[level][i][j] = newValue;
              }
            }
          }
        }
        data = newData;
        var status = "Generation " + generation.toString() + "\n";
        status = status + "Bug count " + bugCount.toString() + "\n";
        // for (level = minLevel; level != maxLevel + 1; ++level) {
        //   status = status + "Level " + level.toString() + "\n";
        //   for (i = 0; i !== 5; ++i) {
        //     for (j = 0; j !== 5; ++j) {
        //       status = status + (data[level][i][j] === 1 ? "#" : " ");
        //     }
        //     status = status + "\n";
        //   }
        // }
        statusText.nodeValue = status;
      }

      var middleRadius = 0.5;
      var deltaRadius = middleRadius / (time + 3);
      var side = Math.min(canvas.width, canvas.height);

      function paintRing(k) {
        context.fillStyle = (k == 0 || k == 1 ? "red" : "green");
        context.setTransform(0.5 * side, 0, 0, 0.5 * side, 0.5 * canvas.width, 0.5 * canvas.height);

        var innerRadius = middleRadius + (0.5 + k) * deltaRadius;
        var outerRadius = middleRadius + (0.5 + k + 1) * deltaRadius;

        var level = Math.floor(k / 2) | 0;
        var i, j, n, value, cellsInRing;
        cellsInRing = (k % 2 === 0 ? 8 : 16);
        for (n = 0; n !== cellsInRing; ++n) {
          value = data[-level][CELLSPEC[k % 2 === 0 ? 0 : 1][n][0]][CELLSPEC[k % 2 === 0 ? 0 : 1][n][1]];
          if (value === 1) {
            var theta0 = (2 * n + 1) * Math.PI / cellsInRing;
            var theta1 = (2 * n + 3) * Math.PI / cellsInRing;
            context.beginPath();
            context.moveTo(outerRadius * Math.cos(theta0), outerRadius * Math.sin(theta0));
            context.arc(0, 0, outerRadius, theta0, theta1);
            context.lineTo(innerRadius * Math.cos(theta1), innerRadius * Math.sin(theta1));
            context.arc(0, 0, innerRadius, theta1, theta0, true);
            context.closePath();
            context.fill();
          }
        }
      }

      var context = canvas.getContext("2d");
      context.setTransform(1, 0, 0, 1, 0, 0);
      context.clearRect(0, 0, canvas.width, canvas.height);

      paintRing(0);
      var k;
      for (k = -generation; k !== generation + 2; ++k) {
        paintRing(k);
      }

      if (generation < GENERATION_COUNT) {
        requestAnimationFrame(paint);
      }
    }

    function startClick() {
      startTime = undefined;
      statusText.nodeValue = "Generation 0";
      requestAnimationFrame(paint);
    }

    function windowResize() {
      var rect = canvas.getBoundingClientRect();
      var dpr = window.devicePixelRatio;
      var width = Math.round(dpr * rect.right) - Math.round(dpr * rect.left);
      var height = Math.round(dpr * rect.bottom) - Math.round(dpr * rect.top);
      canvas.width = width;
      canvas.height = height;
    }

    document.getElementById("input").value = "....#\n#..#.\n#..##\n..#..\n#....\n";
    document.getElementById("start").addEventListener("click", startClick);
    window.addEventListener("resize", windowResize);
    window.dispatchEvent(new Event("resize"));
  }, true);
}());
