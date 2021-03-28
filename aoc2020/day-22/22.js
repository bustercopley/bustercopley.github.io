// -*- coding: utf-8-unix; -*-

(function () {
  "use strict";
  const rate = 125;
  const levelRateExponent = 0.250;
  const roundRateExponent = 0.05;

  let startTime;
  let finished = true;
  let levelText = [];
  let gameTime = 0.0;
  let nextGame = 1;
  let hands, history, round, game, stack;

  function compareHands(hands, ohands) {
    for (let i = 0; i != 2; ++i) {
      if (hands[i].length != ohands[i].length) {
        return false;
      }
      for (let j = 0; j != hands[i].length; ++j) {
        if (hands[i][j] != ohands[i][j]) {
          return false;
        }
      }
    }
    return true;
  }

  function findDuplicate(hands, history) {
    for (let i = 0; i != history.length; ++i) {
      if (compareHands(hands, history[i])) {
        return true;
      }
    }
    return false;
  }

  function playTurn(verbose) {
    let text;
    let lines = 0;
    let level = stack.length;
    if (!hands[0].length || !hands[1].length) {
      throw new Error("logic error: empty hand(s)");
    }
    verbose = true;
    if (verbose) {
      text = "" +
        "=== Game " + game
        + " (level " + stack.length + ") ===\n\n" +
        "-- Round " + round + " --\n";
      for (let i = 0; i != 2; ++i) {
        text = text + "Player " + (i + 1) + "'s deck";
        let sep = ":";
        for (let card of hands[i]) {
          text = text + sep + " " + card
          sep = ",";
        }
        text = text + "\n";
      }
      lines += 5;
    }
    let roundWinner, gameWinner;
    if (findDuplicate(hands, history)) {
      text = text +
        "Duplicate hand detected!\n\n";
      lines += 2;
      gameWinner = 0;
    } else {
      history.push([hands[0].slice(), hands[1].slice()]);
      if (verbose) {
        for (let i = 0; i != 2; ++i) {
          text = text + "Player " + (i + 1) + " plays: " + hands[i][0] + "\n";
        }
        lines += 2;
      }
      if (hands[0].length - 1 < hands[0][0] ||
          hands[1].length - 1 < hands[1][0]) {
        roundWinner = hands[1][0] > hands[0][0] ? 1 : 0;
      } else {
        if (verbose) {
          text = text + "Playing a sub-game to determine the winner...\n";
          ++lines;
        }
        stack.push([hands, history, round, game]);
        let ohands = hands;
        hands = [
          ohands[0].slice(1, ohands[0][0] + 1),
          ohands[1].slice(1, ohands[1][0] + 1)
        ];
        history = [];
        round = 1;
        game = ++nextGame;
      }
    }

    do {
      if (roundWinner !== undefined) {
        gameWinner = undefined;
        hands[roundWinner].push(hands[roundWinner].shift());
        hands[roundWinner].push(hands[roundWinner ^ 1].shift());
        ++round;
        if (!hands[roundWinner ^ 1].length) {
          gameWinner = roundWinner;
        }
        else {
          if (verbose) {
            text = text +
              "Player " + (roundWinner + 1) + " wins round " + round +
              " of game " + game + "!\n";
            ++lines;
          }
        }
        roundWinner = undefined;
      }
      if (gameWinner !== undefined) {
        if (verbose && !stack.length) {
          text = text +
            "The winner of game " + game +
            " is player " + (gameWinner + 1) + "!\n";
          ++lines;
        }
        if (stack.length) {
          roundWinner = gameWinner;
          [hands, history, round, game] = stack.pop();
        } else {
          finished = true;
        }
      }
    } while (roundWinner !== undefined);
    if (verbose) {
      if (lines < 12) {
        text = text + "\n\n\n\n\n\n\n\n\n\n".substr(lines - 10);
      }
      levelText[level] = text;
    }
  }

  function reset() {
    nextGame = 0;
    hands = [[9, 2, 6, 3, 1], [5, 8, 4, 7, 10]];
    let anchor = location.href.match(/#(\d{1,3}(?:,\d{1,3})*);(\d{1,3}(?:,\d{1,3})*)$/);
    if (anchor) {
      hands = [
        anchor[1].split(",").map(parseFloat),
        anchor[2].split(",").map(parseFloat),
      ];
    }
    const ids = ["hand1", "hand2"];
    for (let i = 0; i != 2; ++i) {
      document.getElementById(ids[i]).value = hands[i].join(" ");
    }

    history = [];
    round = 1;
    game = ++nextGame;
    stack = [];
    startTime = undefined;
    gameTime = 0.0;
    levelText = [];
    if (finished) {
      finished = false;
      requestAnimationFrame(update);
    }
  }

  function update(animationTime) {
    let modified;
    if (startTime === undefined ) {
      reset();
      modified = true;
      startTime = animationTime;
    }
    let targetGameTime = animationTime - startTime;
    // don't try to catch up more than one second
    gameTime = Math.max(gameTime, targetGameTime - 1000);
    while (gameTime < targetGameTime) {
      let turnTime = rate *
          Math.pow(levelRateExponent, stack.length) *
          Math.pow(1 - stack.length * roundRateExponent, round);
      // don't try to advance more than 1000 turns per frame
      turnTime = Math.max(turnTime, 0.001);
      gameTime += turnTime;
      playTurn(gameTime >= animationTime - startTime);
      modified = true;
    }
    if (modified) {
      let text = "";
      for (let t of levelText) {
        text = text + t;
      }
      let statusText = document.getElementById("status").firstChild;
      statusText.nodeValue = text;
      modified = false;
    }
    if (!finished) {
      requestAnimationFrame(update);
    }
  }

  function click() {
    const ids = ["hand1", "hand2"];
    hands = new Array(2);
    for (let i = 0; i != 2; ++i) {
      hands[i] = new Array();
      for (let item of document.getElementById(ids[i]).value.split(/[ ,]+/)) {
        if (item.match(/^\d{0,3}$/)) {
          hands[i].push(item);
        }
      }
    }
    location.href = "#" + hands[0].join(",") + ";" + hands[1].join(",");
    reset();
  }

  document.addEventListener("DOMContentLoaded", () => {
    document.getElementById("reset").addEventListener("click", click);
    reset();
  });
}());
