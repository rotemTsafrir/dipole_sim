// Zoom update: wheel over the field; 0 resets magnification. Range 50%-400%.
/*
 * EM simulator — CPU optimization, stage 2 (2026-10-02)
 * Replace your complete p5.js sketch.js with this file; keep your existing p5 setup.
 * Changes: primitive segment integration; real/imaginary phasors; double-precision
 * visible arrays; shared frame sine/cosine; one-pixel-per-cell rendering; avoid
 * redundant segment rebuilds; render the last row/column; restore 60 FPS target.
 * Stage 2: per-antenna unit-drive vector-potential grids, overlap reuse on pan,
 * numeric derivative stencils and immediate recombination on slider release.
 * Double-precision caches are bounded to the last viewport per antenna.
 * Original physical model, derivatives, softening and display transfer functions
 * are retained. The energy-flux heatmap is the original proxy, not calibrated |S|.
 */

// Converted Java abstract class to JavaScript class structure using p5.js-compatible syntax

class Antenna {
  constructor(wavelength, amp, phase, dl = 0.12) {
    this.wavelength = wavelength;
    this.amp = amp;
    this.phase = phase;
    this.I0 = ComplexNum.cis(phase).product(amp);
    this.dl = dl;
    this.segments = [];
    this.currentSegments = [];
    this.segFlags = [];

    this.showBox = false;
    this.xBox = 0;
    this.yBox = 0;
  }

  getWavelength() {
    return this.wavelength;
  }

  getAmp() {
    return this.amp;
  }

  getPhase() {
    return this.phase;
  }

  getI0() {
    return this.I0;
  }

  getWirePoints() {
    return this.wirePoints;
  }

  getSegments() {
    return this.segments;
  }

  getCurrentSegments() {
    return this.currentSegments;
  }

  getSegFlags() {
    return this.segFlags;
  }

  getDl() {
    return this.dl;
  }

  setWavelength(new_wavelength) {
    this.wavelength = new_wavelength;
  }

  setAmp(amp) {
    this.amp = amp;
    this.I0 = ComplexNum.cis(this.phase).product(amp);
  }

  setPhase(phase) {
    this.phase = phase;
    this.I0 = ComplexNum.cis(phase).product(this.amp);
  }

  setDl(dl) {
    this.dl = dl;
  }

  setShowBox(b) {
    this.showBox = b;
  }

  getShowBox() {
    return this.showBox;
  }

  toggleShowBox() {
    this.showBox = !this.showBox;
  }
}

////

class Dipole extends Antenna {
  constructor(
    wavelength,
    endPointA,
    endPointB,
    amp,
    phase,
    thickOrDl,
    thick = null
  ) {
    super(wavelength, amp, phase, thick !== null ? thickOrDl : undefined);

    this.length = 0;

    this.wirePoints = [];
    // antenna endpoints
    this.endPointA = [endPointA[0], endPointA[1]];
    this.endPointB = [endPointB[0], endPointB[1]];

    // bounding rectangle points
    this.P1 = [0, 0];
    this.P2 = [0, 0];
    this.P3 = [0, 0];
    this.P4 = [0, 0];

    // feeding separation of dipole
    this.sep = 0.2;

    // thickness (in my scale) of dipole
    this.thickness = thick !== null ? thick : thickOrDl;

    // horizontal and vertical distance (in my scale) of status box from antenna
    this.delBoxAntenna = 0.1;

    // Add wire points

    this.wirePoints.push([...endPointA]);
    this.wirePoints.push([...endPointB]);

    this.length = this.dist2D(endPointA, endPointB);

    this.setCurrentSegments();
    this.calculateBoundingBox();
  }

  calculateBoundingBox() {
    let tempDel = [
      this.endPointB[1] - this.endPointA[1],
      this.endPointA[0] - this.endPointB[0],
    ];
    let len = Math.sqrt(tempDel[0] * tempDel[0] + tempDel[1] * tempDel[1]);

    tempDel[0] *= this.thickness / len;
    tempDel[1] *= this.thickness / len;

    this.P1[0] = this.endPointA[0] + 0.5 * tempDel[0];
    this.P1[1] = this.endPointA[1] + 0.5 * tempDel[1];

    this.P2[0] = this.endPointB[0] + 0.5 * tempDel[0];
    this.P2[1] = this.endPointB[1] + 0.5 * tempDel[1];

    this.P3[0] = this.endPointB[0] - 0.5 * tempDel[0];
    this.P3[1] = this.endPointB[1] - 0.5 * tempDel[1];

    this.P4[0] = this.endPointA[0] - 0.5 * tempDel[0];
    this.P4[1] = this.endPointA[1] - 0.5 * tempDel[1];

    let maxX = Math.max(this.P1[0], this.P2[0]);
    let Y = this.P4[1];

    maxX = Math.max(maxX, this.P3[0]);
    maxX = Math.max(maxX, this.P4[0]);

    if (maxX === this.P1[0]) {
      Y = this.P1[1];
    } else if (maxX === this.P2[0]) {
      Y = this.P2[1];
    } else if (maxX === this.P3[0]) {
      Y = this.P3[1];
    }

    this.xBox = maxX + this.delBoxAntenna;
    this.yBox = Y - this.delBoxAntenna;
  }

  getSep() {
    return this.sep;
  }

  dist2D(p1, p2) {
    let temp =
      (p2[0] - p1[0]) * (p2[0] - p1[0]) + (p2[1] - p1[1]) * (p2[1] - p1[1]);
    return Math.sqrt(temp);
  }

  setCurrentSegments() {
    this.segments = [];
    this.currentSegments = [];
    this.segFlags = [];

    let l = 0;
    let sideLen;
    let pA;
    let pB;
    let delP = [0, 0];
    let p = [0, 0];
    let k = (2 * Math.PI) / this.wavelength;

    for (let i = 0; i < this.wirePoints.length - 1; i++) {
      pA = this.wirePoints[i];
      pB = this.wirePoints[i + 1];

      sideLen = this.dist2D(pA, pB);

      delP[0] = pB[0] - pA[0];
      delP[1] = pB[1] - pA[1];

      l = 0;

      while (l <= sideLen) {
        p[0] = pA[0] + (l / sideLen) * delP[0];
        p[1] = pA[1] + (l / sideLen) * delP[1];

        this.segments.push([p[0], p[1]]);

        let Ivec = [0, 0];
        let z = l + this.dl / 2 - sideLen / 2;

        if (Math.abs(z) > this.sep / 2) {
          this.segFlags.push(true);

          Ivec[0] = delP[0] / sideLen;
          Ivec[1] = delP[1] / sideLen;
          Ivec[0] *= Math.sin(k * (sideLen / 2 - Math.abs(z)));
          Ivec[1] *= Math.sin(k * (sideLen / 2 - Math.abs(z)));
        } else {
          this.segFlags.push(false);
        }

        this.currentSegments.push(Ivec);

        l += this.dl;
        l = Math.round(10000 * l) / 10000.0;
      }
    }
  }

  setSep(newSep) {
    this.sep = newSep;
    this.setCurrentSegments();
  }

  setDl(new_dl) {
    if (this.dl === new_dl) return;
    this.dl = new_dl;
    this.setCurrentSegments();
  }

  setXBox(newX) {
    this.xBox = newX;
  }

  setYBox(newY) {
    this.yBox = newY;
  }

  setDelBoxAntenna(del) {
    this.delBoxAntenna = del;
  }

  setAmp(newAmp) {
    this.amp = newAmp;
    this.I0 = ComplexNum.cis(this.phase).product(newAmp);
  }

  setPhase(newPhase) {
    this.phase = newPhase;
    this.I0 = ComplexNum.cis(newPhase).product(this.amp);
  }

  setWavelength(newWavelength) {
    this.wavelength = newWavelength;
    this.setCurrentSegments();
  }

  getXBox() {
    return this.xBox;
  }

  getYBox() {
    return this.yBox;
  }

  getLength() {
    return this.length;
  }

  getType() {
    return "Dipole";
  }

  getP1() {
    return this.P1;
  }

  getP2() {
    return this.P2;
  }

  getP3() {
    return this.P3;
  }

  getP4() {
    return this.P4;
  }
}

//

// Converted ComplexNum class from Java to JavaScript

class ComplexNum {
  constructor(a, b) {
    this.a = a;
    this.b = b;
    this.absVal = Math.sqrt(a * a + b * b);

    this.phase = Math.atan2(b, a);
  }

  getA() {
    return this.a;
  }

  getB() {
    return this.b;
  }

  getAbsVal() {
    return this.absVal;
  }

  getPhase() {
    return this.phase;
  }

  toString() {
    return `${this.getA()} + ${this.getB()}i`;
  }

  add(num) {
    return new ComplexNum(this.getA() + num.getA(), this.getB() + num.getB());
  }

  subtract(num) {
    return new ComplexNum(this.getA() - num.getA(), this.getB() - num.getB());
  }

  product(num) {
    if (typeof num === "number") {
      return new ComplexNum(this.a * num, this.b * num);
    } else {
      const a = this.a;
      const b = this.b;
      const x = num.getA();
      const y = num.getB();
      return new ComplexNum(a * x - b * y, a * y + b * x);
    }
  }

  static cis(phase) {
    return new ComplexNum(Math.cos(phase), Math.sin(phase));
  }
}
////

//grid square side length
let sLength = 5;
const dl_hr = 0.04;
const dl_mr = 0.1;
const dl_lr = 0.2;

//spacial scaling factor
const Scale = 50; // Fixed physics scale; camera magnification is independent.
let zoom = 1;
const MIN_ZOOM = 0.5, MAX_ZOOM = 4;
let zoomRebuildAfter = 0;

function gridView() {
  const cell = sLength * zoom;
  const x0 = Math.floor(-orig[0] / cell);
  const y0 = Math.ceil(orig[1] / cell);
  const drawX = orig[0] + x0 * cell;
  const drawY = orig[1] - y0 * cell;
  return {x0, y0, drawX, drawY, cell,
    cols: Math.ceil((width - drawX) / cell),
    rows: Math.ceil((height - drawY) / cell)};
}

function setViewZoom(value, screenX, screenY) {
  const next = Math.max(MIN_ZOOM, Math.min(MAX_ZOOM, value));
  if (!Number.isFinite(next) || next === zoom || one_point) return;
  // Preserve the exact world coordinate beneath the cursor (without snapping).
  const ratio = next / zoom;
  orig[0] = screenX - (screenX - orig[0]) * ratio;
  orig[1] = screenY - (screenY - orig[1]) * ratio;
  zoom = next;
  waitProcess = true;
  processingScheduled = false;
  zoomRebuildAfter = millis() + 100;
}

function mouseWheel(event) {
  if (!inRect(0, 0, width, height, mouseX, mouseY)) return;
  const overBox = antennas.some(a => a.getShowBox() && inRect(
    conScreenX(a.getXBox()), conScreenY(a.getYBox()),
    windowWidth * 400 / 1400, windowHeight * 370 / 1000, mouseX, mouseY));
  if (overBox || one_point || mousePressedFlag) return false;
  let delta = event.deltaY === undefined ? event.delta : event.deltaY;
  if (event.deltaMode === 1) delta *= 16;
  else if (event.deltaMode === 2) delta *= height;
  if (Number.isFinite(delta)) setViewZoom(zoom * Math.exp(-Math.max(-500, Math.min(500, delta)) * 0.0015), mouseX, mouseY);
  return false;
}

function keyPressed(event) {
  if (key === '0' && !(event && (event.ctrlKey || event.metaKey || event.altKey))) {
    setViewZoom(1, width / 2, height / 2);
    return false;
  }
}

let width = 1000;
let height = 800;

let prevWidth;
let prevHeight;

let N = Math.ceil(height / sLength);
let M = Math.ceil(width / sLength);

let extra_Width = 400;
let extra_height = 200;

let freq_button_offset = 0;
let speed_button_offset = 0;

let amp_button_offset = [];
let phase_button_offset = [];
let sep_button_offset = [];
let flags_amp_button = [];
let flags_phase_button = [];
let flags_sep_button = [];
let delButtonPressed = [];

let orig = [500.1, 400.1];

let resolution = 2;

//paramaters used for squiz functions

const k1 = 0.5;
const k1B = 2.5;

const k2 = 0.4;

const k3 = 0.2;
const k4 = 0.05;

const k5 = 0.05;
const k6 = 0.01;

const k7 = 0.6;

const k8 = 0.25;

//total elapsed time
let time = 0;

let timeSim = 0;
//temporal scaling factor
const timeScale = 8000;

//list of antenna objects
let antennas = [];

//electric field phasor
let E_phasor = [];

//magnetic field phasor
let B_phasor;

//EM field phasor
// World-coordinate cache for setup/panning only; values are Cartesian phasors.
let EM_phase_amp_map = new Map();

let A_map = new Map();

//mouse states
let mouseRelease = false;
let mousePressedFlag = false;

//states
let addNew_dipole_tx = false;
let change_dipole_tx = false;
let waitProcess = true;
let simulate = false;
let pause = false;

//field GUI states
let show_BField = false;
let show_EField = true;
let show_EnergyFlux = false;

let dipole_antenna_pressed = false;
let one_point = false;

let p1 = new Float32Array(2);
let p2 = new Float32Array(2);

const minFreq = 0.2;
const maxFreq = 0.8;
let freq = 0.5 * (minFreq + maxFreq);


const minSpeed = 1;
const maxSpeed = 7;
let c = 0.8 * minSpeed + 0.2 * maxSpeed;
let k = (2 * Math.PI * freq) / c;

//maximum current amplitude
const maxAmp = 25;

//default current magnitude
const defAmp = 10;

//
const maxSep = 0.4;
const defaultDipoleSep = 0.2;

//frequency slider
let freq_slider_on = false;
let freq_slider_process = false;

//speed of EM waves slider flags

let speed_slider_on = false;
let speed_slider_process = false;

// amplitude slider flags
let amp_slider_on = false;
let amp_slider_process = false;

// phase slider flags
let phase_slider_on = false;
let phase_slider_process = false;

let arrow_spacing = 5;
const maxArrowLen = 14;
const thickDipole = 12;

//storage
segments = [];
currentSegments = [];
currents = [];
current = [];
r = [];
currentSegement = [];

tempE1 = [];
tempE2 = [];

tempE = [];
tempB1 = new ComplexNum(0, 0);
tempB2 = new ComplexNum(0, 0);

let zero = new ComplexNum(0, 0);

let isMouseInStatBox = false;
let processingScheduled = false

let p1F = new Float32Array(2);
let p2F = new Float32Array(2);
let p3F = new Float32Array(2);
let p4F = new Float32Array(2);

let lastMouseX = 0;
let lastMouseY = 0;

let const_rFactor1 = (Scale * Scale) / (sLength * sLength);
let const_rFactor2 = Scale / (2 * sLength);

function squiz(l, k, levels = 256) {
  if (l <= 0) {
    return 0;
  }

  // Step 1: compute the squiz output (same as original)
  let ans = (k * l) / (1 + k * l);

  return ans;
}

function length2D(p1, p2) {
  return Math.sqrt(
    (p1[0] - p2[0]) * (p1[0] - p2[0]) + (p1[1] - p2[1]) * (p1[1] - p2[1])
  );
}

function conMyX(x) {
  return (sLength / Scale) * Math.round((x - orig[0]) / (sLength * zoom));
}

function conMyY(y) {
  return (sLength / Scale) * Math.round((orig[1] - y) / (sLength * zoom));
}

function conScreenX(x) {
  return x * Scale * zoom + orig[0];
}

function conScreenY(y) {
  return orig[1] - y * Scale * zoom;
}

function inRect(xRect, yRect, W, H, xPoint, yPoint) {
  let ans = xRect <= xPoint && xPoint <= xRect + W;
  ans = ans && yRect <= yPoint && yPoint <= yRect + H;
  return ans;
}

// is mouse point inside a general rectangle. Order of vertices should be clockwise or counter clockwise
function inRectGen(p1, p2, p3, p4, xMouse, yMouse) {
  let ans = false;

  //vectors that are two perpendicular sides of the rectangle

  sideVec1 = [p2[0] - p1[0], p2[1] - p1[1]];
  sideVec2 = [p3[0] - p2[0], p3[1] - p2[1]];

  let len1 = Math.sqrt(sideVec1[0] * sideVec1[0] + sideVec1[1] * sideVec1[1]);
  let len2 = Math.sqrt(sideVec2[0] * sideVec2[0] + sideVec2[1] * sideVec2[1]);

  //normalise side vectors

  sideVec1[0] /= len1;
  sideVec1[1] /= len1;

  sideVec2[0] /= len2;
  sideVec2[1] /= len2;

  let center = [
    0.25 * (p1[0] + p2[0] + p3[0] + p4[0]),
    0.25 * (p1[1] + p2[1] + p3[1] + p4[1]),
  ];

  let xTag = xMouse - center[0];
  let yTag = yMouse - center[1];

  //calculate projection of (xTag,yTag) vector on side vectors

  let proj1 = sideVec1[0] * xTag + sideVec1[1] * yTag;
  let proj2 = sideVec2[0] * xTag + sideVec2[1] * yTag;

  ans = Math.abs(proj1) <= 0.5 * len1 && Math.abs(proj2) <= 0.5 * len2;

  return ans;
}

function thickLine(p1, p2, thickness, c) {
  let dN = [];

  let A = [];
  let B = [];
  let C = [];
  let D = [];

  dN[0] = p2[1] - p1[1];
  dN[1] = p1[0] - p2[0];

  let len = Math.sqrt(dN[0] * dN[0] + dN[1] * dN[1]);

  dN[0] *= thickness / len;
  dN[1] *= thickness / len;

  A[0] = Math.round(p1[0] + 0.5 * dN[0]);
  A[1] = Math.round(p1[1] + 0.5 * dN[1]);

  B[0] = Math.round(p2[0] + 0.5 * dN[0]);
  B[1] = Math.round(p2[1] + 0.5 * dN[1]);

  C[0] = Math.round(p2[0] - 0.5 * dN[0]);
  C[1] = Math.round(p2[1] - 0.5 * dN[1]);

  D[0] = Math.round(p1[0] - 0.5 * dN[0]);
  D[1] = Math.round(p1[1] - 0.5 * dN[1]);

  stroke(c[0], c[1], c[2]);
  fill(c[0], c[1], c[2]);

  beginShape();

  vertex(A[0], A[1]);
  vertex(B[0], B[1]);
  vertex(C[0], C[1]);
  vertex(D[0], D[1]);

  endShape();
}

// Same softened Green function and segment weights as the original.
// Allocate only the two result objects, never inside the segment loop.
function calcA(myX, myY) {
  let axRe = 0, axIm = 0, ayRe = 0, ayIm = 0;
  for (let i = 0; i < antennas.length; i++) {
    const antenna = antennas[i];
    const segments = antenna.getSegments();
    const currents = antenna.getCurrentSegments();
    const dl = antenna.getDl();
    const drive = antenna.getI0();
    const driveRe = drive.a, driveIm = drive.b;
    for (let j = 0; j < segments.length; j++) {
      const dx = myX - segments[j][0];
      const dy = myY - segments[j][1];
      const distance = Math.sqrt(dx * dx + dy * dy) + 0.1;
      const invR = 1 / distance;
      const gRe = Math.cos(-k * distance) * invR;
      const gIm = Math.sin(-k * distance) * invR;
      const jx = currents[j][0] * dl;
      const jy = currents[j][1] * dl;
      const jxRe = jx * driveRe, jxIm = jx * driveIm;
      const jyRe = jy * driveRe, jyIm = jy * driveIm;
      axRe += jxRe * gRe - jxIm * gIm;
      axIm += jxRe * gIm + jxIm * gRe;
      ayRe += jyRe * gRe - jyIm * gIm;
      ayIm += jyRe * gIm + jyIm * gRe;
    }
  }
  return [new ComplexNum(axRe, axIm), new ComplexNum(ayRe, ayIm)];
}

// Each antenna owns one double-precision unit-drive A grid with a one-cell halo.
// WeakMap entries disappear with deleted antennas; panning retains only the last
// viewport per antenna, reusing its overlap instead of growing an unbounded cache.
let antennaBasisCache = new WeakMap();
let fieldWork = null;
let fieldCacheSignature = null;
let lastProcessingStats = null;

function unitPotentialGrid(antenna, cols, rows, x0, y0, stats) {
  const old = antennaBasisCache.get(antenna);
  const step = sLength / Scale;
  const reusable = old && old.segments === antenna.segments &&
    old.currents === antenna.currentSegments && old.dl === antenna.dl &&
    old.k === k && old.step === step;
  if (reusable && old.cols === cols && old.rows === rows && old.x0 === x0 && old.y0 === y0) {
    stats.reusedSamples += cols * rows;
    return old;
  }
  const size = cols * rows;
  const grid = { cols, rows, x0, y0, step, k, dl: antenna.dl,
    segments: antenna.segments, currents: antenna.currentSegments,
    axRe: new Float64Array(size), axIm: new Float64Array(size),
    ayRe: new Float64Array(size), ayIm: new Float64Array(size) };
  const segments = antenna.segments, currents = antenna.currentSegments;
  for (let j = 0; j < rows; j++) {
    const gy = y0 - j;
    const oldRow = reusable ? old.y0 - gy : -1;
    for (let i = 0; i < cols; i++) {
      const gx = x0 + i, idx = j * cols + i;
      const oldCol = reusable ? gx - old.x0 : -1;
      if (reusable && oldCol >= 0 && oldCol < old.cols && oldRow >= 0 && oldRow < old.rows) {
        const src = oldRow * old.cols + oldCol;
        grid.axRe[idx] = old.axRe[src]; grid.axIm[idx] = old.axIm[src];
        grid.ayRe[idx] = old.ayRe[src]; grid.ayIm[idx] = old.ayIm[src];
        stats.reusedSamples++;
        continue;
      }
      const x = gx * step, y = gy * step;
      let xr = 0, xi = 0, yr = 0, yi = 0;
      for (let n = 0; n < segments.length; n++) {
        const dx = x - segments[n][0], dy = y - segments[n][1];
        const distance = Math.sqrt(dx * dx + dy * dy) + 0.1;
        const gr = Math.cos(-k * distance) / distance;
        const gi = Math.sin(-k * distance) / distance;
        const jx = currents[n][0] * antenna.dl, jy = currents[n][1] * antenna.dl;
        xr += jx * gr; xi += jx * gi;
        yr += jy * gr; yi += jy * gi;
      }
      grid.axRe[idx] = xr; grid.axIm[idx] = xi;
      grid.ayRe[idx] = yr; grid.ayIm[idx] = yi;
      stats.integratedSamples++;
    }
  }
  antennaBasisCache.set(antenna, grid);
  return grid;
}

function startProcessingNewSetup() {
  const begin = performance.now();
  const dl = resolution === 1 ? dl_lr : resolution === 2 ? dl_mr : dl_hr;
  for (const a of antennas) a.setDl(dl);
  k = 2 * Math.PI * freq / c;
  const view = gridView();
  N = view.rows; M = view.cols;
  const cols = M + 2, rows = N + 2, size = cols * rows;
  const x0 = view.x0 - 1;
  const y0 = view.y0 + 1;
  const stats = {integratedSamples: 0, reusedSamples: 0, antennaCount: antennas.length};

  // Detect changes independently of UI invalidation, including deletion/reordering.
  const signature = {freq, c, sLength, antennas: antennas.map(a => ({
    antenna: a, segments: a.segments, currents: a.currentSegments, dl: a.dl,
    re: a.I0.a, im: a.I0.b
  }))};
  const prev = fieldCacheSignature;
  const same = prev && prev.freq === freq && prev.c === c && prev.sLength === sLength &&
    prev.antennas.length === antennas.length && signature.antennas.every((a,i) => {
      const b = prev.antennas[i];
      return a.antenna === b.antenna && a.segments === b.segments &&
        a.currents === b.currents && a.dl === b.dl && a.re === b.re && a.im === b.im;
    });
  // Limit the aggregate panning cache to roughly four current viewports.
  if (!same || EM_phase_amp_map.size > 4 * M * N) EM_phase_amp_map = new Map();
  fieldCacheSignature = signature;

  if (!fieldWork || fieldWork.axRe.length !== size) {
    fieldWork = {axRe:new Float64Array(size),axIm:new Float64Array(size),
                 ayRe:new Float64Array(size),ayIm:new Float64Array(size)};
  } else {
    fieldWork.axRe.fill(0); fieldWork.axIm.fill(0);
    fieldWork.ayRe.fill(0); fieldWork.ayIm.fill(0);
  }
  const {axRe,axIm,ayRe,ayIm} = fieldWork;
  for (const a of antennas) {
    const g = unitPotentialGrid(a, cols, rows, x0, y0, stats);
    const re = a.I0.a, im = a.I0.b;
    for (let n = 0; n < size; n++) {
      axRe[n] += re * g.axRe[n] - im * g.axIm[n];
      axIm[n] += re * g.axIm[n] + im * g.axRe[n];
      ayRe[n] += re * g.ayRe[n] - im * g.ayIm[n];
      ayIm[n] += re * g.ayIm[n] + im * g.ayRe[n];
    }
  }
  // Derivatives of the combined A are equivalent to combining per-antenna E/B.
  // Do the finite-difference stencil just once, using numeric array neighbors.
  const d2 = Scale * Scale / (sLength * sLength);
  const d1 = Scale / (2 * sLength), omega = 2 * Math.PI * freq, ck = c / k;
  for (let j = 0; j < N; j++) {
    const y = (view.y0 - j) * (sLength / Scale);
    for (let i = 0; i < M; i++) {
      const n = (j + 1) * cols + i + 1;
      const l = n - 1, r = n + 1, u = n - cols, d = n + cols;
      const xxr = (axRe[r] + axRe[l] - 2 * axRe[n]) * d2;
      const xxi = (axIm[r] + axIm[l] - 2 * axIm[n]) * d2;
      const yyr = (ayRe[u] + ayRe[d] - 2 * ayRe[n]) * d2;
      const yyi = (ayIm[u] + ayIm[d] - 2 * ayIm[n]) * d2;
      const xyyr = (ayRe[d-1] + ayRe[u+1] - (ayRe[d+1] + ayRe[u-1])) * 0.25 * d2;
      const xyyi = (ayIm[d-1] + ayIm[u+1] - (ayIm[d+1] + ayIm[u-1])) * 0.25 * d2;
      const xyxr = (axRe[d-1] + axRe[u+1] - (axRe[d+1] + axRe[u-1])) * 0.25 * d2;
      const xyxi = (axIm[d-1] + axIm[u+1] - (axIm[d+1] + axIm[u-1])) * 0.25 * d2;
      EM_phase_amp_map.set(`${(view.x0 + i) * (sLength / Scale)},${y}`, {
        ExRe: omega * axIm[n] + ck * (xxi + xyyi),
        ExIm: -omega * axRe[n] - ck * (xxr + xyyr),
        EyRe: omega * ayIm[n] + ck * (yyi + xyxi),
        EyIm: -omega * ayRe[n] - ck * (yyr + xyxr),
        BRe: (ayRe[r] - ayRe[l]) * d1 - (axRe[u] - axRe[d]) * d1,
        BIm: (ayIm[r] - ayIm[l]) * d1 - (axIm[u] - axIm[d]) * d1
      });
    }
  }
  refreshVisibleFields(true);
  waitProcess = false; simulate = true; processingScheduled = false;
  stats.elapsedMs = performance.now() - begin;
  lastProcessingStats = stats;
}


// Stage 1: keep double precision to avoid introducing quantization changes.
// The string-keyed world cache is read only when the view or setup changes.
let visibleFields = null;
let fieldImage = null;
let visibleCache = null;
let visibleOriginX = NaN, visibleOriginY = NaN, visibleStep = NaN;
let frameCos = 1, frameSin = 0, lastFramePhase = NaN;

function updateFramePhase() {
  const theta = timeSim * 2 * Math.PI * freq;
  if (theta !== lastFramePhase) {
    frameCos = Math.cos(theta);
    frameSin = Math.sin(theta);
    lastFramePhase = theta;
  }
}

function refreshVisibleFields(force = false) {
  const view = gridView();
  const cols = view.cols, rows = view.rows;
  const ox = view.x0, oy = view.y0;
  if (visibleFields) {
    visibleFields.drawX = view.drawX; visibleFields.drawY = view.drawY;
    visibleFields.cell = view.cell;
  }
  const resized = !visibleFields || visibleFields.cols !== cols ||
                  visibleFields.rows !== rows;
  if (!force && !resized && visibleCache === EM_phase_amp_map &&
      visibleOriginX === ox && visibleOriginY === oy && visibleStep === sLength) return;
  if (resized) {
    const size = cols * rows;
    visibleFields = {
      cols, rows, drawX: view.drawX, drawY: view.drawY, cell: view.cell,
      ExRe: new Float64Array(size), ExIm: new Float64Array(size),
      EyRe: new Float64Array(size), EyIm: new Float64Array(size),
      BRe: new Float64Array(size), BIm: new Float64Array(size),
      valid: new Uint8Array(size)
    };
    fieldImage = createImage(cols, rows);
    fieldImage.loadPixels();
  }
  const f = visibleFields;
  for (let j = 0; j < rows; j++) {
    const y = (oy - j) * (sLength / Scale);
    for (let i = 0; i < cols; i++) {
      const idx = j * cols + i;
      const data = EM_phase_amp_map.get(`${(ox + i) * (sLength / Scale)},${y}`);
      f.valid[idx] = data ? 1 : 0;
      f.ExRe[idx] = data ? data.ExRe : 0;
      f.ExIm[idx] = data ? data.ExIm : 0;
      f.EyRe[idx] = data ? data.EyRe : 0;
      f.EyIm[idx] = data ? data.EyIm : 0;
      f.BRe[idx] = data ? data.BRe : 0;
      f.BIm[idx] = data ? data.BIm : 0;
    }
  }
  visibleCache = EM_phase_amp_map;
  visibleOriginX = ox;
  visibleOriginY = oy;
  visibleStep = sLength;
}

function setup() {
  createCanvas(windowWidth, windowHeight);
  prevWidth = windowWidth
  prevHeight = windowHeight
  
  pixelDensity(1);
  width = Math.round((1 / 1.4) * windowWidth);
  height = Math.round((8 / 10) * windowHeight);

  current = [];

  one_point = false;

  freq_button_offset =
    ((freq - minFreq) / (maxFreq - minFreq)) *
    ((200 - 30) / 1400) *
    windowWidth;
  speed_button_offset =
    ((c - minSpeed) / (maxSpeed - minSpeed)) *
    ((200 - 30) / 1400) *
    windowWidth;

  background(0);
  frameRate(60);

  let myP1 = [conMyX(0.05 + width / 2), conMyY((2 * height) / 3)];
  let myP2 = [conMyX(width / 2), conMyY(height / 3)];

  let amp = defAmp;
  let phase = 0;

  dipole = new Dipole(c / freq, myP1, myP2, amp, phase, thickDipole / Scale);

  if (resolution == 1) {
    dipole.setDl(dl_lr);
  } else if (resolution == 2) {
    dipole.setDl(dl_mr);
  } else {
    dipole.setDl(dl_hr);
  }

  antennas.push(dipole);

  amp_button_offset.push((defAmp / maxAmp) * ((300 - 40) / 1400) * windowWidth);
  phase_button_offset.push(0);
  sep_button_offset.push(
    (defaultDipoleSep * ((300 - 40) / 1400) * windowWidth) / maxSep
  );

  flags_amp_button.push(false);
  flags_phase_button.push(false);
  flags_sep_button.push(false);
  delButtonPressed.push(false);
  waitProcess = true;
  simulate = false;
}


let flagA = true;
let flagB = false;
function draw() {
  


  
  if(mousePressedFlag&&!isMouseInStatBox&& inRect(0, 0, width, height, mouseX, mouseY)
){
    
  
    
  
   if(flagA){
    lastMouseX = mouseX
    lastMouseY = mouseY
    flagA = false
   }
    
    if(mouseX!=lastMouseX||mouseY!=lastMouseY||one_point){
    
    flagB = true
      
    }   
    
   
    orig[0] += mouseX-lastMouseX
    orig[1] += mouseY-lastMouseY
    
   
    
    
    lastMouseX = mouseX
    lastMouseY = mouseY

  }
  
 else if(!mousePressedFlag){
    
     lastMouseX = mouseX
    lastMouseY = mouseY
   flagA = true
    
  }
  
  
  
   if(flagB && mouseRelease&&!one_point){
    
    
    flagB = false 
    flagA = true
   simulate = false
  waitProcess = true
  processingScheduled = false
    mouseRelease = false 
    
    
        
 
  }
  
  
  // All resolutions use the 60 FPS target; work per frame determines actual FPS.


  if (windowWidth != prevWidth || windowHeight != prevHeight) {
    freq_button_offset = (freq_button_offset * windowWidth) / prevWidth;
    speed_button_offset = (speed_button_offset * windowWidth) / prevWidth;

    for (let b = 0; b < antennas.length; b++) {
      amp_button_offset[b] *= windowWidth / prevWidth;
      phase_button_offset[b] *= windowWidth / prevWidth;
      sep_button_offset[b] *= windowWidth / prevWidth;
    }

    // Repack the flat visible grid on both growth and shrinkage.
    simulate = false;
    waitProcess = true;
    processingScheduled = false;

    prevWidth = windowWidth;
    prevHeight = windowHeight;

    resizeCanvas(windowWidth, windowHeight);

    width = Math.floor((1 / 1.4) * windowWidth);
    height = Math.floor((8 / 10) * windowHeight);
  }

  background(0, 0, 0);

  if (simulate) {
    if (!pause) {
      timeSim += millis() / timeScale - time;
    }
    time = millis() / timeScale;

    updateFramePhase();
    refreshVisibleFields();
    const f = visibleFields;
    const fieldPixels = fieldImage.pixels;
    const ct = frameCos, st = frameSin;
    for (let j = 0; j < f.rows; j++) {
      for (let i = 0; i < f.cols; i++) {
        const idx = j * f.cols + i;
        const Ex_t = f.ExRe[idx] * ct - f.ExIm[idx] * st;
        const Ey_t = f.EyRe[idx] * ct - f.EyIm[idx] * st;
        const B_t = f.BRe[idx] * ct - f.BIm[idx] * st;

        // Color calculation
        let r, g, b, a;

        if (show_BField) {
          let colorSizeB_Blue = 255 * squiz(-B_t * Math.abs(B_t), k1, 256);
          let colorSizeB_Red = 255 * squiz(B_t * Math.abs(B_t), k1, 256);
          let colorSizeMag = squiz(Math.abs(B_t), k1B, 256);

          r = Math.round(colorSizeB_Red);
          g = 0;
          b = Math.round(colorSizeB_Blue);
          a = colorSizeMag;
        } else if (show_EField) {
          let E_mag_2 = Ex_t * Ex_t + Ey_t * Ey_t;
          let colorSizeEA = squiz(E_mag_2, k3, 256);
          let colorSizeEB = 255 * squiz(E_mag_2, k4, 256);

          // Convert HSB to RGB
          let hue = (255 - colorSizeEB) / 500;
          let sat = 1.0;
          let bright = 1.0;
          let c_val = bright * sat;
          let x_val = c_val * (1 - Math.abs(((hue * 6) % 2) - 1));
          let m_val = bright - c_val;

          let r1, g1, b1;
          let h = hue * 6;
          if (h < 1) {
            r1 = c_val;
            g1 = x_val;
            b1 = 0;
          } else if (h < 2) {
            r1 = x_val;
            g1 = c_val;
            b1 = 0;
          } else if (h < 3) {
            r1 = 0;
            g1 = c_val;
            b1 = x_val;
          } else if (h < 4) {
            r1 = 0;
            g1 = x_val;
            b1 = c_val;
          } else if (h < 5) {
            r1 = x_val;
            g1 = 0;
            b1 = c_val;
          } else {
            r1 = c_val;
            g1 = 0;
            b1 = x_val;
          }

          r = Math.round((r1 + m_val) * 255);
          g = Math.round((g1 + m_val) * 255);
          b = Math.round((b1 + m_val) * 255);
          a = colorSizeEA;
        } else if (show_EnergyFlux) {
          let Energy_flux_mag = (Ex_t * Ex_t + Ey_t * Ey_t) * Math.abs(B_t);
          let colorSizeEnergyA = squiz(Energy_flux_mag, k5, 256);
          let colorSizeEnergyB = 255 * squiz(Energy_flux_mag, k6, 256);

          let hue = (255 - colorSizeEnergyB) / 500;
          let sat = 1.0;
          let bright = 1.0;
          let c_val = bright * sat;
          let x_val = c_val * (1 - Math.abs(((hue * 6) % 2) - 1));
          let m_val = bright - c_val;

          let r1, g1, b1;
          let h_val = hue * 6;
          if (h_val < 1) {
            r1 = c_val;
            g1 = x_val;
            b1 = 0;
          } else if (h_val < 2) {
            r1 = x_val;
            g1 = c_val;
            b1 = 0;
          } else if (h_val < 3) {
            r1 = 0;
            g1 = c_val;
            b1 = x_val;
          } else if (h_val < 4) {
            r1 = 0;
            g1 = x_val;
            b1 = c_val;
          } else if (h_val < 5) {
            r1 = x_val;
            g1 = 0;
            b1 = c_val;
          } else {
            r1 = c_val;
            g1 = 0;
            b1 = x_val;
          }

          r = Math.round((r1 + m_val) * 255);
          g = Math.round((g1 + m_val) * 255);
          b = Math.round((b1 + m_val) * 255);
          a = colorSizeEnergyA;
        }

        // One pixel per field cell. RGB already includes the original brightness.
        const pixelIndex = idx * 4;
        fieldPixels[pixelIndex] = r * a;
        fieldPixels[pixelIndex + 1] = g * a;
        fieldPixels[pixelIndex + 2] = b * a;
        fieldPixels[pixelIndex + 3] = 255;
      }
    }
    fieldImage.updatePixels();
    // Preserve exact cell size and clip partial edge cells rather than stretching.
    drawingContext.save();
    drawingContext.beginPath();
    drawingContext.rect(0, 0, width, height);
    drawingContext.clip();
    drawingContext.imageSmoothingEnabled = false;
    image(fieldImage, f.drawX, f.drawY, f.cols * f.cell, f.rows * f.cell);
    drawingContext.restore();
  }

  
  
  
  
  for (let b = 0; b < antennas.length && !simulate; b++) {
    antenna = antennas[b];

    if (antenna.getType() == "Dipole") {
      segments = antennas[b].getSegments();

      for (let j = 0; j < segments.length - 1; j++) {
        let x = segments[j][0];
        let x_next = segments[j + 1][0];

        let y = segments[j][1];
        let y_next = segments[j + 1][1];

        if (antennas[b].getSegFlags()[j]) {
          thickLine(
            [conScreenX(x), conScreenY(y)],
            [conScreenX(x_next), conScreenY(y_next)],
            thickDipole * zoom,
            [250, 250, 250]
          );
        }
      }
    }
  }

  if (addNew_dipole_tx) {
    let is_mouse_pos_ok = inRect(0, 0, width, height, mouseX, mouseY);

    if (is_mouse_pos_ok) {
      if (!one_point && mousePressedFlag) {
        p1[0] = mouseX;
        p1[1] = mouseY;
        

        one_point = true;

        mousePressedFlag = false;
      }

      if (one_point) {
        if (mousePressedFlag) {
          p2[0] = mouseX;
          p2[1] = mouseY;

          addNew_dipole_tx = false;
          dipole_antenna_pressed = false;
          one_point = false;

          let myP1 = [conMyX(p1[0]), conMyY(p1[1])];
          let myP2 = [conMyX(p2[0]), conMyY(p2[1])];

          let amp = defAmp;
          let phase = 0;
          

          dipole = new Dipole(
            c / freq,
            myP1,
            myP2,
            amp,
            phase,
            thickDipole / Scale
          );

          if (resolution == 1) {
            dipole.setDl(dl_lr);
          } else if (resolution == 2) {
            dipole.setDl(dl_mr);
          } else {
            dipole.setDl(dl_hr);
          }

          antennas.push(dipole);

          amp_button_offset.push(
            (defAmp / maxAmp) * ((300 - 40) / 1400) * windowWidth
          );
          phase_button_offset.push(0);
          sep_button_offset.push(
            (defaultDipoleSep * ((300 - 40) / 1400) * windowWidth) / maxSep
          );

          flags_amp_button.push(false);
          flags_phase_button.push(false);
          flags_sep_button.push(false);
          delButtonPressed.push(false);
          
          
         
          
        A_map = new Map()
        EM_phase_amp_map = new Map()

          mousePressedFlag = false;
        }

        thickLine([mouseX, mouseY], p1, thickDipole * zoom, [250, 250, 250]);

        let lenAntenna = Math.sqrt(
          (mouseX - p1[0]) * (mouseX - p1[0]) +
            (mouseY - p1[1]) * (mouseY - p1[1])
        );

        let tempX1 =
          0.5 * (mouseX + p1[0]) -
          (0.5 * (p1[0] - mouseX) * defaultDipoleSep * Scale * zoom) / lenAntenna;
        let tempX2 =
          0.5 * (mouseX + p1[0]) +
          (0.5 * (p1[0] - mouseX) * defaultDipoleSep * Scale * zoom) / lenAntenna;

        let tempY1 =
          0.5 * (mouseY + p1[1]) -
          (0.5 * (p1[1] - mouseY) * defaultDipoleSep * Scale * zoom) / lenAntenna;
        let tempY2 =
          0.5 * (mouseY + p1[1]) +
          (0.5 * (p1[1] - mouseY) * defaultDipoleSep * Scale * zoom) / lenAntenna;

        thickLine([tempX1, tempY1], [tempX2, tempY2], thickDipole * zoom, [0, 0, 0]);
      }
    }
  }

  //menu background

  stroke(200, 200, 200);
  fill(200, 200, 200);
  rect(width, 0, windowWidth - width, height);
  rect(0, height, windowWidth, height / 2);

  //dipole antenna button

  stroke(0, 0, 0);
  fill(0, 0, 0);

  fill(150, 150, 150);

  if (dipole_antenna_pressed) {
    fill(120, 120, 120);
    stroke(250, 250, 250);
  }

  rect(
    (35 / 1400) * windowWidth,
    (860 / 1000) * windowHeight,
    (175 / 1400) * windowWidth,
    (70 / 1000) * windowHeight
  );

  fill(150, 150, 150);
  stroke(0, 0, 0);

  if (
    !dipole_antenna_pressed &&
    inRect(
      (35 / 1400) * windowWidth,
      (860 / 1000) * windowHeight,
      (175 / 1400) * windowWidth,
      (70 / 1000) * windowHeight,
      mouseX,
      mouseY
    ) &&
    mousePressedFlag
  ) {
    dipole_antenna_pressed = true;
    mousePressedFlag = false;
    
   
  } else if (
    dipole_antenna_pressed &&
    inRect(
      (35 / 1400) * windowWidth,
      (860 / 1000) * windowHeight,
      (175 / 1400) * windowWidth,
      (70 / 1000) * windowHeight,
      mouseX,
      mouseY
    ) &&
    mousePressedFlag
  ) {
    dipole_antenna_pressed = false;
    mousePressedFlag = false;
    addNew_dipole_tx = false;
    one_point = false
    
    
  }

  if (dipole_antenna_pressed) {
    addNew_dipole_tx = true;
  }

  stroke(0, 0, 0);
  fill(0, 0, 0);

  push();

  scale(windowWidth / 1400, windowHeight / 1000);

  textSize(23);
  text("Dipole Antenna", 44, 902);

  pop();

  //pause button

  stroke(0, 0, 0);

  if (!pause) {
    fill(240, 0, 0);
  } else {
    fill(0, 240, 0);
  }

  if (
    inRect(
      (920 / 1400) * windowWidth,
      (860 / 1000) * windowHeight,
      (80 / 1400) * windowWidth,
      (80 / 1000) * windowHeight,
      mouseX,
      mouseY
    )
  ) {
    if (!pause && mousePressedFlag&&!isMouseInStatBox) {
      mousePressedFlag = false;
      mouseRelease = false;
      pause = true;
    } else if (pause && mousePressedFlag&&!isMouseInStatBox) {
      mousePressedFlag = false;
      mouseRelease = false;
      pause = false;
    }
  }

  rect(
    (920 / 1400) * windowWidth,
    (860 / 1000) * windowHeight,
    (80 / 1400) * windowWidth,
    (80 / 1000) * windowHeight
  );

  if (!pause) {
    fill(0, 0, 0);

    push();

    scale(windowWidth / 1400, windowHeight / 1000);
    textSize(24);
    text("Pause", 925, 852);
    pop();
  } else {
    fill(0, 0, 0);

    push();

    scale(windowWidth / 1400, windowHeight / 1000);
    textSize(24);
    text("Resume", 915, 852);
    pop();
  }

  
  fill(150, 150, 150);

  if (
    inRect(
      (1150 / 1400) * windowWidth,
      (615 / 1000) * windowHeight,
      (180 / 1400) * windowWidth,
      (60 / 1000) * windowHeight,
      mouseX,
      mouseY
    )&&!isMouseInStatBox
  ) {
    fill(120, 120, 120);

    if (mousePressedFlag) {
      antennas = [];
      amp_button_offset = [];
      phase_button_offset = [];
      sep_button_offset = [];
      flags_amp_button = [];
      flags_phase_button = [];
      sep_phase_button = [];
      one_point = false;
      simulate = false;
      waitProcess = true
      processingScheduled = false
      A_map = new Map()
      EM_phase_amp_map = new Map()
    }
  }

  rect(
    (1150 / 1400) * windowWidth,
    (615 / 1000) * windowHeight,
    (180 / 1400) * windowWidth,
    (60 / 1000) * windowHeight
  );

  push();
  fill(0, 0, 0);
  scale(windowWidth / 1400, windowHeight / 1000);
  textSize(26);
  text("Clear All", 1190, 653);

  textSize(20);

  pop();

  fill(150, 150, 150);

  //resolution buttons logic

 

  let flagChangeRes = false;
 if(!isMouseInStatBox){
  if (
    resolution != 3 &&
    inRect(
      (1188 / 1400) * windowWidth,
      (770 / 1000) * windowHeight,
      (110 / 1400) * windowWidth,
      (50 / 1000) * windowHeight,
      mouseX,
      mouseY
    ) &&
    mousePressedFlag
  ) {
    resolution = 3;
    mousePressedFlag = false;
    simulate = false;
    waitProcess = true
    
     A_map = new Map()
    EM_phase_amp_map = new Map()
    flagChangeRes = true;
    sLength = 2;
    arrow_spacing = 10;
  }

  if (
    resolution != 2 &&
    inRect(
      (1188 / 1400) * windowWidth,
      (840 / 1000) * windowHeight,
      (windowWidth * 110) / 1400,
      (windowHeight * 50) / 1000,
      mouseX,
      mouseY
    ) &&
    mousePressedFlag
  ) {
    resolution = 2;
    mousePressedFlag = false;
    simulate = false;
    waitProcess = true
     A_map = new Map()
    EM_phase_amp_map = new Map()
    flagChangeRes = true;
    sLength = 4;
    arrow_spacing = 5;
  }

  if (
    resolution != 1 &&
    inRect(
      (windowWidth * 1188) / 1400,
      (windowHeight * 910) / 1000,
      (windowWidth * 110) / 1400,
      (windowHeight * 50) / 1000,
      mouseX,
      mouseY
    ) &&
    mousePressedFlag
  ) {
    resolution = 1;
    mousePressedFlag = false;
    simulate = false;
    waitProcess = true
     A_map = new Map()
    EM_phase_amp_map = new Map()
    flagChangeRes = true;
    sLength = 5;
    arrow_spacing = 4;
  }

 }

  if (flagChangeRes) {
    N = Math.ceil(height / sLength);
    M = Math.ceil(width / sLength);

    const_rFactor1 = (Scale * Scale) / (sLength * sLength);
    const_rFactor2 = Scale / (2 * sLength);
  }

  // resolution buttons graphics

  if (resolution == 3) {
    fill(120, 120, 120);
    stroke(250, 250, 250);
  }

  rect(
    (windowWidth * 1188) / 1400,
    (windowHeight * 770) / 1000,
    (windowWidth * 110) / 1400,
    (windowHeight * 50) / 1000
  );

  stroke(0, 0, 0);
  fill(150, 150, 150);

  if (resolution == 2) {
    fill(120, 120, 120);
    stroke(250, 250, 250);
  }

  rect(
    (windowWidth * 1188) / 1400,
    (windowHeight * 840) / 1000,
    (windowWidth * 110) / 1400,
    (windowHeight * 50) / 1000
  );

  stroke(0, 0, 0);
  fill(150, 150, 150);

  if (resolution == 1) {
    fill(120, 120, 120);
    stroke(250, 250, 250);
  }

  rect(
    (windowWidth * 1188) / 1400,
    (windowHeight * 910) / 1000,
    (windowWidth * 110) / 1400,
    (windowHeight * 50) / 1000
  );

  stroke(0, 0, 0);
  fill(150, 150, 150);

  stroke(0, 0, 0);
  fill(0, 0, 0);
  textSize(25);
  push();
  scale(windowWidth / 1400, windowHeight / 1000);
  text("Resolution:", 1185, 750);

  text("high", 1220, 800);

  text("medium", 1198, 870);

  text("low", 1224, 942);
  pop();
  stroke(0, 0, 0);
  fill(150, 150, 150);

  //fields GUI

  //E field button
  if (show_EField) {
    fill(120, 120, 120);
    stroke(250, 250, 250);
  }

  rect(
    (windowWidth * 1130) / 1400,
    (windowHeight * 50) / 1000,
    (windowWidth * 148) / 1400,
    (windowHeight * 50) / 1000
  );

  push();
  stroke(0, 0, 0);
  fill(0, 0, 0);
  textSize(20);
  scale(windowWidth / 1400, windowHeight / 1000);
  text("Show E Field", 1144, 80);

  pop();

  stroke(0, 0, 0);
  fill(150, 150, 150);

  if (
    !show_EField &&
    inRect(
      (windowWidth * 1130) / 1400,
      (windowHeight * 50) / 1000,
      (windowWidth * 148) / 1400,
      (windowHeight * 50) / 1000,
      mouseX,
      mouseY
    ) &&
    mousePressedFlag
  ) {
    show_EField = true;
    show_BField = false;
    show_EnergyFlux = false;

    mousePressedFlag = false;
  }

  // B field button

  if (show_BField) {
    fill(120, 120, 120);
    stroke(250, 250, 250);
  }

  rect(
    (windowWidth * 1130) / 1400,
    (windowHeight * 150) / 1000,
    (windowWidth * 148) / 1400,
    (windowHeight * 50) / 1000
  );

  stroke(0, 0, 0);
  fill(0, 0, 0);
  push();
  textSize(20);
  scale(windowWidth / 1400, windowHeight / 1000);
  text("Show B Field", 1144, 180);
  pop();
  textSize(20);

  fill(150, 150, 150);

  if (
    !show_BField &&
    inRect(
      (windowWidth * 1130) / 1400,
      (windowHeight * 150) / 1000,
      (windowWidth * 148) / 1400,
      (windowHeight * 50) / 1000,
      mouseX,
      mouseY
    ) &&
    mousePressedFlag
  ) {
    show_EField = false;
    show_BField = true;
    show_EnergyFlux = false;

    mousePressedFlag = false;
  }

  //energy flux button

  if (show_EnergyFlux) {
    fill(120, 120, 120);
    stroke(250, 250, 250);
  }

  rect(
    (windowWidth * 1130) / 1400,
    (windowHeight * 250) / 1000,
    (windowWidth * 190) / 1400,
    (windowHeight * 50) / 1000
  );

  push();
  stroke(0, 0, 0);
  fill(0, 0, 0);
  scale(windowWidth / 1400, windowHeight / 1000);
  text("Show energy flux", 1144, 280);
  pop();
  textSize(20);

  fill(150, 150, 150);

  if (
    !show_EnergyFlux &&
    inRect(
      (windowWidth * 1130) / 1400,
      (windowHeight * 250) / 1000,
      (windowWidth * 190) / 1400,
      (windowHeight * 50) / 1000,
      mouseX,
      mouseY
    ) &&
    mousePressedFlag
  ) {
    show_EField = false;
    show_BField = false;
    show_EnergyFlux = true;

    mousePressedFlag = false;
  }

  // frequency slider

  stroke(0, 0, 0);
  fill(0, 0, 0);

  push();
  scale(windowWidth / 1400, windowHeight / 1000);
  text("Frequency", 1180, 340);
  pop();
  textSize(20);

  fill(250, 250, 250);

  rect(
    (windowWidth * 1130) / 1400,
    (windowHeight * 350) / 1000,
    (windowWidth * 200) / 1400,
    (windowHeight * 30) / 1000
  );

  //frequency slider button

  fill(110, 110, 110);

  if (
    mousePressedFlag &&
    inRect(
      (windowWidth * 1130) / 1400 + freq_button_offset,
      (windowHeight * 350) / 1000,
      (windowWidth * 30) / 1400,
      (windowHeight * 30) / 1000,
      mouseX,
      mouseY
    )
  ) {
    freq_slider_on = true;
  } else if (
    !mousePressedFlag ||
    !inRect(
      (windowWidth * 1130) / 1400,
      (windowHeight * 350) / 1000,
      (windowWidth * 200) / 1400,
      (windowHeight * 30) / 1000,
      mouseX,
      mouseY
    )
  ) {
    if(freq_slider_on){
      
      simulate = false;
    waitProcess = true
     A_map = new Map()
    EM_phase_amp_map = new Map()
      
       freq =
      minFreq +
      (freq_button_offset / ((windowWidth * (200 - 30)) / 1400)) *
        (maxFreq - minFreq);

    freq_slider_process = false;

    for (let i = 0; i < antennas.length; i++) {
      antennas[i].setWavelength(c / freq);
    }
      
    }
    freq_slider_on = false;
  }

  if (freq_slider_on) {
    freq_slider_process = true;
    

    fill(90, 90, 90);

    freq_button_offset = mouseX - (windowWidth * (1130 + 15)) / 1400;
    if (freq_button_offset < 0) {
      freq_button_offset = 0;
    } else if (
      freq_button_offset + (windowWidth * 30) / 1400 >
      (windowWidth * 200) / 1400
    ) {
      freq_button_offset = ((200 - 30) * windowWidth) / 1400;
    }
  }

  rect(
    (windowWidth * 1130) / 1400 + freq_button_offset,
    (windowHeight * 350) / 1000,
    (windowWidth * 30) / 1400,
    (windowHeight * 30) / 1000
  );

 
stroke(0, 0, 0);
fill(0, 0, 0);
push();
scale(windowWidth / 1400, windowHeight / 1000);
text("Speed", 1200, 420);
pop();
textSize(20);

// Slider bar background
stroke(0, 0, 0);
fill(250, 250, 250);
rect(
  (windowWidth * 1130) / 1400,
  (windowHeight * 430) / 1000,
  (windowWidth * 200) / 1400,
  (windowHeight * 30) / 1000
);
  
  fill(110,110,110)

// --- Handle slider button logic ---

// 1. Check if mouse is pressed on the slider button → Activate slider
if (
  mousePressedFlag &&
  inRect(
    (windowWidth * 1130) / 1400 + speed_button_offset,
    (windowHeight * 430) / 1000,
    (windowWidth * 30) / 1400,
    (windowHeight * 30) / 1000,
    mouseX,
    mouseY
  )
) {
  speed_slider_on = true;
}

// 2. If mouse is released OR moved out of slider bar → finalize speed update
else if (
  !mousePressedFlag ||
  !inRect(
    (windowWidth * 1130) / 1400,
    (windowHeight * 430) / 1000,
    (windowWidth * 200) / 1400,
    (windowHeight * 30) / 1000,
    mouseX,
    mouseY
  )
) {
  if (speed_slider_on) {
    // Stop simulation and reset maps for recalculation
    simulate = false;
    waitProcess = true;
    A_map = new Map();
    EM_phase_amp_map = new Map();

    // Update speed based on button position
    c =
      minSpeed +
      (speed_button_offset / ((windowWidth * (200 - 30)) / 1400)) *
      (maxSpeed - minSpeed);

    speed_slider_process = false;

    // Update wavelength for all antennas
    for (let i = 0; i < antennas.length; i++) {
      antennas[i].setWavelength(c / freq);
    }
  }
  speed_slider_on = false;
}

// 3. If slider is active → update button offset as mouse drags
if (speed_slider_on) {
  speed_slider_process = true;
  fill(90, 90, 90);

  speed_button_offset = mouseX - (windowWidth * (1130 + 15)) / 1400;

  if (speed_button_offset < 0) {
    speed_button_offset = 0;
  } else if (
    speed_button_offset + (windowWidth * 30) / 1400 >
    (windowWidth * 200) / 1400
  ) {
    speed_button_offset = ((200 - 30) * windowWidth) / 1400;
  }
}

// Draw the slider button
rect(
  (windowWidth * 1130) / 1400 + speed_button_offset,
  (windowHeight * 430) / 1000,
  (windowWidth * 30) / 1400,
  (windowHeight * 30) / 1000
);
  
  
   //simulation status button

  if (simulate) {
    fill(0, 250, 0);


    rect(
      (1150 / 1400) * windowWidth,
      (515 / 1000) * windowHeight,
      (180 / 1400) * windowWidth,
      (60 / 1000) * windowHeight
    );

    fill(0, 0, 0);
    push();

    scale(windowWidth / 1400, windowHeight / 1000);

    textSize(28);
    text("Simulating", 1170, 553);

    pop();
  }
  
  else{
    
    
     fill(220, 0, 0);

 
    rect(
      (1150 / 1400) * windowWidth,
      (515 / 1000) * windowHeight,
      (180 / 1400) * windowWidth,
      (60 / 1000) * windowHeight
    );

    fill(0, 0, 0);

    push();

    scale(windowWidth / 1400, windowHeight / 1000);

    textSize(28);
    text("Loading", 1185, 553);

    pop();
    

    
  }
  
  
  // process changes

  if (waitProcess && millis() >= zoomRebuildAfter) {
     
    
    if( !processingScheduled ){
   
     processingScheduled = true  
    }
    
    else{
      
       startProcessingNewSetup()
    }
    
    
  }


  
  //after processing

  if (simulate) {
    push();
    drawingContext.beginPath();
    drawingContext.rect(0, 0, width, height);
    drawingContext.clip();
    updateFramePhase();
    refreshVisibleFields();
    const f = visibleFields;
    const ct = frameCos, st = frameSin;
    for (let i = 0; i < f.cols; i += Math.max(1, Math.round(arrow_spacing / zoom))) {
      for (let j = 0; j < f.rows; j += Math.max(1, Math.round(arrow_spacing / zoom))) {
        const screenX = f.drawX + f.cell * i;
        const screenY = f.drawY + f.cell * j;
        if (screenX < 0 || screenY < 0 || screenX >= width || screenY >= height) continue;
        const idx = j * f.cols + i;
        if (!f.valid[idx]) continue;
        const Ex_t = f.ExRe[idx] * ct - f.ExIm[idx] * st;
        const Ey_t = f.EyRe[idx] * ct - f.EyIm[idx] * st;
        const B_t = f.BRe[idx] * ct - f.BIm[idx] * st;

        // === Draw E field arrows ===
        if (show_EField) {
          let E_mag = Math.sqrt(Ex_t * Ex_t + Ey_t * Ey_t);
          if (E_mag > 1e-8) {
            let arrow_l = maxArrowLen * squiz(E_mag, k2, 256);

            let dx = Ex_t * (arrow_l / E_mag);
            let dy = Ey_t * (arrow_l / E_mag);

            stroke(0, 0, 0);
            fill(0, 0, 0);

            circle(screenX, screenY, 2);
            line(screenX, screenY, screenX + dx, screenY - dy);
            noStroke();
          }
        }

        // === Draw Energy Flux arrows ===
        if (show_EnergyFlux) {
          let E_flux_mag = Math.sqrt(Ex_t * Ex_t + Ey_t * Ey_t) * Math.abs(B_t);
          if (E_flux_mag > 1e-8) {
            let arrow_l = maxArrowLen * squiz(E_flux_mag, k7, 256);

            // Poynting vector direction (S ~ E × B)
            let dx = Ey_t * B_t;
            let dy = -(Ex_t * B_t);

            dx *= arrow_l / E_flux_mag;
            dy *= arrow_l / E_flux_mag;

            fill(0, 0, 0);
            circle(screenX, screenY, 2);

            stroke(0, 0, 0);
            line(screenX, screenY, screenX + dx, screenY - dy);
            noStroke();
          }
        }
      }
    }

    for (let b = 0; b < antennas.length; b++) {
      antenna = antennas[b];

      segments = antenna.getSegments();
      currentSegments = antenna.getCurrentSegments();

      for (let j = 0; j < segments.length - 1; j++) {
        let x = segments[j][0];
        let x_next = segments[j + 1][0];

        let y = segments[j][1];
        let y_next = segments[j + 1][1];

        currentSegement = currentSegments[j];

        let currentAmp = Math.sqrt(
          currentSegement[0] * currentSegement[0] +
            currentSegement[1] * currentSegement[1]
        );

        const drive = antenna.getI0();
        const currentMag = Math.abs(currentAmp * (drive.a * ct - drive.b * st));
        const currentColorMag = 255 * squiz(currentMag, k8, 256);

        if (
          antenna.getSegFlags()[j] &&
          conScreenX(Math.max(x, x_next)) < width &&
          conScreenY(Math.max(y, y_next)) < height
        ) {
          thickLine(
            [conScreenX(x), conScreenY(y)],
            [conScreenX(x_next), conScreenY(y_next)],
            thickDipole * zoom,
            [currentColorMag, currentColorMag, 0]
          );
        }
      }
    }
    pop();
  }

  isMouseInStatBox = false;

  for (let b = 0; b < antennas.length; b++) {
    antenna = antennas[b];

    let xBox = conScreenX(antenna.getXBox());
    let yBox = conScreenY(antenna.getYBox());
    if (
      antenna.getShowBox() &&
      inRect(
        xBox,
        yBox,
        (windowWidth * 400) / 1400,
        (windowHeight * 370) / 1000,
        mouseX,
        mouseY
      )
    ) {
      isMouseInStatBox = true;
    }
  }

  for (let b = 0; b < antennas.length; b++) {
    antenna = antennas[b];

    let xBox = conScreenX(antenna.getXBox());
    let yBox = conScreenY(antenna.getYBox());

    if (antenna.getType() == "Dipole") {
      p1F[0] = conScreenX(antenna.getP1()[0]);
      p1F[1] = conScreenY(antenna.getP1()[1]);

      p2F[0] = conScreenX(antenna.getP2()[0]);
      p2F[1] = conScreenY(antenna.getP2()[1]);

      p3F[0] = conScreenX(antenna.getP3()[0]);
      p3F[1] = conScreenY(antenna.getP3()[1]);

      p4F[0] = conScreenX(antenna.getP4()[0]);
      p4F[1] = conScreenY(antenna.getP4()[1]);

      stroke(200, 0, 0);
      fill(200, 0, 0);

      if (antenna.getShowBox()) {
        stroke(250, 250, 250);
        fill(140, 140, 140);

        rect(
          xBox,
          yBox,
          (windowWidth * 400) / 1400,
          (windowHeight * 370) / 1000
        );

        stroke(0, 0, 0);
        fill(250, 250, 250);

        //amplitude slide

        rect(
          xBox + (windowWidth * 50) / 1400,
          yBox + (windowHeight * 55) / 1000,
          (windowWidth * 300) / 1400,
          (windowHeight * 30) / 1000
        );

        stroke(0, 0, 0);
        fill(0, 0, 0);
        push();
        scale(windowWidth / 1400, windowHeight / 1000);
        textSize(25);
        text(
          "Amplitude",
          (xBox * 1400) / windowWidth + 145,
          (yBox * 1000) / windowHeight + 48
        );
        pop();
        stroke(0, 0, 0);
        fill(110, 110, 110);

        //amplitude slide button

        if (
          inRect(
            xBox + (windowWidth * 50) / 1400 + amp_button_offset[b],
            yBox + (windowHeight * 55) / 1000,
            (windowWidth * 40) / 1400,
            (windowHeight * 30) / 1000,
            mouseX,
            mouseY
          ) &&
          mousePressedFlag
        ) {
          flags_amp_button[b] = true;
        }

        else if (
          !mousePressedFlag ||
          !inRect(
            xBox + (windowWidth * 50) / 1400,
            yBox + (windowHeight * 55) / 1000,
            (windowWidth * 300) / 1400,
            (windowHeight * 30) / 1000,
            mouseX,
            mouseY
          )
        ) {
          
          
          if(flags_amp_button[b]){
          flags_amp_button[b] = false;
           amp_slider_process = true;
          simulate = false;
           waitProcess = true
          processingScheduled = false
           A_map = new Map()
          EM_phase_amp_map = new Map()
            
           antenna.setAmp(
            (maxAmp * amp_button_offset[b]) /
              ((windowWidth * (300 - 40)) / 1400)
          );
          startProcessingNewSetup();

            
          }
        }

        if (flags_amp_button[b]) {
         
         
          fill(90, 90, 90);

          amp_button_offset[b] =
            mouseX - xBox - (windowWidth * (50 + 20)) / 1400;
          if (amp_button_offset[b] < 0) {
            amp_button_offset[b] = 0;
          } else if (
            amp_button_offset[b] + (windowWidth * 40) / 1400 >
            (windowWidth * 300) / 1400
          ) {
            amp_button_offset[b] = ((300 - 40) * windowWidth) / 1400;
          }
        }

        rect(
          xBox + (windowWidth * 50) / 1400 + amp_button_offset[b],
          yBox + (windowHeight * 55) / 1000,
          (windowWidth * 40) / 1400,
          (windowHeight * 30) / 1000
        );

        stroke(0, 0, 0);
        fill(250, 250, 250);
        
        

        //phase slide
        rect(
          xBox + (windowWidth * 50) / 1400,
          yBox + (windowHeight * 155) / 1000,
          (windowWidth * 300) / 1400,
          (windowHeight * 30) / 1000
        );

        stroke(0, 0, 0);
        fill(0, 0, 0);
        push();
        scale(windowWidth / 1400, windowHeight / 1000);
        textSize(25);
        text(
          "Phase",
          (1400 * xBox) / windowWidth + 170,
          (yBox * 1000) / windowHeight + 148
        );

        text(
          "0",
          (1400 * xBox) / windowWidth + 22,
          (yBox * 1000) / windowHeight + 178
        );
        text(
          "2π",
          (1400 * xBox) / windowWidth + 357,
          (yBox * 1000) / windowHeight + 178
        );
        pop();

        fill(110, 110, 110);

        if (
          inRect(
            xBox + (windowWidth * 50) / 1400 + phase_button_offset[b],
            yBox + (155 * windowHeight) / 1000,
            (windowWidth * 40) / 1400,
            (windowHeight * 30) / 1000,
            mouseX,
            mouseY
          ) &&
          mousePressedFlag
        ) {
          flags_phase_button[b] = true;
        }

        if (
          !mousePressedFlag ||
          !inRect(
            xBox + (windowWidth * 50) / 1400,
            yBox + (windowHeight * 155) / 1000,
            (windowWidth * 300) / 1400,
            (windowHeight * 30) / 1000,
            mouseX,
            mouseY
          )
        ) {
          
          if(flags_phase_button[b]){
          
          flags_phase_button[b] = false;
          phase_slider_process = true;
          simulate = false;
           waitProcess = true
          processingScheduled = false
           A_map = new Map()
          EM_phase_amp_map = new Map()
          antenna.setPhase(
            2 *
              Math.PI *
              (1 / 12) *
              Math.round(
                (12 * phase_button_offset[b]) /
                  ((windowWidth * (300 - 40)) / 1400)
              )
          );
          startProcessingNewSetup();  
            
          }
        }

        if (flags_phase_button[b]) {
          

          fill(90, 90, 90);

          phase_button_offset[b] =
            mouseX - xBox - (windowWidth * (50 + 20)) / 1400;
          if (phase_button_offset[b] < 0) {
            phase_button_offset[b] = 0;
          } else if (
            phase_button_offset[b] + (windowWidth * 40) / 1400 >
            (windowWidth * 300) / 1400
          ) {
            phase_button_offset[b] = (windowWidth * (300 - 40)) / 1400;
          }
        }

        rect(
          xBox + (windowWidth * 50) / 1400 + phase_button_offset[b],
          yBox + (windowHeight * 155) / 1000,
          (windowWidth * 40) / 1400,
          (windowHeight * 30) / 1000
        );

        stroke(0, 0, 0);
        fill(250, 250, 250);

        // seperation slide

        rect(
          xBox + (windowWidth * 50) / 1400,
          yBox + (windowHeight * 255) / 1000,
          (windowWidth * 300) / 1400,
          (windowHeight * 30) / 1000
        );

        push();
        stroke(0, 0, 0);
        fill(0, 0, 0);
        scale(windowWidth / 1400, windowHeight / 1000);
        textSize(25);
        text(
          "Seperation",
          (1400 * xBox) / windowWidth + 150,
          (1000 * yBox) / windowHeight + 248
        );
        pop();

        fill(110, 110, 110);

        if (
          inRect(
            xBox + (windowWidth * 50) / 1400 + sep_button_offset[b],
            yBox + (windowHeight * 255) / 1000,
            (windowWidth * 40) / 1400,
            (windowHeight * 30) / 1000,
            mouseX,
            mouseY
          ) &&
          mousePressedFlag
        ) {
          flags_sep_button[b] = true;
        }

        if (
          !mousePressedFlag ||
          !inRect(
            xBox + (windowWidth * 50) / 1400,
            yBox + (windowHeight * 255) / 1000,
            (windowWidth * 300) / 1400,
            (windowHeight * 30) / 1000,
            mouseX,
            mouseY
          )
        ) {
          
          if(flags_sep_button[b]){
          
          phase_slider_process = true;
          simulate = false;
           waitProcess = true
          processingScheduled = false
          A_map = new Map()
          EM_phase_amp_map = new Map()
          antenna.setSep(
            (maxSep * sep_button_offset[b]) /
              ((windowWidth * (300 - 40)) / 1400)
          );  
          flags_sep_button[b] = false;
            
          }
        }

        if (flags_sep_button[b]) {
          

          fill(90, 90, 90);

          sep_button_offset[b] =
            mouseX - xBox - (windowWidth * (50 + 20)) / 1400;
          if (sep_button_offset[b] < 0) {
            sep_button_offset[b] = 0;
          } else if (
            sep_button_offset[b] + (windowWidth * 40) / 1400 >
            (windowWidth * 300) / 1400
          ) {
            sep_button_offset[b] = (windowWidth * (300 - 40)) / 1400;
          }
        }

        rect(
          xBox + (windowWidth * 50) / 1400 + sep_button_offset[b],
          yBox + (windowHeight * 255) / 1000,
          (windowWidth * 40) / 1400,
          (windowHeight * 30) / 1000
        );

        //delete antenna button

        fill(240, 0, 0);

        if (
          inRect(
            xBox + (windowWidth * 128) / 1400,
            yBox + (windowHeight * 308) / 1000,
            (windowWidth * 150) / 1400,
            (windowHeight * 40) / 1000,
            mouseX,
            mouseY
          ) &&
          mousePressedFlag
        ) {
          amp_button_offset.splice(b, 1);
          phase_button_offset.splice(b, 1);
          sep_button_offset.splice(b, 1);
          flags_amp_button.splice(b, 1);
          flags_phase_button.splice(b, 1);
          flags_sep_button.splice(b, 1);
          delButtonPressed.splice(b, 1);
          antennas.splice(b, 1);
          mousePressedFlag = false;
          simulate = false;
           waitProcess = true
           A_map = new Map()
          EM_phase_amp_map = new Map()
        }

        rect(
          xBox + (windowWidth * 128) / 1400,
          yBox + (windowHeight * 308) / 1000,
          (windowWidth * 150) / 1400,
          (windowHeight * 40) / 1000
        );

        stroke(0, 0, 0);
        fill(0, 0, 0);
        push();
        scale(windowWidth / 1400, windowHeight / 1000);
        textSize(26);
        text(
          "delete",
          (xBox * 1400) / windowWidth + 170,
          (yBox * 1000) / windowHeight + 335
        );
        pop();
      }

      if (
        inRectGen(p1F, p2F, p3F, p4F, mouseX, mouseY) &&
        mousePressedFlag &&
        !isMouseInStatBox
      ) {
        mousePressedFlag = false;

        antenna.toggleShowBox();
      }
    }
  }
  

  
  
  push();
  noStroke(); fill(0, 0, 0, 180); rect(8, 8, 275, 25);
  fill(255); textSize(13); textAlign(LEFT, BASELINE);
  text('Zoom ' + Math.round(zoom * 100) + '%  |  wheel: zoom  |  0: reset', 15, 25);
  pop();
}

function mousePressed() {
  mouseRelease = false;
  mousePressedFlag = true;
}

function mouseReleased() {
  mouseRelease = true;
  mousePressedFlag = false;
}



// Mobile
function touchStarted() {
  mousePressedFlag = true;
  mouseRelease = false;
  return false; // prevent default scroll
}

function touchEnded() {
  mousePressedFlag = false;
  mouseRelease = true;
  return false;
}



