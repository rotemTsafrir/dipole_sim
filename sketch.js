
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
    this.segmentLengths = new Float64Array(0);
    this.segFlags = [];

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

  // Stable, normalized source data. Every rebuild replaces these arrays so that
  // geometry/current-profile edits invalidate the basis cache. Drive edits do not.
  // currentSegments[n] = Float64Array([JxRe, JyRe, JxIm, JyIm]).
  // segmentLengths is the integration weight; for a point current moment it is 1.
  getCurrentElements() {
    return {positions: this.segments, currents: this.currentSegments, lengths: this.segmentLengths};
  }

  getProperties() {
    return [
      {key:'amp', label:'Current amplitude', min:0, max:maxAmp, step:0.1, get:()=>this.amp, set:v=>this.setAmp(v)},
      {key:'phase', label:'Phase (°)', min:0, max:360, step:1, get:()=>this.phase*180/Math.PI, set:v=>this.setPhase(v*Math.PI/180)}
    ];
  }
  getReadouts() { return []; }
  getNotes() { return []; }
  getDisplayEdges() { return []; }
  getHitEdges() { return this.getDisplayEdges(); }
  getSelectionPoints() { return []; }

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

        let Ivec = new Float64Array(4);
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
    this.segmentLengths = new Float64Array(this.segments.length).fill(this.dl);
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

  getLength() {
    return this.length;
  }

  getType() {
    return "Dipole";
  }

  getProperties() {
    return [...super.getProperties(),
      {key:'gap', label:'Feed gap', min:0, max:Math.min(maxSep,Math.max(0,this.length-.01)), step:.01, get:()=>this.sep, set:v=>this.setSep(v)}];
  }
  getReadouts() { return [['Length',this.length.toFixed(3)],['Electrical length · L/λ',(this.length/this.wavelength).toFixed(3)]]; }
  getDisplayEdges() {
    const edges=[];
    for(let j=0;j<this.segments.length-1;j++) if(this.segFlags[j]) edges.push({a:screenPoint(this.segments[j]),b:screenPoint(this.segments[j+1]),current:j,width:thickDipole*zoom});
    return edges;
  }
  getHitEdges() { return [{a:screenPoint(this.endPointA),b:screenPoint(this.endPointB),width:thickDipole*zoom}]; }
  getSelectionPoints() { return [screenPoint(this.endPointA),screenPoint(this.endPointB)]; }

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

// Point current element: amp means |I*l|, not current through the display glyph.
// Weight 1 keeps this moment independent of wire discretization and zoom.
class HertzianDipole extends Antenna {
  constructor(wavelength,position,angle=0,moment=1,phase=0) {
    super(wavelength,moment,phase);
    this.position=[...position]; this.angle=angle;
    this.setCurrentSegments();
  }
  setCurrentSegments() {
    this.segments=[[...this.position]];
    this.currentSegments=[new Float64Array([Math.cos(this.angle),Math.sin(this.angle),0,0])];
    this.segmentLengths=new Float64Array([1]); this.segFlags=[true];
  }
  setAngle(angle) { this.angle=angle; this.setCurrentSegments(); }
  getProperties() {
    const props=super.getProperties(); props[0].label='Current moment |Iℓ|'; props[0].step=.01;
    return [...props,{key:'angle',label:'Orientation (°)',min:0,max:360,step:1,
      get:()=>((this.angle*180/Math.PI)%360+360)%360,set:v=>this.setAngle(v*Math.PI/180)}];
  }
  getReadouts() { return [['Source','Point current element'],['Orientation','From +x, counterclockwise']]; }
  getNotes() { return [{text:'The arrow is a fixed-size symbol, not a physical wire. Strength is set by Iℓ. Fields use the existing 0.1-unit source smoothing.'}]; }
  getDisplayEdges() {
    const p=screenPoint(this.position),dx=Math.cos(this.angle),dy=-Math.sin(this.angle);
    const start=[p[0]-18*dx,p[1]-18*dy],end=[p[0]+18*dx,p[1]+18*dy];
    return [
      {a:start,b:end,current:0,width:4},
      {a:end,b:[end[0]-8*dx+5*dy,end[1]-8*dy-5*dx],current:0,width:3},
      {a:end,b:[end[0]-8*dx-5*dy,end[1]-8*dy+5*dx],current:0,width:3}];
  }
  getSelectionPoints() { return [screenPoint(this.position)]; }
}

// Closed circular wire in the XY plane. Midpoint quadrature uses exact arc
// weights R*dphi, no duplicate endpoint and no feed gap. Positive current is CCW.
// At least 64 samples are used even when the entire loop is smaller than dl.
class SmallLoop extends Antenna {
  constructor(wavelength,center,radius,amp=10,phase=0) {
    super(wavelength,amp,phase);
    this.center=[...center]; this.radius=radius;
    this.setCurrentSegments();
  }
  setCurrentSegments() {
    const circumference=2*Math.PI*this.radius;
    const count=Math.max(64,Math.ceil(circumference/this.dl));
    const dphi=2*Math.PI/count;
    this.segments=[]; this.currentSegments=[]; this.segFlags=[];
    this.segmentLengths=new Float64Array(count).fill(circumference/count);
    for(let n=0;n<count;n++) {
      const phi=(n+.5)*dphi;
      this.segments.push([this.center[0]+this.radius*Math.cos(phi),this.center[1]+this.radius*Math.sin(phi)]);
      this.currentSegments.push(new Float64Array([-Math.sin(phi),Math.cos(phi),0,0]));
      this.segFlags.push(true);
    }
  }
  setRadius(radius) { if(!Number.isFinite(radius)||radius<=0) return; this.radius=radius; this.setCurrentSegments(); }
  setDl(dl) { if(this.dl===dl) return; this.dl=dl; this.setCurrentSegments(); }
  getProperties() {
    const props=super.getProperties();
    props[0].max=10000; props[0].scale='log';
    return [...props,{key:'radius',label:'Loop radius',min:.001,max:2,step:.001,get:()=>this.radius,set:v=>this.setRadius(v)}];
  }
  getReadouts() {
    const circumference=2*Math.PI*this.radius;
    return [['Circumference',circumference.toFixed(4)],['Electrical size · C/λ',(circumference/this.wavelength).toFixed(4)],['Magnetic moment |I·area|',(this.amp*Math.PI*this.radius*this.radius).toPrecision(4)],['Current profile','Uniform · positive CCW']];
  }
  getNotes() {
    const notes=[{text:'Current slider uses a logarithmic scale from 0 to 10,000. Small loops may need much higher current for a visible field; the numeric value is the actual source current.'},{text:'Loop lies in the XY plane (normal +z). Current is prescribed uniformly around the closed wire; there is no feed gap.'}];
    if(2*Math.PI*this.radius/this.wavelength>.1+1e-12) notes.push({warning:true,text:'C/λ exceeds 0.1. The uniform current is still prescribed, but is outside the small-loop approximation; no real feed response is being solved.'});
    if(this.radius*Scale*zoom<14) notes.push({text:'The ring symbol is enlarged for visibility. Only the radius value sets the field geometry.'});
    if(this.radius<Math.max(.1,sLength/Scale)) notes.push({text:'This loop is smaller than the smoothing/grid scale. Its near-source field is approximate; inspect the field away from the ring.'});
    return notes;
  }
  getDisplayEdges() {
    const p=screenPoint(this.center),r=Math.max(14,this.radius*Scale*zoom),edges=[];
    for(let n=0;n<64;n++) {
      const a=n*2*Math.PI/64,b=(n+1)*2*Math.PI/64;
      edges.push({a:[p[0]+r*Math.cos(a),p[1]-r*Math.sin(a)],b:[p[0]+r*Math.cos(b),p[1]-r*Math.sin(b)],current:Math.min(this.currentSegments.length-1,Math.floor(n*this.currentSegments.length/64)),width:3});
    }
    return edges;
  }
  getSelectionPoints() {
    const p=screenPoint(this.center),r=Math.max(14,this.radius*Scale*zoom);
    return [[p[0]+r,p[1]],[p[0]-r,p[1]]];
  }
}

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
  if (!Number.isFinite(next) || next === zoom || placementStart) return;
  // Preserve the exact world coordinate beneath the cursor (without snapping).
  const ratio = next / zoom;
  orig[0] = screenX - (screenX - orig[0]) * ratio;
  orig[1] = screenY - (screenY - orig[1]) * ratio;
  zoom = next;
  waitProcess = true;
  processingScheduled = false;
  zoomRebuildAfter = millis() + 100;
}

let width = 1000;
let height = 800;

let N = Math.ceil(height / sLength);
let M = Math.ceil(width / sLength);

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

let waitProcess = true;
let simulate = false;
let pause = false;

//field GUI states
let show_BField = false;
let show_EField = true;
let show_EnergyFlux = false;

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

let processingScheduled = false;

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
    const {positions: segments, currents, lengths} = antenna.getCurrentElements();
    const drive = antenna.getI0();
    const driveRe = drive.a, driveIm = drive.b;
    for (let j = 0; j < segments.length; j++) {
      const dx = myX - segments[j][0];
      const dy = myY - segments[j][1];
      const distance = Math.sqrt(dx * dx + dy * dy) + 0.1;
      const invR = 1 / distance;
      const gRe = Math.cos(-k * distance) * invR;
      const gIm = Math.sin(-k * distance) * invR;
      const jx = currents[j][0] * lengths[j], jy = currents[j][1] * lengths[j];
      const ix = currents[j][2] * lengths[j], iy = currents[j][3] * lengths[j];
      const jxRe = jx * driveRe - ix * driveIm, jxIm = jx * driveIm + ix * driveRe;
      const jyRe = jy * driveRe - iy * driveIm, jyIm = jy * driveIm + iy * driveRe;
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
  const {positions: segments, currents, lengths} = antenna.getCurrentElements();
  const step = sLength / Scale;
  const reusable = old && old.segments === segments &&
    old.currents === currents && old.lengths === lengths &&
    old.k === k && old.step === step;
  if (reusable && old.cols === cols && old.rows === rows && old.x0 === x0 && old.y0 === y0) {
    stats.reusedSamples += cols * rows;
    return old;
  }
  const size = cols * rows;
  const grid = { cols, rows, x0, y0, step, k, dl: antenna.dl,
    segments, currents, lengths,
    axRe: new Float64Array(size), axIm: new Float64Array(size),
    ayRe: new Float64Array(size), ayIm: new Float64Array(size) };
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
        const jx = currents[n][0] * lengths[n], jy = currents[n][1] * lengths[n];
        const ix = currents[n][2] * lengths[n], iy = currents[n][3] * lengths[n];
        xr += jx * gr - ix * gi; xi += jx * gi + ix * gr;
        yr += jy * gr - iy * gi; yi += jy * gi + iy * gr;
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
    antenna: a, segments: a.segments, currents: a.currentSegments, lengths: a.segmentLengths, dl: a.dl,
    re: a.I0.a, im: a.I0.b
  }))};
  const prev = fieldCacheSignature;
  const same = prev && prev.freq === freq && prev.c === c && prev.sLength === sLength &&
    prev.antennas.length === antennas.length && signature.antennas.every((a,i) => {
      const b = prev.antennas[i];
      return a.antenna === b.antenna && a.segments === b.segments &&
        a.currents === b.currents && a.lengths === b.lengths && a.dl === b.dl && a.re === b.re && a.im === b.im;
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

function renderFields() {
  if (simulate) {

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

  
  
  
  
  // Original field vectors and animated wire currents.

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

    pop();
  }

}
// Component palette and inspector share model-supplied geometry/properties.
// Self-contained drop-in sketch: native HTML/CSS are mounted by setup().
let selectedComponent = null;
let activeTool = 'select';
let placementType = null;
let placementStart = null;
let pointerGesture = null;
let hoverPoint = null;
let componentSerial = 0;
let ui = {};

// Descriptors own placement; the solver sees only complex current elements.
const componentTypes = {
  dipole: {
    name:'Center-fed dipole', description:'Two endpoints · sinusoidal current',
    firstHint:'Click the first endpoint', secondHint:'Click the second endpoint',
    validate:(a,b)=>length2D(a,b)>defaultDipoleSep+.01,
    invalidHint:'Choose a second endpoint farther than the feed gap',
    create:(a,b)=>configureSource(new Dipole(c/freq,a,b,defAmp,0,thickDipole/Scale),'dipole','Dipole')
  },
  hertzian: {
    name:'Hertzian dipole', description:'Position + direction · ideal current element',
    firstHint:'Click the source position', secondHint:'Click to set direction (symbol length is not physical)',
    validate:(a,b)=>length2D(a,b)>.001,
    invalidHint:'Choose a different point to set the direction',
    create:(a,b)=>configureSource(new HertzianDipole(c/freq,a,Math.atan2(b[1]-a[1],b[0]-a[0]),10,0),'hertzian','Hertzian')
  },
  smallLoop: {
    name:'Small loop', description:'Click a center · uniform circulating current', oneClick:true,
    firstHint:'Click the loop center; adjust radius in Properties',
    create:a=>configureSource(new SmallLoop(c/freq,a,.08*(c/freq)/(2*Math.PI),defAmp*100,0),'smallLoop','Loop')
  }
};
function configureSource(source,type,name) {
  source.setDl(resolution===1?dl_lr:resolution===2?dl_mr:dl_hr);
  source.componentType=type; source.label=`${name} ${++componentSerial}`;
  return source;
}

const interfaceCSS = `
html,body {margin:0!important;padding:0!important;width:100%;height:100%;overflow:hidden;background:#090d13;}
#em-app {--panel:#121923;--line:#283340;--muted:#93a3b6;--accent:#77e2c3;position:fixed;inset:0;z-index:10;display:grid;grid-template-rows:auto minmax(0,1fr) auto;color:#e6edf5;background:#090d13;font:13px/1.45 system-ui,-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;color-scheme:dark;}
#em-app * {box-sizing:border-box;}
#em-app [hidden] {display:none!important;}
#em-app button,#em-app input,#em-app select {font:inherit;}
#em-app button,#em-app select,#em-app input[type=number] {color:inherit;background:#1a2431;border:1px solid #344252;border-radius:7px;min-height:34px;padding:6px 10px;}
#em-app button {cursor:pointer;white-space:nowrap;}
#em-app button:hover {background:#293849;border-color:#62798f;}
#em-app button[aria-pressed=true],#em-app .primary {color:#0b2520;background:var(--accent);border-color:var(--accent);}
#em-app button:focus-visible,#em-app input:focus-visible,#em-app select:focus-visible {outline:2px solid var(--accent);outline-offset:3px;}
#em-app button:disabled {opacity:.4;cursor:default;}
#em-app .topbar {display:flex;align-items:center;flex-wrap:wrap;gap:12px;padding:12px 18px;background:var(--panel);border-bottom:1px solid var(--line);}
#em-app .brand {display:flex;align-items:center;gap:10px;margin-right:auto;min-width:170px;}
#em-app .brand-mark {color:var(--accent);font-size:27px;line-height:1;}
#em-app .brand strong {display:block;font-size:16px;letter-spacing:-.3px;}
#em-app .brand small {display:block;color:var(--muted);font-size:10px;letter-spacing:1.6px;text-transform:uppercase;}
#em-app .top-control {display:flex;align-items:center;gap:7px;color:var(--muted);font-size:12px;}
#em-app .top-control input {width:76px;}
#em-app .body {display:grid;grid-template-columns:66px minmax(0,1fr) 280px;min-height:0;}
#em-app .tools {display:flex;flex-direction:column;align-items:stretch;gap:9px;padding:14px 7px;background:var(--panel);border-right:1px solid var(--line);}
#em-app .tool {display:flex;flex-direction:column;align-items:center;gap:2px;padding:7px 2px;font-size:10px;}
#em-app .tool svg {width:21px;height:21px;fill:none;stroke:currentColor;stroke-width:1.7;stroke-linecap:round;stroke-linejoin:round;}
#em-app .divider {height:1px;background:var(--line);margin:5px 0;}
#em-app .workspace {position:relative;min-width:0;min-height:0;overflow:hidden;background:#000;}
#em-app .workspace canvas {display:block;touch-action:none;outline:none;}
#em-app .workspace canvas:focus-visible {outline:1px solid var(--accent);outline-offset:-2px;}
#em-app .canvas-tag {position:absolute;left:18px;top:15px;pointer-events:none;color:#c4d2dc;font-size:10px;text-transform:uppercase;letter-spacing:1.4px;background:#091019c9;padding:5px 8px;border-radius:4px;}
#em-app .hint {position:absolute;bottom:16px;left:50%;transform:translateX(-50%);max-width:90%;padding:8px 13px;border:1px solid #334354;border-radius:8px;background:#111b27e8;color:#c2d0de;font-size:12px;text-align:center;pointer-events:none;}
#em-app .palette {position:absolute;top:14px;left:14px;width:min(280px,calc(100% - 28px));z-index:2;background:#151e2a;border:1px solid #405368;border-radius:11px;padding:14px;box-shadow:0 12px 40px #0008;}
#em-app .eyebrow {text-transform:uppercase;letter-spacing:1.5px;font-size:10px;font-weight:650;color:var(--muted);margin:0 0 12px;}
#em-app .palette button {width:100%;text-align:left;padding:12px;white-space:normal;}
#em-app .palette button + button {margin-top:8px;}
#em-app .model-note {font-size:11px;line-height:1.6;color:var(--muted);margin:10px 0;}
#em-app .model-note.warning {color:#f5c781;}
#em-app .palette small {display:block;color:var(--muted);margin-top:4px;}
#em-app .inspector {min-height:0;overflow:auto;background:var(--panel);border-left:1px solid var(--line);padding:20px 18px;}
#em-app h2 {font-size:18px;margin:0 0 4px;letter-spacing:-.4px;}
#em-app .muted {color:var(--muted);font-size:12px;}
#em-app .section {padding-top:18px;margin-top:18px;border-top:1px solid var(--line);}
#em-app .property {display:block;margin:0 0 16px;}
#em-app .property-head {display:flex;align-items:center;justify-content:space-between;gap:10px;margin-bottom:8px;}
#em-app .property-head input {width:90px;text-align:right;font-variant-numeric:tabular-nums;}
#em-app input[type=range] {display:block;width:100%;margin:0;accent-color:var(--accent);cursor:pointer;}
#em-app .readout {display:flex;justify-content:space-between;color:var(--muted);margin:8px 0;font-size:12px;}
#em-app .readout output {color:#e6edf5;font-variant-numeric:tabular-nums;}
#em-app .danger {color:#ffaba7;background:transparent;border-color:#654046;width:100%;}
#em-app .scene-list {display:grid;gap:6px;margin-top:10px;}
#em-app .scene-list button {text-align:left;font-size:12px;overflow:hidden;text-overflow:ellipsis;}
#em-app .footer {display:flex;align-items:center;flex-wrap:wrap;gap:8px 18px;padding:8px 16px;color:var(--muted);background:var(--panel);border-top:1px solid var(--line);font-size:11px;}
#em-app .footer .status {margin-right:auto;}
#em-app .status-dot {display:inline-block;width:6px;height:6px;border-radius:50%;background:var(--accent);margin-right:7px;}
#em-app .zoom-controls {display:flex;gap:5px;align-items:center;}
#em-app .zoom-controls button {min-height:26px;padding:2px 8px;font-size:11px;}
#em-app .empty {padding:18px 0;line-height:1.7;color:var(--muted);}
@media(max-width:1000px) {#em-app .topbar {gap:8px;padding:9px 12px;}#em-app .body {grid-template-columns:58px minmax(0,1fr) 240px;}#em-app .brand {min-width:145px;}#em-app .brand strong{font-size:14px;}#em-app .top-control{font-size:11px;}}
@media(max-width:680px) {#em-app .body{grid-template-columns:52px minmax(0,1fr);grid-template-rows:minmax(180px,1fr) minmax(120px,34%);}#em-app .tools{grid-row:1 / 3;padding:10px 4px;}#em-app .inspector{grid-column:2;grid-row:2;border-left:0;border-top:1px solid var(--line);padding:14px;}#em-app .topbar{gap:7px;}#em-app .brand{min-width:135px;}#em-app .brand small{display:none;}#em-app .top-control input{width:65px;}#em-app .top-control select{max-width:125px;}#em-app .footer{padding:6px 10px;gap:6px 10px;}#em-app .hint{font-size:11px;bottom:9px;}#em-app .property{margin-bottom:12px;}}
`;

function icon(paths) { return `<svg viewBox="0 0 24 24" aria-hidden="true">${paths}</svg>`; }
function mountInterface() {
  const style = document.createElement('style');
  style.textContent = interfaceCSS;
  document.head.appendChild(style);
  const root = document.createElement('main');
  root.id = 'em-app';
  root.innerHTML = `
    <header class="topbar">
      <div class="brand"><span class="brand-mark" aria-hidden="true">∿</span><div><strong>EM Simulator</strong><small>Electromagnetic workspace</small></div></div>
      <button id="em-play" title="Pause / resume (Space)">Ⅱ Pause</button>
      <label class="top-control">Field <select id="em-field"><option value="E">Electric field</option><option value="B">Magnetic field</option><option value="S">Energy flux proxy</option></select></label>
      <label class="top-control">Frequency <input id="em-frequency" type="number" min="${minFreq}" max="${maxFreq}" step="0.01" value="${freq}" title="Simulation frequency"></label>
      <label class="top-control">Wave speed <input id="em-speed" type="number" min="${minSpeed}" max="${maxSpeed}" step="0.1" value="${c}" title="Simulation wave speed"></label>
      <label class="top-control">Quality <select id="em-quality"><option value="1">Low</option><option value="2" selected>Medium</option><option value="3">High</option></select></label>
    </header>
    <div class="body">
      <nav class="tools" aria-label="Workspace tools">
        <button class="tool" id="em-select" title="Select component (V)" aria-pressed="true">${icon('<path d="M5 3l14 9-7 1-3 7z"/>')}Select</button>
        <button class="tool" id="em-pan" title="Pan (H); middle-drag also pans" aria-pressed="false">${icon('<path d="M12 3v18M3 12h18M9 6l3-3 3 3M9 18l3 3 3-3M6 9l-3 3 3 3M18 9l3 3-3 3"/>')}Pan</button>
        <div class="divider"></div>
        <button class="tool" id="em-add" title="Add component (A)" aria-expanded="false" aria-controls="em-palette">${icon('<path d="M12 5v14M5 12h14"/>')}Add</button>
      </nav>
      <section class="workspace" id="em-workspace" aria-label="Simulation workspace">
        <div class="canvas-tag" id="em-field-tag">Electric field · XY plane</div>
        <div class="palette" id="em-palette" hidden><p class="eyebrow">Add component · Antennas</p><div id="em-palette-items"></div></div>
        <div class="hint" id="em-hint"></div>
      </section>
      <aside class="inspector" aria-label="Component properties">
        <p class="eyebrow">Properties</p>
        <div id="em-properties"></div>
        <div class="section"><p class="eyebrow">Scene <span id="em-count"></span></p><div class="scene-list" id="em-scene"></div></div>
        <div class="section"><button class="danger" id="em-clear">Clear scene</button></div>
        <div class="section muted">Wheel to zoom · 0 to reset<br>V select · H pan · A add<br>Esc cancel · Delete selected<br><br>Values use the original simulation units.</div>
      </aside>
    </div>
    <footer class="footer"><span class="status" role="status" id="em-status"></span><span id="em-grid"></span><span id="em-wavelength"></span><div class="zoom-controls"><button id="em-zoom-out" aria-label="Zoom out">−</button><button id="em-zoom-reset" title="Reset zoom (0)">100%</button><button id="em-zoom-in" aria-label="Zoom in">+</button></div></footer>`;
  document.body.appendChild(root);
  ui.root = root;
  const find = id => root.querySelector(`#em-${id}`);
  for (const id of ['workspace','properties','scene','palette','hint','status','count','grid','wavelength']) ui[id] = find(id);
  ui.play = find('play');
  ui.select = find('select'); ui.pan = find('pan'); ui.add = find('add');
  ui.select.onclick = () => setTool('select');
  ui.pan.onclick = () => setTool('pan');
  ui.add.onclick = togglePalette;
  for (const [type, descriptor] of Object.entries(componentTypes)) {
    const button = document.createElement('button');
    button.innerHTML = `<strong>${descriptor.name}</strong><small>${descriptor.description}</small>`;
    button.onclick = () => { placementType = type; setTool('add'); ui.canvas.focus(); };
    find('palette-items').appendChild(button);
  }
  ui.play.onclick = togglePause;
  find('field').onchange = e => {
    show_EField = e.target.value === 'E'; show_BField = e.target.value === 'B'; show_EnergyFlux = e.target.value === 'S';
    find('field-tag').textContent = e.target.selectedOptions[0].textContent + ' · XY plane';
  };
  function bindNumber(id, read, write) {
    const el = find(id);
    el.onchange = () => {
      const value = el.valueAsNumber;
      if (!Number.isFinite(value)) { el.value = read(); return; }
      const next = Math.max(+el.min, Math.min(+el.max, value));
      write(next); el.value = read();
      for (const a of antennas) a.setWavelength(c / freq);
      renderInspector(); requestFieldUpdate();
    };
  }
  bindNumber('frequency', () => freq, v => { freq = v; });
  bindNumber('speed', () => c, v => { c = v; });
  find('quality').onchange = e => {
    resolution = +e.target.value;
    sLength = resolution === 3 ? 2 : resolution === 2 ? 4 : 5;
    arrow_spacing = resolution === 3 ? 10 : resolution === 2 ? 5 : 4;
    renderInspector(); requestFieldUpdate();
  };
  find('clear').onclick = () => {
    antennas = []; selectedComponent = null; setTool('select');
    renderInspector(); renderScene(); requestFieldUpdate();
  };
  find('zoom-out').onclick = () => setViewZoom(zoom / 1.25, width / 2, height / 2);
  find('zoom-in').onclick = () => setViewZoom(zoom * 1.25, width / 2, height / 2);
  find('zoom-reset').onclick = () => setViewZoom(1, width / 2, height / 2);
  document.addEventListener('keydown', handleKeyboard);
  root.addEventListener('pointerdown', e => {
    if (!ui.palette.hidden && !ui.palette.contains(e.target) && !ui.add.contains(e.target)) closePalette();
  });
}

function requestFieldUpdate(delay = 0) {
  simulate = false;
  waitProcess = true;
  processingScheduled = false;
  zoomRebuildAfter = millis() + delay;
}
function closePalette() { ui.palette.hidden = true; ui.add.setAttribute('aria-expanded','false'); }
function togglePalette() {
  const open = ui.palette.hidden;
  setTool('select');
  ui.palette.hidden = !open;
  ui.add.setAttribute('aria-expanded', String(open));
  if (open) ui.palette.querySelector('button').focus();
}
function setTool(tool) {
  activeTool = tool; placementStart = null;
  closePalette();
  ui.select.setAttribute('aria-pressed',String(tool === 'select'));
  ui.pan.setAttribute('aria-pressed',String(tool === 'pan'));
  ui.add.setAttribute('aria-pressed',String(tool === 'add'));
  if (ui.canvas) ui.canvas.style.cursor = tool === 'pan' ? 'grab' : tool === 'add' ? 'crosshair' : 'default';
  updateHint();
}
function updateHint(message) {
  const descriptor=componentTypes[placementType];
  ui.hint.textContent = message || (activeTool==='add' && descriptor ? `${placementStart?descriptor.secondHint:descriptor.firstHint} · Esc to cancel` : activeTool==='pan' ? 'Drag to pan · Wheel to zoom' : 'Select a source to edit its properties · Add to place an antenna');
}
function selectComponent(component) {
  selectedComponent = component;
  renderInspector(); renderScene();
}
function renderScene() {
  ui.scene.replaceChildren();
  ui.count.textContent = `(${antennas.length})`;
  ui.root.querySelector('#em-clear').disabled = !antennas.length;
  for (const a of antennas) {
    const b = document.createElement('button');
    b.textContent = a.label;
    b.setAttribute('aria-pressed',String(a === selectedComponent));
    b.onclick = () => { setTool('select'); selectComponent(a); };
    ui.scene.appendChild(b);
  }
  if (!antennas.length) ui.scene.innerHTML = '<span class="muted">No components. Choose Add to begin.</span>';
}
function renderInspector() {
  const a = selectedComponent;
  ui.properties.replaceChildren();
  if (!a) { ui.properties.innerHTML = '<h2>No selection</h2><div class="empty">Select a component in the workspace or the scene list to edit its properties.</div>'; return; }
  const title = document.createElement('h2'); title.textContent = a.label;
  ui.properties.appendChild(title);
  const subtitle = document.createElement('div'); subtitle.className = 'muted'; subtitle.textContent = componentTypes[a.componentType].name;
  ui.properties.appendChild(subtitle);
  const section = document.createElement('div'); section.className = 'section'; ui.properties.appendChild(section);
  // A shared inspector consumes property descriptors; no per-source floating boxes.
  const properties = a.getProperties();
  for (const property of properties) {
    const row = document.createElement('div'); row.className = 'property';
    row.innerHTML = `<div class="property-head"><label for="em-prop-${property.key}">${property.label}</label><input id="em-prop-${property.key}" type="number" min="${property.min}" max="${property.max}" step="${property.step}"></div><input type="range" min="${property.min}" max="${property.max}" step="${property.step}" aria-label="${property.label} slider">`;
    const [number, range] = row.querySelectorAll('input');
    const logarithmic=property.scale==='log';
    if(logarithmic) { range.min=0; range.max=1000; range.step=1; range.setAttribute('aria-label',property.label+' logarithmic slider'); }
    const fromRange=()=>logarithmic?Math.expm1(range.valueAsNumber/1000*Math.log1p(property.max)):range.valueAsNumber;
    const syncInputs=value=>{
      number.value=Number(value.toFixed(4));
      range.value=logarithmic?1000*Math.log1p(value)/Math.log1p(property.max):value;
      range.setAttribute('aria-valuetext',String(Number(value.toFixed(4))));
    };
    syncInputs(property.get());
    range.oninput = () => { number.value=Number(fromRange().toFixed(4)); range.setAttribute('aria-valuetext',number.value); };
    const commit = input => {
      if (selectedComponent !== a || !antennas.includes(a)) return;
      let value = input===range?fromRange():input.valueAsNumber;
      if (!Number.isFinite(value)) { syncInputs(property.get()); return; }
      value = Math.max(property.min, Math.min(property.max,value));
      property.set(value); syncInputs(value);
      renderModelReadouts(a,readout);
      requestFieldUpdate();
    };
    number.onchange = () => commit(number); range.onchange = () => commit(range);
    section.appendChild(row);
  }
  const readout = document.createElement('div');
  renderModelReadouts(a,readout);
  section.appendChild(readout);
  const del = document.createElement('button'); del.className = 'danger'; del.textContent = 'Delete component';
  del.onclick = deleteSelected; section.appendChild(del);
}
function renderModelReadouts(source,container) {
  container.replaceChildren();
  for(const [label,value] of source.getReadouts()) {
    const row=document.createElement('div'); row.className='readout';
    const caption=document.createElement('span'); caption.textContent=label;
    const output=document.createElement('output'); output.textContent=value;
    row.append(caption,output); container.appendChild(row);
  }
  for(const note of source.getNotes()) {
    const el=document.createElement('p'); el.className='model-note'+(note.warning?' warning':''); el.textContent=note.text; container.appendChild(el);
  }
}
function deleteSelected() {
  if (!selectedComponent) return;
  antennas = antennas.filter(a => a !== selectedComponent);
  selectComponent(null); requestFieldUpdate();
}
function togglePause() { pause = !pause; time = millis()/timeScale; ui.play.textContent = pause ? '▶ Resume' : 'Ⅱ Pause'; }
function handleKeyboard(e) {
  if (/^(INPUT|SELECT|TEXTAREA)$/.test(e.target.tagName) || e.target.isContentEditable || e.ctrlKey || e.metaKey || e.altKey) return;
  switch (e.key.toLowerCase()) {
    case 'escape': setTool('select'); break;
    case 'v': setTool('select'); break;
    case 'h': setTool('pan'); break;
    case 'a': togglePalette(); break;
    case '0': setViewZoom(1,width/2,height/2); break;
    case 'delete': case 'backspace': deleteSelected(); break;
    case ' ': if (e.target.tagName === 'BUTTON') return; togglePause(); break;
    default: return;
  }
  e.preventDefault();
}
function pointerLocation(e) {
  const rect = ui.canvas.getBoundingClientRect();
  return [(e.clientX-rect.left)*width/rect.width,(e.clientY-rect.top)*height/rect.height];
}
function hitComponent(point) {
  for (let i=antennas.length-1;i>=0;i--) {
    const a=antennas[i];
    for(const edge of a.getHitEdges()) {
      const p=edge.a,q=edge.b,dx=q[0]-p[0],dy=q[1]-p[1],d=dx*dx+dy*dy;
      const t=d?Math.max(0,Math.min(1,((point[0]-p[0])*dx+(point[1]-p[1])*dy)/d)):0;
      if(Math.hypot(point[0]-p[0]-t*dx,point[1]-p[1]-t*dy)<=Math.max(9,edge.width/2+4)) return a;
    }
  }
  return null;
}
function bindCanvasEvents() {
  const canvas=ui.canvas;
  canvas.tabIndex=0; canvas.setAttribute('aria-label','Electromagnetic field. Choose a source from Add and follow the placement hint.');
  canvas.addEventListener('pointerdown',e=>{
    if (pointerGesture || (e.button!==0 && e.button!==1)) return;
    e.preventDefault(); canvas.focus(); closePalette();
    const point=pointerLocation(e);
    pointerGesture={id:e.pointerId, start:point, last:point, pan:activeTool==='pan'||e.button===1, moved:false};
    canvas.setPointerCapture(e.pointerId);
    if (pointerGesture.pan) canvas.style.cursor='grabbing';
  });
  canvas.addEventListener('pointermove',e=>{
    const point=pointerLocation(e); hoverPoint=point;
    if (!pointerGesture || pointerGesture.id!==e.pointerId) return;
    const g=pointerGesture;
    if (Math.hypot(point[0]-g.start[0],point[1]-g.start[1])>4) g.moved=true;
    if (g.pan) { orig[0]+=point[0]-g.last[0]; orig[1]+=point[1]-g.last[1]; waitProcess=true; processingScheduled=false; }
    g.last=point;
  });
  canvas.addEventListener('pointerup',e=>{
    if (!pointerGesture || pointerGesture.id!==e.pointerId) return;
    const g=pointerGesture, point=pointerLocation(e); pointerGesture=null;
    if (canvas.hasPointerCapture(e.pointerId)) canvas.releasePointerCapture(e.pointerId);
    canvas.style.cursor=activeTool==='pan'?'grab':activeTool==='add'?'crosshair':'default';
    if (g.pan) { requestFieldUpdate(); return; }
    if (g.moved || point[0]<0 || point[0]>width || point[1]<0 || point[1]>height) return;
    if (activeTool==='add') {
      const world=[conMyX(point[0]),conMyY(point[1])];
      const descriptor=componentTypes[placementType];
      if (!descriptor.oneClick && !placementStart) { placementStart=world; updateHint(); return; }
      if (!descriptor.oneClick && !descriptor.validate(placementStart,world)) { updateHint(descriptor.invalidHint); return; }
      const a=descriptor.oneClick?descriptor.create(world):descriptor.create(placementStart,world);
      antennas.push(a); setTool('select'); selectComponent(a); requestFieldUpdate();
    } else selectComponent(hitComponent(point));
  });
  const cancel=e=>{
    if (!pointerGesture || pointerGesture.id!==e.pointerId) return;
    const wasPan=pointerGesture.pan; pointerGesture=null;
    canvas.style.cursor=activeTool==='pan'?'grab':activeTool==='add'?'crosshair':'default';
    if (wasPan) requestFieldUpdate();
  };
  canvas.addEventListener('pointercancel',cancel); canvas.addEventListener('lostpointercapture',cancel);
  canvas.addEventListener('pointerleave',()=>{ if (!pointerGesture) hoverPoint=null; });
  canvas.addEventListener('wheel',e=>{
    e.preventDefault(); if (pointerGesture || placementStart) return;
    const point=pointerLocation(e);
    const delta=e.deltaY*(e.deltaMode===1?16:e.deltaMode===2?height:1);
    setViewZoom(zoom*Math.exp(-Math.max(-500,Math.min(500,delta))*.0015),...point);
  },{passive:false});
  canvas.addEventListener('contextmenu',e=>e.preventDefault());
}

function setup() {
  mountInterface();
  width=Math.max(1,ui.workspace.clientWidth); height=Math.max(1,ui.workspace.clientHeight);
  const canvas=createCanvas(width,height); canvas.parent(ui.workspace); ui.canvas=canvas.elt;
  pixelDensity(1); frameRate(60);
  orig=[width/2+.1,height/2+.1];
  // Same default source construction as before, now centered in the workspace.
  const a=componentTypes.dipole.create([conMyX(.05+width/2),conMyY(2*height/3)],[conMyX(width/2),conMyY(height/3)]);
  antennas.push(a); selectedComponent=a;
  bindCanvasEvents(); setTool('select'); renderInspector(); renderScene();
  time=millis()/timeScale;
}
function draw() {
  const w=Math.max(1,ui.workspace.clientWidth), h=Math.max(1,ui.workspace.clientHeight);
  if (w!==width || h!==height) {
    orig[0]+=(w-width)/2; orig[1]+=(h-height)/2;
    resizeCanvas(w,h); width=w; height=h; requestFieldUpdate(80);
  }
  const now=millis()/timeScale;
  if (!pause && simulate) timeSim+=now-time;
  time=now;
  background(0);
  renderFields();
  drawComponentOverlay();
  updateStatus();
  // Allow a painted "Updating" state before the synchronous, unchanged solver.
  if (waitProcess && !pointerGesture && millis()>=zoomRebuildAfter) {
    if (!processingScheduled) processingScheduled=true;
    else startProcessingNewSetup();
  }
}
function screenPoint(p) { return [conScreenX(p[0]),conScreenY(p[1])]; }
function currentColor(source,index) {
  if(!simulate) return [200,213,226];
  const j=source.currentSegments[index], drive=source.I0;
  const re=drive.a*frameCos-drive.b*frameSin, im=drive.a*frameSin+drive.b*frameCos;
  const magnitude=Math.hypot(j[0]*re-j[2]*im,j[1]*re-j[3]*im);
  const bright=255*squiz(magnitude,k8,256);
  return [bright,bright,0];
}
function drawComponentOverlay() {
  push();
  for(const a of antennas) {
    for(const edge of a.getDisplayEdges()) {
      if(length2D(edge.a,edge.b)>1e-10) thickLine(edge.a,edge.b,edge.width,currentColor(a,edge.current));
    }
  }
  if (selectedComponent) {
    stroke(119,226,195); strokeWeight(1.5); noFill();
    for (const p of selectedComponent.getSelectionPoints()) circle(p[0],p[1],10);
  }
  if(activeTool==='add' && placementStart && hoverPoint) {
    const end=[conMyX(hoverPoint[0]),conMyY(hoverPoint[1])];
    stroke(119,226,195); strokeWeight(2); drawingContext.setLineDash([6,5]);
    if(placementType==='hertzian') {
      const center=screenPoint(placementStart),angle=Math.atan2(end[1]-placementStart[1],end[0]-placementStart[0]);
      line(center[0]-18*Math.cos(angle),center[1]+18*Math.sin(angle),center[0]+18*Math.cos(angle),center[1]-18*Math.sin(angle));
    } else line(conScreenX(placementStart[0]),conScreenY(placementStart[1]),conScreenX(end[0]),conScreenY(end[1]));
    drawingContext.setLineDash([]); noFill(); circle(conScreenX(placementStart[0]),conScreenY(placementStart[1]),10);
  }
  pop();
}
function updateStatus() {
  const state=waitProcess?'Updating fields…':pause?'Paused':'Running';
  const status=`${state} · ${antennas.length} ${antennas.length===1?'source':'sources'}`;
  if(ui.status.dataset.text!==status) { ui.status.innerHTML='<span class="status-dot"></span>'+status; ui.status.dataset.text=status; }
  ui.grid.textContent=`Grid ${(sLength/Scale).toFixed(2)}`;
  ui.wavelength.textContent=`λ ${(c/freq).toFixed(2)}`;
  ui.root.querySelector('#em-zoom-reset').textContent=`${Math.round(zoom*100)}%`;
}
