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
    this.setCurrentSegments();
    this.I0 = ComplexNum.cis(this.phase).product(newAmp);
  }

  setPhase(newPhase) {
    this.phase = newPhase;
    this.I0 = ComplexNum.cis(newPhase).product(this.amp);
    this.setCurrentSegments();
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
const Scale = 50;

let width = 1000;
let height = 800;

let prevWidth = width;
let prevHeight = height;

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
  return (sLength/(Scale))*Math.round((x - orig[0])/sLength );
}

function conMyY(y) {
  return (sLength/(Scale))*Math.round((orig[1] - y)/sLength );
}

function conScreenX(x) {
  return x * Scale + orig[0];
}

function conScreenY(y) {
  return orig[1] - y * Scale;
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

function calcA(myX, myY) {
  let A = [];

  A[0] = zero;
  A[1] = zero;

  for (let i = 0; i < antennas.length; i++) {
    antenna = antennas[i];

    segments = antenna.getSegments();

    currents = antenna.getCurrentSegments();

    for (let j = 0; j < segments.length; j++) {
      let r_tag = [segments[j][0], segments[j][1]];
       
      let distance = length2D([myX,myY], r_tag);

      //for numerical stability near current sources
      distance = distance + 0.1;

      current[0] = new ComplexNum(currents[j][0] * antenna.getDl(), 0).product(
        antenna.getI0()
      );
      current[1] = new ComplexNum(currents[j][1] * antenna.getDl(), 0).product(
        antenna.getI0()
      );

      phase = ComplexNum.cis(-k * distance);

      factor = phase.product(1 / distance);

      A[0] = A[0].add(current[0].product(factor));
      A[1] = A[1].add(current[1].product(factor));
           
      
    }
  }
  
  return A
}

function startProcessingNewSetup(){
  for (let i = 0; i < antennas.length; i++) {
          antenna = antennas[i];

          if (resolution == 1) {
            antenna.setDl(dl_lr);
          } else if (resolution == 2) {
            antenna.setDl(dl_mr);
          } else {
            antenna.setDl(dl_hr);
          }
        }
    
    
    
    
     k = (2 * Math.PI * freq) / c;

    cFACTOR1 = new ComplexNum(0, -2 * Math.PI * freq);
    cFACTOR2 = new ComplexNum(0, -c / k);
    N = Math.ceil(height / sLength);
    M = Math.ceil(width / sLength);
    

    for (let y = -sLength; y < height+sLength; y += sLength) {
      for (let x = -sLength; x < width+sLength; x += sLength) {
      
        r[0] = conMyX(x);
        r[1] = conMyY(y);

        let key = `${r[0]},${r[1]}`;
        
             
        
        if(!A_map.get(key)){
           let A_vec = calcA(r[0],r[1]) 
        
        A_map.set(key,{x:A_vec[0], y:A_vec[1]})
          
        }
        
        
  
        
      }
    }
    
    
for (let i = 0; i + 1 < width / sLength; i++) {
  for (let j = 0; j + 1 < height / sLength; j++) {
   // Compute coordinates
let x_coord_c = conMyX(i * sLength);

let x_coord_l = conMyX((i - 1) * sLength);
let x_coord_r = conMyX((i + 1) * sLength);

let y_coord_c = conMyY(j * sLength);
let y_coord_u = conMyY((j - 1) * sLength);
let y_coord_d = conMyY((j + 1) * sLength);
    
    
    if(! EM_phase_amp_map.get(`${x_coord_c},${y_coord_c}`)){

// Neighbor keys (using corrected variables)
let key_center = `${x_coord_c},${y_coord_c}`;
let key_up = `${x_coord_c},${y_coord_u}`;
let key_down = `${x_coord_c},${y_coord_d}`;
let key_left = `${x_coord_l},${y_coord_c}`;
let key_right = `${x_coord_r},${y_coord_c}`;
let key_up_left = `${x_coord_l},${y_coord_u}`;
let key_up_right = `${x_coord_r},${y_coord_u}`;
let key_down_left = `${x_coord_l},${y_coord_d}`;
let key_down_right = `${x_coord_r},${y_coord_d}`;


    // Fetch from A_map
    let A_center = A_map.get(key_center);
    let A_up = A_map.get(key_up);
    let A_down = A_map.get(key_down);
    let A_left = A_map.get(key_left);
    let A_right = A_map.get(key_right);
    let A_up_left = A_map.get(key_up_left);
    let A_up_right = A_map.get(key_up_right);
    let A_down_left = A_map.get(key_down_left);
    let A_down_right = A_map.get(key_down_right);
    
    

    if (!A_center || !A_up || !A_down || !A_left || !A_right) continue; // Safety
      
    // --- Calculate E phasor ---
    let E_phasor = [];
    E_phasor[0] = A_center.x.product(cFACTOR1);
    E_phasor[1] = A_center.y.product(cFACTOR1);

    // Compute tempE1 (second derivatives along x and y)
    let tempE1 = [];
    tempE1[0] = A_right.x.add(A_left.x).add(A_center.x.product(-2));
    tempE1[1] = A_up.y.add(A_down.y).add(A_center.y.product(-2));

    tempE1[0] = tempE1[0].product(const_rFactor1);
    tempE1[1] = tempE1[1].product(const_rFactor1);

    // Compute tempE2 (cross derivatives)
    let tempE2 = [];
    if (A_up_right && A_down_left && A_down_right && A_up_left) {
      tempE2[0] = A_down_left.y.add(A_up_right.y).subtract(A_down_right.y.add(A_up_left.y));
      tempE2[1] = A_down_left.x.add(A_up_right.x).subtract(A_down_right.x.add(A_up_left.x));
    } else {
      tempE2[0] = new ComplexNum(0, 0);
      tempE2[1] = new ComplexNum(0, 0);
    }

    tempE2[0] = tempE2[0].product(0.25 * const_rFactor1);
    tempE2[1] = tempE2[1].product(0.25 * const_rFactor1);

    // Combine
    let tempE = [];
    tempE[0] = tempE1[0].add(tempE2[0]).product(cFACTOR2);
    tempE[1] = tempE1[1].add(tempE2[1]).product(cFACTOR2);

    E_phasor[0] = E_phasor[0].add(tempE[0]);
    E_phasor[1] = E_phasor[1].add(tempE[1]);

    // --- Calculate B phasor ---
    let tempB1 = A_up.x.subtract(A_down.x).product(const_rFactor2);
    let tempB2 = A_right.y.subtract(A_left.y).product(const_rFactor2);
    let B_phasor = tempB2.subtract(tempB1);

    // Store in EM_phase_amp_map
    EM_phase_amp_map.set(key_center, {
      Ex_amp: E_phasor[0].getAbsVal(),
      Ex_phase: E_phasor[0].getPhase(),
      Ey_amp: E_phasor[1].getAbsVal(),
      Ey_phase: E_phasor[1].getPhase(),
      B_amp: B_phasor.getAbsVal(),
      B_phase: B_phasor.getPhase()
    });
  }
}
}
    waitProcess = false;
    simulate = true;
  
  processingScheduled =  false
}


function setup() {
  createCanvas(windowWidth, windowHeight);
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

  background(0); // draw black background
  loadPixels(); // fill the pixel buffer from canvas

  updatePixels(); // commit black pixels to canvas

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
  
  
  if (resolution == 3) {
    frameRate(20);
  }


  if (windowWidth != prevWidth || windowHeight != prevHeight) {
    freq_button_offset = (freq_button_offset * windowWidth) / prevWidth;
    speed_button_offset = (speed_button_offset * windowWidth) / prevWidth;

    for (let b = 0; b < antennas.length; b++) {
      amp_button_offset[b] *= windowWidth / prevWidth;
      phase_button_offset[b] *= windowWidth / prevWidth;
      sep_button_offset[b] *= windowWidth / prevWidth;
    }

    if (simulate) {
      if (prevWidth < windowWidth || prevHeight < windowHeight) {
        simulate = false;
        waitProcess = true;
      }
    }

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

    timePhase = ComplexNum.cis(timeSim * w);

    let pixelIndex = 0;

    // Load the pixel buffer once at the start
    // Load the pixel buffer once at the start

    loadPixels();

    for (let i = 0; i  < width / sLength; i++) {
      for (let j = 0; j < height / sLength; j++) {
        let screenX = sLength * i;
        let screenY = sLength * j;

        // Convert to coordinate key
        let x_coord = conMyX(i * sLength);
        let y_coord = conMyY(j * sLength);
        let key = `${x_coord},${y_coord}`;

        // Get EM data from the map
        let data = EM_phase_amp_map.get(key);
        if (!data) continue; // Safety check

        // Extract precomputed amplitude and phase
        let Ex_amp = data.Ex_amp;
        let Ex_phase = data.Ex_phase;
        let Ey_amp = data.Ey_amp;
        let Ey_phase = data.Ey_phase;
        let B_amp = data.B_amp;
        let B_phase = data.B_phase;

        // Calculate real-time values for E and B
        let Ex_t = Ex_amp * Math.cos(Ex_phase + timeSim * w);
        let Ey_t = Ey_amp * Math.cos(Ey_phase + timeSim * w);
        let B_t = B_amp * Math.cos(B_phase + timeSim * w);

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

        // Draw block
        let maxX = width;
        let maxY = height;
        let endX = Math.min(screenX + sLength, maxX);
        let endY = Math.min(screenY + sLength, maxY);

        if (screenX < maxX && screenY < maxY) {
          for (let di = 0; di < sLength; di++) {
            for (let dj = 0; dj < sLength; dj++) {
              let pixelIndex =
                ((j * sLength + dj) * windowWidth + i * sLength + di) * 4;
              pixels[pixelIndex] = r * a;
              pixels[pixelIndex + 1] = g * a;
              pixels[pixelIndex + 2] = b * a;
            }
          }
        }
      }
    }

    updatePixels();
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
            thickDipole,
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

        thickLine([mouseX, mouseY], p1, thickDipole, [250, 250, 250]);

        let lenAntenna = Math.sqrt(
          (mouseX - p1[0]) * (mouseX - p1[0]) +
            (mouseY - p1[1]) * (mouseY - p1[1])
        );

        let tempX1 =
          0.5 * (mouseX + p1[0]) -
          (0.5 * (p1[0] - mouseX) * defaultDipoleSep * Scale) / lenAntenna;
        let tempX2 =
          0.5 * (mouseX + p1[0]) +
          (0.5 * (p1[0] - mouseX) * defaultDipoleSep * Scale) / lenAntenna;

        let tempY1 =
          0.5 * (mouseY + p1[1]) -
          (0.5 * (p1[1] - mouseY) * defaultDipoleSep * Scale) / lenAntenna;
        let tempY2 =
          0.5 * (mouseY + p1[1]) +
          (0.5 * (p1[1] - mouseY) * defaultDipoleSep * Scale) / lenAntenna;

        thickLine([tempX1, tempY1], [tempX2, tempY2], thickDipole, [0, 0, 0]);
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

  if (waitProcess) {
     
    
    if( !processingScheduled ){
   
     processingScheduled = true  
    }
    
    else{
      
       startProcessingNewSetup()
    }
    
    
  }


  
  //after processing

  if (simulate) {
    w = 2 * Math.PI * freq;

    timePhase = ComplexNum.cis(timeSim * w);

    for (let i = 0; i < width / sLength; i += arrow_spacing) {
      for (let j = 0; j < height / sLength; j += arrow_spacing) {
        let screenX = sLength * i;
        let screenY = sLength * j;

        // Convert to coordinate key
        let x_coord = conMyX(i * sLength);
        let y_coord = conMyY(j * sLength);
        let key = `${x_coord},${y_coord}`;

        // Retrieve data from the map
        let data = EM_phase_amp_map.get(key);
        if (!data) continue;

        let Ex_amp = data.Ex_amp;
        let Ex_phase = data.Ex_phase;
        let Ey_amp = data.Ey_amp;
        let Ey_phase = data.Ey_phase;
        let B_amp = data.B_amp;
        let B_phase = data.B_phase;

        // Compute real-time values
        let Ex_t = Ex_amp * Math.cos(Ex_phase + timeSim * w);
        let Ey_t = Ey_amp * Math.cos(Ey_phase + timeSim * w);
        let B_t = B_amp * Math.cos(B_phase + timeSim * w);

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

        currentPhasor = antenna.getI0().product(currentAmp);

        let currentMag = Math.abs(currentPhasor.product(timePhase).getA());

        let currentColorMag = 255 * squiz(currentMag, k8, (levels = 256));

        if (
          antenna.getSegFlags()[j] &&
          conScreenX(Math.max(x, x_next)) < width &&
          conScreenY(Math.max(y, y_next)) < height
        ) {
          thickLine(
            [conScreenX(x), conScreenY(y)],
            [conScreenX(x_next), conScreenY(y_next)],
            thickDipole,
            [currentColorMag, currentColorMag, 0]
          );
        }
      }
    }
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
  

  
  
}

function mousePressed() {
  mouseRelease = false;
  mousePressedFlag = true;
}

function mouseReleased() {
  mouseRelease = true;
  mousePressedFlag = false;
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
