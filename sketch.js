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

      
    this.wirePoints = []
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
    let len = Math.fround(
      Math.sqrt(tempDel[0] * tempDel[0] + tempDel[1] * tempDel[1])
    );

    tempDel[0] *= Math.fround(this.thickness / len);
    tempDel[1] *= Math.fround(this.thickness / len);

    this.P1[0] = Math.fround(this.endPointA[0] + 0.5 * tempDel[0]);
    this.P1[1] = Math.fround(this.endPointA[1] + 0.5 * tempDel[1]);

    this.P2[0] = Math.fround(this.endPointB[0] + 0.5 * tempDel[0]);
    this.P2[1] = Math.fround(this.endPointB[1] + 0.5 * tempDel[1]);

    this.P3[0] = Math.fround(this.endPointB[0] - 0.5 * tempDel[0]);
    this.P3[1] = Math.fround(this.endPointB[1] - 0.5 * tempDel[1]);

    this.P4[0] = Math.fround(this.endPointA[0] - 0.5 * tempDel[0]);
    this.P4[1] = Math.fround(this.endPointA[1] - 0.5 * tempDel[1]);

    let maxX = Math.fround(Math.max(this.P1[0], this.P2[0]));
    let Y = this.P4[1];

    maxX = Math.fround(Math.max(maxX, this.P3[0]));
    maxX = Math.fround(Math.max(maxX, this.P4[0]));

    if (maxX === this.P1[0]) {
      Y = this.P1[1];
    } else if (maxX === this.P2[0]) {
      Y = this.P2[1];
    } else if (maxX === this.P3[0]) {
      Y = this.P3[1];
    }

    this.xBox = Math.fround(maxX + this.delBoxAntenna);
    this.yBox = Math.fround(Y - this.delBoxAntenna);
  }

  getSep() {
    return this.sep;
  }

  dist2D(p1, p2) {
    let temp = Math.fround(
      (p2[0] - p1[0]) * (p2[0] - p1[0]) + (p2[1] - p1[1]) * (p2[1] - p1[1])
    );
    return Math.fround(Math.sqrt(temp));
  }

  setCurrentSegments() {
    this.segments = [];
    this.currentSegments = [];
    this.segFlags = []

    let l = 0;
    let sideLen;
    let pA;
    let pB;
    let delP = [0, 0];
    let p = [0, 0];
    let k = Math.fround((2 * Math.PI) / this.wavelength);

    for (let i = 0; i < this.wirePoints.length - 1; i++) {
      pA = this.wirePoints[i];
      pB = this.wirePoints[i + 1];

      sideLen = this.dist2D(pA, pB);

      delP[0] = pB[0] - pA[0];
      delP[1] = pB[1] - pA[1];

      l = 0;

      while (l <= sideLen) {
        p[0] = Math.fround(pA[0] + (l / sideLen) * delP[0]);
        p[1] = Math.fround(pA[1] + (l / sideLen) * delP[1]);

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
        l = Math.fround(Math.round(10000 * l) / 10000.0);
      }
    }
  }

  setSep(newSep) {
    this.sep = newSep;
    this.setCurrentSegments();
    
  }
  
   setDl(new_dl){
    
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
    this.absVal = Math.fround(Math.sqrt(a * a + b * b));
    if (a > 0 && b > 0) this.phase = Math.fround(Math.atan(b / a));
    else if (a > 0 && b < 0) this.phase = Math.fround(Math.atan(b / a));
    else if (a == 0 && b > 0) this.phase = Math.fround(Math.PI / 2);
    else if (a == 0 && b < 0) this.phase = -Math.fround(Math.PI / 2);
    else if (a > 0 && b == 0) this.phase = 0;
    else if (a < 0 && b == 0) this.phase = Math.fround(Math.PI);
    else if (a < 0 && b > 0)
      this.phase = Math.fround(Math.atan(b / a) + Math.PI);
    else this.phase = Math.fround(Math.atan(b / a) - Math.PI);
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
    return new ComplexNum(
      Math.fround(Math.cos(phase)),
      Math.fround(Math.sin(phase))
    );
  }
}
////

//grid square side length
let sLength = 5;
const dl_hr = 0.04;
const dl_mr = 0.1
const dl_lr = 0.2

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




const orig = [500.1, 400.1];

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
let B_phasor = [];

//electric field phasor
let E_phase_amp = [];

//magnetic field phasor
let B_phase_amp =  [];

//electric field at each simulated pixel
let E =  [];

//magnetic field  at each simulated pixel
let B =  [];

let A = [];

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


//maximum current amplitude
const maxAmp = 25;

//default current magnitude
const defAmp = 10;


//
const maxSep = 0.4;
const defaultDipoleSep = 0.2

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

let lastMouseXB = 0;
let lastMouseYB = 0;

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

let p1F = new Float32Array(2);
let p2F = new Float32Array(2);
let p3F = new Float32Array(2);
let p4F = new Float32Array(2);

let lastMouseX;
let lastMouseY;

let const_rFactor1 = Math.fround((Scale * Scale) / (sLength * sLength));
let const_rFactor2 = Math.fround(Scale / (2 * sLength));

function squiz(l, k, levels = 256) {
  
  if(l<=0){
    
    return 0;
  }
  
  // Step 1: compute the squiz output (same as original)
  let ans = k*l / (1 + k*l)

  return ans;
}

function length2D(p1, p2) {
  return Math.fround(
    Math.sqrt(
      (p1[0] - p2[0]) * (p1[0] - p2[0]) + (p1[1] - p2[1]) * (p1[1] - p2[1])
    )
  );
}

function conMyX(x) {
  return (x - orig[0]) / Scale;
}

function conMyY(y) {
  return (orig[1] - y) / Scale;
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

function setup() {
  createCanvas(windowWidth, windowHeight);
  
  width = Math.round((1/1.4)*windowWidth)
  height = Math.round((8/10)*windowHeight)

  

  background(0, 0, 0);

  pixelDensity(1);

  current = [];

  A = [];

  one_point = false;

  freq_button_offset = ((freq - minFreq) / (maxFreq - minFreq)) * ((200 - 30)/1400)*windowWidth;
  speed_button_offset = ((c - minSpeed) / (maxSpeed - minSpeed)) * ((200 - 30)/1400)*windowWidth;

  for (let i = 0; i < width / sLength; i++) {
    for (let j = 0; j < height / sLength; j++) {
      let tempInd = M * j + i;

      E[2 * tempInd] = 0;
      E[2 * tempInd + 1] = 0;
      B[tempInd] = 0;
    }
  }

  background(0); // draw black background
  loadPixels(); // fill the pixel buffer from canvas

  updatePixels(); // commit black pixels to canvas
  
  
   let myP1 = [conMyX(0.05+width/2), conMyY(2*height/3)];
    let myP2 = [conMyX(width/2), conMyY(height/3)];

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
          
          if(resolution==1){
            
            dipole.setDl(dl_lr)
            
          }
          
          
          else if(resolution==2){
            
            dipole.setDl(dl_mr)
            
          }
          
          
          else {
            
            dipole.setDl(dl_hr)
            
          }
          

          antennas.push(dipole);

          amp_button_offset.push((defAmp / maxAmp) * ((300 - 40)/1400)*windowWidth);
          phase_button_offset.push(0);
          sep_button_offset.push(defaultDipoleSep*((300 - 40)/1400)*windowWidth/maxSep);

          flags_amp_button.push(false);
          flags_phase_button.push(false);
          flags_sep_button.push(false);
          delButtonPressed.push(false)
          wairProcess = true;
          simulate = false;
        
  
}

function draw() {
  
  
  
 if(windowWidth!=prevWidth||windowHeight!=prevHeight){
   
   freq_button_offset = freq_button_offset*windowWidth/prevWidth
   speed_button_offset = speed_button_offset*windowWidth/prevWidth
   
   
   
    for (let b = 0; b < antennas.length; b++) {
 
    amp_button_offset[b]*=windowWidth/prevWidth;
    phase_button_offset[b] *=windowWidth/prevWidth;
    sep_button_offset[b] *=windowWidth/prevWidth;     
      
    }
   
   
  
   
   
     if(simulate){
 
   if(prevWidth<windowWidth||prevHeight<windowHeight){
        
    simulate  = false;
    waitProcess = true;
                                
    }
      
      
  }
   
    
   prevWidth = windowWidth;
   prevHeight = windowHeight;
   
   resizeCanvas(windowWidth, windowHeight);
   
  width = Math.floor((1/1.4)*windowWidth)
  height = Math.floor((8/10)*windowHeight)

    
    }
     
 

  
 
  
  background(0, 0, 0);
 
  

  if (simulate) {
    
     if(!pause){
    timeSim += (millis() / timeScale)-time
      
    }
    time = millis() / timeScale;
   
    
    
    timePhase = ComplexNum.cis(timeSim * w);

    let pixelIndex = 0; 

    // Load the pixel buffer once at the start
    // Load the pixel buffer once at the start

    loadPixels();
    

    for (let i = 2; i + 2 < width / sLength; i++) {
      for (let j = 2; j + 2 < height / sLength; j++) {
        let tempInd = j * M + i;

        if (!show_BField) {
          E[2 * tempInd] = Math.fround(
            E_phase_amp[4 * tempInd] *
              Math.cos(E_phase_amp[4 * tempInd + 1] + timeSim * w)
          );
          E[2 * tempInd + 1] = Math.fround(
            E_phase_amp[4 * tempInd + 2] *
              Math.cos(E_phase_amp[4 * tempInd + 3] + timeSim * w)
          );
        }

        if (show_BField || show_EnergyFlux) {
          B[tempInd] = Math.fround(
            B_phase_amp[2 * tempInd] *
              Math.cos(B_phase_amp[2 * tempInd + 1] + timeSim * w)
          );
        }

        let screenX = sLength * i;
        let screenY = sLength * j;

        // Calculate color values
        let r, g, b, a;

        if (show_BField) {
          let B_temp = B[tempInd];
          let colorSizeB_Blue =
            255 *
            squiz(
              -B[tempInd] * Math.abs(B_temp),
              
              k1,
              (level = 256)
            );
          let colorSizeB_Red =
            255 *
            squiz(B_temp * Math.abs(B_temp), k1, (colorRes = 256));
          let colorSizeMag = squiz(
            Math.abs(B_temp),
          
            k1B,
            (levels = 256)
          );

          r = Math.round(colorSizeB_Red);
          g = 0;
          b = Math.round(colorSizeB_Blue);
          a = colorSizeMag;
        } else if (show_EField) {
          let E_mag_2 =
            E[2 * tempInd] * E[2 * tempInd] +
            E[2 * tempInd + 1] * E[2 * tempInd + 1];
          let colorSizeEA = squiz(E_mag_2, k3, (levels = 256));
          let colorSizeEB =
            255 * squiz(E_mag_2, k4, (levels = 256));

          // Simplified HSB to RGB - assuming your original colorMode(HSB, 500, 100, 100)
          let hue = (255 - colorSizeEB) / 500; // 0-1 range
          let sat = 1.0;
          let bright = 1.0;

          // Simple HSB to RGB conversion
          let c = bright * sat;
          let x = c * (1 - Math.abs(((hue * 6) % 2) - 1));
          let m = bright - c;

          let r1, g1, b1;
          let h = hue * 6;
          if (h < 1) {
            r1 = c;
            g1 = x;
            b1 = 0;
          } else if (h < 2) {
            r1 = x;
            g1 = c;
            b1 = 0;
          } else if (h < 3) {
            r1 = 0;
            g1 = c;
            b1 = x;
          } else if (h < 4) {
            r1 = 0;
            g1 = x;
            b1 = c;
          } else if (h < 5) {
            r1 = x;
            g1 = 0;
            b1 = c;
          } else {
            r1 = c;
            g1 = 0;
            b1 = x;
          }

          r = Math.round((r1 + m) * 255);
          g = Math.round((g1 + m) * 255);
          b = Math.round((b1 + m) * 255);
          a = colorSizeEA;
        } else if (show_EnergyFlux) {
          let Energy_flux_mag =
           (
              E[2 * tempInd] * E[2 * tempInd] +
                E[2 * tempInd + 1] * E[2 * tempInd + 1]
            ) * Math.abs(B[tempInd]);
          let colorSizeEnergyA = squiz(
            Energy_flux_mag,
            
            k5,
            (levels = 256)
          );
          let colorSizeEnergyB =
            255 * squiz(Energy_flux_mag, k6, (levels = 256));

          // Similar HSB conversion
          let hue = (255 - colorSizeEnergyB) / 500;
          let sat = 1.0;
          let bright = 1.0;

          let c_val = bright * sat;
          let x = c_val * (1 - Math.abs(((hue * 6) % 2) - 1));
          let m = bright - c_val;

          let r1, g1, b1;
          let h = hue * 6;
          if (h < 1) {
            r1 = c_val;
            g1 = x;
            b1 = 0;
          } else if (h < 2) {
            r1 = x;
            g1 = c_val;
            b1 = 0;
          } else if (h < 3) {
            r1 = 0;
            g1 = c_val;
            b1 = x;
          } else if (h < 4) {
            r1 = 0;
            g1 = x;
            b1 = c_val;
          } else if (h < 5) {
            r1 = x;
            g1 = 0;
            b1 = c_val;
          } else {
            r1 = c_val;
            g1 = 0;
            b1 = x;
          }

          r = Math.round((r1 + m) * 255);
          g = Math.round((g1 + m) * 255);
          b = Math.round((b1 + m) * 255);
          a = colorSizeEnergyA;
        }

        // Ensure values are within valid range
        r = Math.max(0, Math.min(255, r));
        g = Math.max(0, Math.min(255, g));
        b = Math.max(0, Math.min(255, b));
        a = Math.max(0, Math.min(255, a));

        // Fill the rectangular block - make sure we're only in the 1000x800 region
        let maxX = width; // Limit to your actual drawing area
        let maxY = height;

        let endX = Math.min(screenX + sLength, maxX);
        let endY = Math.min(screenY + sLength, maxY);

        // Only draw if we're within the intended region
        if (screenX < maxX && screenY < maxY) {
          for (let di = 0; di < sLength; di++) {
            for (let dj = 0; dj < sLength; dj++) {
              pixelIndex = ((j * sLength + dj) * windowWidth + i * sLength + di) * 4;

              pixels[pixelIndex] = r * a; // Red
              pixels[pixelIndex + 1] = g * a; // Green
              pixels[pixelIndex + 2] = b * a; // Blue
              pixels[pixelIndex + 3] = 255; // Alpha
            }
          }
        }
      }
    }

    // Update pixels once at the end
    updatePixels();
  }

    for (let b = 0; b < antennas.length && !simulate; b++) {
    antenna = antennas[b];
    
    if(antenna.getType()=="Dipole"){

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
          
          if(resolution==1){
            
            dipole.setDl(dl_lr)
            
          }
          
          
          else if(resolution==2){
            
            dipole.setDl(dl_mr)
            
          }
          
          
          else {
            
            dipole.setDl(dl_hr)
            
          }
          

          antennas.push(dipole);

          amp_button_offset.push((defAmp / maxAmp) * ((300 - 40)/1400)*windowWidth);
          phase_button_offset.push(0);
          sep_button_offset.push(defaultDipoleSep*((300 - 40)/1400)*windowWidth/maxSep);

          flags_amp_button.push(false);
          flags_phase_button.push(false);
          flags_sep_button.push(false);
          delButtonPressed.push(false)
            

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
  rect(width, 0, windowWidth-width, height);
  rect(0, height, windowWidth, height/2);

  //dipole antenna button

  stroke(0, 0, 0);
  fill(0, 0, 0);
  

  fill(150, 150, 150);

  if (dipole_antenna_pressed) {
    fill(120, 120, 120);
    stroke(250, 250, 250);
  }

  rect((35/1400)*windowWidth, (860/1000)*windowHeight, (175/1400)*windowWidth, (70/1000)*windowHeight);

  fill(150, 150, 150);
  stroke(0, 0, 0);

  if (
    !dipole_antenna_pressed &&
    inRect((35/1400)*windowWidth,(860/1000)*windowHeight , (175/1400)*windowWidth, (70/1000)*windowHeight, mouseX, mouseY) &&
    mousePressedFlag
  ) {
    dipole_antenna_pressed = true;
    mousePressedFlag = false;
    simulate = false;
  } else if (
    dipole_antenna_pressed &&
   inRect((35/1400)*windowWidth,(860/1000)*windowHeight , (175/1400)*windowWidth, (70/1000)*windowHeight, mouseX, mouseY) &&
    mousePressedFlag
  ) {
    dipole_antenna_pressed = false;
    mousePressedFlag = false;
    addNew_dipole_tx = false;
  }

  if (dipole_antenna_pressed) {
    addNew_dipole_tx = true;
  }
  
  
   stroke(0, 0, 0);
  fill(0, 0, 0);
  
  push();

  scale(windowWidth/1400, windowHeight/1000);


  textSize(23);
  text("Dipole Antenna", (44), (902));
  
  pop();
  
  
  
  //pause button
  
  stroke(0,0,0)
  
  if(!pause){
  fill(240,0,0)
    
  }
  
  else{
    
    fill(0,240,0)
    
  }
  
  
   if (inRect((920/1400)*windowWidth, (860/1000)*windowHeight, (80/1400)*windowWidth, (80/1000)*windowHeight, mouseX, mouseY)) {
  
    
     if((!pause)&&mousePressedFlag){
       
       mousePressedFlag = false;
       mouseRelease = false;
       pause = true;
     }
     
     
     else if(pause&&mousePressedFlag){
       
       mousePressedFlag = false;
       mouseRelease = false ;
       pause = false;
       
     }
     
   }
  
  rect((920/1400)*windowWidth, (860/1000)*windowHeight, (80/1400)*windowWidth, (80/1000)*windowHeight);
  
  
   if(!pause){
  
  fill(0, 0, 0);
     
   push();

  scale(windowWidth/1400, windowHeight/1000);   
  textSize((24));
  text("Pause", (925), (852));
    pop();
  }
  
  else{
    
  fill(0, 0, 0);
    
  push();

  scale(windowWidth/1400, windowHeight/1000);
  textSize(24);
  text("Resume", (915), (852));
  pop()
    
  }
  
  //simulation status button

  if (simulate) {
    fill(0, 250, 0);
    

    if (inRect((1150/1400)*windowWidth, (515/1000)*windowHeight, (180/1400)*windowWidth, (60/1000)*windowHeight, mouseX, mouseY)) {
      fill(0, 150, 0);
       
      if (mousePressedFlag) {
        waitProcess = false;
        simulate = false;
        mousePressedFlag = false;
      }
    }

    rect((1150/1400)*windowWidth, (515/1000)*windowHeight, (180/1400)*windowWidth, (60/1000)*windowHeight);

    fill(0, 0, 0);
     push();

  scale(windowWidth/1400, windowHeight/1000);
    
    textSize(28);
    text("Simulating", (1170), (553));
    
    pop()
  } else {
    fill(220, 0, 0);

    if (inRect((1150/1400)*windowWidth, (515/1000)*windowHeight, (180/1400)*windowWidth, (60/1000)*windowHeight, mouseX, mouseY)) {
      fill(255, 50, 0);

      if (mousePressedFlag) {
        waitProcess = true;
        simulate = false;
        pause = false;
        mousePressedFlag = false;
      }
    }

    rect((1150/1400)*windowWidth, (515/1000)*windowHeight, (180/1400)*windowWidth, (60/1000)*windowHeight);

    fill(0, 0, 0);
    
     push();

  scale(windowWidth/1400, windowHeight/1000);
    
    textSize(28);
    text("Simulate", 1185, 553);
    
    pop()
  }

  fill(150, 150, 150);

  if (inRect((1150/1400)*windowWidth, (615/1000)*windowHeight, (180/1400)*windowWidth, (60/1000)*windowHeight, mouseX, mouseY)) {
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
    }
  }

  rect((1150/1400)*windowWidth, (615/1000)*windowHeight, (180/1400)*windowWidth, (60/1000)*windowHeight);

 
  push()
  fill(0, 0, 0);
  scale(windowWidth/1400, windowHeight/1000);
  textSize(26);
  text("Clear All", 1190, 653);

  textSize(20);

  
  
  pop()
  
  fill(150, 150, 150);
  
  
  //resolution buttons logic
  
  let flagChangeRes = false;
  
   if (
    resolution!=3 &&
    inRect((1188/1400)*windowWidth, (770/1000)*windowHeight, (110/1400)*windowWidth, (50/1000)*windowHeight, mouseX, mouseY) &&
    mousePressedFlag
  ) {
   
      resolution=3;
    mousePressedFlag = false;
      simulate = false;
      flagChangeRes = true;
     sLength = 2    
      arrow_spacing = 10
  }

  
  
   if (
    resolution!=2 &&
    inRect((1188/1400)*windowWidth, (840/1000)*windowHeight, windowWidth*110/1400, windowHeight*50/1000, mouseX, mouseY) &&
    mousePressedFlag
  ) {
   
      resolution=2;
    mousePressedFlag = false;
      simulate = false;
     flagChangeRes = true;
      sLength = 4
      arrow_spacing = 5
     
      
  }

 
   if (
    resolution!=1 &&
    inRect(windowWidth*1188/1400, windowHeight*910/1000, windowWidth*110/1400, windowHeight*50/1000, mouseX, mouseY) &&
    mousePressedFlag
  ) {
   
      resolution=1;
    mousePressedFlag = false;
     simulate = false;
     flagChangeRes= true;
      sLength = 5
      arrow_spacing = 4
      
       
  }

  
  if(flagChangeRes){
    
    N = Math.ceil(height / sLength);
    M = Math.ceil(width / sLength);
    
     const_rFactor1 = Math.fround((Scale * Scale) / (sLength * sLength));
     const_rFactor2 = Math.fround(Scale / (2 * sLength));


    
     for (let i = 0; i < width / sLength; i++) {
    for (let j = 0; j < height / sLength; j++) {
      let tempInd = M * j + i;

      E[2 * tempInd] = 0;
      E[2 * tempInd + 1] = 0;
      B[tempInd] = 0;
    }
  }
    
  }
  
  
// resolution buttons graphics 
  
  
  if(resolution==3){
    
     fill(120, 120, 120);
    stroke(250, 250, 250);     
  }
  
  rect(windowWidth*1188/1400, windowHeight*770/1000, windowWidth*110/1400, windowHeight*50/1000);
  
   stroke(0, 0, 0);
  fill(150, 150, 150);
  
    if(resolution==2){
    
     fill(120, 120, 120);
    stroke(250, 250, 250);     
  }
  
 
  
   rect(windowWidth*1188/1400, windowHeight*840/1000, windowWidth*110/1400, windowHeight*50/1000);
  
   stroke(0, 0, 0);
  fill(150, 150, 150);


  
    if(resolution==1){
    
     fill(120, 120, 120);
    stroke(250, 250, 250);     
  }
  
  
   rect(windowWidth*1188/1400, windowHeight*910/1000, windowWidth*110/1400, windowHeight*50/1000);
  
   stroke(0, 0, 0);
  fill(150, 150, 150);

  
  
  
   stroke(0, 0, 0);
  fill(0, 0, 0);
  textSize(25);
  push()
  scale(windowWidth/1400, windowHeight/1000);
  text("Resolution:", 1185, 750);

  text("high",1220, 800)
  
   text("medium",1198, 870)
  
   text("low",1224, 942)
  pop()
   stroke(0, 0, 0);
  fill(150, 150, 150);

  
  
  

  //fields GUI

  //E field button
  if (show_EField) {
    fill(120, 120, 120);
    stroke(250, 250, 250);
  }

  rect(windowWidth*1130/1400, windowHeight*50/1000, windowWidth*148/1400, windowHeight*50/1000);

  push()
  stroke(0, 0, 0);
  fill(0, 0, 0);
  textSize(20);
  scale(windowWidth/1400, windowHeight/1000);
  text("Show E Field", 1144, 80);
  
  pop()
 

  stroke(0, 0, 0);
  fill(150, 150, 150);

  if (
    !show_EField &&
    inRect(windowWidth*1130/1400, windowHeight*50/1000, windowWidth*148/1400, windowHeight*50/1000, mouseX, mouseY) &&
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

  rect(windowWidth*1130/1400, windowHeight*150/1000, windowWidth*148/1400, windowHeight*50/1000);

  stroke(0, 0, 0);
  fill(0, 0, 0);
  push()
  textSize(20)
  scale(windowWidth/1400, windowHeight/1000);
  text("Show B Field", 1144, 180);
  pop()
  textSize(20);

  fill(150, 150, 150);

  if (
    !show_BField &&
    inRect(windowWidth*1130/1400,windowHeight* 150/1000, windowWidth*148/1400, windowHeight*50/1000, mouseX, mouseY) &&
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

  rect(windowWidth*1130/1400, windowHeight*250/1000, windowWidth*190/1400, windowHeight*50/1000);

  push()
  stroke(0, 0, 0);
  fill(0, 0, 0);
  scale(windowWidth/1400, windowHeight/1000);
  text("Show energy flux", 1144, 280);
  pop()
  textSize(20);

  fill(150, 150, 150);

  if (
    !show_EnergyFlux &&
    inRect(windowWidth*1130/1400, windowHeight*250/1000, windowWidth*190/1400, windowHeight*50/1000, mouseX, mouseY) &&
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
  
  push()
  scale(windowWidth/1400, windowHeight/1000);
  text("Frequency", 1180, 340);
  pop()
  textSize(20);

  fill(250, 250, 250);

  rect(windowWidth*1130/1400, windowHeight*350/1000, windowWidth*200/1400, windowHeight*30/1000);
  


  //frequency slider button

  fill(110, 110, 110);

  if (
    mousePressedFlag &&
    inRect(windowWidth*1130/1400 + freq_button_offset, windowHeight*350/1000, windowWidth*30/1400, windowHeight*30/1000, mouseX, mouseY)
  ) {
    freq_slider_on = true;
  } else if (!mousePressedFlag ||  !inRect(windowWidth*1130/1400, windowHeight*350/1000, windowWidth*200/1400, windowHeight*30/1000,mouseX,mouseY)) {
    freq_slider_on = false;
  }

  if (freq_slider_on) {
    freq_slider_process = true;
    simulate = false;

    fill(90, 90, 90);

    freq_button_offset =mouseX - windowWidth*(1130 + 15)/1400;
    if (freq_button_offset < 0) {
      freq_button_offset = 0;
    } else if (freq_button_offset + windowWidth*30/1400 > windowWidth*200/1400) {
      freq_button_offset = (200 - 30)*windowWidth/1400;
    }
  }

  rect(windowWidth*1130/1400 + freq_button_offset, windowHeight*350/1000, windowWidth*30/1400, windowHeight*30/1000);
  
  

  if (freq_slider_process && mouseRelease) {
    freq = minFreq + (freq_button_offset / (windowWidth*(200 - 30)/1400)) * (maxFreq - minFreq);
    

    freq_slider_process = false;

    for (let i = 0; i < antennas.length; i++) {
      antennas[i].setWavelength(c / freq);
    }
  }

  //speed slider
  stroke(0, 0, 0);
  fill(0, 0, 0);
  push()
  scale(windowWidth/1400, windowHeight/1000);
  text("Speed", 1200, 420);
  pop()
  textSize(20);

  stroke(0, 0, 0);
  fill(250, 250, 250);

  rect(windowWidth*1130/1400, windowHeight*430/1000, windowWidth*200/1400, windowHeight*30/1000);

  //speed of EM waves slider button

  stroke(0, 0, 0);
  fill(110, 110, 110);

  if (
    mousePressedFlag &&
    inRect(windowWidth*1130/1400 + speed_button_offset, windowHeight*430/1000, windowWidth*30/1400, windowHeight*30/1000, mouseX, mouseY)
  ) {
    speed_slider_on = true;
  } else if (!mousePressedFlag || !inRect(windowWidth*1130/1400, windowHeight*430/1000, windowWidth*200/1400, windowHeight*30/1000,mouseX,mouseY)) {
    speed_slider_on = false;
  }

  if (speed_slider_on) {
    speed_slider_process = true;
    simulate = false;

    fill(90, 90, 90);

    speed_button_offset = mouseX - windowWidth*(1130 + 15)/1400;
    if (speed_button_offset < 0) {
      speed_button_offset = 0;
    } else if (speed_button_offset + windowWidth*30/1400 > windowWidth*200/1400) {
      speed_button_offset = windowWidth*(200 - 30)/1400;
    }
  }

  rect(windowWidth*1130/1400 + speed_button_offset, windowHeight*430/1000, windowWidth*30/1400, windowHeight*30/1000);

  if (speed_slider_process && mouseRelease) {
    c = minSpeed + (speed_button_offset /(windowWidth*(200 - 30)/1400)) * (maxSpeed - minSpeed);

    speed_slider_process = false;

    for (let i = 0; i < antennas.length; i++) {
      antennas[i].setWavelength(c / freq);
    }
  }
  
  
  // process changes 

  if (waitProcess) {
    
    
    let k = (2 * Math.PI * freq) / c;
   
    cFACTOR1 = new ComplexNum(0, -2 * Math.PI * freq);
    cFACTOR2 = new ComplexNum(0, -c / k);
    N = Math.ceil(height / sLength);
    M = Math.ceil(width / sLength);
    
    for (let y = 0; y < height; y += sLength) {
      for (let x = 0; x < width; x += sLength) {
        let x_ind = int(x / sLength);
        let y_ind = int(y / sLength);

        r[0] = conMyX(x);
        r[1] = conMyY(y);

        let tempInd = y_ind * M + x_ind;

        A[2 * tempInd] = zero;
        A[2 * tempInd + 1] = zero;

        for (let i = 0; i < antennas.length; i++) {
          
          
          antenna = antennas[i];
          
          
          
          if(resolution==1){
            
            antenna.setDl(dl_lr)
            
          }
          
          
          else if(resolution==2){
            
            antenna.setDl(dl_mr)
            
          }
          
          
          else {
            
            antenna.setDl(dl_hr)
            
          }

          segments = antenna.getSegments();

          currents = antenna.getCurrentSegments();

          for (let j = 0; j < segments.length; j++) {
            let r_tag = [segments[j][0], segments[j][1]];
            let distance = length2D(r, r_tag);
            
            
            //for numerical stability near current sources
            distance = distance + 0.1;

            current[0] = new ComplexNum(
              currents[j][0] * antenna.getDl(),
              0
            ).product(antenna.getI0());
            current[1] = new ComplexNum(
              currents[j][1] * antenna.getDl(),
              0
            ).product(antenna.getI0());

            phase = ComplexNum.cis(-k * distance);

            factor = phase.product(1 / distance);

            A[2 * tempInd] = A[2 * tempInd].add(current[0].product(factor));
            A[2 * tempInd + 1] = A[2 * tempInd + 1].add(
              current[1].product(factor)
            );
          }
        }
      }
    }
    

    for (let i = 1; i + 2 < width / sLength; i++) {
      for (let j = 1; j + 2 < height / sLength; j++) {
        let tempInd = M * j + i;

        E_phasor[2 * tempInd] = A[2 * tempInd].product(cFACTOR1);
        E_phasor[2 * tempInd + 1] = A[2 * tempInd + 1].product(cFACTOR1);

        //calculate E phasor

        tempE1[0] = A[2 * tempInd + 2]
          .add(A[2 * tempInd - 2])
          .add(A[2 * tempInd].product(-2));

        tempE1[1] = A[2 * tempInd - 2 * M + 1]
          .add(A[2 * tempInd + 2 * M + 1])
          .add(A[2 * tempInd + 1].product(-2));

        tempE1[0] = tempE1[0].product(const_rFactor1);
        tempE1[1] = tempE1[1].product(const_rFactor1);

        tempE2[0] = A[2 * tempInd + 2 * M - 1].add(A[2 * tempInd - 2 * M + 3]);

        tempE2[0] = tempE2[0].subtract(
          A[2 * tempInd + 2 * M + 3].add(A[2 * tempInd - 2 * M - 1])
        );

        tempE2[1] = A[2 * tempInd + 2 * M - 2].add(A[2 * tempInd - 2 * M + 2]);

        tempE2[1] = tempE2[1].subtract(
          A[2 * tempInd + 2 * M + 2].add(A[2 * tempInd - 2 * M - 2])
        );

        tempE2[0] = tempE2[0].product(0.25 * const_rFactor1);
        tempE2[1] = tempE2[1].product(0.25 * const_rFactor1);

        tempE[0] = tempE1[0].add(tempE2[0]);
        tempE[1] = tempE1[1].add(tempE2[1]);

        tempE[0] = tempE[0].product(cFACTOR2);
        tempE[1] = tempE[1].product(cFACTOR2);

        E_phasor[2 * tempInd] = E_phasor[2 * tempInd].add(tempE[0]);
        E_phasor[2 * tempInd + 1] = E_phasor[2 * tempInd + 1].add(tempE[1]);

        E_phase_amp[4 * tempInd] = E_phasor[2 * tempInd].getAbsVal();
        E_phase_amp[4 * tempInd + 1] = E_phasor[2 * tempInd].getPhase();
        E_phase_amp[4 * tempInd + 2] = E_phasor[2 * tempInd + 1].getAbsVal();
        E_phase_amp[4 * tempInd + 3] = E_phasor[2 * tempInd + 1].getPhase();

        //calculate B phasor
        tempB1 = A[2 * tempInd - 2 * M].subtract(A[2 * tempInd + 2 * M]);
        tempB2 = A[2 * tempInd + 3].subtract(A[2 * tempInd - 1]);

        tempB1 = tempB1.product(const_rFactor2);
        tempB2 = tempB2.product(const_rFactor2);

        B_phasor[tempInd] = tempB2.subtract(tempB1);

        B_phase_amp[2 * tempInd] = B_phasor[tempInd].getAbsVal();
        B_phase_amp[2 * tempInd + 1] = B_phasor[tempInd].getPhase();
      }
    }

  
    waitProcess = false;
    simulate = true;
  }

  //after processing

  if (simulate) {
   

    w  = 2 * Math.PI * freq;

    timePhase = ComplexNum.cis(timeSim * w);

    for (let i = 0; i < width / sLength; i += arrow_spacing) {
      for (let j = 0; j < height / sLength; j += arrow_spacing) {
        let tempInd = M * j + i;

        if (show_EField) {
          let E_mag = Math.fround(
            Math.sqrt(
              E[2 * tempInd] * E[2 * tempInd] +
                E[2 * tempInd + 1] * E[2 * tempInd + 1]
            )
          );

          let arrow_l =
            maxArrowLen * squiz(E_mag, k2, (levels = 256));

          let screenX = sLength * i;
          let screenY = sLength * j;

          let dx = E[2 * tempInd];
          let dy = E[2 * tempInd + 1];

          dx *= arrow_l / E_mag;
          dy *= arrow_l / E_mag;

          stroke(0, 0, 0);
          fill(0, 0, 0);

          circle(screenX, screenY, 2);

          line(screenX, screenY, screenX + dx, screenY - dy);

          noStroke();
        }

        if (show_EnergyFlux) {
          let E_flux_mag =
            Math.sqrt(
              E[2 * tempInd] * E[2 * tempInd] +
                E[2 * tempInd + 1] * E[2 * tempInd + 1]
            ) * Math.abs(B[tempInd]);

          let arrow_l =
            maxArrowLen * squiz(E_flux_mag, k7, (levels = 256));

          let screenX = sLength * i;
          let screenY = sLength * j;

          let dx = E[2 * tempInd + 1] * B[tempInd];
          let dy = -(E[2 * tempInd] * B[tempInd]);

          dx *= arrow_l / E_flux_mag;
          dy *= arrow_l / E_flux_mag;

          fill(0, 0, 0);

          circle(screenX, screenY, 2);

          stroke(0, 0, 0);

          line(screenX, screenY, screenX + dx, screenY - dy);
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

        let currentColorMag =
          255 * squiz(currentMag, k8, (levels = 256));

        if (antenna.getSegFlags()[j]&&conScreenX(Math.max(x,x_next))<width&&conScreenY(Math.max(y,y_next))<height) {
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
     if(antenna.getShowBox()&&inRect(xBox , yBox,  windowWidth*400/1400,  windowHeight*370/1000, mouseX, mouseY)){
        
        isMouseInStatBox = true;
        
      }
  
  
   }
  
  
  
 for (let b = 0; b < antennas.length; b++) {
    antenna = antennas[b];

    let xBox = conScreenX(antenna.getXBox());
    let yBox = conScreenY(antenna.getYBox());
   
   
    
    if(antenna.getType()=="Dipole"){

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

      rect(xBox, yBox, windowWidth*400/1400, windowHeight*370/1000);
      
    

      stroke(0, 0, 0);
      fill(250, 250, 250);

      //amplitude slide

      rect(xBox + windowWidth*50/1400, yBox + windowHeight*55/1000, windowWidth*300/1400, windowHeight*30/1000);

      stroke(0, 0, 0);
      fill(0, 0, 0);
      push()
      scale(windowWidth/1400, windowHeight/1000);
      textSize(25);
      text("Amplitude", xBox*1400/windowWidth +145, yBox*1000/windowHeight + 48);
      pop()
      stroke(0, 0, 0);
      fill(110, 110, 110);

      //amplitude slide button

      if (
        inRect(
          xBox + windowWidth*50/1400 + amp_button_offset[b],
          yBox + windowHeight*55/1000,
          windowWidth*40/1400,
          windowHeight*30/1000,
          mouseX,
          mouseY
        ) &&
        mousePressedFlag
      ) {
        flags_amp_button[b] = true;
      }

      if (
        !mousePressedFlag ||
        !inRect(xBox + windowWidth*50/1400, yBox + windowHeight*55/1000, windowWidth*300/1400, windowHeight*30/1000, mouseX, mouseY)
      ) {
        flags_amp_button[b] = false;
      }

      if (flags_amp_button[b]) {
        amp_slider_process = true;
        simulate = false;
        antenna.setAmp((maxAmp * amp_button_offset[b]) / (windowWidth*(300 - 40)/1400));

        fill(90, 90, 90);

        amp_button_offset[b] = mouseX - xBox - windowWidth*(50 +20)/1400;
        if (amp_button_offset[b] < 0) {
          amp_button_offset[b] = 0;
        } else if (amp_button_offset[b] + windowWidth*40/1400 > windowWidth*300/1400) {
          amp_button_offset[b] = (300 - 40)*windowWidth/1400;
        }
      }

      rect(xBox + windowWidth*50/1400 + amp_button_offset[b], yBox + windowHeight*55/1000, windowWidth*40/1400,windowHeight* 30/1000);

      stroke(0, 0, 0);
      fill(250, 250, 250);

      //phase slide
      rect(xBox + windowWidth*50/1400, yBox + windowHeight*155/1000, windowWidth*300/1400, windowHeight*30/1000);

      stroke(0, 0, 0);
      fill(0, 0, 0);
      push()
        scale(windowWidth/1400, windowHeight/1000);
      textSize(25);
      text("Phase", 1400*xBox/windowWidth + 170, yBox*1000/windowHeight + 148);
      
      
       text("0", 1400*xBox/windowWidth +22, yBox*1000/windowHeight+178);
       text("2π", 1400*xBox/windowWidth +357, yBox*1000/windowHeight+178);
       pop()

      fill(110, 110, 110);

      if (
        inRect(
          xBox + windowWidth*50/1400 + phase_button_offset[b],
          yBox + 155*windowHeight/1000,
          windowWidth*40/1400,
          windowHeight*30/1000,
          mouseX,
          mouseY
        ) &&
        mousePressedFlag
      ) {
        flags_phase_button[b] = true;
      }

      if (
        !mousePressedFlag ||
        !inRect(xBox + windowWidth*50/1400, yBox + windowHeight*155/1000, windowWidth*300/1400, windowHeight*30/1000, mouseX, mouseY)
      ) {
        flags_phase_button[b] = false;
      }

      if (flags_phase_button[b]) {
        phase_slider_process = true;
        simulate = false;
        antenna.setPhase((2 * Math.PI * phase_button_offset[b]) / (windowWidth*(300 - 40)/1400));

        fill(90, 90, 90);

        phase_button_offset[b] = mouseX - xBox - windowWidth*(50 +20)/1400;
        if (phase_button_offset[b] < 0) {
          phase_button_offset[b] = 0;
        } else if (phase_button_offset[b] + windowWidth*40/1400 > windowWidth*300/1400) {
          phase_button_offset[b] = windowWidth*(300 - 40)/1400;
        }
      }

      rect(xBox + windowWidth*50/1400 + phase_button_offset[b], yBox + windowHeight*155/1000, windowWidth*40/1400, windowHeight*30/1000);
      
       stroke(0, 0, 0);
      fill(250, 250, 250);
      
      // seperation slide
      
      
     
      rect(xBox + windowWidth*50/1400, yBox + windowHeight*255/1000, windowWidth*300/1400, windowHeight*30/1000);

      push()
      stroke(0, 0, 0);
      fill(0, 0, 0);
      scale(windowWidth/1400, windowHeight/1000);
      textSize(25);
      text("Seperation", 1400*xBox/windowWidth + 150, 1000*yBox/windowHeight + 248);
      pop()

      fill(110, 110, 110);

      if (
        inRect(
          xBox + windowWidth*50/1400 + sep_button_offset[b],
          yBox + windowHeight*255/1000,
          windowWidth*40/1400,
          windowHeight*30/1000,
          mouseX,
          mouseY
        ) &&
        mousePressedFlag
      ) {
        flags_sep_button[b] = true;
      }

      if (
        !mousePressedFlag ||
        ! inRect(
          xBox + windowWidth*50/1400 + sep_button_offset[b],
          yBox + windowHeight*255/1000,
          windowWidth*40/1400,
          windowHeight*30/1000,
          mouseX,
          mouseY
        ) 
      ) {
        flags_sep_button[b] = false;
      }

      if (flags_sep_button[b]) {
        phase_slider_process = true;
        simulate = false;
        antenna.setSep((maxSep* sep_button_offset[b]) / (windowWidth*(300 - 40)/1400));

        fill(90, 90, 90);

        sep_button_offset[b] = mouseX - xBox - windowWidth*(50 + 20)/1400;
        if (sep_button_offset[b] < 0) {
          sep_button_offset[b] = 0;
        } else if (sep_button_offset[b] + windowWidth*40/1400 > windowWidth*300/1400) {
          sep_button_offset[b] = windowWidth*(300 - 40)/1400;
        }
      }

      rect(xBox + windowWidth*50/1400 + sep_button_offset[b], yBox +windowHeight* 255/1000, windowWidth*40/1400, windowHeight*30/1000);
      

      //delete antenna button
      
      fill(240,0,0)
      
      if(inRect(xBox + windowWidth*128/1400, yBox + windowHeight*308/1000, windowWidth*150/1400, windowHeight*40/1000, mouseX, mouseY)&&mousePressedFlag){
         amp_button_offset.splice(b,1)
         phase_button_offset.splice(b,1);
         sep_button_offset.splice(b,1);
         flags_amp_button.splice(b,1);
         flags_phase_button.splice(b,1);
         flags_sep_button.splice(b,1);
        delButtonPressed.splice(b,1)
        antennas.splice(b,1)
        mousePressedFlag = false;
        simulate = false;
      }
     
      
    rect(xBox + windowWidth*128/1400, yBox + windowHeight*308/1000, windowWidth*150/1400, windowHeight*40/1000);
      
    stroke(0, 0, 0);
    fill(0, 0, 0);
    push()
    scale(windowWidth/1400, windowHeight/1000);
    textSize(26);
    text("delete", xBox*1400/windowWidth + 170, yBox*1000/windowHeight + 335);
    pop()
      


    }
      
      
    if (inRectGen(p1F, p2F, p3F, p4F, mouseX, mouseY) && mousePressedFlag&&!isMouseInStatBox) {
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
