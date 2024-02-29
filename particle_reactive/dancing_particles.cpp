// Yuehao Gao | MAT201B
// 2022-02-23 | Final Project
// Audio-reactive visualizer
//

#include "al/app/al_App.hpp"
#include "al/app/al_GUIDomain.hpp"
//#include "al/ui/al_ControlGUI.hpp"
#include "al/graphics/al_Shapes.hpp"
#include "al/math/al_Functions.hpp"    // al::abs
#include "al/math/al_Random.hpp"
#include "al/ui/al_Parameter.hpp"

#include "Gamma/Analysis.h"
#include "Gamma/DFT.h"
#include "Gamma/Effects.h"
#include "Gamma/SamplePlayer.h"
#include <cmath>
#include <fstream>
#include <vector>

using namespace al;
using namespace std;
#define FFT_SIZE 4048

const float sphereSize = 0.8;

const int pattern1CylinderHeight = 70;
const float pattern1CylinderHalfLength = 0.6;
const int pattern1CylinderNumParticlePerlayer = 60;
const float pattern1AngleIncrement = 2.0 * M_PI / pattern1CylinderNumParticlePerlayer;

const int pattern2BallNumParticle = 1000;

// Common-used statements or parameters
struct CommonState {
  Vec3f pattern1Position[17000];
  Vec3f pattern2Position[17000];
  int size;
  Nav nav;
  float pointSize;
};

// To generate a random 3D vector for position, color, ...
Vec3f randomVec3f(float scale) {
  return Vec3f(rnd::uniformS(), rnd::uniformS(), rnd::uniformS()) * scale;
}

// To slurp a file (initialized)
string slurp(string fileName);  // forward declaration



// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

struct MyApp : App {
  // ----------- All parameters Start -----------
  Parameter valueL{"value_left", 0, 0, 1};
  Parameter valueR{"value_right", 0, 0, 1};
  Parameter pattern{"visual_pattern", 0, 0, 3};     // This will be limited to int only
  Parameter musicPower{"/musicPower", "", 2.5, 0.5, 6.0};
  Parameter pointSize{"/pointSize", "", 0.5, 0.1, 1.5};
  Parameter timeStep{"/timeStep", "", 1.0, 0.0, 1.5};
  Parameter dragFactor{"/dragFactor", "", 0.15, 0.0, 0.6};
  Parameter springConstant{"/springConstant", "", 1.0, 0.1, 5.0};
  Parameter pattern1CylinderRadius{"/cylinderRadius", "", 0.25, 0.1, 1.0};
  Parameter sphereRadius{"/sphereRadius", "", 2.0, 0.5, 4.0};
  Parameter repellingConstant{"/repellingConstant", "", 0.0, 0.0, 0.01};

  // The two spectrum tables for both channels
  vector<float> spectrumL;
  vector<float> spectrumR;
  vector<float> spectrum;

  // STFT variables
  // Format of frequency samples: COMPLEX, MAG_PHASE, or MAG_FREQ
  gam::STFT stft = gam::STFT(FFT_SIZE, FFT_SIZE / 4, 0, gam::HANN, gam::MAG_FREQ);
  // --------------------------------------------

  // The sample player for 1 sound track only
  // The envelop follower for the left and right channel
  gam::SamplePlayer<float, gam::ipl::Linear, gam::phsInc::Loop> player;
  gam::EnvFollow<> followLeft;
  gam::EnvFollow<> followRight;

  // The Point-Shader Program
  ShaderProgram pointShader;

  // ----------- All Meshes -------------
  // * Pattern 0: two spheres representing L and R envelop values
  Mesh pattern0SphereL;
  Mesh pattern0SphereR;

  // * Pattern 1: a cylinder of particles dancing according to the spectrum (enveloped)
  Mesh pattern1ParticleCylinder;
  vector<Vec3f> pattern1Velocity;
  vector<Vec3f> pattern1Force;
  vector<Vec3f> pattern1OriginalPosition;   // Archive of original positions of particles
  vector<float> pattern1Mass;

  // * Pattern 2:
  //Mesh pattern2ParticleBall;
  //Vec3f pattern2SphereCenter = {0.0, 0.0, 0.0};

  // * Pattern 3:

  // ----------- All Meshes -------------

  // onInit
  // When the app is first initialized
  //
  void onInit() override {
    // Load the sound file
    player.load("../tv.mp3");
    
    // set up GUI
    auto GUIdomain = GUIDomain::enableGUI(defaultWindowDomain());
    auto& gui = GUIdomain->newGUI();
    gui.add(valueL);
    parameterServer() << valueL;
    gui.add(valueR);
    parameterServer() << valueR;
    gui.add(pattern);
    parameterServer() << pattern;
    gui.add(musicPower);
    parameterServer() << musicPower;
    gui.add(pointSize);
    parameterServer() << pointSize;
    gui.add(timeStep);
    parameterServer() << timeStep;
    gui.add(dragFactor);
    parameterServer() << dragFactor;
    gui.add(springConstant); 
    parameterServer() << springConstant;
    gui.add(pattern1CylinderRadius);
    parameterServer() << pattern1CylinderRadius;
    gui.add(sphereRadius);
    parameterServer() << sphereRadius;
    gui.add(repellingConstant);
    parameterServer() << repellingConstant;

    // Declare the size of the spectrums
    spectrumL.resize(FFT_SIZE / 2 + 1);
    spectrumR.resize(FFT_SIZE / 2 + 1);
    spectrum.resize(FFT_SIZE / 2 + 1);

  }

  // onCreate
  // Create all the initial visual patterns
  // Slurp corresponding files for the point shader program
  // Set envelop followers
  // Initialize the position of camera
  //
  void onCreate() override { 
    pointShader.compile(slurp("../point-vertex.glsl"),
                        slurp("../point-fragment.glsl"),
                        slurp("../point-geometry.glsl"));

    // ----------- Initialize Parameters for All Meshes Begins -------------

    // Pattern 0: move the two spheres seperately to the left and right
    addSphere(pattern0SphereL);
    addSphere(pattern0SphereR);

    // Pattern 1: create all the particles with their initial positions, colors, velocity, and forces
    pattern1ParticleCylinder.primitive(Mesh::POINTS);
  
    for (int layerIndex = 0; layerIndex < pattern1CylinderHeight; layerIndex++) {
      float y = -1.0 * pattern1CylinderHalfLength + (2.0 * pattern1CylinderHalfLength / pattern1CylinderHeight) * layerIndex;
      
      for (int particleIndex = 0; particleIndex < pattern1CylinderNumParticlePerlayer; particleIndex++) {
        float angle = particleIndex * pattern1AngleIncrement;
        float x = pattern1CylinderRadius * cos(angle);
        float z = pattern1CylinderRadius * sin(angle);
        pattern1ParticleCylinder.vertex(Vec3f(x, y, z));
        pattern1OriginalPosition.push_back(Vec3f(x, y, z));  // Archive of original positions of particles

        float hue = 0.7 / pattern1CylinderHeight * layerIndex;
        pattern1ParticleCylinder.color(HSV(hue, 1.0f, 1.0f));

        pattern1Mass.push_back(3.0);
        pattern1ParticleCylinder.texCoord(pow(3.0, 1.0f / 3), 0);
        pattern1Velocity.push_back(Vec3f(0.0, 0.0, 0.0));
        pattern1Force.push_back(Vec3f(0.0, 0.0, 0.0));
      
      }
    }
    //const vector<Vec3f> &pattern1OriginalPositionVec(pattern1ParticleCylinder.vertices());

    // Pattern 2: ...
    // Pattern 3: ...
    
    // ----------- Initialize Parameters for All Meshes Ends -------------
  
    // Set the magnitude of the envelop followers
    followLeft.lag(0.5); 
    followRight.lag(0.5);
  
    // Initialize the position of "camera"
    nav().pos(0, 0, 2.5);
  }

  // onSound
  // Read in the stereo sound file
  //
  void onSound(AudioIOData& io) override {
    while (io()) {
      // Read in the value of left and right channel
      // Seperately envelop them
      // And feed into the "valueL" and "valueR" parameters
      player();
      float sLeft = player.read(0);
      float sRight = player.read(1);

      io.out(0) = sLeft;
      io.out(1) = sRight;

      valueL.set(followLeft(sLeft));
      valueR.set(followRight(sRight));

      // STFT
      ///*
      if (stft(io.out(0))) { 
        // Loop through all the frequency bins
        //int testing = 0;
        for (unsigned k = 0; k < stft.numBins(); ++k) {
          spectrum[k] = 8.0 * tanh(pow(stft.bin(k).real(), 1.5));
          //testing++;
        }
        //cout << testing;
      }
      /*
      if (stft(io.out(1))) { 
        // Loop through all the frequency bins
        for (unsigned k = 0; k < stft.numBins(); ++k) {
          spectrumR[k] = 10.0 * tanh(pow(stft.bin(k).real(), 1.5) );
        }
      }
      */
      //*/

      // https://github.com/adamstark/Gist
    }
  }

  // onDraw
  // Decide how each frame should be drawn
  // The current frame shows which mesh(es) is depending on the "pattern" parameter
  //
  void onDraw(Graphics& g) override {
    // The background changes brightness according to the volume
    g.clear(0.6 * (valueL.get() + valueR.get()) * (valueL.get() + valueR.get()) * 2.0 * musicPower);
    float redColorChange = 5.0 * valueL.get();
    if (redColorChange > 2.0) {
      redColorChange = 2.0;
    }
    float blueColorChange = 5.0 * valueR.get();
    if (blueColorChange > 2.0) {
      blueColorChange = 2.0;
    } 

    // Show different pattern according to the "pattern" value
    switch((int)pattern.get()) {
      case 0:
        g.pushMatrix();
        g.translate(-0.5, 0, 0);
        g.scale(sphereSize * valueL.get() * musicPower);
        //g.lighting(true);
        g.color(RGB(1.0, 1.0 - redColorChange * redColorChange, 1.0 - redColorChange * redColorChange));
        g.draw(pattern0SphereL);
        g.popMatrix();

        g.pushMatrix();
        g.translate(0.5, 0, 0);
        g.scale(sphereSize * valueR.get() * musicPower);
        //g.lighting(true);
        g.color(RGB(1.0 - blueColorChange * blueColorChange, 1.0, 1.0 - blueColorChange * blueColorChange));
        g.draw(pattern0SphereR);
        g.popMatrix();
        break;
      case 1:
        g.shader(pointShader);
        g.shader().uniform("pointSize", pointSize / 100);
        g.blending(true);
        g.blendTrans();
        g.depthTesting(true);
        g.draw(pattern1ParticleCylinder);
        break;
      case 2:
        break;
      case 3:
        break;
    }

  }

  // onAnimate
  // Decides the logic between the current frame to the next frame
  // For pattern 2, 3: calculate the new force excerted on each particle
  //
  void onAnimate(double dt) override {
    // Set the parameter of "pattern" to its floored-down int value
    int flooredPatternIndex = (int) (std::floor(pattern));
    pattern = flooredPatternIndex;

    // Pattern 1:
    // Deal with all the particle forces
    vector<Vec3f> &pattern1PositionVec(pattern1ParticleCylinder.vertices());
    for (int i = 0; i < pattern1PositionVec.size(); i++) {
      // Add spring force to particles' original positions
      Vec3f centerAtItsLayer = {0.0, pattern1PositionVec[i].y, 0.0};
      Vec3f particleToCenter = centerAtItsLayer - pattern1PositionVec[i];
      float distanceToSurface = particleToCenter.mag() - pattern1CylinderRadius;
      Vec3f springForce = particleToCenter.normalize() * springConstant * distanceToSurface;
      pattern1Force[i] += springForce;
    }
    
    // Exert the force and accelerations
    for (int i = 0; i < pattern1Velocity.size(); i++) {
      pattern1Force[i] += - pattern1Velocity[i] * dragFactor;

      // Important: add the force excerted by audio
      float musicForce = 0.0;
      
      int positionInSpectrum = std::floor(i / 50);
      //int positionInSpectrum = std::floor(((pattern1PositionVec[i].y + 0.6) / 1.2) * 2025);
      //int positionInSpectrum = (int) (FFT_SIZE * 0.5 * ((pattern1PositionVec[i].y + pattern1CylinderHalfLength) / (2 * pattern1CylinderHalfLength))) - 0.5 * FFT_SIZE;
      /*
      if (pattern1OriginalPosition[i].x < 0) {
        musicForce = spectrumL[positionInSpectrum];
      } else {
        musicForce = spectrumR[positionInSpectrum];
      }
      */
      musicForce = spectrum[positionInSpectrum];
      pattern1Force[i] += Vec3f(pattern1OriginalPosition[i].x, 0.0, pattern1OriginalPosition[i].z) * musicForce * musicPower;

      pattern1Velocity[i] += pattern1Force[i] / pattern1Mass[i] * timeStep;
      pattern1PositionVec[i] += pattern1Velocity[i] * timeStep;
    }


    // Clear all accelerations
    for (auto &a : pattern1Force) a.set(0);
    

  }

  void onMessage(osc::Message& m) override { m.print(); }
  bool onKeyDown(const Keyboard& k) override { 
    if (k.key() == '0') {
      pattern = 0;
    }
    if (k.key() == '1') {
      pattern = 1;
    }
    if (k.key() == '2') {
      pattern = 2;
    }
    if (k.key() == '3') {
      pattern = 3;
    }
    if (k.key() == 'u') {
      player.load("../tv.mp3");
    }
    if (k.key() == 'i') {
      player.load("../cd.mp3");
    }
    if (k.key() == 'o') {
      player.load("../lf.mp3");
    }
    return true; 
  }
  
};

// slurp
// To slurp from a file
string slurp(string fileName) {
  fstream file(fileName);
  string returnValue = "";
  while (file.good()) {
    string line;
    getline(file, line);
    returnValue += line + "\n";
  }
  return returnValue;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

int main() {
  MyApp app;
  app.configureAudio(48000, 512, 2, 2);
  app.start();
}

