// 03-10-2024

// Yuehao Gao | MAT201B
// 2022-02-23 | Final Project
// Audio-reactive visualizer
//

#include "al/app/al_DistributedApp.hpp"
#include "al/app/al_GUIDomain.hpp"
#include "al_ext/statedistribution/al_CuttleboneDomain.hpp"
#include "al_ext/statedistribution/al_CuttleboneStateSimulationDomain.hpp"
#include "al/graphics/al_Shapes.hpp"
#include "al/math/al_Functions.hpp"
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


const float dragFactor = 0.15;
const float particleMass = 3.0;
const float sphereSize = 0.8;
const float timeStep = 1.0;

const int pattern1CylinderHeight = 56;
const int pattern1CylinderNumParticlePerlayer = 80;
const float pattern1AngleIncrement = 2.0 * M_PI / pattern1CylinderNumParticlePerlayer;
const float pattern1CylinderHalfLength = 0.5;
vector<Vec3f> pattern1OriginalPosition;

/*
const int pattern2MostInnerLayerNumParticle = 20;
const int pattern2NumLayers = 30;
const float pattern2Displacement = 0.5;
const float pattern2LayerRadiusIncrement = 1.0 / 3.0;
*/

// Common-used statements or parameters
struct CommonState {
  // Public parameters for all patterns
  int pattern;
  float valueL;
  float valueR;
  float musicPower;
  float pointSize;
  float springConstant;
  float spectrum[FFT_SIZE / 2 + 500];

  // Pattern 1 specific parameters
  float pattern1CylinderRadius;
  Vec3f pattern1Velocity[5000];
  Vec3f pattern1Force[5000];
  //Vec3f pattern1OriginalPosition[5000];   // Archive of original positions of particles
  
  
  // Pattern 2 specific parameters
  //Vec3f pattern2Position[5000];
  //float sphereRadius;
  //float repellingConstant;

};

// To generate a random 3D vector for position, color, ...
Vec3f randomVec3f(float scale) {
  return Vec3f(rnd::uniformS(), rnd::uniformS(), rnd::uniformS()) * scale;
}

// To slurp a file (initialized)
string slurp(string fileName);  // forward declaration



// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

struct MyApp : DistributedAppWithState<CommonState> {
  // ----------- All parameters Start -----------
  Parameter valueL{"value_left", 0, 0, 1};
  Parameter valueR{"value_right", 0, 0, 1};
  Parameter pattern{"visual_pattern", 1, 0, 3};     // This will be limited to int only
  Parameter musicPower{"/musicPower", "", 3.0, 0.1, 8.0};
  Parameter pointSize{"/pointSize", "", 0.5, 0.1, 1.5};
  Parameter springConstant{"/springConstant", "", 0.3, 0.05, 2.0};
  Parameter pattern1CylinderRadius{"/cylinderRadius", "", 0.75, 0.2, 2.0};
  // Parameter sphereRadius{"/sphereRadius", "", 2.0, 0.5, 4.0};
  // Parameter repellingConstant{"/repellingConstant", "", 0.0, 0.0, 0.01};


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

  // * Pattern 2:
  // Mesh pattern2ParticleBall;
  // Vec3f pattern2SphereCenter = {0.0, 0.0, 0.0};

  // * Pattern 3:
  // ------------------------------------


  // onInit
  // When the app is first initialized
  //
  void onInit() override {

    // Try starting the program
    auto cuttleboneDomain =
        CuttleboneStateSimulationDomain<CommonState>::enableCuttlebone(this);
    if (!cuttleboneDomain) {
      std::cerr << "ERROR: Could not start Cuttlebone. Quitting." << std::endl;
      quit();
    }
    

    if (isPrimary()) {
      // Load the sound file
      player.load("../tv.mp3");
    
      // set up GUI
      auto GUIdomain = GUIDomain::enableGUI(defaultWindowDomain());
      auto& gui = GUIdomain->newGUI();
      gui.add(valueL);
      state().valueL = valueL;
      gui.add(valueR);
      state().valueR = valueR;
      gui.add(pattern);
      state().pattern = pattern;
      gui.add(musicPower);
      state().musicPower = musicPower;
      gui.add(pointSize);
      state().pointSize = pointSize;
      gui.add(springConstant);
      state().springConstant = springConstant;
      gui.add(pattern1CylinderRadius);
      state().pattern1CylinderRadius = pattern1CylinderRadius;
      // gui.add(sphereRadius);
      // gui.add(repellingConstant);

      // Declare the size of the spectrums
      //state().spectrum.resize(FFT_SIZE / 2 + 1);
    }

  }

  // onCreate
  // Create all the initial visual patterns
  // Slurp corresponding files for the point shader program
  // Set envelop followers
  // Initialize the position of camera
  //
  void onCreate() override { 
    if (isPrimary()) {
      pointShader.compile(slurp("../point-vertex.glsl"),
                          slurp("../point-fragment.glsl"),
                          slurp("../point-geometry.glsl"));

      // ----------- Initialize Parameters for All Meshes Begins -------------

      // Pattern 0: move the two spheres seperately to the left and right
      addSphere(pattern0SphereL);
      addSphere(pattern0SphereR);

      // Pattern 1: create all the particles with their initial positions, colors, velocity, and forces
      pattern1ParticleCylinder.primitive(Mesh::POINTS);

      // Here, "i" means the index of the article, from 0 to the last one
      int i = 0;
  
      for (int layerIndex = 0; layerIndex < pattern1CylinderHeight; layerIndex++) {
        float y = -1.0 * pattern1CylinderHalfLength + (2.0 * pattern1CylinderHalfLength / pattern1CylinderHeight) * layerIndex;
      
        for (int particleIndex = 0; particleIndex < pattern1CylinderNumParticlePerlayer; particleIndex++) {
          float angle = particleIndex * pattern1AngleIncrement;
          float x = pattern1CylinderRadius * cos(angle);
          float z = pattern1CylinderRadius * sin(angle);

          // Place the particle's position
          pattern1ParticleCylinder.vertex(Vec3f(x, y, z));
          //state().pattern1OriginalPosition[i] = Vec3f(x, y, z);  // Archive of original positions of particles
          pattern1OriginalPosition.push_back(Vec3f(x, y, z));

          // Color the particle
          float hue = 0.7 / pattern1CylinderHeight * layerIndex;
          pattern1ParticleCylinder.color(HSV(hue, 1.0f, 1.0f));

          // Set the particle's physical force system
          pattern1ParticleCylinder.texCoord(pow(particleMass, 1.0f / 3), 0);
          state().pattern1Velocity[i] = Vec3f(0.0, 0.0, 0.0);
          state().pattern1Force[i] = Vec3f(0.0, 0.0, 0.0);
        }
      }
    
      // Pattern 2: ...
      // Pattern 3: ...
    
      // ----------- Initialize Parameters for All Meshes Ends -------------
  
      // Set the magnitude of the envelop followers
      followLeft.lag(0.5); 
      followRight.lag(0.5);
  
      // Initialize the position of "camera"
      nav().pos(0, 0, 5.0);
    }
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

      valueL = (followLeft(sLeft));
      valueR = (followRight(sRight));

      // STFT
      if (stft(io.out(0))) { 
        // Loop through all the frequency bins
        for (unsigned k = 0; k < stft.numBins(); ++k) {
          state().spectrum[k] = 8.0 * tanh(pow(stft.bin(k).real(), 1.5));
        }
      }

      // https://github.com/adamstark/Gist
    }
  }

  // onDraw
  // Decide how each frame should be drawn
  // The current frame shows which mesh(es) is depending on the "pattern" parameter
  //
  void onDraw(Graphics& g) override {
    // The background changes brightness according to the volume
    float valueLAndR = state().valueL + state().valueR;
    g.clear(1.2 * pow(valueLAndR, 2.0) * state().musicPower);
    float redColorChange = 5.0 * state().valueL;
    if (redColorChange > 2.0) {
      redColorChange = 2.0;
    }
    float greenColorChange = 5.0 * state().valueR;
    if (greenColorChange > 2.0) {
      greenColorChange = 2.0;
    } 

    // Show different pattern according to the "pattern" value
    switch((int)(state().pattern)) {
      case 0:
        g.pushMatrix();
        g.translate(-1.0, 0, 0);
        g.scale(sphereSize * state().valueL * state().musicPower * 1.5);
        //g.lighting(true);
        g.color(RGB(1.0, 1.0 - pow(redColorChange, 2.0), 1.0 - pow(redColorChange, 2.0)));
        g.draw(pattern0SphereL);
        g.popMatrix();

        g.pushMatrix();
        g.translate(1.0, 0, 0);
        g.scale(sphereSize * state().valueR * state().musicPower * 1.5);
        //g.lighting(true);
        g.color(RGB(1.0 - pow(greenColorChange, 2.0), 1.0, 1.0 - pow(greenColorChange, 2.0)));
        g.draw(pattern0SphereR);
        g.popMatrix();
        break;
      case 1:
        g.shader(pointShader);
        g.shader().uniform("pointSize", state().pointSize / 100);
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
    if (isPrimary()) {
      // Set the parameter of "pattern" to its floored-down int value
      int flooredPatternIndex = (int) (std::floor(pattern));
      pattern = flooredPatternIndex;

      // Unify local data with "state" data (only pattern)
      state().pattern = pattern;
      
      // Pattern 1:
      // Deal with all the particle forces
      vector<Vec3f> &pattern1PositionVec(pattern1ParticleCylinder.vertices());
      for (int i = 0; i < pattern1PositionVec.size(); i++) {
        // Add spring force to particles' original positions
        Vec3f centerAtItsLayer = {0.0, pattern1PositionVec[i].y, 0.0};
        Vec3f particleToCenter = centerAtItsLayer - pattern1PositionVec[i];
        float distanceToSurface = particleToCenter.mag() - pattern1CylinderRadius;
        Vec3f springForce = particleToCenter.normalize() * springConstant * distanceToSurface;
        state().pattern1Force[i] += springForce;
      
    
        // Exert the force and accelerations
      
        state().pattern1Force[i] -= state().pattern1Velocity[i] * dragFactor;

        // Important: add the force excerted by audio
        float musicForce = 0.0;
      
        int positionInSpectrum = std::floor(pow(i, 1.8) / 10000);
        //int positionInSpectrum = std::floor(((pattern1PositionVec[i].y + 0.6) / 1.2) * 2025);
        //int positionInSpectrum = (int) (FFT_SIZE * 0.5 * ((pattern1PositionVec[i].y + pattern1CylinderHalfLength) / (2 * pattern1CylinderHalfLength))) - 0.5 * FFT_SIZE;

        musicForce = state().spectrum[positionInSpectrum];
        float indexForceBoostLimit = 3.0;
        float fftForce = ((pow(i, 1.5)) / 100000.0) * musicForce;
        if (fftForce > indexForceBoostLimit ) {
          fftForce = indexForceBoostLimit;
        }
        state().pattern1Force[i] -= Vec3f(pattern1OriginalPosition[i].x, 0.0, pattern1OriginalPosition[i].z) * fftForce * state().musicPower;

        state().pattern1Velocity[i] += state().pattern1Force[i] / particleMass * timeStep;
        pattern1PositionVec[i] += state().pattern1Velocity[i] * timeStep;
      }

      pointSize = musicPower * 0.5 * (valueL + valueR);

      // Clear all accelerations
      for (auto &a : state().pattern1Force) a.set(0);

      // Unify other parameters between local and "state"
      state().valueL = valueL;
      state().valueR = valueR;
      state().pointSize = pointSize;
      state().musicPower = musicPower;
      state().springConstant = springConstant;
      state().pattern1CylinderRadius = pattern1CylinderRadius;

    }

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


