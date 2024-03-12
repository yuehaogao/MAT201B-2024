// Yuehao Gao | MAT201B
// 2022-03-11 | Final Project
// Audio-reactive 3D visualizer

// GitHub: https://github.com/yuehaogao/MAT201B-2024_Yuehao_Gao
// Full Introduction: https://github.com/yuehaogao/MAT201B-2024_Yuehao_Gao/blob/main/README.md


// Introduction
/*
  This is my project implementing a 3D visualizer for input music files. 
  The project brings 2 (future: 4) distinct visual patterns that moves along the input audio signal.
  
  --- Express Keys ---
  [0], [1], [2]: changing visual patterns
  [u], [i], [o], [k], [l]: switch between music files

  --- Controllable Parameters ---
  VisualPattern (range 0 - 3, acquiscent: 1): 
    controlls the current visual pattern shown to follow the music. 
    This parameter will automatically be floored to the lower integer by the "onAnimate" fuction.

  MusicPower (range 0 - 8, acquiscent: 3): 
    how strong the music affects the visual patterns,
    you may understand this parameter as "music volume".

  SpringConstant (range 0.05 - 2, acquiscent: 0.3): 
    how stiff the spring is, 
    used for pattern {1, 2, 3} for calculating according to the Hooke's Law.

  Pattern1CylinderRadius (range 0.2 - 2, acquiscent: 0.75): 
    the radius of the cylinder, 
    used for pattern {1}.

  --- Visuals ---
  Background: 
    The background is acquiscently all-black, and turns white according to the enveloped music volume.
    It mimics a flashing-lighting environment of dancing halls or party venues.

  Pattern 0 (Express key [0]): 
    Two spheres seperately representing the envelopsed volume (NOT the spectrum of FFT) of the left and right channel separately. 
    When the input signal increases, the left sphere turns larger and "redder" from white,
    while the right sphere turns larger and "greener" from white. 
    The amount the two spheres change color and size are controlled by the "musicPower" parameter as well. 
    Each of the two spheres are seperately drawn with a distinct "Mesh" object in the Distributed App structure,
    respectively named "pattern0SphereL" and "pattern0SphereR".

  Pattern 1 (Express key [1]): 
    A SINGULAR hollow cylinder made by particles that dances according spectrum of the FFT of the input audio. 
    Specifically, from the bottom part to the top part of the particle cylinder, 
    the particles are gradually colored from red to blue (fixed color), 
    and dances according to the dB value from the lower frequencies to the higher frequencies on the STFT table. 
    The position that each particle from the bottom to the top that locate themselves on the spectrum is non-linear, 
    nor does the value of "musicForce" as higher frequencies have lower powers. 
    How ardently the particles dance could also be controlled by the "musicPower" parameter. 
    The sizes of the particles are also determined by the enveloped signal value of the music (average of left and right channel", 
    just like how the background works. 
    The cylinder is drawn by a SINGULAR Mesh called "pattern1ParticleCylinder".
*/



// ---- TO DO -----
// fix the minor bug for pattern 1: the particles should not float horizontally
// fix the 167 broken numbers
// implement pattern 2
// implement pattern 3
// make the music selection system by "dropdown list" instead of triple express keys



// All files, libraries, functions imported:
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

// Name spaces and fixed parameters
using namespace al;
using namespace std;
#define FFT_SIZE 4048

const float dragFactor = 0.15;
const float particleMass = 3.0;
const float sphereSize = 0.8;
const float timeStep = 1.0;

const int pattern1CylinderHeight = 56;
const int pattern1CylinderNumParticlePerlayer = 80;
const int pattern1NumParticle = pattern1CylinderHeight * pattern1CylinderNumParticlePerlayer;
const float pattern1AngleIncrement = 2.0 * M_PI / pattern1CylinderNumParticlePerlayer;
const float pattern1CylinderHalfLength = 0.5;
vector<Vec3f> pattern1OriginalPosition;



// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

// Information shared with the "Distributed App"
struct CommonState {
  // Public parameters for all patterns
  int pattern;                             // The visual pattern, either being 0, 1, 2, or 3
  Pose pose;                               // The position and angel of the "camera"
  float valueL;                            // Envelopsed right channel value
  float valueR;                            // Enveloped left channel value
  float musicPower;                        // The "volume" of music: how much the music affect visual patterns
  float pointSize;                         // The size of the particles
  float springConstant;                    // The spring constant for Hooke's Law
  float spectrum[FFT_SIZE / 2 + 100];      // The added number is arbitrarily decided

  // Pattern 1 specific parameters
  float pattern1CylinderRadius;            // The radius of the cylinder
  Vec3f pattern1RealTimePosition[5000];    // The updated (at this moment) position of each particle on the cylinder
  HSV pattern1FixedColors[5000];           // The fixed color (HSV) of each particle on the cylinder
  
  
  // Pattern 2 specific parameters
  // Pattern 3 specific parameters
  
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
  Parameter pointSize{"/pointSize", "", 0.5, 0.1, 1.5};
  Parameter pattern{"visual_pattern", 1, 0, 3};     // This will be limited to int only
  Parameter musicPower{"/musicPower", "", 0.0, 0.0, 8.0};
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
  Vec3f pattern1Velocity[pattern1NumParticle];
  Vec3f pattern1Force[pattern1NumParticle];

  // * Pattern 2:

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

    }

  }


  // onCreate
  // Create all the initial visual patterns
  // Slurp corresponding files for the point shader program
  // Set envelop followers
  // Initialize the position of camera
  //
  void onCreate() override { 
    bool createPointShaderSuccess = pointShader.compile(slurp("../point-vertex.glsl"),
                          slurp("../point-fragment.glsl"),
                          slurp("../point-geometry.glsl"));

    if (!createPointShaderSuccess) {
      exit(1);
    }


    // ----------- Initialize Parameters for All Meshes Begins -------------

    // Pattern 0: move the two spheres seperately to the left and right


    // Pattern 1: create all the particles with their initial positions, colors, velocity, and forces
      

    // Here, "i" means the index of the article, from 0 to the last one
    int i = 0;

    for (int layerIndex = 0; layerIndex < pattern1CylinderHeight; layerIndex++) {
      float y = -1.0 * pattern1CylinderHalfLength + (2.0 * pattern1CylinderHalfLength / pattern1CylinderHeight) * layerIndex;
      
      for (int particleIndex = 0; particleIndex < pattern1CylinderNumParticlePerlayer; particleIndex++) {
        float angle = particleIndex * pattern1AngleIncrement;
        float x = pattern1CylinderRadius * cos(angle);
        float z = pattern1CylinderRadius * sin(angle);

        // Place the particle's position
        pattern1ParticleCylinder.vertex(Vec3f(x, y, z));        // Push the position into the mesh
        pattern1OriginalPosition.push_back(Vec3f(x, y, z));     // Local archive
        state().pattern1RealTimePosition[i] = Vec3f(x, y, z);   // Common real-time position

        // Color the particle
        float hue = 0.7 / pattern1CylinderHeight * layerIndex;
        pattern1ParticleCylinder.color(HSV(hue, 1.0f, 1.0f));
        state().pattern1FixedColors[i] = HSV(hue, 1.0f, 1.0f);       

        // Set the particle's physical force system
        pattern1ParticleCylinder.texCoord(pow(particleMass, 1.0f / 3), 0);
        pattern1Velocity[i] = Vec3f(0.0, 0.0, 0.0);
        pattern1Force[i] = Vec3f(0.0, 0.0, 0.0);
        i++;


      }
    }

    //cout << "Invalid Value Initialized: " << nanInitialized << endl;

    // Pattern 2: ...
    // Pattern 3: ...
    
    // ----------- Initialize Parameters for All Meshes Ends -------------


    // Add the shapes and points
    // These should be sent to BOTH primary and secondary sites
    addSphere(pattern0SphereL);
    addSphere(pattern0SphereR);
    pattern1ParticleCylinder.primitive(Mesh::POINTS);

    // LOCAL-ONLY initialization process
    if (isPrimary()) {
      // Set the magnitude of the envelop followers
      followLeft.lag(0.5); 
      followRight.lag(0.5);
  
      // Initialize the position of "camera"
      nav().pos(0, 0, 0.0);
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
    g.clear(0.5 * pow(valueLAndR, 1.25) * state().musicPower);
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
      // Update the camera position
      state().pose = nav();

      // Set the parameter of "pattern" to its floored-down int value
      // Unify local data with "state" data (only pattern)
      int flooredPatternIndex = (int) (std::floor(pattern));
      pattern = flooredPatternIndex;
      state().pattern = pattern;
      
      // Pattern 1: --------------------------------------------------------------------------
      // Deal with all the particle forces


      vector<Vec3f> &pattern1PositionVec(pattern1ParticleCylinder.vertices());
      for (int i = 0; i < pattern1PositionVec.size(); i++) {

        // Add spring force to particles
        // According to their original positions and current positions
        Vec3f centerAtItsLayer = {0.0, pattern1PositionVec[i].y, 0.0};
        Vec3f particleToCenter = centerAtItsLayer - pattern1PositionVec[i];
        float distanceToSurface = particleToCenter.mag() - pattern1CylinderRadius;
        Vec3f springForce = particleToCenter.normalize() * springConstant * distanceToSurface;
        pattern1Force[i] += springForce;

        // Add the reverse the force caused by acceleration
        pattern1Force[i] -= pattern1Velocity[i] * dragFactor;
      

        // Important: add the force excerted by audio
        int positionInSpectrum = std::floor(pow(i, 1.8) / 10000);
        //int positionInSpectrum = std::floor(((pattern1PositionVec[i].y + 0.6) / 1.2) * 2025);
        //int positionInSpectrum = (int) (FFT_SIZE * 0.5 * ((pattern1PositionVec[i].y + pattern1CylinderHalfLength) / (2 * pattern1CylinderHalfLength))) - 0.5 * FFT_SIZE;

        float musicForce = state().spectrum[positionInSpectrum];
        float indexForceBoostLimit = 2.5;
        float fftForce = ((pow(i, 1.5)) / 100000.0) * musicForce;
        if (fftForce > indexForceBoostLimit ) {
          fftForce = indexForceBoostLimit;
        }
        pattern1Force[i] -= Vec3f(pattern1OriginalPosition[i].x, 0.0, pattern1OriginalPosition[i].z) * fftForce * state().musicPower;
        pattern1Velocity[i] += pattern1Force[i] / particleMass * timeStep;

        // Exert the force and accelerations
        pattern1PositionVec[i] += pattern1Velocity[i] * timeStep;        // Move the particle locally
        
        

        // if (isnan(musicForce)) {
        //     cout << "NAN MUSIC FORCE is at particle indexed: " << i << endl;
        //   }
         // if (isnan(positionInSpectrum)) {
            //cout << "NAN POS IN SPEC is at particle indexed: " << i << endl;
         // }
        
          // if (isnan(pattern1PositionVec[i].x)) {
          //   cout << "NAN X is at particle indexed: " << i << endl;
          // }
          // if (isnan(pattern1PositionVec[i].y)) {
          //   cout << "NAN Y is at particle indexed: " << i << endl;
          // }
          // if (isnan(pattern1PositionVec[i].z)) {
          //   cout << "NAN Z is at particle indexed: " << i << endl;
          // }
        
      }

      

      pointSize = musicPower * 0.5 * (valueL + valueR);

      // Clear all accelerations
      for (auto &a : pattern1Force) a.set(0);

      // Unify other parameters between local and "state"
      state().valueL = valueL;
      state().valueR = valueR;
      state().pointSize = pointSize;
      state().musicPower = musicPower;
      state().springConstant = springConstant;
      state().pattern1CylinderRadius = pattern1CylinderRadius;

      for (int i = 0; i < pattern1NumParticle; i++) {
        state().pattern1RealTimePosition[i] = pattern1PositionVec[i];    // Tell the distributed app about it
      }
   
    } else {
      nav().set(state().pose);
      // Unify the position and hue of distributed app
      // Clear the mesh from previous frame
      pattern1ParticleCylinder.vertices().clear();
      pattern1ParticleCylinder.colors().clear();

      int nanupdated = 0;  // USED FOR DEBUGGING

      for (int i = 0; i < pattern1NumParticle; i++) {
        // Update each particle's position and 
        pattern1ParticleCylinder.vertex(state().pattern1RealTimePosition[i]);
        pattern1ParticleCylinder.color(state().pattern1FixedColors[i]);

        if (isnan(state().pattern1RealTimePosition[i].x)) {
          nanupdated++; // USED FOR DEBUGGING
        }
        
        
      }
      //cout << nanupdated << endl;
    }

  }

  // onMessage
  // When the system receives an input message
  //
  void onMessage(osc::Message& m) override { m.print(); }

  
  // onKeyDown
  // React to any key pressed by the user
  //
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
    if (k.key() == 'k') {
      player.load("../dc.mp3");
    }
    if (k.key() == 'l') {
      player.load("../it.mp3");
    }
    return true; 
  }
  
};


// slurp
// To slurp from a file
//
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


