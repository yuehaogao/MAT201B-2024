// Yuehao Gao | MAT201B
// 2022-03-17 | Final Project
// Audio-reactive 3D visualizer

// GitHub: https://github.com/yuehaogao/MAT201B-2024_Yuehao_Gao
// Full Introduction: https://github.com/yuehaogao/MAT201B-2024_Yuehao_Gao/blob/main/README.md


// Introduction
/*
  This is my project implementing a 3D visualizer for input music files. 
  The project brings 2 (future: 4) distinct visual patterns that moves along the input audio signal.
  
  --- Express Shortcut Keys ---
  [0], [1], [2]: changing visual patterns
  [u], [i], [o], [k], [l]: switch between music files
  [p]: proceed to the position (0.0, 0.0, 5.0)

  --- Controllable Parameters ---
  VisualPattern (range 0 - 3, acquiscent: 1): 
    controlls the current visual pattern shown to follow the music. 
    This parameter will automatically be floored to the lower integer by the "onAnimate" fuction.

  MusicPower (range 0 - 6, acquiscent: 3): 
    how strong the music affects the visual patterns,
    you may understand this parameter as "music volume".
    used for pattern {0, 1, 2, 3}.

  SpringConstant (range 0.05 - 2, acquiscent: 0.6): 
    how stiff the spring is, 
    used for pattern {1, 2, 3} for calculating according to the Hooke's Law.

  Radius (range 0.2 - 2, acquiscent: 1): 
    the radius of all the patterns, 
    used for pattern {0, 1, 2, 3}.

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
// FIX: spectrum L R issue
// enable auto-changing pattern
// express p: rotate to original as well
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

const int pattern1CylinderHeight = 40;
const int pattern1CylinderNumParticlePerlayer = 112;
const int pattern1NumParticle = pattern1CylinderHeight * pattern1CylinderNumParticlePerlayer;
const float pattern1AngleIncrement = 2.0 * M_PI / pattern1CylinderNumParticlePerlayer;
const float pattern1CylinderHalfLength = 1.0;
const float pattern1RadiusIncrement = 1.2;

const int pattern2NumParticle = 1000;
const float pattern2RadiusIncrement = 1.1;

const int pattern3NumParticle = 300;          // This is the amount of particle for EACH mesh
const float pattern3RadiusIncrement = 0.5;

int frameCount = 0;                           // Count the frame
const int frameChangeThreshold = 100;        // How many frame to change a pattern automatically
bool randomPattern = false;                   // Initially, the app should not display random pattern
                                              // Press [p] or [4] to switch between


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

// Information shared with the "Distributed App"
struct CommonState {
  // Public parameters for all patterns
  Pose pose;                               // The position and angel of the "camera"
  int pattern;                             // The visual pattern, either being 0, 1, 2, or 3
  float valueL;                            // Envelopsed right channel value
  float valueR;                            // Enveloped left channel value
  float musicPower;                        // The "volume" of music: how much the music affect visual patterns
  float pointSize;                         // The size of the particles
  float timeStep;                          // How fast the particles proceed
  float springConstant;                    // The spring constant for Hooke's Law
  float spectrum[FFT_SIZE / 2 + 100];      // The added number is arbitrarily decided (The spectrum for both channels)
  float spectrumL[FFT_SIZE / 2 + 100];     // The spectrum for the left channel only
  float spectrumR[FFT_SIZE / 2 + 100];     // The spectrum for the right channel only
  float spectrumWidth;                     // The range of spectrum that the visual patterns grab data from
  float radius;                            // The radius of the cylinder / single sphere / double sphere
  float distance;                          // The distance between the mesh for left and right channel

  // Pattern 1 specific parameters
  Vec3f pattern1RealTimePosition[5000];    // The updated (at this moment) position of each particle on the cylinder
  HSV pattern1FixedColors[5000];           // The fixed color (HSV) of each particle on the cylinder
  
  // Pattern 2 specific parameters
  float pattern2RealTimePointSize[1050];   // The updated (at this moment) particle size on the sphere
  Vec3f pattern2RealTimePosition[1050];    // The updated (at this moment) position of each particle on the sphere
  HSV pattern2RealTimeColors[1050];        // The updated (at this moment) color of each particle on the sphere

  // Pattern 3 specific parameters
  float pattern3LRealTimePointSize[350];
  float pattern3RRealTimePointSize[350];
  Vec3f pattern3LRealTimePosition[350];
  Vec3f pattern3RRealTimePosition[350];
  HSV pattern3LRealTimeColors[350];
  HSV pattern3RRealTimeColors[350];
  
};


// To generate a random 3D vector for position, color, ...
Vec3f randomVec3f(float scale) {
  return Vec3f(rnd::uniformS(), rnd::uniformS(), rnd::uniformS()) * scale;
}

Vec3f randomVec3fWithCenter(float scale, float x, float y, float z) {
  return Vec3f(rnd::uniformS() * scale + x, rnd::uniformS() * scale+ y, rnd::uniformS() * scale + z);
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
  Parameter pattern{"visual_pattern", 1.0, 0.0, 3.9};      // This will be limited to int only
  Parameter musicPower{"/musicPower", "", 2.5, 0.0, 6.0};
  Parameter spectrumWidth{"/spectrumWidth", "", 1.0, 0.25, 1.5};
  Parameter timeStep{"/timeStep", "", 1.0, 0.25, 2.0};
  Parameter springConstant{"/springConstant", "", 0.75, 0.05, 2.0};
  Parameter radius{"/radius", "", 1.0, 0.2, 2.0};
  Parameter distance{"/distance", "", 2.0, 0.5, 4.0};


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
  float pattern1CylinderRadius;

  // * Pattern 2:
  Mesh pattern2ParticleSphere;
  Vec3f pattern2Velocity[pattern2NumParticle];
  Vec3f pattern2Force[pattern2NumParticle];
  float pattern2SphereRadius;

  // * Pattern 3:
  Mesh pattern3LParticleSphere;
  Vec3f pattern3LVelocity[pattern3NumParticle];
  Vec3f pattern3LForce[pattern3NumParticle];
  Mesh pattern3RParticleSphere;
  Vec3f pattern3RVelocity[pattern3NumParticle];
  Vec3f pattern3RForce[pattern3NumParticle];
  float pattern3SphereRadius;


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
      player.load("../wav_files/Dancin.wav");
    
      // set up GUI
      auto GUIdomain = GUIDomain::enableGUI(defaultWindowDomain());
      auto& gui = GUIdomain->newGUI();
      gui.add(valueL);
      state().valueL = valueL;
      gui.add(valueR);
      state().valueR = valueR;
      gui.add(pointSize);
      state().pointSize = pointSize;
      gui.add(pattern);
      state().pattern = pattern;
      gui.add(musicPower);
      state().musicPower = musicPower;
      gui.add(spectrumWidth);
      state().spectrumWidth = spectrumWidth;
      gui.add(timeStep);
      state().timeStep = timeStep;
      gui.add(springConstant);
      state().springConstant = springConstant;
      gui.add(radius);
      state().radius = radius;
      gui.add(distance);
      state().distance = distance;

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

    // * Pattern 0: move the two spheres seperately to the left and right, executed in the "onDraw" function

    // * Pattern 1: create all the particles with their initial positions, colors, velocity, and forces
      
    // Here, "p1i" means the index of the particle, from 0 to the last one
    pattern1CylinderRadius = radius * pattern1RadiusIncrement;
    int p1i = 0;

    for (int layerIndex = 0; layerIndex < pattern1CylinderHeight; layerIndex++) {
      float y = -1.0 * pattern1CylinderHalfLength + (2.0 * pattern1CylinderHalfLength / pattern1CylinderHeight) * layerIndex;
      
      for (int particleIndex = 0; particleIndex < pattern1CylinderNumParticlePerlayer; particleIndex++) {
        float angle = particleIndex * pattern1AngleIncrement;
        float x = pattern1CylinderRadius * cos(angle);
        float z = pattern1CylinderRadius * sin(angle);

        // Place the particle's position
        pattern1ParticleCylinder.vertex(Vec3f(x, y, z));          // Push the position into the mesh
        state().pattern1RealTimePosition[p1i] = Vec3f(x, y, z);   // Common real-time position

        // Color the particle
        float hue = 0.7 / pattern1CylinderHeight * layerIndex;
        pattern1ParticleCylinder.color(HSV(hue, 1.0f, 1.0f));
        state().pattern1FixedColors[p1i] = HSV(hue, 1.0f, 1.0f);       

        // Set the particle's physical force system
        pattern1ParticleCylinder.texCoord(pow(particleMass, 1.0f / 3), 0);
        pattern1Velocity[p1i] = Vec3f(0.0, 0.0, 0.0);
        pattern1Force[p1i] = Vec3f(0.0, 0.0, 0.0);

        // Finally, increment p1i by 1, proceed to the next particle in the mesh
        p1i++;

      }
    }


    // * Pattern 2: 

    // Here, "p2i" means the index of the particle, from 0 to the last one
    // Initialize the position and color of each particle
    pattern2SphereRadius = radius * pattern2RadiusIncrement;
    for (int p2i = 0; p2i < pattern2NumParticle; p2i++) {

      // Place the particle's position (initialized with random position)
      Vec3f randomInitPositionOfThisArticle = randomVec3f(pattern2SphereRadius);
      pattern2ParticleSphere.vertex(randomInitPositionOfThisArticle);
      state().pattern2RealTimePosition[p2i] = randomInitPositionOfThisArticle;

      // Color the particle
      float p2hue = 0.7 * ((randomInitPositionOfThisArticle.y + pattern2SphereRadius) / (2.0 * pattern2SphereRadius));
      pattern2ParticleSphere.color(HSV(p2hue, 1.0f, 1.0f));
      state().pattern2RealTimeColors[p2i] = HSV(p2hue, 1.0f, 1.0f);

      // Set the particle's physical force system
      pattern2ParticleSphere.texCoord(pow(particleMass, 1.0f / 3), 0);
      pattern2Velocity[p2i] = Vec3f(0.0, 0.0, 0.0);
      pattern2Force[p2i] = Vec3f(0.0, 0.0, 0.0);

    }


    // Pattern 3:
    pattern3SphereRadius = radius * pattern3RadiusIncrement;


    for (int p3li = 0; p3li < pattern3NumParticle; p3li++) {

      // Place the particle's position (initialized with random position)
      Vec3f randomInitPositionOfThisParticle = randomVec3fWithCenter(pattern3SphereRadius, distance * -0.5, 0.0, -2.0);
      pattern3LParticleSphere.vertex(randomInitPositionOfThisParticle);
      state().pattern3LRealTimePosition[p3li] = randomInitPositionOfThisParticle;

      // Color the particle
      float p3hue = 0.7 * ((randomInitPositionOfThisParticle.y + pattern3SphereRadius) / (2.0 * pattern3SphereRadius));
      pattern3LParticleSphere.color(HSV(p3hue, 1.0f, 1.0f));
      state().pattern3LRealTimeColors[p3li] = HSV(p3hue, 1.0f, 1.0f);

      // Set the particle's physical force system
      pattern3LParticleSphere.texCoord(pow(particleMass, 1.0f / 3), 0);
      pattern3LVelocity[p3li] = Vec3f(0.0, 0.0, 0.0);
      pattern3LForce[p3li] = Vec3f(0.0, 0.0, 0.0);

    }

    for (int p3ri = 0; p3ri < pattern3NumParticle; p3ri++) {

      // Place the particle's position (initialized with random position)
      Vec3f randomInitPositionOfThisParticle = randomVec3fWithCenter(pattern3SphereRadius, distance * 0.5, 0.0, -2.0);
      pattern3RParticleSphere.vertex(randomInitPositionOfThisParticle);
      state().pattern3RRealTimePosition[p3ri] = randomInitPositionOfThisParticle;

      // Color the particle
      float p3hue = 0.7 * ((randomInitPositionOfThisParticle.y + pattern3SphereRadius) / (2.0 * pattern3SphereRadius));
      pattern3RParticleSphere.color(HSV(p3hue, 1.0f, 1.0f));
      state().pattern3RRealTimeColors[p3ri] = HSV(p3hue, 1.0f, 1.0f);

      // Set the particle's physical force system
      pattern3RParticleSphere.texCoord(pow(particleMass, 1.0f / 3), 0);
      pattern3RVelocity[p3ri] = Vec3f(0.0, 0.0, 0.0);
      pattern3RForce[p3ri] = Vec3f(0.0, 0.0, 0.0);

    }


    
    // ----------- Initialize Parameters for All Meshes Ends -------------


    // Add the shapes and points
    // These should be sent to BOTH primary and secondary sites
    addSphere(pattern0SphereL);
    addSphere(pattern0SphereR);
    pattern1ParticleCylinder.primitive(Mesh::POINTS);
    pattern2ParticleSphere.primitive(Mesh::POINTS);
    pattern3LParticleSphere.primitive(Mesh::POINTS);
    pattern3RParticleSphere.primitive(Mesh::POINTS);

    // LOCAL-ONLY initialization process
    if (isPrimary()) {
      // Set the magnitude of the envelop followers
      followLeft.lag(0.5); 
      followRight.lag(0.5);
  
      // Initialize the position of "camera"
      nav().pos(0.0, 0.0, 0.0);
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
      float sRight = player.read(0);
      float sLeft = player.read(1);
      
      io.out(0) = sLeft;
      io.out(1) = sRight;

      valueL = (followLeft(sLeft));
      valueR = (followRight(sRight));



      // STFT
      // Spectrum for right channel
      if (stft(io.out(0))) { 
        // Loop through all the frequency bins
        for (unsigned k = 0; k < stft.numBins(); ++k) {
          state().spectrum[k] = 8.0 * tanh(pow(stft.bin(k).real(), 1.5));
          state().spectrumR[k] = 8.0 * tanh(pow(stft.bin(k).real(), 1.5));
        }
      }

      // Spectrum for left channel
      if (stft(io.out(1))) { 
        // Loop through all the frequency bins
        for (unsigned k = 0; k < stft.numBins(); ++k) {
          state().spectrum[k] = 8.0 * tanh(pow(stft.bin(k).real(), 1.5));
          state().spectrumL[k] = 8.0 * tanh(pow(stft.bin(k).real(), 1.5));
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
    g.clear(0.3 * pow(valueLAndR, 0.9) * state().musicPower);
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
        g.translate(-0.5 * distance, 0, -2.0);
        g.scale(sphereSize * radius * state().valueL * state().musicPower * 1.5);
        //g.lighting(true);
        g.color(RGB(1.0, 1.0 - pow(redColorChange, 2.0), 1.0 - pow(redColorChange, 2.0)));
        g.draw(pattern0SphereL);
        g.popMatrix();

        g.pushMatrix();
        g.translate(0.5 * distance, 0, -2.0);
        g.scale(sphereSize * radius * state().valueR * state().musicPower * 1.5);
        //g.lighting(true);
        g.color(RGB(1.0 - pow(greenColorChange, 2.0), 1.0, 1.0 - pow(greenColorChange, 2.0)));
        g.draw(pattern0SphereR);
        g.popMatrix();
        break;

      case 1:
        g.shader(pointShader);
        g.shader().uniform("pointSize", state().pointSize / 100);
        g.blending(true);
        //g.lighting(true);
        g.blendTrans();
        g.depthTesting(true);
        g.draw(pattern1ParticleCylinder);
        break;

      case 2:
        g.shader(pointShader);
        g.shader().uniform("pointSize", state().pointSize / 100);
        g.blending(true);
        //g.lighting(true);
        g.blendTrans();
        g.depthTesting(true);
        g.draw(pattern2ParticleSphere);
        break;

      case 3:
        //g.pushMatrix();
        g.shader(pointShader);
        g.shader().uniform("pointSize", state().pointSize / 100);
        g.blending(true);
        //g.lighting(true);
        g.blendTrans();
        g.depthTesting(true);
        g.draw(pattern3LParticleSphere);
        g.draw(pattern3RParticleSphere);

        break;
    }

  }


  // onAnimate
  // Decides the logic between the current frame to the next frame
  // For pattern 2, 3: calculate the new force excerted on each particle
  //
  void onAnimate(double dt) override { 

    frameCount++;  // used for debugging only

    if (isPrimary()) {
      // Update the camera position
      state().pose = nav();

      // Set the parameter of "pattern" to its floored-down int value
      // Unify local data with "state" data (only pattern)
      int flooredPatternIndex = (int) (std::floor(pattern));
      pattern = flooredPatternIndex;
      

      // If currently in random pattern, change
      int previousPattern = -1;
      if (randomPattern) {
        if (frameCount < frameChangeThreshold) {
          frameCount++;
        } else {
          frameCount = 0;
          int randNextPattern = rand() % 4;
          while(randNextPattern == previousPattern) {
            int randNextPattern = rand() % 4;
          }
          pattern = randNextPattern;
        }
      } else {
        frameCount = 0;
      }

      // Unify the visual pattern index
      state().pattern = pattern;


      // Unify the radius
      pattern1CylinderRadius = pattern1RadiusIncrement * radius;
      pattern2SphereRadius = pattern2RadiusIncrement * radius;
      pattern3SphereRadius = pattern3RadiusIncrement * radius;


      // Unify other parameters between local and "state"
      state().valueL = valueL;
      state().valueR = valueR;
      state().pointSize = pointSize;
      state().musicPower = musicPower;
      state().spectrumWidth = spectrumWidth;
      state().timeStep = timeStep;
      state().springConstant = springConstant;
      state().radius = radius;
      state().distance = distance;

      
      // Pattern 1: --------------------------------------------------------------------------
      // Deal with all the particle forces
      vector<Vec3f> &pattern1PositionVec(pattern1ParticleCylinder.vertices());
      
      for (int i = 0; i < pattern1PositionVec.size(); i++) {
        // Add spring force to particles
        // According to their original positions and current positions
        Vec3f centerAtItsLayer = {0.0, pattern1PositionVec[i].y, 0.0};
        Vec3f particleToCenter = centerAtItsLayer - pattern1PositionVec[i];
        float distanceToSurface = particleToCenter.mag() - pattern1CylinderRadius;
        Vec3f inNOutspringForce = particleToCenter.normalize() * springConstant * distanceToSurface;
        pattern1Force[i] += inNOutspringForce;

        // Add the reverse the force caused by acceleration
        pattern1Force[i] -= pattern1Velocity[i] * dragFactor;
      
        // Important: add the force excerted by audio
        int positionInSpectrum = 3 + std::floor(pow(i, 1.85) / std::floor(28000.0 / spectrumWidth));

        float musicForce = state().spectrum[positionInSpectrum];
        float indexForceBoostLimit = 5.0;
        float fftForce = ((pow(i, 1.5)) / 100000.0) * musicForce;
        if (fftForce > indexForceBoostLimit ) {
          fftForce = indexForceBoostLimit;
        }

        float bassBoost = 3.0 - 2.0 * ((pattern1PositionVec.size() - i) / pattern1PositionVec.size());

        
        pattern1Force[i] -= Vec3f(pattern1PositionVec[i].x, 0.0, pattern1PositionVec[i].z) * fftForce * state().musicPower * bassBoost;
        pattern1Velocity[i] += pattern1Force[i] / particleMass * timeStep;

        // Exert the force and accelerations
        // Move the particle locally
        pattern1PositionVec[i] += pattern1Velocity[i] * timeStep; 

      }


      // Pattern 2: --------------------------------------------------------------------------
      // Deal with all the particle forces
      vector<Vec3f> &pattern2PositionVec(pattern2ParticleSphere.vertices());

      // First clear the sphere particles' color
      // They will be updated in real-time on every frame
      pattern2ParticleSphere.colors().clear();

      for (int i = 0; i < pattern2NumParticle; i++) {
        // Add spring force to particles
        // According to their original positions and current positions
        Vec3f particleToCenter = Vec3f(0.0, 0.0, 0.0) - pattern2PositionVec[i];
        float distanceToSurface = particleToCenter.mag() - pattern2SphereRadius;
        Vec3f inNOutspringForce = particleToCenter.normalize() * (springConstant) * distanceToSurface;
        pattern2Force[i] += inNOutspringForce;

        // For every other particle in the mesh system
        // Calculating the repelling force excerted within each other
        for (int j = i + 1; j < pattern2NumParticle; j++) {
          Vec3f displacement = pattern2PositionVec[j] - pattern2PositionVec[i];
          float distanceSquared = displacement.magSqr();
          Vec3f repellingForce = displacement.normalize() * (0.00025 / distanceSquared);
          pattern2Force[i] -= repellingForce;
          pattern2Force[j] += repellingForce;
        }

        pattern2Force[i] -= pattern2Velocity[i] * dragFactor;
        
        // Important: add the force excerted by audio
        // "particleHeightIndex" is a float number that is SUPPOSED to range from 0 to 10
        // which indicates the relavant "height" the particle is in the entire sphere
        // if the particle exceeds the sphere range, the number will be limited
        float particleHeightIndex = (pattern2PositionVec[i].y + pattern2SphereRadius) * 10.0 / (2.0 * pattern2SphereRadius);
        if (particleHeightIndex < 0.0) {
          particleHeightIndex = 0.0;
        }
        if (particleHeightIndex > 10.0) {
          particleHeightIndex = 10.0;
        }

        int positionInSpectrum = 5 + std::floor(pow(particleHeightIndex, 4.0) / std::floor(15.0 / spectrumWidth)); 
        float musicForce = state().spectrum[positionInSpectrum];
        float fftForce = musicForce * 0.1;
    
        float indexForceBoostLimit = 0.01;
        if (fftForce > indexForceBoostLimit) {
          fftForce = indexForceBoostLimit;
        }

        float tweetBoost = 1.5 + 8.5 * (particleHeightIndex / 10);
        pattern2Force[i] += Vec3f(pattern2PositionVec[i].x, pattern2PositionVec[i].y, pattern2PositionVec[i].z) * fftForce * tweetBoost * state().musicPower;

        pattern2Velocity[i] += pattern2Force[i] / particleMass * (timeStep * 0.66);
        pattern2PositionVec[i] += pattern2Velocity[i] * timeStep;
        state().pattern2RealTimePosition[i] = pattern2PositionVec[i];

        // Meanwhile, change its color according to its new position
        float yIndexForColor = state().pattern2RealTimePosition[i].y + pattern2SphereRadius;
        if (yIndexForColor < 0.0) {
          yIndexForColor = 0.0;
        }
        if (yIndexForColor > 2.0 * pattern2SphereRadius) {
          yIndexForColor = 2.0 * pattern2SphereRadius;
        }

        float newHue = 0.7 * (yIndexForColor / (2.0 * pattern2SphereRadius));
        state().pattern2RealTimeColors[i] = HSV(newHue, 1.0, 1.0);
        pattern2ParticleSphere.color(HSV(newHue, 1.0, 1.0));


      }

      // Pattern 3L&R: --------------------------------------------------------------------------
      // Deal with all the particle forces
      vector<Vec3f> &pattern3LPositionVec(pattern3LParticleSphere.vertices());
      vector<Vec3f> &pattern3RPositionVec(pattern3RParticleSphere.vertices());

      // First clear the sphere particles' color
      // They will be updated in real-time on every frame
      pattern3LParticleSphere.colors().clear();
      pattern3RParticleSphere.colors().clear();


      // Deal with each particle in the left mesh
      for (int i = 0; i < pattern3NumParticle; i++) {
        // Add spring force to particles
        // According to their original positions and current positions
        Vec3f particleToCenter = Vec3f(-0.5 * distance, 0.0, -2.0) - pattern3LPositionVec[i];
        //Vec3f particleToCenter = Vec3f(0.0, 0.0, 0.0) - pattern3LPositionVec[i];
        float distanceToSurface = particleToCenter.mag() - pattern3SphereRadius;
        Vec3f inNOutspringForce = particleToCenter.normalize() * (springConstant) * distanceToSurface;
        pattern3LForce[i] += inNOutspringForce * 0.75;

        // For every other particle in the mesh system
        // Calculating the repelling force excerted within each other
        for (int j = i + 1; j < pattern3NumParticle; j++) {
          Vec3f displacement = pattern3LPositionVec[j] - pattern3LPositionVec[i];
          float distanceSquared = displacement.magSqr();
          Vec3f repellingForce = displacement.normalize() * (0.0002 / distanceSquared);
          pattern3LForce[i] -= repellingForce;
          pattern3LForce[j] += repellingForce;
        }

        pattern3LForce[i] -= pattern3LVelocity[i] * dragFactor;
        
        // Important: add the force excerted by audio
        // "particleHeightIndex" is a float number that is SUPPOSED to range from 0 to 10
        // which indicates the relavant "height" the particle is in the entire sphere
        // if the particle exceeds the sphere range, the number will be limited
        float particleHeightIndex = (pattern3LPositionVec[i].y + pattern3SphereRadius) * 10.0 / (2.0 * pattern3SphereRadius);
        if (particleHeightIndex < 0.0) {
          particleHeightIndex = 0.0;
        }
        if (particleHeightIndex > 10.0) {
          particleHeightIndex = 10.0;
        }

        int positionInSpectrum = 4 + std::floor(pow(particleHeightIndex, 4.0) / std::floor(12.5 / spectrumWidth));
        float musicForce = state().spectrumL[positionInSpectrum];
        float fftForce = musicForce * 0.2;
    
        float indexForceBoostLimit = 0.01;
        if (fftForce > indexForceBoostLimit) {
          fftForce = indexForceBoostLimit;
        }

        float tweetBoost = 1.5 + 8.5 * (particleHeightIndex / 10);

        // CHANGE THIS
        pattern3LForce[i] += Vec3f(pattern3LPositionVec[i].x + (0.5 * distance), pattern3LPositionVec[i].y, pattern3LPositionVec[i].z + 2.0) * fftForce * tweetBoost * state().musicPower;

        pattern3LVelocity[i] += pattern3LForce[i] / particleMass * (timeStep * 0.66);
        pattern3LPositionVec[i] += pattern3LVelocity[i] * timeStep;
        state().pattern3LRealTimePosition[i] = pattern3LPositionVec[i];

        // Meanwhile, change its color according to its new position
        float yIndexForColor = state().pattern3LRealTimePosition[i].y + pattern3SphereRadius;
        if (yIndexForColor < 0.0) {
          yIndexForColor = 0.0;
        }
        if (yIndexForColor > 2.0 * pattern3SphereRadius) {
          yIndexForColor = 2.0 * pattern3SphereRadius;
        }

        float newHue = 0.7 * (yIndexForColor / (2.0 * pattern3SphereRadius));
        state().pattern3LRealTimeColors[i] = HSV(newHue, 1.0, 1.0);
        pattern3LParticleSphere.color(HSV(newHue, 1.0, 1.0));

      }

      
      // Deal with each particle in the right mesh
      for (int i = 0; i < pattern3NumParticle; i++) {
        // Add spring force to particles
        // According to their original positions and current positions
        Vec3f particleToCenter = Vec3f(0.5 * distance, 0.0, -2.0) - pattern3RPositionVec[i];
        float distanceToSurface = particleToCenter.mag() - pattern3SphereRadius;
        Vec3f inNOutspringForce = particleToCenter.normalize() * (springConstant) * distanceToSurface;
        pattern3RForce[i] += inNOutspringForce * 0.75;

        // For every other particle in the mesh system
        // Calculating the repelling force excerted within each other
        for (int j = i + 1; j < pattern3NumParticle; j++) {
          Vec3f displacement = pattern3RPositionVec[j] - pattern3RPositionVec[i];
          float distanceSquared = displacement.magSqr();
          Vec3f repellingForce = displacement.normalize() * (0.0002 / distanceSquared);
          pattern3RForce[i] -= repellingForce;
          pattern3RForce[j] += repellingForce;
        }

        pattern3RForce[i] -= pattern3RVelocity[i] * dragFactor;
        
        // Important: add the force excerted by audio
        // "particleHeightIndex" is a float number that is SUPPOSED to range from 0 to 10
        // which indicates the relavant "height" the particle is in the entire sphere
        // if the particle exceeds the sphere range, the number will be limited
        float particleHeightIndex = (pattern3RPositionVec[i].y + pattern3SphereRadius) * 10.0 / (2.0 * pattern3SphereRadius);
        if (particleHeightIndex < 0.0) {
          particleHeightIndex = 0.0;
        }
        if (particleHeightIndex > 10.0) {
          particleHeightIndex = 10.0;
        }

        int positionInSpectrum = 4 + std::floor(pow(particleHeightIndex, 4.0) / std::floor(12.5 / spectrumWidth));
        float musicForce = state().spectrum[positionInSpectrum];
        float fftForce = musicForce * 0.2;
    
        float indexForceBoostLimit = 0.01;
        if (fftForce > indexForceBoostLimit) {
          fftForce = indexForceBoostLimit;
        }

        float tweetBoost = 1.5 + 8.5 * (particleHeightIndex / 10);

        pattern3RForce[i] += Vec3f(pattern3RPositionVec[i].x + (-0.5 * distance), pattern3RPositionVec[i].y, pattern3RPositionVec[i].z + 2.0) * fftForce * tweetBoost * state().musicPower;

        pattern3RVelocity[i] += pattern3RForce[i] / particleMass * (timeStep * 0.66);
        pattern3RPositionVec[i] += pattern3RVelocity[i] * timeStep;
        state().pattern3RRealTimePosition[i] = pattern3RPositionVec[i];

        // Meanwhile, change its color according to its new position
        float yIndexForColor = state().pattern3RRealTimePosition[i].y + pattern3SphereRadius;
        if (yIndexForColor < 0.0) {
          yIndexForColor = 0.0;
        }
        if (yIndexForColor > 2.0 * pattern3SphereRadius) {
          yIndexForColor = 2.0 * pattern3SphereRadius;
        }

        float newHue = 0.7 * (yIndexForColor / (2.0 * pattern3SphereRadius));
        state().pattern3RRealTimeColors[i] = HSV(newHue, 1.0, 1.0);
        pattern3RParticleSphere.color(HSV(newHue, 1.0, 1.0));


      }


      pointSize = 0.1 + musicPower * 0.3 * pow((valueL + valueR), 0.4);

      // Clear all forces for pattern 1, 2, 3
      for (auto &a : pattern1Force) a.set(0);
      for (auto &a : pattern2Force) a.set(0);
      for (auto &a : pattern3LForce) a.set(0);
      for (auto &a : pattern3RForce) a.set(0);


      // Tell the distributed app about the refresh of every particle
      for (int i = 0; i < pattern1NumParticle; i++) {
        state().pattern1RealTimePosition[i] = pattern1PositionVec[i];    
      }
      for (int i = 0; i < pattern2NumParticle; i++) {
        state().pattern2RealTimePosition[i] = pattern2PositionVec[i];    
      }
      for (int i = 0; i < pattern3NumParticle; i++) {
        state().pattern3LRealTimePosition[i] = pattern3LPositionVec[i]; 
        //state().pattern3RRealTimePosition[i] = pattern3RPositionVec[i];   
      }
   
    } else {
      // Unify the position and hue of distributed app
      // Clear the mesh from previous frame
      nav().set(state().pose);

      pattern1ParticleCylinder.vertices().clear();
      pattern1ParticleCylinder.colors().clear();
      for (int i = 0; i < pattern1NumParticle; i++) {
        // Update each particle's position and color
        pattern1ParticleCylinder.vertex(state().pattern1RealTimePosition[i]);
        pattern1ParticleCylinder.color(state().pattern1FixedColors[i]);
      }

      pattern2ParticleSphere.vertices().clear();
      pattern2ParticleSphere.colors().clear();
      for (int i = 0; i < pattern2NumParticle; i++) {
        // Update each particle's position and color
        pattern2ParticleSphere.vertex(state().pattern2RealTimePosition[i]);
        pattern2ParticleSphere.color(state().pattern2RealTimeColors[i]);
      }

      pattern3LParticleSphere.vertices().clear();
      pattern3LParticleSphere.colors().clear();
      pattern3RParticleSphere.vertices().clear();
      pattern3RParticleSphere.colors().clear();
      for (int i = 0; i < pattern3NumParticle; i++) {
        // Update each particle's position and color
        pattern3LParticleSphere.vertex(state().pattern3LRealTimePosition[i]);
        pattern3LParticleSphere.color(state().pattern3LRealTimeColors[i]);
        pattern3RParticleSphere.vertex(state().pattern3RRealTimePosition[i]);
        pattern3RParticleSphere.color(state().pattern3RRealTimeColors[i]);
      }

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
    if (k.key() == 'p') {
      nav().pos(0.0, 0.0, 5.0);
      // QUESTION: HOW TO MAKE THE ROTATION STRAIGHT?
    }
    if (k.key() == 'r' || k.key() == '4') {
      randomPattern = !randomPattern;
    }

    return true; 
  }

  // onExit
  // Clear all meshes by removing all vertices and colors
  //
  void onExit() override {
    pattern1ParticleCylinder.vertices().clear();
    pattern2ParticleSphere.vertices().clear();
    pattern3LParticleSphere.vertices().clear();
    pattern3RParticleSphere.vertices().clear();

    pattern1ParticleCylinder.colors().clear();
    pattern2ParticleSphere.colors().clear();
    pattern3LParticleSphere.colors().clear();
    pattern3RParticleSphere.colors().clear();
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
  // MyApp app;
  // app.configureAudio(48000, 512, 2, 2);
  // app.start();

  MyApp app;
  if (al::Socket::hostName() == "ar01.1g") {
    AudioDevice device = AudioDevice("ECHO X5");
    app.configureAudio(device, 44100, 128, device.channelsOutMax(), 2);
  } else {
    app.configureAudio(48000, 512, 2, 2);
  }
  app.start();  
}



