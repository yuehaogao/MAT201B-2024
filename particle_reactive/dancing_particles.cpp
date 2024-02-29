// Yuehao Gao | MAT201B
// 2022-02-23 | Final Project
// Audio-reactive visualizer
//

#include "al/app/al_App.hpp"
#include "al/app/al_GUIDomain.hpp"
//#include "al/ui/al_ControlGUI.hpp"
#include "al/graphics/al_Shapes.hpp"
#include "al/math/al_Functions.hpp"    // al::abs
#include "al/ui/al_Parameter.hpp"
#include "al/math/al_Random.hpp"

#include "Gamma/SamplePlayer.h"
#include "Gamma/Analysis.h"
#include "Gamma/Effects.h"
#include <cmath>
#include <fstream>
#include <vector>

using namespace std;
using namespace al;

const float sphereSize = 0.8;

const int pattern1CylinderHeight = 45;
const float pattern1CylinderHalfLength = 0.12;
const int pattern1CylinderNumParticlePerlayer = 40;
const float pattern1CylinderRadius = 0.06;
const float pattern1AngleIncrement = 2.0 * M_PI / pattern1CylinderNumParticlePerlayer;

const int pattern2BallNumParticle = 1000;

struct CommonState {
  Vec3f pattern1Position[17000];
  int size;
  Nav nav;
  float pointSize;
};


Vec3f randomVec3f(float scale) {
  return Vec3f(rnd::uniformS(), rnd::uniformS(), rnd::uniformS()) * scale;
}
string slurp(string fileName);  // forward declaration

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

struct MyApp : App {
  // ----------- All parameters -----------
  Parameter valueL{"value_left", 0, 0, 1};
  Parameter valueR{"value_right", 0, 0, 1};
  Parameter pattern{"visual_pattern", 0, 0, 3};
  Parameter pointSize{"/pointSize", "", 0.15, 0.05, 1.0};
  Parameter timeStep{"/timeStep", "", 0.0, 0.0, 1.0};
  Parameter dragFactor{"/dragFactor", "", 0.1, 0.0, 0.9};
  Parameter springConstant{"/springConstant", "", 0.3, 0.1, 1.0};
  Parameter repellingConstant{"/repellingConstant", "", 0.0, 0.0, 0.01};
  // ----------- All parameters -----------

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
  vector<float> pattern1Mass;

  // * Pattern 2:
  //Mesh pattern2ParticleBall;

  // * Pattern 3:

  // ----------- All Meshes -------------

  // onInit
  // When the app is first initialized
  //
  void onInit() override {
    player.load("../tv.mp3");

    auto GUIdomain = GUIDomain::enableGUI(defaultWindowDomain());
    auto& gui = GUIdomain->newGUI();
    gui.add(valueL);
    parameterServer() << valueL;
    gui.add(valueR);
    parameterServer() << valueR;
    gui.add(pattern);
    parameterServer() << pattern;
    gui.add(pointSize);
    parameterServer() << pointSize;  // Have this?
    gui.add(timeStep);                    // add parameter to GUI
    parameterServer() << timeStep;
    gui.add(dragFactor);                  // add parameter to GUI
    parameterServer() << dragFactor;
    gui.add(springConstant);              // add parameter to GUI
    parameterServer() << springConstant;
    //gui.add(sphereRadius);                // add parameter to GUI
    //parameterServer() << sphereRadius;
    gui.add(repellingConstant);           // add parameter to GUI
    parameterServer() << repellingConstant;
    
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

    // Pattern 1: 
    
    pattern1ParticleCylinder.primitive(Mesh::POINTS);
  
    for (int layerIndex = 0; layerIndex < pattern1CylinderHeight; layerIndex++) {
      float y = -1.0 * pattern1CylinderHalfLength + (2.0 * pattern1CylinderHalfLength / pattern1CylinderHeight) * layerIndex;
      for (int particleIndex = 0; particleIndex < pattern1CylinderNumParticlePerlayer; particleIndex++) {
        float angle = particleIndex * pattern1AngleIncrement;
        float x = pattern1CylinderRadius * cos(angle);
        float z = pattern1CylinderRadius * sin(angle);

        float hue = 0.7 / pattern1CylinderHeight * layerIndex;

        pattern1ParticleCylinder.vertex(Vec3f(x, y, z));
        pattern1ParticleCylinder.color(HSV(hue, 1.0f, 1.0f));    // THIS LINE IS NOT WORKING FOR SOME REASON
        pattern1Mass.push_back(3.0);
        pattern1ParticleCylinder.texCoord(pow(3.0, 1.0f / 3), 0);
        pattern1Velocity.push_back(Vec3f(0.0, 0.0, 0.0));
        pattern1Force.push_back(Vec3f(0.0, 0.0, 0.0));
      
      }
    }
    
    // ----------- Initialize Parameters for All Meshes Ends -------------
  
    // Set the magnitude of the envelop followers
    followLeft.lag(0.5); 
    followRight.lag(0.5);
  
    // Initialize the position of "camera"
    nav().pos(0, 0, 0.6);
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

      // https://github.com/adamstark/Gist
    }
  }

  // onDraw
  // Decide how each frame should be drawn
  // The current frame shows which mesh(es) is depending on the "pattern" parameter
  //
  void onDraw(Graphics& g) override {
    g.clear(0.6 * (valueL.get() + valueR.get()));
    float redColorChange = 5.0 * valueL.get();
    if (redColorChange > 2.0) {
      redColorChange = 2.0;
    }
    float blueColorChange = 5.0 * valueR.get();
    if (blueColorChange > 2.0) {
      blueColorChange = 2.0;
    } 
    switch((int)pattern.get()) {
      case 0:
        g.pushMatrix();
        g.translate(-0.1, 0, 0);
        g.scale(sphereSize * valueL.get() * 0.5);
        
        g.color(RGB(1.0, 1.0 - redColorChange * redColorChange, 1.0 - redColorChange * redColorChange));
        g.draw(pattern0SphereL);
        g.popMatrix();

        g.pushMatrix();
        g.translate(0.1, 0, 0);
        g.scale(sphereSize * valueR.get() * 0.5);
        
        g.color(RGB(1.0 - blueColorChange * blueColorChange, 1.0, 1.0 - blueColorChange * blueColorChange));
        g.draw(pattern0SphereR);
        g.popMatrix();
        break;
      case 1:
        g.shader(pointShader);   // THIS LINE IS HAVING SOME ISSUE
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
      //
    }
    

  }

  void onMessage(osc::Message& m) override { m.print(); }
  bool onKeyDown(const Keyboard& k) override { 
    if (k.key() == '1') {
      // introduce some "random" forces
      for (int i = 0; i < pattern1Velocity.size(); i++) {
        // F = ma
        pattern1Force[i] += randomVec3f(1);
      }

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

