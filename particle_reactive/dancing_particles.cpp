// Yuehao Gao | MAT201B
// 2022-02-23 | Final Project
// Audio-reactive visualizer
//

#include "al/app/al_App.hpp"
#include "al/app/al_GUIDomain.hpp"
//#include "al/ui/al_ControlGUI.hpp"
#include "al/graphics/al_Shapes.hpp"
#include "al/math/al_Functions.hpp"  // al::abs
#include "al/ui/al_Parameter.hpp"

#include "Gamma/SamplePlayer.h"  // XXX
#include "Gamma/Analysis.h"
#include "Gamma/Effects.h"
#include <cmath>
#include <fstream>
#include <vector>

using namespace std;
using namespace al;

const float sphereSize = 0.8;
const int pattern1CylinderHeight = 55;
const float pattern1CylinderHalfLength = 0.12;
const int pattern1CylinderNumParticlePerlayer = 30;
const float pattern1CylinderRadius = 0.06;
const float pattern1AngleIncrement = 2.0 * M_PI / pattern1CylinderNumParticlePerlayer;
const int pattern2BallNumParticle = 1000;

struct CommonState {
  //Vec3f position[17000];
  //int size;
  Nav nav;
  float pointSize;
};

string slurp(string fileName);  // forward declaration

struct MyApp : App {
  // All parameters
  Parameter valueL{"value_left", 0, 0, 1};
  Parameter valueR{"value_right", 0, 0, 1};
  Parameter pattern{"visual_pattern", 0, 0, 3};
  Parameter pointSize{"/pointSize", "", 1.0, 0.0, 2.0};

  // The sample player for 1 sound track
  // The envelop follower for the left and right channel
  gam::SamplePlayer<float, gam::ipl::Linear, gam::phsInc::Loop> player;
  gam::EnvFollow<> followLeft;
  gam::EnvFollow<> followRight;

  ShaderProgram pointShader;

  // ----------- All Meshes -------------
  // * Pattern 0: two spheres representing L and R envelop values
  Mesh pattern0SphereL;
  Mesh pattern0SphereR;

  // * Pattern 1:
  Mesh pattern1ParticleCylinder;

  // * Pattern 2:
  //Mesh pattern2ParticleBall;

  // * Pattern 3:

  // ------------------------------------

  void onCreate() override { 
    pointShader.compile(slurp("../point-vertex.glsl"),
                        slurp("../point-fragment.glsl"),
                        slurp("../point-geometry.glsl"));

    // ----------- Initialize Parameters for All Meshes -------------

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

        float hue = 120.0 / pattern1CylinderHeight * layerIndex;

        pattern1ParticleCylinder.vertex(Vec3f(x, y, z));
        pattern1ParticleCylinder.color(HSV(hue, 1.0f, 1.0f));    // THIS LINE IS NOT WORKING FOR SOME REASON
      }
    }
    

    // --------------------------------------------------------------
  
    // Set the magnitude of the envelop followers
    followLeft.lag(0.5); 
    followRight.lag(0.5);
  
    // Initialize the position of "camera"
    nav().pos(0, 0, 0.6);
  }

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
    
  }

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


  void onDraw(Graphics& g) override {
    g.clear(0.2);
    switch((int)pattern.get()) {
      case 0:
        g.pushMatrix();
        g.translate(-0.1, 0, 0);
        g.scale(sphereSize * valueL.get());
        g.draw(pattern0SphereL);
        g.popMatrix();

        g.pushMatrix();
        g.translate(0.1, 0, 0);
        g.scale(sphereSize * valueR.get()); 
        g.draw(pattern0SphereR);
        g.popMatrix();
        break;
      case 1:
        //g.shader(pointShader);   // THIS LINE IS HAVING SOME ISSUE
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

  void onAnimate(double dt) override {
    // Set the parameter of "pattern" to its floored-down int value
    int flooredPatternIndex = (int) (std::floor(pattern));
    pattern = flooredPatternIndex;

    vector<Vec3f> &positionVec(pattern1ParticleCylinder.vertices());
    

  }

  void onMessage(osc::Message& m) override { m.print(); }
  bool onKeyDown(const Keyboard& k) override { return false; }
  
};

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

int main() {
  MyApp app;
  app.configureAudio(48000, 512, 2, 2);
  app.start();
}

