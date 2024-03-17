
## Dancing Particles - Yuehao Gao, University of California, Santa Barbara

## MAT201B, Professor Karl Yerkes, March 2024

This is the repository of my project implementing a 3D visualizer for input music files. The project brings 4 distinct visual patterns that moves along the input audio signal. There is a controllable sliding parameter named "pattern" that allows changing the visual parameter from one to another. You may also use express keys (number keys: [0], [1], [2], [3]) to change visual patterns. Additionally, express keys ([4] and [r]) will start or pause randomly-changing visual patterns.

You may also press on [p] to locate yourself at (0, 0, 5) outside the mesh and press on [`] to get back to the center.

The current version allows switching between input pieces with the [u], [i], [o], [k], [l] express key. In the future, it will allow a dropdown list to choose between multiple input music files.

The major c++ file dealing with all the algorithms is the "dancing_particles.cpp", and it takes in the one of the audio files (limited to be located in the same folder) at a time using the "player()" object in the Distributed App structure. The input file will be played in loops, and limited to ".wav" only. Specifically, the input audio is acquiescently to be in stereo, meaning having two channels. In the function named "onSound", two sound samples are grabbed each time, and sent to the left and right channels separately. If not, then each sample will be sent to both the left and right channels at the same time.

------------------------------ Introduction of Controllable Parameters ------------------------------
- VisualPattern (range 0-3, acquiscent: 1): controlls the current visual pattern shown to follow the music. This parameter will automatically be floored to the lower integer by the "onAnimate" fuction.

- MusicPower (range 0-6, acquiscent: 2.5): how strong the music affects the visual patterns, you may understand this parameter as "music volume".

- SpectrumWidth(range 0.25-1.5, acquiscent: 1.0): from the fixed lowest frequency value on the spectrum, how large is the range of frequency-amplitude values this pattern grabs from the spectrum. The larger the number is, the narrower the "dancing bass" will be, but you shall see more particles on the top dancing to higher frequencies.

- TimeStep(range 0.25-2.0, acquiscent: 1.0): how fast the particles move.

- SpringConstant (range 0.05-2, acquiscent: 0.75): how stiff the spring is, used for pattern {1, 2, 3} for calculating according to the Hook's Law.

- Radius (range 0.2-2, acquiscent: 1): the radius of all the patterns. Changing this will affect every pattern of {0, 1, 2, 3}.

- Distance (range 0.5-4, acquiscent: 2): the distance between the mesh for left and right channels. Changing this will only affect pattern {0, 3}.


------------------------------ Introduction of Uncontrollable Parameters ------------------------------
- ValueL: the enveloped value of the left input channel
  
- ValueR: the enveloped value of the right input channel
  
- PointSize: the size of the particles, which are automatically controlled by the magnitude of input signal and the parameter "musicPower".


------------------------------ Introduction of Background and Visual Patterns ------------------------------

- Background: the background is acquiscently all-black, and turns white according to the enveloped music volume, which mimics a flashing-lighting environment of dancing halls or party venues.

- Pattern 0 (Express key [0]): two spheres seperately representing the envelopsed volume (NOT the spectrum of FFT) of the left and right channel separately. When the input signal increases, the left sphere turns larger and "redder" from white, while the right sphere turns larger and "greener" from white. The amount the two spheres change color and size are controlled by the "musicPower" parameter as well. Each of the two spheres are seperately drawn with a distinct "Mesh" object in the Distributed App structure, respectively named "pattern0SphereL" and "pattern0SphereR". 

- Pattern 1 (Express key [1]): a SINGULAR hollow cylinder made by particles that dances according to the spectrum of the FFT of the input audio. Specifically, from the bottom part to the top part of the particle cylinder, the particles are gradually colored from red to blue (fixed color), and dances according to the dB value from the lower frequencies to the higher frequencies on the STFT table. The position that each particle from the bottom to the top that locate themselves on the spectrum is non-linear, nor does the value of "musicForce" as higher frequencies have lower powers. How ardently the particles dance could also be controlled by the "musicPower" parameter. The sizes of the particles are also determined by the enveloped signal value of the music (average of left and right channel", just like how the background works. The cylinder is drawn by a SINGULAR Mesh called "pattern1ParticleCylinder".

- Pattern 2 (Express key [2]): a SINGULAR hollow sphere made by particles that dances according to the spectrum of the FFT of the input audio. The lower to higher part of this mesh also dance differently according to different values on the spectrum table. The particles will be pushed towards outward by the music force. It may also be controlled by musicPower. The sphere is drawn by a SINGULAR Mesh called "pattern2SphereCylinder".
  
- Pattern 3 (Express key [3]): a PAIR OF hollow spheres made by particles that dances separately according to the left and right channel spectrum of the FFT of the input audio. The lower to higher part of these two mesh also dance differently according to different values on the left and right spectrum table. The particles will be pushed towards outward by the music force. They may also be controlled by musicPower. The spheres are drawn by a PAIR OF Meshes separately called "pattern3LSphereCylinder" and "pattern3RSphereCylinder".



Thank you very much!
