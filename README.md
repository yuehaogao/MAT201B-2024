
Dancing Particles - Yuehao Gao, University of California, Santa Barbara

MAT201B, Professor Karl Yerkes, March 2024

This is the repository of my project implementing a 3D visualizer for input music files. The project brings 2 (future: 4) distinct visual patterns that moves along the input audio signal. There is a controllable sliding parameter named "pattern" that allows changing the visual parameter from one to another. You may also use express keys (number keys: [0], [1] (future: [2], [3])) to change visual patterns.

The current version allows switching between the 3 input pieces with the [u], [i], [o] express key. In the future, it will allow a dropdown list to choose between multiple input music files.

The major c++ file dealing with all the algorithms is the "dancing_particles.cpp", and it takes in the one of the audio files (limited to be located in the same folder) at a time using the "player()" object in the Distributed App structure. The input file will be played in loops, and limited to ".mp3" or ".wav" only. Specifically, the input audio is acquiescently to be in stereo, meaning having two channels. In the function named "onSound", two sound samples are grabbed each time, and sent to the left and right channels separately. If not, then each sample will be sent to both the left and right channels at the same time.

------------------------------ Introduction of Background and Visual Patterns ------------------------------

- Background: the background is acquiscently all-black, and turns white according to the enveloped music volume, which mimics a flashing-lighting environment of dancing halls or party venues.

- Pattern 0 (Express key [0]): two spheres seperately representing the envelopsed volume (NOT the spectrum of FFT) of the left and right channel separately. When the input signal increases, the left sphere turns larger and "redder" from white, while the right sphere turns larger and "greener" from white. The amount the two spheres change color and size are controlled by the "musicPower" parameter as well. Each of the two spheres are seperately drawn with a distinct "Mesh" object in the Distributed App structure, respectively named "pattern0SphereL" and "pattern0SphereR". 

- Pattern 1 (Express key [1]): a SINGULAR hollow cylinder made by particles that dances according spectrum of the FFT of the input audio. Specifically, from the bottom part to the top part of the particle cylinder, the particles are gradually colored from red to blue (fixed color), and dances according to the dB value from the lower frequencies to the higher frequencies on the STFT table. The position that each particle from the bottom to the top that locate themselves on the spectrum is non-linear, nor does the value of "musicForce" as higher frequencies have lower powers. How ardently the particles dance could also be controlled by the "musicPower" parameter. The sizes of the particles are also determined by the enveloped signal value of the music (average of left and right channel", just like how the background works. The cylinder is drawn by a SINGULAR Mesh called "pattern1ParticleCylinder".


--------- Future patterns ---------
- Pattern 2 (Express key [2]): two pans of particles that dances according to the FFT of the left and right channel signal *(hypothesized, to be continued...)*

- Pattern 3 (Express key [3]): two spheres of particles that dances according to the FFT of the left and right channel signal *(hypothesized, to be continued...)*


------------------------------ Introduction of Controllable Parameters ------------------------------
- VisualPattern (range 0-3, acquiscent: 1): controlls the current visual pattern shown to follow the music. This parameter will automatically be floored to the lower integer by the "onAnimate" fuction.

- MusicPower (range 0-8, acquiscent: 3): how strong the music affects the visual patterns, you may understand this parameter as "music volume".

- SpringConstant (range 0.05-2, acquiscent: 0.3): how stiff the spring is, used for pattern {1, 2, 3} for calculating according to the Hook's Law.

- Pattern1CylinderRadius (range 0.2-2, acquiscent: 0.75): the radius of the cylinder, used for pattern {1}.


--------- Future controllable parameters ---------
- Pattern2Size (TBD): the size of the two particle disks, used for pattern {2}.

- Pattern3Radius (TBD): the radius of the two particle spheres, used for pattern {3}.


Thank you very much!
