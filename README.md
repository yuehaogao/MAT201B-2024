
Dancing Particles - Yuehao Gao (高悦淏)

University of California, Santa Barbara

MAT201B (2024), Professor Karl Yerkes

This is the repository of my project implementing a 3D visualizer of input music files. The project brings [2-4] distinct visual patterns that moves along the input audio signal. There is a sliding parameter controller named "pattern" that allows changing the visual parameter from one to another. You may also use express keys (number keys: [0[, [1] (future: [2], [3])) to change visual patterns.

The major file dealing with all the algorithms is the "dancing_particles.cpp", and it takes in the one of the audio files (limited to be located in the same folder) using the "player()" object in the Distributed App structure. The input file will be played in loops, and limited to ".mp3" or ".wav" only. Specifically, the input audio is acquiescently having two channels, and two sound samples are grabbed each time, and sent to the left and right channels separately. If not, then each sample will be sent to both the left and right channels at the same time. The change of pattern can be controlled by a sliding
parameter, such that the rounded-down integer will tell the application about which mesh should
be displayed at this time.

The first pattern is two spherical balls, separated with four times their initial diameters regarding
their centers, that respectively visualize the current value of each channel. Specifically, this
pattern does not utilize the spectrum analysis (FFT) of the input audio.
The following patterns are all based on “repelled particles” and dance according to the frequency
spectrum. 

The second pattern is a hollow cylinder of particles: those toward the lower part mean lower frequencies,
and vice versa. They are exerted with pushing forces that distort them from their original stringed
positions according to the magnitude on different parts of the spectrum. So, if a DJ music with
strong bass is imported, one can see the bottom part dancing much more vigorously than the
upper part.

The third pattern is heavily based on Assignment 3, in which we have hollow
particle sphere. The upper part of the sphere (particles with 0 or positive y values) represent the
left channel, and vice versa. From the left to right, the particles are moved by forces with the
magnitude according to the “dB” value respective to the frequency spectrum, from lower
frequencies to higher frequency. Hence, if a DJ music with strong bass is imported, one can see
the left part dancing much more vigorously than the right part.


With this being said, the fourth one is also similar, but this time the particles are aligned in a 2D
plate, with the outer particles representing lower frequencies, and inner particles representing
higher frequencies on the spectrum. An up-and-down distortion force will be exerted accordingly
to the spectrum. Music with stronger beats will make the outer rings dance more vigorously than
the inner rings.
In order to acquire the FFT analysis and the spectrum data to feed into the second, third and
fourth patterns, I am going to utilize the functions inside the “integrated.cpp” file inside the
“tutorial” folder.
