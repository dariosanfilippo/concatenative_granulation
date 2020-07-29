# zero-crossing_granulator
Faust implementation of zero-crossing-synchronous rectangular windowing for continuous streams of non-overlapping sonic fragments. The main characteristic of the study is to provide a multi-modal system that, depending on the parameters, can operate as a click-free looper, wavetable oscillator, or granulator.

The "Looped/Live" button switches between live and looped inputs. The granulator is looped is the button is unchecked, live otherwise.

The "Grain rate" parameter determines the number of sonic fragments per second or, conversely, the size of each fragment. Please note that the grain rate or size is only desirable: in most cases, the length of each segment is slightly longer as the beginning of the next grain is dependent on the zero-crossing occurrence in the previous grain. Depending on the grain rate, whether below or above the perceptual lower frequency threshold, the algorithm operates in the looping or wavetable oscillatiion regions. This also requires the time-stretching parameter to be set to zero. (See below.)

The "Pitch factor" parameter can be both positive or negative. Negative pitch factors result in reversed grains. It determines the playback speed of each grain.

The "Pitch modulation" parameter creates an exponential curve in the delay shift, which in turn results in a pitch modulation for each grain. A constant pitch is given by a delay shift whose second derivative is zero. In this case, the parameter ranges from -1 to 1, mapped over the range fro 1/16 to 16, which is the exponent that shapes the line segment performing the delay shift. Positive values of the parameter result in a positive second derivative, or a function curving up creating a raising pitch. Negative values of the parameter behave in the opposite way.

The "Time factor" parameter determines the speed and direction for the exploration of the buffer. This parameter, provided that the position modulation parameters are inactive, is the time-stretching parameter. Particularly, when the parameter is zero, the algorithm operates as a looper or wavetable oscillator.

The "Buffer region" parameter is can be used to explore different regions of the buffer when the time-stretching factor is set to zero, or it can be used for manual time-stretching.

The "Position self-modulation depth" sets the magnitude of the perturbation over the curve that determines the time factor. Specifically, the output amplitude of the system is sent back, expanded through this parameter, wrapped around, and summed to the time factor function. The nonlinear iteration process can be used to achieve chaotic beahviours.

Lastly, the "Position self-modulation rate" lowpasses the output of the granulator before being processed in the position self-modulation depth parameter. By changing the rate of this signal, we can smoothly transition between correlated to uncorrelated positions in the buffer to move from time-stratching-like behaviours to sonic dusts.
