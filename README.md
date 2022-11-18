# Link to the DAFx paper

https://dafx2020.mdw.ac.at/proceedings/papers/DAFx20in21_paper_38.pdf.

# Live concatenative granular processing

This algorithm addresses signal discontinuity and concatenation artefacts in real-time granular processing with rectangular windowing. By combining zero-crossing synchronicity, first-order derivative analysis, and Lagrange polynomials, we can generate streams of uncorrelated and non-overlapping sonic fragments with minimal low-order derivatives discontinuities. The resulting open-source algorithm, implemented in the Faust language, provides a versatile real-time software for dynamical looping, wavetable oscillation, and granulation with reduced artefacts due to rectangular windowing and no artefacts from overlap-add-to-one techniques commonly deployed in granular processing.

# Parameters

Interpolation length: length of the reconstructed segment at the junction through 5th-order Lagrange. Short lenghts may work best for high-frequency signals and vice versa.

Grain length: determines the approximate length of each grain in seconds, consequently setting the grain rate (1 / grain_length). Note that the grain length is approximate as the triggering of each grain is dependent on the zero-crossing occurrences in the output of the system.

Buffer position: offset to move the granulators reading head along the buffer, where 0 is the left-most area, and 1 is the right-most area. Positions 0 and 1 are equivalent as the buffer is circular. This parameter is particularly useful when using the granulator with a zero-factor time transposition, that is, as a wavetable oscillator to explore different waveforms.

Time transposition: time transposition factor; negative factors correspond to reversed buffer indexing.

Time async degree: degree of asynchronicity in the recursive chaotic factors variations.

Time async depth: depth of the oscillations in the recursive chaotic factors variations.

Pitch transposition: pitch transposition factor; negative factors correspond to reversed grain playback.

Pitch async degree: degree of asynchronicity in the recursive chaotic factors variations.

Pitch async depth: depth of the oscillations in the recursive chaotic factors variations.

Freeze buffer: check this box to prevent the buffer from being updated with new input signals.

Volume: linear scaling factor.

# Compilation

Please compile using double-precision.

# Optimality

For best results, run the software at high sample rates, i.e., 192 kHz.

# Acknowledgements

I want to thank Julian Parker for his valuable help in improving this algorithm.
