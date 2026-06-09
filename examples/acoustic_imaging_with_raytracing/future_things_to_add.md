What you're currently extracting:
    - Per-cell speed of sound — from transducer-to-transducer travel times via tomographic inversion
    - Per-cell 3D velocity — from self-to-self pulse-echo cross-correlation + multi-beam least-squares
    - Per-cell attenuation — from transducer-to-transducer return volume → gas_holdup × bubble_radii
    - Per-cell bubble acceleration — from the time derivative of velocity → constrains drag_coefficient

What you're leaving on the table:

Low-hanging fruit:
    - Backscatter amplitude profile. You're cross-correlating the pulse-echo pairs to get velocity, but the envelope amplitude at each depth window tells you the local scatterer reflectivity, which is proportional to gas_holdup * bubble_cross_section. You already have the windowed data in extract_velocities_from_data — extracting abs(hilbert(window)) or just rms(window) alongside the velocity would give you per-cell backscatter intensity essentially for free. This is an independent constraint on gas_holdup that doesn't come through the transducer-to-transducer path.
    
    - Spectral broadening / turbulence intensity. The width of the cross-correlation peak (not just the peak location) contains information about velocity variance within the measurement volume. A sharp peak means uniform flow; a broad peak means turbulent mixing. You could fit a Gaussian to the peak region and extract the standard deviation as a turbulence metric per cell.

Medium effort:
    - Frequency-dependent attenuation → bubble size distribution. Right now your attenuation model treats bubbles as having a single bubble_radii per cell. But acoustic attenuation vs. frequency follows different scaling laws depending on bubble size relative to wavelength (ka << 1 thermal/viscous regime vs. ka ~ 1 resonance). If you can measure attenuation at multiple frequencies (either by using multiple transducer frequencies or by spectral analysis of broadband pulses), you can decompose the bubble population into a size distribution rather than a single radius.

    - Temperature field from speed of sound. You're computing speed of sound and using it for velocity reconstruction, but with a known EOS (which you have from the drift flux model), the speed of sound field also constrains the temperature field. For water, c(T) is well-characterized, and since you already have c per cell, you could invert for T per cell — giving you an independent check on your thermal solver.

Higher effort but high value:
    - Bubble passage frequency / number density. The temporal fluctuation rate of backscatter intensity at a given cell tells you how often bubbles pass through. Combined with bubble_radii and gas_velocity, this gives you bubble number density n = gas_holdup / bubble_volume, which would directly close the constraint you've been worried about in calculate_bubble_velocity! between gas_holdup and bubble_radii.

Bottom line
The two biggest gaps are backscatter amplitude (it's right there in your data, you're just not using it) and frequency-dependent attenuation (needs more signal processing but gives you bubble size distribution). The backscatter one in particular would give you an independent per-cell constraint on void fraction that doesn't go through the lossy tomographic inversion, since it comes from the self-to-self echo path which has per-beam spatial resolution.