# FourierDemo

Three interactive pages about spectra, sampling, filtering and Fourier series. Plain HTML, CSS and JavaScript: no build step, no dependencies.

## Run it

```
python -m http.server 8000 --bind 127.0.0.1
```

Run it from this folder, then open <http://127.0.0.1:8000/>. Opening `index.html` straight from disk also works.

## Pages

- [Sampling & aliasing](#sampling--aliasing) — `index.html`
- [Spectrum × Filter](#spectrum--filter) — `filter.html`
- [Epicycles](#epicycles) — `epicycles.html`

## Sampling & aliasing

A single-page interactive demo of sampling, the Nyquist frequency and aliasing, using a square wave.

- **Original:** a square wave of frequency f0, with its Fourier series as amplitude and phase stems
  (odd harmonics n·f0, amplitude 4/(πn), phase −90° against a cosine reference). The frequency axis runs to
  2.5·fs, with fs/2 (the Nyquist frequency), fs and 2·fs marked.
- **Sampled:** the samples, the ideal reconstruction from them (band-limited to fs/2), and the spectrum of what
  the samples contain over 0…fs/2. Every harmonic above fs/2 folds back to an alias frequency in that band.
- **Reconstruction switch:** Ideal (sinc, band-limited to fs/2, the default), Linear (straight lines between the
  samples) or Zero-order hold (a staircase, like a simple DAC). Linear and ZOH leave images of the spectrum around
  fs, 2·fs, … and bend the band down, so the error stays even with the anti-alias filter on; their spectrum axis
  always runs to 2.5·fs.
- **Show difference from original in red:** the error between reconstruction and original in the time plot,
  and in the spectra the aliases, the changed components and (with "wide axis") the harmonics lost above fs/2.
- **Anti-alias filter:** removes the harmonics at or above a cutoff (default fs/2) before sampling. The
  reconstruction then matches the filtered signal to rounding error.

The point: the sampling theorem holds for every signal. A square wave is not band-limited, so it always aliases;
band-limit it first and sampling above twice its highest remaining harmonic gives it back exactly.

Hover (or tap) a stem to see its frequency, amplitude and phase, and which harmonics a sampled component is made of.

### How it is computed

Analytically, from the Fourier series. A harmonic at f with phase φ, sampled at t_k = (k + offset)/fs, gives the
same samples as a component at f − m·fs (m = round(f/fs)) with phase φ + 2π·m·offset; a negative result is
mirrored and its phase sign flipped. Components on the same frequency are added as complex phasors; at DC and at
exactly fs/2 only the part the samples can see is kept. The reconstruction is the sum of these components, which
equals ideal sinc interpolation of the samples. Linear interpolation is the sample train convolved with a triangle
2·Ts wide, centred on the samples: each baseband phasor P at f gives lines P·e^{−i2πm·offset} at f + m·fs and
conj(P)·e^{−i2πm·offset} at m·fs − f, each times sinc²(ν/fs) (sinc x = sin πx/πx). A zero-order hold uses a
rectangle Ts wide instead: sinc(ν/fs)·e^{−iπν/fs}. The time plot draws these reconstructions straight from the
samples. The harmonic sum runs to 40·fs (at least 200, at most 8000 odd
harmonics). Ranges: f0 1–200 Hz, fs 10 Hz–5 kHz.

## Spectrum × Filter

Draw a spectrum and a filter, and see their product.

- **Input:** draw the amplitude (0…1) and phase (±180°) of 65 harmonic bins, k = 0…64, and see the waveform
  x(t) = Σ A_k cos(2πk t/T + φ_k) over two periods. Presets: sine, square, sawtooth, triangle, impulse, random.
- **Filter:** draw a gain H_k (0…1.5) and a phase on the same bins. Presets: all-pass, low-pass, high-pass,
  band-pass, notch, comb and a T/8 delay, each ideal or smooth, with cutoff and width sliders.
- **Interpolation kernels:** Linear, Cubic (Catmull-Rom), Gaussian, Sinc (ideal) and Lanczos-3 as filters. Each is
  the kernel's Fourier transform Ĥ (Ĥ(0) = 1) for samples Ts = T/(2k₀) apart, so H_k = Ĥ(k/(2k₀)) and the k₀ slider
  is their Nyquist bin. Linear is sinc², Gaussian (σ = Ts/2) is a Gaussian, Sinc is a rect (½ at k₀); Catmull-Rom
  and Lanczos-3 are transformed numerically and their negative lobes show as 180° phase. h(t) shows the kernel.
- **Result:** Y_k = X_k · H_k (amplitudes multiply, phases add), as spectra and as the waveform y(t), with x(t)
  drawn faintly for comparison. The filter's impulse response h(t) is shown too: multiplying spectra is the same
  as circular convolution in time.

## Epicycles

Load an image (file picker, drag and drop, or the built-in sample) and watch it become two chains of rotating circles. `epicycles.html`, `epicycles.js`, `epicycles.css`.

- **Vectorize:** blur, then Canny edges or a dark-line threshold, Zhang–Suen thinning, tracing into polylines, Ramer–Douglas–Peucker simplification. Controls: mode, blur, threshold, min stroke length, tolerance.
- **Connect:** strokes are ordered by nearest neighbour, then improved with 2-opt and Or-opt, each entered from its closer end and joined into one closed path. The bridges are shown dashed and can be hidden in the final drawing (pen up).
- **Signals:** the path is resampled by arc length to N = 256…4096 points (default 1024), giving x(t) and y(t).
- **Fourier:** a radix-2 FFT of x and of y separately; amplitude and phase stems for each. Harmonic k is a circle of radius 2|X_k|/N turning k times per lap.
- **Epicycles:** the X chain above the picture sets the pen's x, the Y chain beside it sets the pen's y. Controls: **X circles** and **Y circles** (two independent log sliders, 1…N/2, largest circles first, with 1/5/20/100/all presets and an optional "link X and Y" checkbox), speed, pause, reset trail, show circles. Fewer X circles make the drawing coarse horizontally while Y stays exact; the unused harmonics are dimmed in each spectrum.
- **Cam mechanism:** the last section swaps the circles for one cam per axis (same circle counts). The cam outline is r(θ) = R_base + x(θ) − min x, so a follower sliding on it moves exactly like the X (or Y) series; the two followers drive a cross-slide whose crossing is the pen. A knife-edge follower is exact; a roller follower (radius slider) rides the inward offset of the outline.
- **Cam checks on the page:** the maximum pressure angle (tan φ = r′/r) and any undercut (convex pitch curve with radius of curvature below the roller radius, drawn in red). With many circles the cam is very wiggly and cannot be made for a roller; real harmonic machines (Kelvin's tide predictor) use one simple cam or crank per harmonic instead. That one-cam-per-harmonic view is not implemented.

Files are read with `createImageBitmap`/`FileReader`, so it works from disk. Photos with a lot of texture give many strokes and a busy path; use a higher threshold or min stroke length.
