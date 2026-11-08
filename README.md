# FourierDemo

Interactive pages about Fourier series, spectra, sampling and filtering, for seeing in the browser what the
formulas say. Each page shows a signal in time and its amplitude and phase spectrum side by side; everything
updates live as you move a slider or draw.

Plain HTML, CSS and JavaScript on a canvas: no build step, no libraries, no network access. The spectra are
computed analytically from Fourier series, not estimated with a windowed FFT, so what you see is exact.

## Contents

- [Run it](#run-it)
- [Pages](#pages)
  - [Sampling & aliasing](#sampling--aliasing) (`index.html`)
  - [Spectrum × Filter](#spectrum--filter) (`filter.html`)
- [Project layout](#project-layout)
- [Conventions](#conventions)
- [Author](#author)

## Run it

From this folder:

```
python -m http.server 8000 --bind 127.0.0.1
```

Then open <http://127.0.0.1:8000/>. Any static web server works, and so does opening `index.html` or
`filter.html` straight from disk. The pages link to each other from the bar at the top.

After an update, reload with Ctrl+F5 so the browser does not keep an old copy of the scripts.

## Pages

### Sampling & aliasing

`index.html` — sampling, the Nyquist frequency and aliasing, using a square wave.

- **Original:** a square wave of frequency f0, with its Fourier series as amplitude and phase stems
  (odd harmonics n·f0, amplitude 4/(πn), phase −90° against a cosine reference). The frequency axis runs to
  2.5·fs, with fs/2 (the Nyquist frequency), fs and 2·fs marked.
- **Sampling frequency:** a slider from 10 Hz to 5 kHz with presets, including ones tied to f0 (2·f0, 10·f0,
  40·f0), and a sampling-offset slider that moves the samples within a period.
- **Sampled:** the samples, the reconstruction from them, and the spectrum of what the samples contain. Every
  harmonic above fs/2 folds back to an alias frequency in 0…fs/2.
- **Reconstruction switch:** Ideal (sinc, band-limited to fs/2, the default), Linear (straight lines between the
  samples) or Zero-order hold (a staircase, like a simple DAC). Linear and ZOH leave images of the spectrum around
  fs, 2·fs, … and bend the band down, so the error stays even with the anti-alias filter on; their spectrum axis
  always runs to 2.5·fs.
- **Show difference from original in red:** the error between reconstruction and original in the time plot,
  and in the spectra the aliases, the changed components and (with "wide axis") the harmonics lost above fs/2.
- **Anti-alias filter:** removes the harmonics at or above a cutoff (default fs/2) before sampling. With Ideal
  reconstruction the result then matches the filtered signal to rounding error.

The point: the sampling theorem holds for every signal. A square wave is not band-limited, so it always aliases;
band-limit it first and sampling above twice its highest remaining harmonic gives it back exactly.

Hover (or tap) a stem to see its frequency, amplitude and phase, and which harmonics a sampled component is made
of. A status line shows how many harmonics alias and what share of the power they carry.

#### How it is computed

Analytically, from the Fourier series. A harmonic at f with phase φ, sampled at t_k = (k + offset)/fs, gives the
same samples as a component at f − m·fs (m = round(f/fs)) with phase φ + 2π·m·offset; a negative result is
mirrored and its phase sign flipped. Components on the same frequency are added as complex phasors; at DC and at
exactly fs/2 only the part the samples can see is kept. The ideal reconstruction is the sum of these components,
which equals sinc interpolation of the samples.

Linear interpolation is the sample train convolved with a triangle 2·Ts wide, centred on the samples: each
baseband phasor P at f gives lines P·e^{−i2πm·offset} at f + m·fs and conj(P)·e^{−i2πm·offset} at m·fs − f, each
times sinc²(ν/fs) (sinc x = sin πx/πx). A zero-order hold uses a rectangle Ts wide instead:
sinc(ν/fs)·e^{−iπν/fs}. The time plot draws these reconstructions straight from the samples.

The harmonic sum runs to 40·fs (at least 200, at most 8000 odd harmonics), so Gibbs ringing is visible, as it
should be. Ranges: f0 1–200 Hz, fs 10 Hz–5 kHz.

### Spectrum × Filter

`filter.html` — draw a spectrum and a filter, and see their product.

- **Input:** draw the amplitude (0…1) and phase (±180°) of 65 harmonic bins, k = 0…64, with the mouse or a
  finger, and see the waveform x(t) = Σ A_k cos(2πk t/T + φ_k) over two periods. Presets: sine, square,
  sawtooth, triangle, impulse, random.
- **Filter:** draw a gain H_k (0…1.5, dashed line at 1) and a phase on the same bins. Presets: all-pass,
  low-pass, high-pass, band-pass, notch, comb and a T/8 delay, each ideal or smooth, with cutoff and width sliders.
- **Interpolation kernels:** Linear, Cubic (Catmull-Rom), Gaussian, Sinc (ideal) and Lanczos-3 as filters. Each is
  the kernel's Fourier transform Ĥ (Ĥ(0) = 1) for samples Ts = T/(2k₀) apart, so H_k = Ĥ(k/(2k₀)) and the k₀ slider
  is their Nyquist bin. Linear is sinc², Gaussian (σ = Ts/2) is a Gaussian, Sinc is a rect (½ at k₀); Catmull-Rom
  and Lanczos-3 are transformed numerically and their negative lobes show as 180° phase.
- **Result:** Y_k = X_k · H_k (amplitudes multiply, phases add), as spectra and as the waveform y(t), with x(t)
  drawn faintly for comparison.
- **Impulse response:** h(t) = H_0 + Σ 2|H_k| cos(2πk t/T + θ_k), plotted over −T/2…T/2. Multiplying spectra is
  the same as circular convolution in time: y(t) = (1/T)·∫₀ᵀ x(τ)·h(t − τ) dτ.

#### Trying a filter elsewhere

The impulse response of any filter drawn on the page can be typed into a plotter such as Desmos, with the gains
A_k and phases θ_k (in radians) as lists:

```
h(x)=A_0+2\sum_{k=1}^{K}A_k\cos\left(\frac{2\pi kx}{T}+\theta_k\right)
y(x)=\frac{1}{T}\int_{0}^{T}f(t)h(x-t)dt
```

For example the T/8 delay preset (all gains 1, θ_k = −2πk/8) gives
h(x) = 1 + 2·Σ_{k=1}^{64} cos(2πk(x − T/8)/T) = sin((2K+1)π(x − T/8)/T) / sin(π(x − T/8)/T),
a band-limited spike at T/8, so y(x) is f(x − T/8) for any f made of harmonics up to 64.

## Project layout

| File | What it is |
|------|------------|
| `index.html`, `app.js`, `style.css` | Sampling & aliasing page |
| `filter.html`, `filter.js`, `filter.css` | Spectrum × Filter page |
| `README.md` | This file |

Each page is self-contained: its own HTML, script and stylesheet, sharing only the link bar at the top.

## Conventions

- One-sided cosine amplitudes: a real signal is written Σ A_k cos(2πk t/T + φ_k); phases are in degrees on
  screen and wrapped to ±180°.
- Red is kept for differences and errors only.
- Light and dark themes follow the system setting; the layout works down to phone width, and drawing works by
  touch.

## Author

Cristian Mori.
