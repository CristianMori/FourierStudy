# FourierDemo

Interactive pages about Fourier series, spectra, sampling and filtering, for seeing in the browser what the
formulas say. Each page shows a signal in time and its amplitude and phase spectrum side by side; everything
updates live as you move a slider or draw.

Plain HTML, CSS and JavaScript on a canvas: no build step, no libraries, no network access. The spectra are
computed analytically from Fourier series, not estimated with a windowed FFT, so what you see is exact (the
Epicycles page, whose signal comes from an image, uses an exact FFT of the resampled path).

## Contents

- [Run it](#run-it)
- [Pages](#pages)
  - [Sampling & aliasing](#sampling--aliasing) (`index.html`)
  - [Spectrum × Filter](#spectrum--filter) (`filter.html`)
  - [Epicycles](#epicycles) (`epicycles.html`)
- [Project layout](#project-layout)
- [Conventions](#conventions)
- [Author](#author)

## Run it

From this folder:

```
python -m http.server 8000 --bind 127.0.0.1
```

Then open <http://127.0.0.1:8000/>. Any static web server works (the repo is also served by GitHub Pages at
<https://cristianmori.github.io/FourierStudy/>), and so does opening the HTML files straight from disk. The pages link to each other from the bar at the top.

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

### Epicycles

`epicycles.html` — load an image (file picker, drag and drop, or the built-in sample) and watch it become a
single line drawn by two chains of rotating circles, and then by two cams.

- **Vectorize:** blur, then Canny edges or a dark-line threshold, Zhang–Suen thinning, tracing into polylines,
  Ramer–Douglas–Peucker simplification. Controls: mode, blur, threshold, min stroke length, tolerance.
- **Connect:** strokes are ordered by nearest neighbour, then improved with 2-opt and Or-opt, each entered from
  its closer end and joined into one closed path. The bridges are shown dashed and can be hidden in the final
  drawing (pen up).
- **Signals:** the path is resampled by arc length to N = 256…4096 points (default 1024), giving x(t) and y(t), and a
  third signal z(t): 0 (pen down) on stroke samples, 1 (pen up) on bridge samples. z is a control signal, not part of
  the drawing, so it is not Fourier-filtered and has no circles or spectrum. The trail is inked only while the pen is
  down; the pen tip is drawn hollow while it is up, and "show pen-up moves" draws the lifted moves faintly.
- **Fourier:** a radix-2 FFT of x and of y separately; amplitude and phase stems for each. Harmonic k is a circle
  of radius 2|X_k|/N turning k times per lap.
- **Epicycles:** the X chain above the picture sets the pen's x, the Y chain beside it sets the pen's y.
  **X circles** and **Y circles** are two independent log sliders (1…N/2, largest circles first, presets
  1/5/20/100/all, optional "link X and Y"). Fewer X circles make the drawing coarse horizontally while Y stays
  exact; the unused harmonics are dimmed in each spectrum. Also: speed, pause, reset trail, show circles.
- **Cam mechanism:** the last section swaps the circles for one cam per axis, with the same circle counts. The cam
  outline is r(θ) = R_base + x(θ) − min x, so a follower sliding on it moves exactly like the X (or Y) series. There
  is no lever: each follower drives its bar (the X bar vertical, the Y bar horizontal) directly, so bar travel =
  follower travel = pen travel, 1:1. The whole section (cams, followers, bars, drawing) is drawn at one scale, so the
  drawing is smaller than in the circles section; that is the real mechanism. A knife-edge follower is exact; a roller
  follower (radius slider, mm) rides the inward offset of the outline. The page shows the maximum pressure angle
  (tan φ = r′/r, with the true radius) and marks undercut in red (convex pitch curve with radius of curvature below the
  roller radius).
- **Z cam (pen lift):** a third cam on the same shaft timing (same θ, same D-bore flat at θ = 0, same direction) stores
  z(t). Its outline is r(θ) = R_z + h·z̃(θ), with lift h (default 1 mm), base radius R_z (default auto) and z̃ the
  signal with its steps replaced by cycloidal ramps, z̃(u) = u − sin(2πu)/2π, because a cam cannot step. The ramps lie
  inside the bridges: a rise starts on the last stroke sample before a bridge and the fall ends on the first stroke sample
  after it, so every stroke sample is at z̃ = 0 (pen down). Ramp width is the larger of the "ramp" control (default 10°)
  and the angle that keeps the pressure angle at or below 30° at that radius (tan φ = h(1 − cos 2πu)/(β(R_z + h z̃))). A bridge
  too short for two such ramps gets all of its span and the largest lift that keeps the pressure angle in range; a lift
  under 10 % of h is dropped and the pen draws that short bridge. The page reports how many bridges were lifted fully,
  reduced or left down, and the Z cam's pressure angle and undercut, like X and Y. A side-view gauge beside the Z cam shows
  the pen tip over the paper moving with the Z follower (a scale of 0…h mm). **Auto base radius** (the default; the "auto"
  button returns to it after moving the slider) is the smallest whole mm, from 8 up to a cap of 80 mm, at which even the
  shortest bridge gets the full lift within the 30° limit. With a roller follower the ramp width also keeps the pitch
  curve's radius of curvature at or above the roller radius, so the Z cam is not undercut whenever that can be met. If the
  cap is reached the page says how many bridges are still reduced or left down (on the sample: 3 of 16 with a knife-edge
  follower and 11 of 16 with the 4 mm roller, whose short bridges would need a 114 or 163 mm base radius).
- **Switch for servo (Z cam mode):** instead of lifting the pen, the Z cam only trips a switch, and a hobby servo lifts the
  pen. Lift is then the switch travel (default 0.8 mm), the switch trips at half of it, and the pressure-angle limit is its
  own control (default 45°, since a switch needs little force); auto base radius and the roller curvature rule use that limit.
  The servo is modelled as a pure dead time `delay` (default 120 ms) after the switch changes state, and the shaft speed is
  given as seconds per lap (default 60 s), so the lead angle is λ = 360° · delay / lap time (0.72° = 2.0 samples at the
  defaults). The switch command is advanced by λ: the half-lift crossings are put at (bridge start − λ) and (bridge end − λ),
  using the exact bridge ends, so the servo pen leaves the paper exactly at the start of a bridge and lands exactly at the
  next stroke. Each ramp is centred on its crossing and may be as wide as the distance to the neighbouring crossings, so
  ramps never overlap and no stroke with samples is lifted; a bridge that cannot get ramps meeting the pressure-angle (and
  roller) rule is skipped and drawn with the pen down, and the page counts them. The inked trail, the hollow pen tip and the
  side-view gauge follow the delayed servo pen; the gauge also shows the follower on its travel scale with the trip level, and
  a switch lamp, and a z_switch(t) trace appears next to z(t). The servo's own travel time is not modelled.
- **Physical size:** "drawing width" (20–200 mm, default 60) sets the mm per image pixel; "base radius" (default 25 mm)
  is R_base. The page shows the follower lift and both cam outer diameters in mm. Both cams use the same scale.
- **Download STL:** one binary STL per cam (`cam-x.stl`, `cam-y.stl`, `cam-z.stl`, or all three), in mm, z up, a closed watertight solid:
  the outline (the pitch curve, or the roller profile) extruded by the thickness (default 5 mm) around a centre bore
  (default 5 mm), built as an annulus at matching angle samples (N samples, so 8N triangles). Refused with a message when
  the roller profile has undercut (it self-intersects) or when the wall around the bore would be under 1.5 mm. The file
  is made in the browser with Blob and an `<a download>` link; it works from disk and from GitHub Pages.
- **Timing:** with "D-bore" on (default) the bore has a flat at θ = 0. At the start of the lap the flat faces that cam's
  follower (X and Z followers on the right of their cams, Y follower below its cam), so mount the cams with the flat in
  the same relation to its follower and they are in phase. Seen from the top face (+z) each cam turns anticlockwise.
  With a larger lift the cams get big (outer diameter of about 2 × (R_base + lift)); real harmonic machines (Kelvin's
  tide predictor) use one simple cam or crank per harmonic instead. That one-cam-per-harmonic view is not
  implemented.

Unlike the other two pages this one uses a numeric FFT, because the signal comes from an image; with all circles
the pen reproduces the path to about 1e-12 px. Images are read with `createImageBitmap`/`FileReader`, so it works
from disk. Photos with a lot of texture give many strokes and a busy path; line drawings, logos and text work
best, or raise the threshold or min stroke length.

## Project layout

| File | What it is |
|------|------------|
| `index.html`, `app.js`, `style.css` | Sampling & aliasing page |
| `filter.html`, `filter.js`, `filter.css` | Spectrum × Filter page |
| `epicycles.html`, `epicycles.js`, `epicycles.css` | Epicycles page (image → circles → cams) |
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
