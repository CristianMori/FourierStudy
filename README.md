# Square wave aliasing demo

A single-page interactive demo of sampling, the Nyquist frequency and aliasing, using a square wave.

## What it shows

- **Original:** a square wave of frequency f0, with its Fourier series as amplitude and phase stems
  (odd harmonics n·f0, amplitude 4/(πn), phase −90° against a cosine reference). The frequency axis runs to
  2.5·fs, with fs/2 (the Nyquist frequency), fs and 2·fs marked.
- **Sampled:** the samples, the ideal reconstruction from them (band-limited to fs/2), and the spectrum of what
  the samples contain over 0…fs/2. Every harmonic above fs/2 folds back to an alias frequency in that band.
- **Show difference from original in red:** the error between reconstruction and original in the time plot,
  and in the spectra the aliases, the changed components and (with "wide axis") the harmonics lost above fs/2.
- **Anti-alias filter:** removes the harmonics at or above a cutoff (default fs/2) before sampling. The
  reconstruction then matches the filtered signal to rounding error.

The point: the sampling theorem holds for every signal. A square wave is not band-limited, so it always aliases;
band-limit it first and sampling above twice its highest remaining harmonic gives it back exactly.

Hover (or tap) a stem to see its frequency, amplitude and phase, and which harmonics a sampled component is made of.

## How to open it

Open `index.html` in a browser, straight from disk. No build step, no server, no dependencies.

## How it is computed

Analytically, from the Fourier series. A harmonic at f with phase φ, sampled at t_k = (k + offset)/fs, gives the
same samples as a component at f − m·fs (m = round(f/fs)) with phase φ + 2π·m·offset; a negative result is
mirrored and its phase sign flipped. Components on the same frequency are added as complex phasors; at DC and at
exactly fs/2 only the part the samples can see is kept. The reconstruction is the sum of these components, which
equals ideal sinc interpolation of the samples. The harmonic sum runs to 40·fs (at least 200, at most 8000 odd
harmonics). Ranges: f0 1–200 Hz, fs 10 Hz–5 kHz.
