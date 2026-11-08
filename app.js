'use strict';

(function () {
  const TAU = 2 * Math.PI;
  const F0_MIN = 1, F0_MAX = 200;
  const FS_MIN = 10, FS_MAX = 5000;
  const VIEW = 2.5;         // spectrum axis extends to VIEW * fs
  const HARM_TOP = 40;      // harmonic sum runs up to HARM_TOP * fs ...
  const HARM_MIN = 200;     // ... with at least this many odd harmonics
  const HARM_MAX = 8000;    // ... and at most this many
  const DIFF_EPS = 1e-9;

  const $ = (id) => document.getElementById(id);

  const state = { f0: 60, fs: 1000, delta: 0.13, filter: false, cutFrac: 1, diff: true, wide: false, method: 'ideal' };
  const METHODS = ['ideal', 'linear', 'zoh'];

  // ---------------------------------------------------------------- numerics

  // Ideal square wave sign(sin(2π f0 t)); 0 exactly on a jump, like its Fourier series.
  function square(t, f0) {
    let u = (t * f0) % 1;
    if (u < 0) u += 1;
    if (u === 0 || u === 0.5) return 0;
    return u < 0.5 ? 1 : -1;
  }

  // out[j] = Σ Re{(re + i·im) · e^{i2π f (t0 + j·dt)}}, evaluated by rotating each phasor.
  function evalSum(comps, t0, dt, N) {
    const out = new Float64Array(N);
    for (const c of comps) {
      if (c.re === 0 && c.im === 0) continue;
      const a0 = TAU * c.f * t0, da = TAU * c.f * dt;
      const ca = Math.cos(a0), sa = Math.sin(a0);
      let zr = c.re * ca - c.im * sa;
      let zi = c.re * sa + c.im * ca;
      const wr = Math.cos(da), wi = Math.sin(da);
      for (let j = 0; j < N; j++) {
        out[j] += zr;
        const nr = zr * wr - zi * wi;
        zi = zr * wi + zi * wr;
        zr = nr;
      }
    }
    return out;
  }

  // Spectrum of the reconstruction kernel, relative to the sample instants, as [re, im].
  // linear: triangle of width 2·Ts → sinc²(f/fs); zoh: rectangle [0, Ts) → sinc(f/fs)·e^{−iπ f/fs}.
  function kernel(method, f, fs) {
    const x = f / fs;
    const sc = x === 0 ? 1 : Math.sin(Math.PI * x) / (Math.PI * x);
    if (method === 'linear') return [sc * sc, 0];
    return [sc * Math.cos(Math.PI * x), -sc * Math.sin(Math.PI * x)];
  }

  // Line spectrum of (sample train ⊛ kernel) up to fMax. A baseband phasor P at f has the same samples
  // as P·e^{−i2πmδ} at f + m·fs and conj(P)·e^{−i2πmδ} at m·fs − f; each is weighted by the kernel there.
  // The original's components are merged in so every line can be compared with what was really there.
  function imageLines(R, fMax) {
    const { s, bins, kept, nyq, tol } = R;
    const fs = s.fs;
    const items = [];
    const add = (b, m, mir, f) => {
      const pr = b.re, pi = mir ? -b.im : b.im;
      const a = -TAU * m * s.delta, ca = Math.cos(a), sa = Math.sin(a);
      const qr = pr * ca - pi * sa, qi = pr * sa + pi * ca;
      const [wr, wi] = kernel(s.method, f, fs);
      items.push({ f, b, m, mir, w: Math.hypot(wr, wi), re: qr * wr - qi * wi, im: qr * wi + qi * wr });
    };
    for (const b of bins) {
      if (b.re === 0 && b.im === 0) continue;
      for (let m = 0; b.f + m * fs <= fMax; m++) add(b, m, false, b.f + m * fs);
      for (let m = 1; m * fs - b.f <= fMax; m++) add(b, m, true, m * fs - b.f);
    }
    for (const h of kept) if (h.f <= fMax) items.push({ f: h.f, h });
    items.sort((a, c) => a.f - c.f);
    const lines = [];
    for (const it of items) {
      let l = lines[lines.length - 1];
      if (!l || it.f - l.f > tol) {
        l = { f: it.f, re: 0, im: 0, ore: 0, oim: 0, contrib: [] };
        lines.push(l);
      }
      if (it.h) {
        l.ore += it.h.A * Math.cos(it.h.phi);
        l.oim += it.h.A * Math.sin(it.h.phi);
      } else {
        l.re += it.re;
        l.im += it.im;
        l.contrib.push(it);
      }
    }
    for (const l of lines) {
      l.amp = Math.hypot(l.re, l.im);
      l.ph = Math.atan2(l.im, l.re);
      l.oamp = Math.hypot(l.ore, l.oim);
      l.oph = Math.atan2(l.oim, l.ore);
      l.isDiff = Math.hypot(l.re - l.ore, l.im - l.oim) > DIFF_EPS;
      l.image = l.f > nyq + tol;
      l.contrib.sort((a, c) => Math.hypot(c.re, c.im) - Math.hypot(a.re, a.im));
    }
    return lines.filter((l) => l.amp > 1e-12 || l.oamp > 0);
  }

  // Linear interpolation or zero-order hold of xs (xs[i] is sample k = i − 1) on t = j·dt.
  function holdRecon(xs, method, fs, delta, dt, N) {
    const out = new Float64Array(N);
    for (let j = 0; j < N; j++) {
      const u = j * dt * fs - delta;
      const k = Math.floor(u + 1e-9), i = clamp(k + 1, 0, xs.length - 2);  // a sample instant takes the new sample
      out[j] = method === 'zoh' ? xs[i] : xs[i] + (u - k) * (xs[i + 1] - xs[i]);
    }
    return out;
  }

  function compute(s) {
    const { f0, fs, delta } = s;
    const nyq = fs / 2;
    const fView = VIEW * fs;
    const tol = 1e-7 * fs;
    const fc = s.filter ? s.cutFrac * nyq : Infinity;

    const fTop = Math.max(HARM_TOP * fs, 1.05 * fView);
    const H = Math.min(HARM_MAX, Math.max(HARM_MIN, Math.ceil((fTop / f0 + 1) / 2)));

    // Every odd harmonic, and where its samples put it in 0..fs/2.
    const harm = new Array(H);
    for (let k = 0; k < H; k++) {
      const n = 2 * k + 1;
      const f = n * f0;
      const A = 4 / (Math.PI * n);
      const phi = -Math.PI / 2;
      const m = Math.round(f / fs);
      let fa = f - m * fs;
      let psi = phi + TAU * m * delta;
      let conj = false;
      if (fa < 0) { fa = -fa; psi = -psi; conj = true; }
      if (fa < tol) fa = 0;
      if (Math.abs(fa - nyq) < tol) fa = nyq;
      harm[k] = {
        n, f, A, phi, m, fa, psi, conj,
        kept: f < fc * (1 - 1e-12),
        inBand: f <= nyq + tol,
      };
    }
    const kept = harm.filter((h) => h.kept);

    // Add components that land on the same frequency as complex phasors.
    const sorted = kept.slice().sort((a, b) => a.fa - b.fa);
    const bins = [];
    for (const h of sorted) {
      let b = bins[bins.length - 1];
      if (!b || h.fa - b.f > tol) {
        b = { f: h.fa, re: 0, im: 0, ore: 0, oim: 0, contrib: [] };
        bins.push(b);
      }
      b.re += h.A * Math.cos(h.psi);
      b.im += h.A * Math.sin(h.psi);
      if (h.inBand) {
        b.ore += h.A * Math.cos(h.phi);
        b.oim += h.A * Math.sin(h.phi);
      }
      b.contrib.push(h);
    }
    for (const b of bins) {
      if (b.f === 0) {
        // DC: the samples see only the real part.
        b.im = 0;
      } else if (b.f === nyq) {
        // At fs/2 the samples are (−1)^k·c with c = Re{P·e^{iπδ}}; sinc
        // interpolation turns that into c·cos(π fs t − πδ).
        const cd = Math.cos(Math.PI * delta), sd = Math.sin(Math.PI * delta);
        const c = b.re * cd - b.im * sd;
        b.re = c * cd;
        b.im = -c * sd;
      }
      b.amp = Math.hypot(b.re, b.im);
      b.ph = Math.atan2(b.im, b.re);
      b.oamp = Math.hypot(b.ore, b.oim);
      b.oph = Math.atan2(b.oim, b.ore);
      b.isDiff = Math.hypot(b.re - b.ore, b.im - b.oim) > DIFF_EPS;
      b.contrib.sort((a, c) => c.A - a.A);
    }

    // Time signals.
    const T = Math.max(3 / f0, 12 / fs);
    const N = Math.min(8000, Math.max(2400, Math.ceil(T * fs * 8)));
    const dt = T / (N - 1);
    const inputComps = kept.map((h) => ({ f: h.f, re: h.A * Math.cos(h.phi), im: h.A * Math.sin(h.phi) }));

    // Samples k = −1 … K, so that a hold or a line covers the whole plot.
    const K = Math.max(0, Math.floor(T * fs - delta) + 1);
    const t0s = delta / fs;
    let xs;
    if (s.filter) {
      xs = evalSum(inputComps, t0s - 1 / fs, 1 / fs, K + 2);
    } else {
      xs = new Float64Array(K + 2);
      for (let i = 0; i < K + 2; i++) xs[i] = square((i - 1 + delta) / fs, f0);
    }
    const sampleVals = xs.subarray(1, K + 1);

    const recon = s.method === 'ideal' ? evalSum(bins, 0, dt, N) : holdRecon(xs, s.method, fs, delta, dt, N);
    let input;
    if (s.filter) {
      input = evalSum(inputComps, 0, dt, N);
    } else {
      input = new Float64Array(N);
      for (let j = 0; j < N; j++) input[j] = square(j * dt, f0);
    }
    let se = 0, si = 0;
    for (let j = 0; j < N; j++) {
      const e = recon[j] - input[j];
      se += e * e;
      si += input[j] * input[j];
    }
    const rmsErr = Math.sqrt(se / N);
    const rmsIn = Math.sqrt(si / N);

    // Bookkeeping for the status line.
    let pBelow = 0, nBelow = 0, nAt = 0, nAboveView = 0, pKept = 0;
    for (const h of harm) {
      if (h.kept) pKept += h.A * h.A / 2;
      if (h.f < nyq - tol) { pBelow += h.A * h.A / 2; nBelow++; }
      else if (h.f <= nyq + tol) nAt++;
      else if (h.f <= fView) nAboveView++;
    }
    const highestKept = kept.length ? kept[kept.length - 1].f : 0;

    const R = {
      s, nyq, fView, tol, fc, H, harm, kept, bins,
      T, N, dt, recon, input, rmsErr, rmsIn,
      K, t0s, xs, sampleVals,
      aliasedPowerFrac: s.filter ? 0 : Math.max(0, 1 - pBelow),
      nBelow, nAt, nAboveView, pKept, highestKept,
    };
    R.lines = s.method === 'ideal' ? null : imageLines(R, fView * 1.0001);
    return R;
  }

  // ---------------------------------------------------------------- formatting

  function trimNum(str) {
    return str.indexOf('.') >= 0 ? str.replace(/0+$/, '').replace(/\.$/, '') : str;
  }
  function fmtHz(v) {
    const a = Math.abs(v);
    if (a >= 1000) return trimNum(v.toFixed(1));
    if (a >= 100) return trimNum(v.toFixed(1));
    if (a >= 10) return trimNum(v.toFixed(2));
    return trimNum(v.toFixed(3));
  }
  function fmtTick(v) {
    const a = Math.abs(v);
    if (a >= 1000) return trimNum((v / 1000).toFixed(3)) + 'k';
    return trimNum(v.toFixed(6));
  }
  function fmtAmp(v) { return v < 1e-3 && v > 0 ? v.toExponential(2) : v.toFixed(4); }
  function fmtDeg(rad) {
    let d = rad * 180 / Math.PI;
    if (d <= -180 + 1e-9) d += 360;
    return (Math.abs(d) < 0.05 ? '0' : d.toFixed(1)) + '°';
  }
  function niceTicks(min, max, count) {
    const raw = (max - min) / Math.max(1, count);
    const p = Math.pow(10, Math.floor(Math.log10(raw)));
    const r = raw / p;
    const step = (r < 1.5 ? 1 : r < 3 ? 2 : r < 7 ? 5 : 10) * p;
    const out = [];
    for (let v = Math.ceil(min / step - 1e-9) * step; v <= max + step * 1e-9; v += step) out.push(+v.toPrecision(12));
    return out;
  }

  // ---------------------------------------------------------------- drawing

  let C = null;
  function readColors() {
    const cs = getComputedStyle(document.documentElement);
    const g = (n) => cs.getPropertyValue(n).trim();
    C = {
      fg: g('--fg'), muted: g('--muted'), grid: g('--grid'), axis: g('--axis'), shade: g('--shade'),
      card: g('--card'), orig: g('--c-orig'), recon: g('--c-recon'), sample: g('--c-sample'),
      ghost: g('--c-ghost'), err: g('--c-err'), errFill: g('--c-err-fill'),
      font: getComputedStyle(document.body).fontFamily,
    };
  }

  function frame(cv, o) {
    const dpr = window.devicePixelRatio || 1;
    const w = cv.clientWidth, h = cv.clientHeight;
    const bw = Math.round(w * dpr), bh = Math.round(h * dpr);
    if (cv.width !== bw || cv.height !== bh) { cv.width = bw; cv.height = bh; }
    const ctx = cv.getContext('2d');
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.clearRect(0, 0, w, h);
    const m = { l: 44, r: 12, t: 8, b: 22 };
    const pw = Math.max(10, w - m.l - m.r), ph = Math.max(10, h - m.t - m.b);
    const X = (x) => m.l + (x - o.x0) / (o.x1 - o.x0) * pw;
    const Y = (y) => m.t + (1 - (y - o.y0) / (o.y1 - o.y0)) * ph;

    ctx.font = '11px ' + C.font;
    ctx.lineWidth = 1;
    ctx.fillStyle = C.muted;
    ctx.textAlign = 'right';
    ctx.textBaseline = 'middle';
    for (const v of o.yTicks) {
      const y = Math.round(Y(v)) + 0.5;
      ctx.strokeStyle = v === 0 ? C.axis : C.grid;
      hline(ctx, m.l, m.l + pw, y);
      ctx.fillText(o.yFmt ? o.yFmt(v) : fmtTick(v), m.l - 5, y);
    }
    ctx.textAlign = 'center';
    ctx.textBaseline = 'top';
    const xt = niceTicks(o.x0, o.x1, Math.max(2, Math.floor(pw / 70)));
    for (const v of xt) {
      const x = Math.round(X(v)) + 0.5;
      ctx.strokeStyle = C.grid;
      vline(ctx, x, m.t, m.t + ph);
      ctx.fillText(o.xFmt ? o.xFmt(v) : fmtTick(v), x, m.t + ph + 5);
    }
    ctx.strokeStyle = C.axis;
    ctx.strokeRect(m.l + 0.5, m.t + 0.5, pw, ph);
    return { ctx, X, Y, m, pw, ph, w, h, o };
  }

  function hline(ctx, x0, x1, y) { ctx.beginPath(); ctx.moveTo(x0, y); ctx.lineTo(x1, y); ctx.stroke(); }
  function vline(ctx, x, y0, y1) { ctx.beginPath(); ctx.moveTo(x, y0); ctx.lineTo(x, y1); ctx.stroke(); }

  function clip(F) {
    F.ctx.save();
    F.ctx.beginPath();
    F.ctx.rect(F.m.l - 4, F.m.t - 4, F.pw + 8, F.ph + 8);
    F.ctx.clip();
  }

  function polyline(F, t0, dt, ys, color, width, dash) {
    const ctx = F.ctx;
    ctx.strokeStyle = color;
    ctx.lineWidth = width;
    ctx.setLineDash(dash || []);
    ctx.beginPath();
    for (let j = 0; j < ys.length; j++) {
      const x = F.X(t0 + j * dt), y = F.Y(ys[j]);
      if (j === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
    }
    ctx.stroke();
    ctx.setLineDash([]);
  }

  // Ideal square wave drawn exactly, edge by edge.
  function drawSquare(F, f0, T, color, width, dash) {
    const ctx = F.ctx;
    ctx.strokeStyle = color;
    ctx.lineWidth = width;
    ctx.setLineDash(dash || []);
    ctx.beginPath();
    let v = 1;
    ctx.moveTo(F.X(0), F.Y(v));
    const half = 1 / (2 * f0);
    for (let k = 1; k * half < T; k++) {
      const x = F.X(k * half);
      ctx.lineTo(x, F.Y(v));
      v = -v;
      ctx.lineTo(x, F.Y(v));
    }
    ctx.lineTo(F.X(T), F.Y(v));
    ctx.stroke();
    ctx.setLineDash([]);
  }

  // Linear or zero-order-hold reconstruction drawn exactly, sample by sample (corners on the sample instants).
  function drawHold(F, R, method, color, width) {
    const ctx = F.ctx, fs = R.s.fs, xs = R.xs;
    const tk = (i) => (i - 1 + R.s.delta) / fs;
    ctx.strokeStyle = color;
    ctx.lineWidth = width;
    ctx.beginPath();
    ctx.moveTo(F.X(tk(0)), F.Y(xs[0]));
    for (let i = 1; i < xs.length; i++) {
      if (method === 'zoh') ctx.lineTo(F.X(tk(i)), F.Y(xs[i - 1]));
      ctx.lineTo(F.X(tk(i)), F.Y(xs[i]));
    }
    if (method === 'zoh') ctx.lineTo(F.X(tk(xs.length)), F.Y(xs[xs.length - 1]));
    ctx.stroke();
  }

  // Nyquist / fs markers on a spectrum plot.
  function markers(F, R, withCut) {
    const ctx = F.ctx, { nyq } = R, fs = R.s.fs;
    const x1 = F.o.x1;
    if (x1 > nyq * 1.001) {
      ctx.fillStyle = C.shade;
      const xa = F.X(nyq);
      ctx.fillRect(xa, F.m.t, F.m.l + F.pw - xa, F.ph);
    }
    const marks = [[nyq, 'fs/2']];
    for (let k = 1; k * fs <= x1 * 1.0001; k++) marks.push([k * fs, k === 1 ? 'fs' : k + 'fs']);
    if (withCut && R.s.filter && R.s.cutFrac < 0.999) marks.push([R.fc, 'cutoff']);
    ctx.font = '11px ' + C.font;
    ctx.textBaseline = 'top';
    for (const [f, label] of marks) {
      const x = Math.round(F.X(f)) + 0.5;
      ctx.strokeStyle = label === 'cutoff' ? C.orig : C.muted;
      ctx.setLineDash([4, 3]);
      vline(ctx, x, F.m.t, F.m.t + F.ph);
      ctx.setLineDash([]);
      ctx.fillStyle = label === 'cutoff' ? C.orig : C.muted;
      const atEdge = x > F.m.l + F.pw - 30;
      ctx.textAlign = atEdge ? 'right' : 'left';
      ctx.fillText(label, x + (atEdge ? -3 : 3), F.m.t + 2 + (label === 'cutoff' ? 12 : 0));
    }
  }

  // Vertical stems from y0 to v; returns hover targets.
  function stem(F, f, v, y0, color, opt) {
    const ctx = F.ctx;
    const x = F.X(f), y = F.Y(v), yb = F.Y(y0);
    ctx.globalAlpha = opt.alpha == null ? 1 : opt.alpha;
    ctx.strokeStyle = color;
    ctx.lineWidth = opt.width || 1.5;
    ctx.setLineDash(opt.dash || []);
    vline(ctx, x, yb, y);
    ctx.setLineDash([]);
    if (opt.dot) {
      ctx.fillStyle = color;
      ctx.beginPath();
      ctx.arc(x, y, 2.6, 0, TAU);
      ctx.fill();
    }
    ctx.globalAlpha = 1;
    return { x, y };
  }
  function ring(F, f, v, color) {
    const ctx = F.ctx;
    ctx.strokeStyle = color;
    ctx.lineWidth = 1.5;
    ctx.beginPath();
    ctx.arc(F.X(f), F.Y(v), 4, 0, TAU);
    ctx.stroke();
  }

  function timeAxis(T) {
    const ms = T < 1;
    return { unit: ms ? 'time (ms)' : 'time (s)', fmt: (v) => fmtTick(ms ? v * 1000 : v) };
  }
  const ampOpts = (x1, yMax) => ({ x0: 0, x1, y0: 0, y1: yMax, yTicks: niceTicks(0, yMax, 3) });
  const phaseOpts = (x1) => ({
    x0: 0, x1, y0: -200, y1: 200, yTicks: [-180, -90, 0, 90, 180], yFmt: (v) => v + '°',
  });
  const deg = (rad) => {
    let d = rad * 180 / Math.PI;
    if (d <= -180 + 1e-9) d += 360;
    return d;
  };

  const cv = {};
  let R = null;

  function drawAll() {
    if (!R) return;
    readColors();
    const s = R.s;
    const ta = timeAxis(R.T);
    $('otUnit').textContent = '· ' + ta.unit;
    $('stUnit').textContent = '· ' + ta.unit;
    const tOpts = { x0: 0, x1: R.T, y0: -1.5, y1: 1.5, yTicks: [-1, 0, 1], xFmt: ta.fmt };

    // Original: time.
    {
      const F = frame(cv.OT, tOpts);
      clip(F);
      if (s.filter) {
        drawSquare(F, s.f0, R.T, C.ghost, 1.2, [5, 4]);
        polyline(F, 0, R.dt, R.input, C.orig, 1.8);
      } else {
        drawSquare(F, s.f0, R.T, C.orig, 1.8);
      }
      F.ctx.restore();
    }

    // Original: spectrum.
    {
      const inView = R.harm.filter((h) => h.f <= R.fView);
      const dot = inView.length <= 120;
      const FA = frame(cv.OA, ampOpts(R.fView, 1.4));
      const FP = frame(cv.OP, phaseOpts(R.fView));
      markers(FA, R, true);
      markers(FP, R, false);
      clip(FA); clip(FP);
      const ta = [], tp = [];
      for (const h of inView) {
        const removed = !h.kept;
        const color = removed ? C.ghost : C.orig;
        const opt = removed ? { dash: [3, 2], dot: false, width: 1.2 } : { dot };
        const a = stem(FA, h.f, h.A, 0, color, opt);
        const p = stem(FP, h.f, deg(h.phi), 0, color, opt);
        let html = `<b>n = ${h.n}</b> · ${fmtHz(h.f)} Hz<br>amplitude ${fmtAmp(h.A)} · phase ${fmtDeg(h.phi)}`;
        if (removed) html += '<br><span class="m">removed by the filter</span>';
        else if (!h.inBand) {
          html += `<br><span class="r">above fs/2: the samples show it at ${fmtHz(h.fa)} Hz` +
            (h.conj ? ' (mirrored, phase sign flipped)' : '') + '</span>';
        } else if (h.f >= R.nyq - R.tol) html += '<br><span class="r">exactly at fs/2: only part of it survives</span>';
        ta.push({ x: a.x, y: a.y, amp: h.A, html });
        tp.push({ x: p.x, y: p.y, amp: h.A, html });
      }
      FA.ctx.restore(); FP.ctx.restore();
      cv.OA.__targets = ta;
      cv.OP.__targets = tp;
    }

    // Sampled: time.
    {
      const F = frame(cv.ST, tOpts);
      clip(F);
      const ctx = F.ctx;
      if (s.filter) polyline(F, 0, R.dt, R.input, C.ghost, 1.2);
      else drawSquare(F, s.f0, R.T, C.ghost, 1.2);
      if (s.diff) {
        ctx.fillStyle = C.errFill;
        ctx.beginPath();
        for (let j = 0; j < R.N; j++) {
          const x = F.X(j * R.dt), y = F.Y(R.recon[j]);
          if (j === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
        }
        for (let j = R.N - 1; j >= 0; j--) ctx.lineTo(F.X(j * R.dt), F.Y(R.input[j]));
        ctx.closePath();
        ctx.fill();
      }
      const K = R.K, y0 = F.Y(0);
      ctx.fillStyle = C.sample;
      ctx.strokeStyle = C.sample;
      ctx.lineWidth = 1;
      for (let k = 0; k < K; k++) {
        const x = F.X(R.t0s + k / s.fs), y = F.Y(R.sampleVals[k]);
        if (K <= 250) {
          ctx.globalAlpha = 0.55;
          vline(ctx, x, y0, y);
          ctx.globalAlpha = 1;
        }
        if (K <= 1500) {
          ctx.beginPath();
          ctx.arc(x, y, K <= 250 ? 3 : 2, 0, TAU);
          ctx.fill();
        } else {
          ctx.fillRect(x - 0.6, y - 0.6, 1.2, 1.2);
        }
      }
      if (s.method === 'ideal') polyline(F, 0, R.dt, R.recon, C.recon, 1.8);
      else drawHold(F, R, s.method, C.recon, 1.8);
      ctx.restore();
    }

    // Sampled: spectrum.
    {
      const ideal = s.method === 'ideal';
      const x1 = (s.wide || !ideal) ? R.fView : R.nyq;
      const shown = ideal ? R.bins.filter((b) => b.f <= x1 * 1.0001) : R.lines;
      let maxAmp = 0;
      for (const b of shown) maxAmp = Math.max(maxAmp, b.amp);
      const yMax = Math.max(1.4, maxAmp * 1.1);
      const lost = (ideal && s.wide && s.diff) ? R.kept.filter((h) => !h.inBand && h.f <= R.fView) : [];
      const dot = shown.length + lost.length <= 120;
      const FA = frame(cv.SA, ampOpts(x1, yMax));
      const FP = frame(cv.SP, phaseOpts(x1));
      markers(FA, R, false);
      markers(FP, R, false);
      clip(FA); clip(FP);
      const ta = [], tp = [];
      for (const h of lost) {
        const opt = { dash: [4, 2], alpha: 0.6, width: 1.5, dot: true };
        const a = stem(FA, h.f, h.A, 0, C.err, opt);
        const p = stem(FP, h.f, deg(h.phi), 0, C.err, opt);
        const html = `<b>n = ${h.n}</b> · ${fmtHz(h.f)} Hz<br><span class="r">lost: above fs/2, ` +
          `reappears at ${fmtHz(h.fa)} Hz</span><br><span class="m">original amplitude ${fmtAmp(h.A)}</span>`;
        ta.push({ x: a.x, y: a.y, amp: h.A * 0.5, html });
        tp.push({ x: p.x, y: p.y, amp: h.A * 0.5, html });
      }
      for (const b of shown) {
        const red = s.diff && b.isDiff;
        const color = red ? C.err : C.recon;
        const a = stem(FA, b.f, b.amp, 0, color, { dot });
        let p = null;
        if (b.amp > 1e-6) p = stem(FP, b.f, deg(b.ph), 0, color, { dot });
        if (red && b.oamp > 0) {
          ring(FA, b.f, b.oamp, C.orig);
          ring(FP, b.f, deg(b.oph), C.orig);
        }
        const html = ideal ? binHtml(b, R) : lineHtml(b, R);
        ta.push({ x: a.x, y: a.y, amp: b.amp, html });
        if (p) tp.push({ x: p.x, y: p.y, amp: b.amp, html });
      }
      FA.ctx.restore(); FP.ctx.restore();
      cv.SA.__targets = ta;
      cv.SP.__targets = tp;
    }
  }

  function binHtml(b, R) {
    let html = `<b>${fmtHz(b.f)} Hz</b><br>amplitude ${fmtAmp(b.amp)} · phase ${fmtDeg(b.ph)}`;
    if (b.isDiff) {
      html += b.oamp > 0
        ? `<br><span class="r">original here: ${fmtAmp(b.oamp)} · ${fmtDeg(b.oph)}</span>`
        : '<br><span class="r">not in the original: pure alias</span>';
    }
    html += '<br><span class="m">made of:</span>';
    const list = b.contrib.slice(0, 5);
    for (const h of list) {
      let what;
      if (h.inBand && h.m === 0) what = 'original';
      else if (h.inBand) what = 'at fs/2, sine part not seen';
      else what = `image ${h.m}·fs` + (h.conj ? ', mirrored, phase flipped' : '');
      html += `<br><span class="m">n=${h.n} (${fmtHz(h.f)} Hz, amp ${fmtAmp(h.A)}): ${what}</span>`;
    }
    if (b.contrib.length > list.length) html += `<br><span class="m">… and ${b.contrib.length - list.length} smaller ones</span>`;
    return html;
  }

  function lineHtml(l, R) {
    let html = `<b>${fmtHz(l.f)} Hz</b>` + (l.image ? ' · <span class="r">image (above fs/2)</span>' : '') +
      `<br>amplitude ${fmtAmp(l.amp)} · phase ${fmtDeg(l.ph)}`;
    if (l.isDiff) {
      html += l.oamp > 0
        ? `<br><span class="r">original here: ${fmtAmp(l.oamp)} · ${fmtDeg(l.oph)}</span>`
        : `<br><span class="r">not in the original: ${l.image ? 'an image left by the reconstruction' : 'pure alias'}</span>`;
    }
    const wName = R.s.method === 'linear' ? 'sinc²' : '|sinc|';
    html += '<br><span class="m">made of:</span>';
    const list = l.contrib.slice(0, 4);
    for (const c of list) {
      const fb = fmtHz(c.b.f);
      const what = c.m === 0 ? `the ${fb} Hz component in 0…fs/2`
        : `image ${c.mir ? `${c.m}·fs − ${fb}` : `${c.m}·fs + ${fb}`} Hz` + (c.mir ? ' (mirrored, phase flipped)' : '');
      html += `<br><span class="m">${what}, × ${wName} = ${c.w.toFixed(4)}</span>`;
    }
    if (l.contrib.length > list.length) html += `<br><span class="m">… and ${l.contrib.length - list.length} more</span>`;
    if (R.s.method === 'zoh' && l.contrib.length) html += '<br><span class="m">plus a half-sample delay (phase −180°·f/fs)</span>';
    return html;
  }

  const RECON_TEXT = {
    ideal: {
      legend: 'ideal reconstruction (band-limited to fs/2)',
      amp: 'of what the samples contain',
      note: 'Ideal reconstruction puts a sinc pulse on every sample: it keeps exactly the band 0…fs/2 and nothing above it.',
    },
    linear: {
      legend: 'linear interpolation (lines between samples)',
      amp: 'of the reconstruction',
      note: 'Joining the samples with straight lines leaves copies of the spectrum around fs, 2fs, … (attenuated by sinc²) ' +
        'and slightly attenuates the band — a real DAC needs a reconstruction filter after it.',
    },
    zoh: {
      legend: 'zero-order hold (staircase)',
      amp: 'of the reconstruction',
      note: 'Holding each sample until the next one, as a simple DAC does, leaves copies of the spectrum around fs, 2fs, … ' +
        '(attenuated by |sinc|), attenuates the band more than lines do, and delays everything by half a sample.',
    },
  };

  // ---------------------------------------------------------------- text

  function updateText() {
    const s = R.s;
    const st = $('status');
    const ratio = s.fs / s.f0;
    let html = `fs/2 = <b>${fmtHz(R.nyq)} Hz</b> · ${trimNum(ratio.toFixed(2))} samples per period. `;
    if (s.filter) {
      if (R.kept.length === 0) {
        html += `The filter (cutoff ${fmtHz(R.fc)} Hz) removes every harmonic, even the fundamental at ${fmtHz(s.f0)} Hz: nothing is left to sample.`;
      } else {
        html += `Anti-alias filter on: it keeps <b>${R.kept.length}</b> harmonic${R.kept.length > 1 ? 's' : ''} below ` +
          `${fmtHz(R.fc)} Hz; the highest is ${fmtHz(R.highestKept)} Hz, below fs/2. ` +
          `Nothing at or above fs/2 is left, so <b>nothing aliases</b>.`;
      }
    } else {
      html += `${R.nBelow} harmonic${R.nBelow === 1 ? '' : 's'} lie${R.nBelow === 1 ? 's' : ''} below fs/2. ` +
        `Every harmonic above fs/2 folds back into 0…fs/2: <b>${R.nAboveView}</b> of them in the plotted range up to 2.5·fs, ` +
        `and infinitely many beyond. They carry <b>${(100 * R.aliasedPowerFrac).toPrecision(3)}%</b> of the square wave's power.`;
      if (R.nAt) html += ' One harmonic sits exactly at fs/2, where the samples see only part of it.';
    }
    st.innerHTML = html;

    const mt = $('metric');
    const rel = R.rmsIn > 0 ? R.rmsErr / R.rmsIn : 0;
    if (R.rmsErr < 1e-9) {
      mt.className = 'metric good';
      mt.textContent = `RMS difference, reconstruction − original: ${R.rmsErr.toExponential(1)}. ` +
        'Zero to rounding error: the samples give the signal back exactly.';
    } else {
      mt.className = 'metric' + (s.diff ? ' bad' : '');
      mt.textContent = `RMS difference, reconstruction − original: ${R.rmsErr.toPrecision(3)} ` +
        `(${(100 * rel).toPrecision(3)}% of the original's RMS).`;
      if (s.filter && s.method !== 'ideal') {
        mt.textContent += ' Not zero even with the anti-alias filter: the images above fs/2 and the droop in the band are the difference.';
      }
    }

    const rt = RECON_TEXT[s.method];
    $('reconNote').textContent = rt.note;
    $('stLegendRecon').textContent = rt.legend;
    $('saWhat').textContent = rt.amp;

    const lostN = R.kept.filter((h) => !h.inBand && h.f <= R.fView).length;
    const aliasBins = R.bins.filter((b) => b.isDiff).length;
    let note = '';
    if (s.diff && s.method !== 'ideal') {
      const img = R.lines.filter((l) => l.isDiff && l.image).length;
      const inb = R.lines.filter((l) => l.isDiff && !l.image).length;
      note = `Up to 2.5·fs, ${img} line${img === 1 ? '' : 's'} above fs/2 and ${inb} in 0…fs/2 differ from the original.`;
    } else if (s.diff && !s.filter) {
      note = `${aliasBins} frequenc${aliasBins === 1 ? 'y' : 'ies'} in 0…fs/2 differ from the original.`;
      if (!s.wide && lostN) note += ` ${lostN} harmonics between fs/2 and 2.5·fs are lost from their true place; tick "wide axis" to see them.`;
    }
    $('lostNote').textContent = note;
  }

  // ---------------------------------------------------------------- controls

  const logPos = (v, lo, hi) => 1000 * Math.log(v / lo) / Math.log(hi / lo);
  const logVal = (p, lo, hi) => lo * Math.pow(hi / lo, p / 1000);
  const clamp = (v, lo, hi) => Math.min(hi, Math.max(lo, v));
  const roundF0 = (v) => (v < 10 ? Math.round(v * 10) / 10 : Math.round(v));
  const roundFs = (v) => (v < 1000 ? Math.round(v) : Math.round(v / 10) * 10);

  function syncControls() {
    $('f0Range').value = logPos(state.f0, F0_MIN, F0_MAX);
    $('fsRange').value = logPos(state.fs, FS_MIN, FS_MAX);
    if (document.activeElement !== $('f0Num')) $('f0Num').value = trimNum(state.f0.toFixed(3));
    if (document.activeElement !== $('fsNum')) $('fsNum').value = trimNum(state.fs.toFixed(3));
    $('offRange').value = state.delta;
    const ts = 1 / state.fs;
    $('offVal').textContent = `${state.delta.toFixed(2)} · Ts (first sample at ${fmtHz(state.delta * ts * 1000)} ms)`;
    $('filter').checked = state.filter;
    $('cutRange').value = Math.round(state.cutFrac * 100);
    $('cutRange').disabled = !state.filter;
    $('cutCtl').classList.toggle('off', !state.filter);
    $('cutVal').textContent = `${Math.round(state.cutFrac * 100)}% of fs/2 = ${fmtHz(state.cutFrac * state.fs / 2)} Hz`;
    $('diff').checked = state.diff;
    $('wide').checked = state.wide || state.method !== 'ideal';
    $('wide').disabled = state.method !== 'ideal';
    for (const r of document.querySelectorAll('input[name=recon]')) r.checked = r.value === state.method;
    document.body.classList.toggle('diff-on', state.diff);
    document.body.classList.toggle('wide-on', state.wide && state.diff && state.method === 'ideal');
    $('otLegendBL').hidden = !state.filter;
    $('otLegendSq').hidden = state.filter;
    $('oaLegendRem').hidden = !state.filter;
  }

  let pending = false;
  function schedule() {
    syncControls();
    if (pending) return;
    pending = true;
    requestAnimationFrame(() => {
      pending = false;
      R = compute(state);
      drawAll();
      updateText();
      window.__demo.last = R;
    });
  }

  function bind() {
    $('f0Range').addEventListener('input', (e) => {
      state.f0 = clamp(roundF0(logVal(+e.target.value, F0_MIN, F0_MAX)), F0_MIN, F0_MAX);
      schedule();
    });
    $('fsRange').addEventListener('input', (e) => {
      state.fs = clamp(roundFs(logVal(+e.target.value, FS_MIN, FS_MAX)), FS_MIN, FS_MAX);
      schedule();
    });
    $('f0Num').addEventListener('change', (e) => {
      const v = parseFloat(e.target.value);
      if (isFinite(v)) state.f0 = clamp(v, F0_MIN, F0_MAX);
      e.target.value = trimNum(state.f0.toFixed(3));
      schedule();
    });
    $('fsNum').addEventListener('change', (e) => {
      const v = parseFloat(e.target.value);
      if (isFinite(v)) state.fs = clamp(v, FS_MIN, FS_MAX);
      e.target.value = trimNum(state.fs.toFixed(3));
      schedule();
    });
    document.querySelectorAll('.presets button').forEach((btn) => {
      btn.addEventListener('click', () => {
        const v = btn.dataset.rel ? state.f0 * +btn.dataset.rel : +btn.dataset.fs;
        state.fs = clamp(v, FS_MIN, FS_MAX);
        $('fsNum').value = trimNum(state.fs.toFixed(3));
        schedule();
      });
    });
    $('offRange').addEventListener('input', (e) => { state.delta = +e.target.value; schedule(); });
    $('filter').addEventListener('change', (e) => { state.filter = e.target.checked; schedule(); });
    $('cutRange').addEventListener('input', (e) => { state.cutFrac = +e.target.value / 100; schedule(); });
    $('diff').addEventListener('change', (e) => { state.diff = e.target.checked; schedule(); });
    $('wide').addEventListener('change', (e) => { state.wide = e.target.checked; schedule(); });
    for (const r of document.querySelectorAll('input[name=recon]')) {
      r.addEventListener('change', (e) => {
        if (e.target.checked && METHODS.includes(e.target.value)) { state.method = e.target.value; schedule(); }
      });
    }

    const tip = $('tip');
    for (const c of [cv.OA, cv.OP, cv.SA, cv.SP]) {
      const show = (e) => {
        const targets = c.__targets || [];
        const r = c.getBoundingClientRect();
        const mx = e.clientX - r.left;
        let best = null, bestD = 9;
        for (const t of targets) {
          const d = Math.abs(t.x - mx);
          if (d < bestD - 0.5 || (Math.abs(d - bestD) <= 0.5 && best && t.amp > best.amp)) { best = t; bestD = d; }
        }
        if (!best) { tip.hidden = true; return; }
        tip.innerHTML = best.html;
        tip.hidden = false;
        const tw = tip.offsetWidth, th = tip.offsetHeight;
        let x = e.clientX + 14, y = e.clientY + 14;
        if (x + tw > window.innerWidth - 8) x = e.clientX - tw - 14;
        if (y + th > window.innerHeight - 8) y = e.clientY - th - 14;
        tip.style.left = Math.max(8, x) + 'px';
        tip.style.top = Math.max(8, y) + 'px';
      };
      c.addEventListener('pointermove', show);
      c.addEventListener('pointerdown', show);
      c.addEventListener('pointerleave', () => { tip.hidden = true; });
    }
    window.addEventListener('scroll', () => { tip.hidden = true; }, { passive: true });

    const ro = new ResizeObserver(() => drawAll());
    ro.observe(document.querySelector('.panels'));
    const mq = window.matchMedia('(prefers-color-scheme: dark)');
    if (mq.addEventListener) mq.addEventListener('change', drawAll);
  }

  // Exposed for testing from the console.
  window.__demo = {
    state, compute, last: null,
    set(patch) { Object.assign(state, patch); schedule(); },
    // Largest gap between the reconstruction and the samples at the sample instants.
    sampleDeviation() {
      const r = window.__demo.last;
      const at = evalSum(r.bins, r.t0s, 1 / r.s.fs, r.K);
      let mx = 0;
      for (let k = 0; k < r.K; k++) mx = Math.max(mx, Math.abs(at[k] - r.sampleVals[k]));
      return mx;
    },
    // The reconstruction at time t, straight from the samples (k = −1 … K).
    reconAt(t) {
      const r = window.__demo.last;
      if (r.s.method === 'ideal') return evalSum(r.bins, t, 0, 1)[0];
      return holdRecon(r.xs, r.s.method, r.s.fs, r.s.delta, t, 2)[1];
    },
    // The analytic line spectrum up to fMax, summed back in time on t = j·dt.
    lineSum(fMax, dt, N) {
      const r = window.__demo.last;
      return evalSum(imageLines(r, fMax), 0, dt, N);
    },
  };

  for (const id of ['OT', 'OA', 'OP', 'ST', 'SA', 'SP']) cv[id] = $('cv' + id);
  bind();
  schedule();
})();
