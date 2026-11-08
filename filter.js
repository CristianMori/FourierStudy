// Spectrum × Filter: draw a line spectrum X and a filter H on harmonic bins
// k = 0..K, and see the product Y = X·H as spectra and waveforms.
//
// Convention: one-sided cosine amplitudes. A spectrum (A, φ) stands for the
// real periodic signal x(t) = Σ_{k=0..K} A_k cos(2π k t/T + φ_k), i.e. the
// complex coefficients c_k = (A_k/2) e^{iφ_k} for k ≥ 1, c_0 = A_0, with
// Hermitian symmetry c_{-k} = conj(c_k) implied.
(() => {
  'use strict';

  const K = 64;
  const NB = K + 1;
  const DEG = Math.PI / 180;
  const WAVE_SAMPLES = 1024; // samples across the plotted time span

  const state = {
    xA: new Float64Array(NB), xP: new Float64Array(NB),
    hA: new Float64Array(NB), hP: new Float64Array(NB),
    yA: new Float64Array(NB), yP: new Float64Array(NB),
  };

  const wrapDeg = (p) => {
    const w = ((p + 180) % 360 + 360) % 360 - 180;
    return w === -180 && p > 0 ? 180 : w;
  };
  const clamp = (v, lo, hi) => Math.min(hi, Math.max(lo, v));

  function computeProduct() {
    for (let k = 0; k < NB; k++) {
      state.yA[k] = state.xA[k] * state.hA[k];
      state.yP[k] = wrapDeg(state.xP[k] + state.hP[k]);
    }
  }

  // Harmonic sum over t in [t0, t1) (t in periods), amplitudes scaled by ampScale(k).
  function synth(A, P, t0, t1, ampScale) {
    const out = new Float64Array(WAVE_SAMPLES);
    for (let k = 0; k < NB; k++) {
      const a = A[k] * (ampScale ? ampScale(k) : 1);
      if (a === 0) continue;
      const ph = P[k] * DEG;
      const w = 2 * Math.PI * k;
      for (let i = 0; i < WAVE_SAMPLES; i++) {
        const t = t0 + (t1 - t0) * i / WAVE_SAMPLES;
        out[i] += a * Math.cos(w * t + ph);
      }
    }
    return out;
  }

  // ---------- presets ----------

  const inputPresets = {
    clear(A, P) { A.fill(0); P.fill(0); },
    sine(A, P) { A.fill(0); P.fill(0); A[1] = 1; P[1] = -90; },
    // Σ_{k odd} sin(kωt)/k: a square wave (π/4 high), normalised so A_1 = 1.
    square(A, P) {
      A.fill(0); P.fill(0);
      for (let k = 1; k <= K; k += 2) { A[k] = 1 / k; P[k] = -90; }
    },
    // Σ (-1)^{k+1} sin(kωt)/k: a rising sawtooth.
    sawtooth(A, P) {
      A.fill(0); P.fill(0);
      for (let k = 1; k <= K; k++) { A[k] = 1 / k; P[k] = k % 2 ? -90 : 90; }
    },
    // Σ_{k odd} cos(kωt)/k²: a triangle with its peak at t = 0.
    triangle(A, P) {
      A.fill(0); P.fill(0);
      for (let k = 1; k <= K; k += 2) A[k] = 1 / (k * k);
    },
    // Flat spectrum: band-limited impulse train (c_k equal for all k, so A_0 = A_k / 2).
    impulse(A, P) { A.fill(1); P.fill(0); A[0] = 0.5; },
    random(A, P) {
      A.fill(0); P.fill(0);
      const top = 8 + Math.floor(Math.random() * 20);
      for (let k = 1; k <= top; k++) {
        A[k] = Math.random() ** 2 / Math.sqrt(k);
        P[k] = Math.round(Math.random() * 360 - 180);
      }
      const m = Math.max(...A);
      if (m > 0) for (let k = 0; k < NB; k++) A[k] /= m;
    },
  };

  const filterCtl = { shape: 'ideal', k0: 10, width: 6, preset: 'lowpass', drawn: false };

  // Each returns a gain for bin k; the phase is set to 0 unless noted.
  const filterGains = {
    allpass: () => 1,
    lowpass: (k, f) => f.shape === 'ideal'
      ? (k <= f.k0 ? 1 : 0)
      : 1 / Math.sqrt(1 + (k / f.k0) ** 8),
    highpass: (k, f) => f.shape === 'ideal'
      ? (k >= f.k0 ? 1 : 0)
      : (k === 0 ? 0 : 1 / Math.sqrt(1 + (f.k0 / k) ** 8)),
    bandpass: (k, f) => f.shape === 'ideal'
      ? (Math.abs(k - f.k0) <= f.width / 2 ? 1 : 0)
      : Math.exp(-(((k - f.k0) / (f.width / 2)) ** 2) * Math.LN2),
    notch: (k, f) => f.shape === 'ideal'
      ? (Math.abs(k - f.k0) <= f.width / 2 ? 0 : 1)
      : 1 - Math.exp(-(((k - f.k0) / (f.width / 2)) ** 2) * Math.LN2),
    comb: (k, f) => f.shape === 'ideal'
      ? (k % f.k0 === 0 ? 1 : 0)
      : (1 + Math.cos(2 * Math.PI * k / f.k0)) / 2,
    delay: () => 1,
  };

  function applyFilterPreset(name) {
    const g = filterGains[name];
    for (let k = 0; k < NB; k++) {
      state.hA[k] = g(k, filterCtl);
      // A delay of T/8 is a linear phase −2π k/8.
      state.hP[k] = name === 'delay' ? wrapDeg(-360 * k / 8) : 0;
    }
    filterCtl.preset = name;
    filterCtl.drawn = false;
  }

  // ---------- drawing helpers ----------

  let colors = {};
  function readColors() {
    const cs = getComputedStyle(document.documentElement);
    for (const n of ['ink', 'muted', 'grid', 'axis', 'panel', 'c-in', 'c-filt', 'c-out', 'c-ref']) {
      colors[n] = cs.getPropertyValue('--' + n).trim();
    }
  }

  const PAD = { l: 38, r: 8, t: 8, b: 18 };

  function setupCanvas(canvas) {
    const dpr = window.devicePixelRatio || 1;
    const w = canvas.clientWidth, h = canvas.clientHeight;
    if (canvas.width !== Math.round(w * dpr) || canvas.height !== Math.round(h * dpr)) {
      canvas.width = Math.round(w * dpr);
      canvas.height = Math.round(h * dpr);
    }
    const ctx = canvas.getContext('2d');
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.clearRect(0, 0, w, h);
    return { ctx, w, h, pw: w - PAD.l - PAD.r, ph: h - PAD.t - PAD.b };
  }

  function fmt(v) {
    if (v === 0) return '0';
    return (+v.toPrecision(3)).toString().replace('-', '−');
  }

  function niceCeil(v) {
    if (!(v > 1e-9)) return 1;
    const e = Math.floor(Math.log10(v));
    for (const m of [1, 1.2, 1.5, 2, 2.5, 3, 4, 5, 6, 8, 10]) {
      const c = m * 10 ** e;
      if (c >= v * 0.999) return c;
    }
    return 10 ** (e + 1);
  }

  function drawYAxis(g, min, max, ticks, refs) {
    const { ctx, pw, ph } = g;
    const yOf = (v) => PAD.t + (1 - (v - min) / (max - min)) * ph;
    ctx.font = '10px system-ui, sans-serif';
    ctx.textAlign = 'right';
    ctx.textBaseline = 'middle';
    for (const v of ticks) {
      const y = Math.round(yOf(v)) + 0.5;
      ctx.strokeStyle = v === 0 ? colors.axis : colors.grid;
      ctx.lineWidth = 1;
      ctx.beginPath(); ctx.moveTo(PAD.l, y); ctx.lineTo(PAD.l + pw, y); ctx.stroke();
      ctx.fillStyle = colors.muted;
      ctx.fillText(fmt(v), PAD.l - 5, y);
    }
    for (const v of refs || []) {
      const y = Math.round(yOf(v)) + 0.5;
      ctx.strokeStyle = colors['c-ref'];
      ctx.setLineDash([4, 4]);
      ctx.beginPath(); ctx.moveTo(PAD.l, y); ctx.lineTo(PAD.l + pw, y); ctx.stroke();
      ctx.setLineDash([]);
    }
    return yOf;
  }

  // ---------- spectrum plot (optionally drawable) ----------

  class SpecPlot {
    // opts: { values: () => array, kind: 'amp'|'phase', range: () => [min,max],
    //         ticks: () => [...], refs, color, dimBy: () => amplitude array,
    //         editable, snap, onEdit }
    constructor(canvas, opts) {
      this.canvas = canvas;
      this.o = opts;
      this.hover = null;
      this.drag = null;
      if (opts.editable) this.bindPointer();
    }

    geom() {
      const w = this.canvas.clientWidth, h = this.canvas.clientHeight;
      const pw = w - PAD.l - PAD.r, ph = h - PAD.t - PAD.b;
      return { pw, ph, bw: pw / NB };
    }

    toPoint(e) {
      const r = this.canvas.getBoundingClientRect();
      const x = e.clientX - r.left, y = e.clientY - r.top;
      const { ph, bw } = this.geom();
      const [min, max] = this.o.range();
      const bin = clamp(Math.floor((x - PAD.l) / bw), 0, K);
      let v = clamp(min + (1 - (y - PAD.t) / ph) * (max - min), min, max);
      if (this.o.snap) v = this.o.snap(v);
      return { bin, v };
    }

    setSegment(a, b) {
      const arr = this.o.values();
      if (a.bin === b.bin) { arr[b.bin] = b.v; return; }
      const step = b.bin > a.bin ? 1 : -1;
      for (let k = a.bin; k !== b.bin + step; k += step) {
        const t = (k - a.bin) / (b.bin - a.bin);
        arr[k] = a.v + t * (b.v - a.v);
      }
    }

    bindPointer() {
      const c = this.canvas;
      c.addEventListener('pointerdown', (e) => {
        if (e.pointerType === 'mouse' && e.button !== 0) return;
        if (this.drag) return;
        e.preventDefault();
        try { c.setPointerCapture(e.pointerId); } catch (_) { /* synthetic events */ }
        const p = this.toPoint(e);
        this.drag = { id: e.pointerId, last: p };
        this.setSegment(p, p);
        this.hover = p;
        this.o.onEdit();
      });
      c.addEventListener('pointermove', (e) => {
        if (this.drag && e.pointerId === this.drag.id) {
          const evs = e.getCoalescedEvents ? e.getCoalescedEvents() : [];
          for (const ev of evs.length ? evs : [e]) {
            const p = this.toPoint(ev);
            this.setSegment(this.drag.last, p);
            this.drag.last = p;
          }
          this.hover = this.drag.last;
          this.o.onEdit();
        } else if (!this.drag) {
          this.hover = this.toPoint(e);
          scheduleRender();
        }
      });
      const end = (e) => {
        if (this.drag && e.pointerId === this.drag.id) this.drag = null;
        if (e.type !== 'pointerup' || e.pointerType !== 'mouse') this.hover = null;
        scheduleRender();
      };
      c.addEventListener('pointerup', end);
      c.addEventListener('pointercancel', end);
      c.addEventListener('pointerleave', (e) => { if (!this.drag) { this.hover = null; scheduleRender(); } });
    }

    draw() {
      const g = setupCanvas(this.canvas);
      const { ctx, pw, ph } = g;
      const [min, max] = this.o.range();
      const yOf = drawYAxis(g, min, max, this.o.ticks(), this.o.refs);
      const arr = this.o.values();
      const dim = this.o.dimBy ? this.o.dimBy() : null;
      const bw = pw / NB;
      const gap = bw > 4 ? 1 : 0;
      const y0 = yOf(clamp(0, min, max));
      ctx.fillStyle = this.o.color();
      for (let k = 0; k < NB; k++) {
        const v = clamp(arr[k], min, max);
        const x = PAD.l + k * bw;
        const y = yOf(v);
        ctx.globalAlpha = dim && dim[k] < 1e-6 ? 0.22 : 1;
        if (Math.abs(y - y0) < 1 && v !== 0 && this.o.kind === 'amp') {
          ctx.fillRect(x + gap / 2, y0 - 1, bw - gap, 1);
        } else {
          ctx.fillRect(x + gap / 2, Math.min(y, y0), bw - gap, Math.abs(y - y0));
        }
        if (arr[k] > max) { // clipped: mark the top
          ctx.fillRect(x + gap / 2, PAD.t - 3, bw - gap, 2);
        }
      }
      ctx.globalAlpha = 1;

      // bin labels
      ctx.fillStyle = colors.muted;
      ctx.font = '10px system-ui, sans-serif';
      ctx.textAlign = 'center';
      ctx.textBaseline = 'top';
      const every = pw < 360 ? 16 : 8;
      for (let k = 0; k <= K; k += every) {
        ctx.fillText(String(k), PAD.l + (k + 0.5) * bw, PAD.t + ph + 4);
      }

      if (this.hover) {
        const { bin } = this.hover;
        const v = arr[bin];
        ctx.strokeStyle = colors.ink;
        ctx.globalAlpha = 0.35;
        ctx.strokeRect(PAD.l + bin * bw + 0.5, PAD.t + 0.5, bw - 1, ph - 1);
        ctx.globalAlpha = 1;
        const txt = this.o.kind === 'phase'
          ? `k ${bin} · ${Math.round(v)}°`
          : `k ${bin} · ${v.toFixed(3)}`;
        ctx.font = '11px system-ui, sans-serif';
        ctx.textAlign = 'right';
        ctx.textBaseline = 'top';
        ctx.fillStyle = colors.ink;
        ctx.fillText(txt, PAD.l + pw - 4, PAD.t + 2);
      }
    }
  }

  // ---------- waveform plot ----------

  function drawWave(canvas, series, t0, t1, scale, xLabels) {
    const g = setupCanvas(canvas);
    const { ctx, pw, ph } = g;
    const s = niceCeil(scale);
    const yOf = drawYAxis(g, -s, s, [-s, -s / 2, 0, s / 2, s]);

    ctx.fillStyle = colors.muted;
    ctx.font = '10px system-ui, sans-serif';
    ctx.textAlign = 'center';
    ctx.textBaseline = 'top';
    for (const [t, label] of xLabels) {
      const x = PAD.l + (t - t0) / (t1 - t0) * pw;
      ctx.strokeStyle = colors.grid;
      ctx.beginPath(); ctx.moveTo(Math.round(x) + 0.5, PAD.t); ctx.lineTo(Math.round(x) + 0.5, PAD.t + ph); ctx.stroke();
      ctx.fillText(label, clamp(x, PAD.l + 6, PAD.l + pw - 6), PAD.t + ph + 4);
    }

    ctx.save();
    ctx.beginPath();
    ctx.rect(PAD.l, PAD.t - 1, pw, ph + 2);
    ctx.clip();
    for (const sr of series) {
      ctx.strokeStyle = sr.color;
      ctx.globalAlpha = sr.alpha ?? 1;
      ctx.lineWidth = sr.width ?? 1.6;
      ctx.setLineDash(sr.dash || []);
      ctx.lineJoin = 'round';
      ctx.beginPath();
      const n = sr.data.length;
      for (let i = 0; i <= n; i++) {
        const x = PAD.l + i / n * pw;
        const y = yOf(sr.data[i % n]); // periodic: close the span with sample 0
        if (i === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
      }
      ctx.stroke();
    }
    ctx.restore();
    ctx.setLineDash([]);
    ctx.globalAlpha = 1;
  }

  const maxAbs = (a) => { let m = 0; for (const v of a) m = Math.max(m, Math.abs(v)); return m; };

  // ---------- wiring ----------

  const $ = (id) => document.getElementById(id);
  const phaseSnap = (v) => (Math.abs(v) < 6 ? 0 : Math.round(v));
  const gainSnap = (v) => (Math.abs(v - 1) < 0.03 ? 1 : v < 0.015 ? 0 : v);
  const ampSnap = (v) => (v < 0.01 ? 0 : v);
  const color = (n) => () => colors[n];
  const phaseRange = () => [-180, 180];
  const phaseTicks = () => [-180, -90, 0, 90, 180];

  const onInputEdit = () => { markInputPreset(null); recompute(); };
  const onFilterEdit = () => { filterCtl.drawn = true; markFilterPreset(); recompute(); };

  const resultMax = () => Math.max(1, niceCeil(maxAbs(state.yA)));

  const plots = [
    new SpecPlot($('in-amp'), {
      values: () => state.xA, kind: 'amp', range: () => [0, 1], ticks: () => [0, 0.5, 1],
      color: color('c-in'), editable: true, snap: ampSnap, onEdit: onInputEdit,
    }),
    new SpecPlot($('in-phase'), {
      values: () => state.xP, kind: 'phase', range: phaseRange, ticks: phaseTicks,
      color: color('c-in'), dimBy: () => state.xA, editable: true, snap: phaseSnap, onEdit: onInputEdit,
    }),
    new SpecPlot($('filt-amp'), {
      values: () => state.hA, kind: 'amp', range: () => [0, 1.5], ticks: () => [0, 0.5, 1, 1.5], refs: [1],
      color: color('c-filt'), editable: true, snap: gainSnap, onEdit: onFilterEdit,
    }),
    new SpecPlot($('filt-phase'), {
      values: () => state.hP, kind: 'phase', range: phaseRange, ticks: phaseTicks,
      color: color('c-filt'), dimBy: () => state.hA, editable: true, snap: phaseSnap, onEdit: onFilterEdit,
    }),
    new SpecPlot($('out-amp'), {
      values: () => state.yA, kind: 'amp', range: () => [0, resultMax()],
      ticks: () => { const m = resultMax(); return [0, m / 2, m]; },
      color: color('c-out'),
    }),
    new SpecPlot($('out-phase'), {
      values: () => state.yP, kind: 'phase', range: phaseRange, ticks: phaseTicks,
      color: color('c-out'), dimBy: () => state.yA,
    }),
  ];

  const periodLabels = [[0, '0'], [0.5, 'T/2'], [1, 'T'], [1.5, '3T/2'], [2, '2T']];
  const centredLabels = [[-0.5, '−T/2'], [-0.25, '−T/4'], [0, '0'], [0.25, 'T/4'], [0.5, 'T/2']];

  const waves = { x: null, y: null, h: null };

  function recompute() {
    computeProduct();
    waves.x = synth(state.xA, state.xP, 0, 2);
    waves.y = synth(state.yA, state.yP, 0, 2);
    // h(t) = Σ_{k=-K..K} H_k e^{ikωt} = H_0 + Σ_{k≥1} 2|H_k| cos(kωt + θ_k)
    waves.h = synth(state.hA, state.hP, -0.5, 0.5, (k) => (k === 0 ? 1 : 2));
    scheduleRender();
  }

  let pending = false;
  function scheduleRender() {
    if (pending) return;
    pending = true;
    requestAnimationFrame(render);
  }

  function render() {
    pending = false;
    for (const p of plots) p.draw();
    const mx = maxAbs(waves.x), my = maxAbs(waves.y);
    drawWave($('in-wave'), [{ data: waves.x, color: colors['c-in'] }], 0, 2, mx, periodLabels);
    drawWave($('out-wave'), [
      { data: waves.x, color: colors['c-in'], alpha: 0.45, width: 1.3, dash: [5, 4] },
      { data: waves.y, color: colors['c-out'], width: 1.8 },
    ], 0, 2, Math.max(mx, my), periodLabels);
    drawWave($('filt-impulse'), [{ data: waves.h, color: colors['c-filt'] }], -0.5, 0.5, maxAbs(waves.h), centredLabels);
  }

  function markInputPreset(name) {
    for (const b of document.querySelectorAll('#input-presets button')) {
      b.classList.toggle('active', b.dataset.preset === name && name !== 'zerophase');
    }
  }

  function markFilterPreset() {
    for (const b of document.querySelectorAll('#filter-presets button')) {
      b.classList.toggle('active', !filterCtl.drawn && b.dataset.preset === filterCtl.preset);
    }
    $('filter-status').textContent = filterCtl.drawn ? 'hand-drawn: pick a preset to use the controls again' : '';
  }

  function setInputPreset(name) {
    if (name === 'zerophase') state.xP.fill(0);
    else inputPresets[name](state.xA, state.xP);
    if (name !== 'zerophase') markInputPreset(name);
    recompute();
  }

  function setFilterPreset(name) {
    if (name === 'zerophase') {
      state.hP.fill(0);
      if (filterCtl.preset === 'delay') filterCtl.drawn = true;
    } else {
      applyFilterPreset(name);
    }
    markFilterPreset();
    recompute();
  }

  $('input-presets').addEventListener('click', (e) => {
    const b = e.target.closest('button[data-preset]');
    if (b) setInputPreset(b.dataset.preset);
  });
  $('filter-presets').addEventListener('click', (e) => {
    const b = e.target.closest('button[data-preset]');
    if (b) setFilterPreset(b.dataset.preset);
  });

  function reapplyFilter() {
    if (!filterCtl.drawn) { applyFilterPreset(filterCtl.preset); markFilterPreset(); recompute(); }
  }
  function setShape(shape) {
    filterCtl.shape = shape;
    $('shape-ideal').setAttribute('aria-pressed', String(shape === 'ideal'));
    $('shape-smooth').setAttribute('aria-pressed', String(shape === 'smooth'));
    reapplyFilter();
  }
  $('shape-ideal').addEventListener('click', () => setShape('ideal'));
  $('shape-smooth').addEventListener('click', () => setShape('smooth'));
  for (const id of ['k0', 'width']) {
    const input = $(id), out = $(id + '-out');
    input.addEventListener('input', () => {
      filterCtl[id] = +input.value;
      out.textContent = input.value;
      reapplyFilter();
    });
  }

  readColors();
  const mq = window.matchMedia('(prefers-color-scheme: dark)');
  (mq.addEventListener ? mq.addEventListener.bind(mq, 'change') : mq.addListener.bind(mq))(() => { readColors(); scheduleRender(); });
  if (window.ResizeObserver) {
    const ro = new ResizeObserver(() => scheduleRender());
    for (const c of document.querySelectorAll('canvas.plot')) ro.observe(c);
  }
  window.addEventListener('resize', scheduleRender);

  // Start with a square wave through an ideal low-pass at k₀ = 10.
  inputPresets.square(state.xA, state.xP);
  markInputPreset('square');
  applyFilterPreset('lowpass');
  markFilterPreset();
  recompute();

  // Handle for poking at the state from the console.
  window.filterDemo = { K, state, waves, setInputPreset, setFilterPreset, setShape, filterCtl };
})();
