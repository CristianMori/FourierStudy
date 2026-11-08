// Epicycles: image -> strokes -> one closed path -> x(t), y(t) -> two Fourier series -> two chains of circles.
// Plain JS and canvas, no libraries. Images are read with createImageBitmap / FileReader, so the canvas is never tainted.
(function () {
'use strict';

const TAU = Math.PI * 2;
const MAX_SIDE = 512;        // working image size
const MAX_STROKES = 1500;    // keep the longest if a texture gives more
const BASE_PERIOD = 12;      // seconds per lap at speed 1
const $ = id => document.getElementById(id);

// ---------- state ----------
const S = {
  work: null, w: 0, h: 0, gray: null, name: '',
  mode: 'lines',
  blurCache: null,
  skeleton: null, strokes: [], path: null, sig: null, ax: null, ay: null,
  N: 1024, fx: 1, fy: 1, Mx: 1, My: 1,
  t: 0, playing: true, speed: 1,
  prev: null, cprev: null, cam: null, follower: 'knife', rho: 10, fps: 0, frameMs: 0, drawMs: 0,
  dirty: 0, theme: null, layout: null, stats: {},
};

// ---------- small helpers ----------
function prep(c) {
  const dpr = window.devicePixelRatio || 1;
  const w = c.clientWidth, h = c.clientHeight;
  const pw = Math.max(1, Math.round(w * dpr)), ph = Math.max(1, Math.round(h * dpr));
  if (c.width !== pw || c.height !== ph) { c.width = pw; c.height = ph; }
  const g = c.getContext('2d');
  g.setTransform(dpr, 0, 0, dpr, 0, 0);
  g.clearRect(0, 0, w, h);
  return { g, w, h };
}

function readTheme() {
  const cs = getComputedStyle(document.documentElement);
  const v = n => cs.getPropertyValue(n).trim();
  S.theme = { ink: v('--ink'), muted: v('--muted'), line: v('--line'), grid: v('--grid'), axis: v('--axis'),
    x: v('--c-x'), y: v('--c-y'), pen: v('--c-pen'), bridge: v('--c-bridge'), bad: v('--c-bad'), ref: v('--c-ref'),
    dark: matchMedia('(prefers-color-scheme: dark)').matches };
}

function fmt(v, d) { return v.toFixed(d); }

// ---------- 2. image processing ----------
function setSource(img, iw, ih, name) {
  const sc = MAX_SIDE / Math.max(iw, ih);
  const w = Math.max(8, Math.round(iw * sc)), h = Math.max(8, Math.round(ih * sc));
  let c, d;
  for (let attempt = 0; attempt < 6; attempt++) {
    c = document.createElement('canvas');
    c.width = w; c.height = h;
    const g = c.getContext('2d', { willReadFrequently: true });
    g.fillStyle = '#fff'; g.fillRect(0, 0, w, h);      // transparent PNGs sit on white
    g.drawImage(img, 0, 0, w, h);
    d = g.getImageData(0, 0, w, h).data;
    if (d[3] === 255 && d[4 * (w * h - 1) + 3] === 255) break;   // a failed readback comes back all zero: draw again
  }
  const gray = new Float32Array(w * h);
  for (let i = 0; i < w * h; i++) gray[i] = 0.299 * d[4 * i] + 0.587 * d[4 * i + 1] + 0.114 * d[4 * i + 2];
  S.work = c; S.w = w; S.h = h; S.gray = gray; S.name = name; S.blurCache = null;
  $('src-name').textContent = name + ' (' + w + '×' + h + ')';
  layoutViews();
  invalidate(3);
}

function gaussBlur(src, w, h, sigma) {
  const r = Math.max(1, Math.ceil(sigma * 3));
  const k = new Float32Array(2 * r + 1);
  let sum = 0;
  for (let i = -r; i <= r; i++) { k[i + r] = Math.exp(-i * i / (2 * sigma * sigma)); sum += k[i + r]; }
  for (let i = 0; i < k.length; i++) k[i] /= sum;
  const tmp = new Float32Array(w * h), out = new Float32Array(w * h);
  for (let y = 0; y < h; y++) for (let x = 0; x < w; x++) {
    let a = 0;
    for (let i = -r; i <= r; i++) a += k[i + r] * src[y * w + Math.min(w - 1, Math.max(0, x + i))];
    tmp[y * w + x] = a;
  }
  for (let y = 0; y < h; y++) for (let x = 0; x < w; x++) {
    let a = 0;
    for (let i = -r; i <= r; i++) a += k[i + r] * tmp[Math.min(h - 1, Math.max(0, y + i)) * w + x];
    out[y * w + x] = a;
  }
  return out;
}

// Canny: Sobel gradient, non-maximum suppression along the gradient, hysteresis. thr in 0..1 of a robust gradient scale.
function cannyEdges(g, w, h, thr) {
  const mag = new Float32Array(w * h), dir = new Uint8Array(w * h);
  for (let y = 1; y < h - 1; y++) for (let x = 1; x < w - 1; x++) {
    const i = y * w + x;
    const gx = (g[i - w + 1] + 2 * g[i + 1] + g[i + w + 1]) - (g[i - w - 1] + 2 * g[i - 1] + g[i + w - 1]);
    const gy = (g[i + w - 1] + 2 * g[i + w] + g[i + w + 1]) - (g[i - w - 1] + 2 * g[i - w] + g[i - w + 1]);
    mag[i] = Math.hypot(gx, gy);
    let a = Math.atan2(gy, gx) * 180 / Math.PI; if (a < 0) a += 180;
    dir[i] = a < 22.5 || a >= 157.5 ? 0 : a < 67.5 ? 1 : a < 112.5 ? 2 : 3;
  }
  // neighbours across the edge for each quantised direction: horizontal gradient, diagonal \, vertical, diagonal /
  const offs = [[1, -1], [w + 1, -w - 1], [w, -w], [w - 1, -w + 1]];
  const nms = new Uint8Array(w * h);
  let mx = 0; const nz = [];
  for (let y = 1; y < h - 1; y++) for (let x = 1; x < w - 1; x++) {
    const i = y * w + x, m = mag[i];
    if (m <= 0) continue;
    const o = offs[dir[i]];
    if (m >= mag[i + o[0]] && m >= mag[i + o[1]]) { nms[i] = 1; if (m > 1) nz.push(m); if (m > mx) mx = m; }
  }
  nz.sort((a, b) => a - b);
  const scale = nz.length ? nz[Math.floor(nz.length * 0.9)] : mx;     // 90th percentile of the ridge strengths
  const hi = Math.max(thr * 2 * scale, 12), lo = 0.4 * hi;
  const out = new Uint8Array(w * h), stack = [];
  for (let i = 0; i < w * h; i++) if (nms[i] && mag[i] >= hi) { out[i] = 1; stack.push(i); }
  const n8 = [-w - 1, -w, -w + 1, -1, 1, w - 1, w, w + 1];
  while (stack.length) {
    const i = stack.pop();
    for (const o of n8) {
      const j = i + o;
      if (!out[j] && nms[j] && mag[j] >= lo) { out[j] = 1; stack.push(j); }
    }
  }
  return out;
}

function darkMask(g, w, h, level) {
  const out = new Uint8Array(w * h);
  for (let i = 0; i < w * h; i++) out[i] = g[i] < level ? 1 : 0;
  return out;
}

// Zhang-Suen thinning: peel boundary pixels until every line is one pixel wide.
function thin(m, w, h) {
  let list = [];
  for (let i = 0; i < w * h; i++) {
    if (!m[i]) continue;
    const x = i % w, y = (i / w) | 0;
    if (x === 0 || y === 0 || x === w - 1 || y === h - 1) m[i] = 0; else list.push(i);
  }
  let changed = true;
  while (changed) {
    changed = false;
    for (let step = 0; step < 2; step++) {
      const del = [];
      for (const i of list) {
        const p2 = m[i - w], p3 = m[i - w + 1], p4 = m[i + 1], p5 = m[i + w + 1];
        const p6 = m[i + w], p7 = m[i + w - 1], p8 = m[i - 1], p9 = m[i - w - 1];
        const B = p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9;
        if (B < 2 || B > 6) continue;
        const A = (p2 === 0 && p3 === 1) + (p3 === 0 && p4 === 1) + (p4 === 0 && p5 === 1) + (p5 === 0 && p6 === 1) +
                  (p6 === 0 && p7 === 1) + (p7 === 0 && p8 === 1) + (p8 === 0 && p9 === 1) + (p9 === 0 && p2 === 1);
        if (A !== 1) continue;
        if (step === 0) { if (p2 && p4 && p6) continue; if (p4 && p6 && p8) continue; }
        else { if (p2 && p4 && p8) continue; if (p2 && p6 && p8) continue; }
        del.push(i);
      }
      if (del.length) {
        changed = true;
        for (const i of del) m[i] = 0;
        list = list.filter(i => m[i]);
      }
    }
  }
  return m;
}

// Drop corner pixels: a pixel whose neighbours already touch each other is redundant for connectivity.
// This removes the false junctions staircases leave behind, so long lines are not chopped into short pieces.
function dropCorners(m, w, h) {
  const OFF = [-w - 1, -w, -w + 1, 1, w + 1, w, w - 1, -1];
  const DX = [-1, 0, 1, 1, 1, 0, -1, -1], DY = [-1, -1, -1, 0, 1, 1, 1, 0];
  for (let pass = 0; pass < 2; pass++) for (let i = w + 1; i < w * h - w - 1; i++) {
    if (!m[i]) continue;
    const nb = [];
    for (let k = 0; k < 8; k++) if (m[i + OFF[k]]) nb.push(k);
    if (nb.length < 2) continue;
    const comp = nb.map((_, a) => a);
    const find = a => { while (comp[a] !== a) a = comp[a]; return a; };
    for (let a = 0; a < nb.length; a++) for (let b = a + 1; b < nb.length; b++)
      if (Math.max(Math.abs(DX[nb[a]] - DX[nb[b]]), Math.abs(DY[nb[a]] - DY[nb[b]])) <= 1) comp[find(a)] = find(b);
    const roots = new Set(nb.map((_, a) => find(a)));
    if (roots.size === 1) m[i] = 0;
  }
  return m;
}

// Walk the one-pixel skeleton from endpoints and junctions (and then around closed loops) into pixel chains.
function traceStrokes(m, w, h) {
  const OFF = [-w - 1, -w, -w + 1, -1, 1, w - 1, w, w + 1];
  const deg = new Uint8Array(w * h), visited = new Uint8Array(w * h);
  const px = [];
  for (let i = 0; i < w * h; i++) if (m[i]) {
    px.push(i);
    let d = 0; for (const o of OFF) d += m[i + o];
    deg[i] = d;
  }
  const chains = [];
  const toPt = i => [(i % w) + 0.5, ((i / w) | 0) + 0.5];
  function walk(s, n) {
    const pts = [s]; let prev = s, cur = n;
    for (;;) {
      if (deg[cur] !== 2) { pts.push(cur); break; }
      if (visited[cur]) break;
      visited[cur] = 1; pts.push(cur);
      let next = -1;
      for (const o of OFF) { const j = cur + o; if (m[j] && j !== prev && (deg[j] !== 2 || !visited[j])) { next = j; break; } }
      if (next < 0) break;
      prev = cur; cur = next;
    }
    return pts;
  }
  for (const s of px) {
    if (deg[s] === 2) continue;
    for (const o of OFF) {
      const n = s + o;
      if (!m[n]) continue;
      if (deg[n] === 2 ? visited[n] : n < s) continue;     // node-to-node steps would be traced twice
      chains.push(walk(s, n));
    }
  }
  for (const s of px) {      // closed loops with no endpoint or junction
    if (deg[s] !== 2 || visited[s]) continue;
    visited[s] = 1;
    const pts = [s]; let prev = s, cur = -1;
    for (const o of OFF) if (m[s + o]) { cur = s + o; break; }
    while (cur >= 0 && cur !== s && !visited[cur]) {
      visited[cur] = 1; pts.push(cur);
      let next = -1;
      for (const o of OFF) { const j = cur + o; if (m[j] && j !== prev && !visited[j]) { next = j; break; } }
      prev = cur; cur = next;
    }
    pts.push(s);
    chains.push(pts);
  }
  return chains.map(c => c.map(toPt));
}

// Ramer-Douglas-Peucker, iterative.
function rdp(pts, eps) {
  const n = pts.length;
  if (n < 3) return pts.slice();
  const keep = new Uint8Array(n); keep[0] = keep[n - 1] = 1;
  const stack = [[0, n - 1]];
  while (stack.length) {
    const [a, b] = stack.pop();
    const ax = pts[a][0], ay = pts[a][1], bx = pts[b][0], by = pts[b][1];
    const dx = bx - ax, dy = by - ay, len = Math.hypot(dx, dy);
    let best = -1, bd = eps;
    for (let i = a + 1; i < b; i++) {
      const d = len > 1e-9 ? Math.abs(dy * (pts[i][0] - ax) - dx * (pts[i][1] - ay)) / len : Math.hypot(pts[i][0] - ax, pts[i][1] - ay);
      if (d > bd) { bd = d; best = i; }
    }
    if (best >= 0) { keep[best] = 1; stack.push([a, best], [best, b]); }
  }
  return pts.filter((_, i) => keep[i]);
}

function polyLen(p) { let l = 0; for (let i = 1; i < p.length; i++) l += Math.hypot(p[i][0] - p[i - 1][0], p[i][1] - p[i - 1][1]); return l; }

// ---------- 3. one closed path ----------
// Order strokes (and pick which end each is entered from) to keep the bridges short:
// greedy nearest neighbour, then 2-opt and Or-opt until nothing improves.
function connect(strokes) {
  const n = strokes.length;
  const A = strokes.map(s => s[0]), B = strokes.map(s => s[s.length - 1]);
  const d = (p, q) => Math.hypot(p[0] - q[0], p[1] - q[1]);
  let first = 0;
  for (let i = 1; i < n; i++) if (strokes[i].length > strokes[first].length) first = i;
  const used = new Uint8Array(n);
  let tour = [{ s: first, r: false }]; used[first] = 1;
  for (let c = 1; c < n; c++) {
    const e = tour[tour.length - 1], tip = e.r ? A[e.s] : B[e.s];
    let bi = -1, br = false, bd = Infinity;
    for (let i = 0; i < n; i++) {
      if (used[i]) continue;
      const da = d(tip, A[i]), db = d(tip, B[i]);
      if (da < bd) { bd = da; bi = i; br = false; }
      if (db < bd) { bd = db; bi = i; br = true; }
    }
    used[bi] = 1; tour.push({ s: bi, r: br });
  }
  const st = e => e.r ? B[e.s] : A[e.s];       // where the pen enters
  const en = e => e.r ? A[e.s] : B[e.s];       // where the pen leaves
  const flip = e => ({ s: e.s, r: !e.r });
  const cost = () => { let c = 0; for (let i = 0; i < tour.length; i++) c += d(en(tour[i]), st(tour[(i + 1) % tour.length])); return c; };
  const greedyCost = cost();
  const t0 = performance.now();
  let improved = true, pass = 0;
  while (improved && n > 2 && pass++ < 40 && performance.now() - t0 < 500) {
    improved = false;
    // 2-opt: reverse tour[i..j], flipping each stroke in it
    for (let i = 0; i < n; i++) for (let j = i; j < n; j++) {
      if (i === 0 && j === n - 1) continue;
      const pe = en(tour[(i - 1 + n) % n]), ns = st(tour[(j + 1) % n]);
      const delta = d(pe, en(tour[j])) + d(st(tour[i]), ns) - d(pe, st(tour[i])) - d(en(tour[j]), ns);
      if (delta < -1e-9) {
        const seg = tour.slice(i, j + 1).reverse().map(flip);
        for (let k = 0; k < seg.length; k++) tour[i + k] = seg[k];
        improved = true;
      }
    }
    // Or-opt: move one stroke elsewhere, either way round
    for (let i = 0; i < n; i++) {
      const pe = en(tour[(i - 1 + n) % n]), ns = st(tour[(i + 1) % n]), e = tour[i];
      const removed = d(pe, st(e)) + d(en(e), ns) - d(pe, ns);
      let bg = 1e-9, bj = -1, bf = false;
      for (let j = 0; j < n; j++) {
        if (j === i || j === (i - 1 + n) % n) continue;
        const a = en(tour[j]), b = st(tour[(j + 1) % n]), base = d(a, b);
        const g1 = removed - (d(a, st(e)) + d(en(e), b) - base);
        const g2 = removed - (d(a, en(e)) + d(st(e), b) - base);
        if (g1 > bg) { bg = g1; bj = j; bf = false; }
        if (g2 > bg) { bg = g2; bj = j; bf = true; }
      }
      if (bj >= 0) {
        const mv = bf ? flip(e) : e, after = tour[bj];
        tour.splice(i, 1);
        tour.splice(tour.indexOf(after) + 1, 0, mv);
        improved = true;
      }
    }
  }
  // flatten; br[i] = 1 when the segment from point i to point i+1 (cyclic) is a bridge
  const xs = [], ys = [], br = [];
  let strokeLen = 0;
  for (const e of tour) {
    const pts = e.r ? strokes[e.s].slice().reverse() : strokes[e.s];
    strokeLen += polyLen(pts);
    for (let k = 0; k < pts.length; k++) { xs.push(pts[k][0]); ys.push(pts[k][1]); br.push(k === pts.length - 1 ? 1 : 0); }
  }
  const P = xs.length;
  let bridgeLen = 0;
  for (let i = 0; i < P; i++) if (br[i]) bridgeLen += Math.hypot(xs[(i + 1) % P] - xs[i], ys[(i + 1) % P] - ys[i]);
  return { x: Float64Array.from(xs), y: Float64Array.from(ys), br: Uint8Array.from(br), P, tour, strokes: tour.length,
           strokeLen, bridgeLen, greedyCost, finalCost: cost() };
}

// ---------- 4. resample uniformly by arc length ----------
function resample(path, N) {
  const { x, y, br, P } = path;
  const cum = new Float64Array(P + 1);
  for (let i = 0; i < P; i++) cum[i + 1] = cum[i] + Math.hypot(x[(i + 1) % P] - x[i], y[(i + 1) % P] - y[i]);
  const L = cum[P];
  const sx = new Float64Array(N), sy = new Float64Array(N), sb = new Uint8Array(N);
  let j = 0;
  for (let n = 0; n < N; n++) {
    const s = n * L / N;
    while (j < P - 1 && cum[j + 1] <= s) j++;
    const seg = cum[j + 1] - cum[j], f = seg > 0 ? (s - cum[j]) / seg : 0;
    const k = (j + 1) % P;
    sx[n] = x[j] + f * (x[k] - x[j]); sy[n] = y[j] + f * (y[k] - y[j]); sb[n] = br[j];
  }
  return { N, x: sx, y: sy, b: sb, L };
}

// ---------- 5. FFT ----------
// Iterative radix-2 Cooley-Tukey, O(N log N), forward transform X_k = sum x_n e^{-i2πkn/N}, in place.
const twCache = {};
function fft(re, im) {
  const n = re.length;
  let tw = twCache[n];
  if (!tw) {
    tw = { c: new Float64Array(n / 2), s: new Float64Array(n / 2) };
    for (let i = 0; i < n / 2; i++) { tw.c[i] = Math.cos(TAU * i / n); tw.s[i] = -Math.sin(TAU * i / n); }
    twCache[n] = tw;
  }
  for (let i = 1, j = 0; i < n; i++) {
    let bit = n >> 1;
    for (; j & bit; bit >>= 1) j ^= bit;
    j ^= bit;
    if (i < j) { let t = re[i]; re[i] = re[j]; re[j] = t; t = im[i]; im[i] = im[j]; im[j] = t; }
  }
  for (let len = 2; len <= n; len <<= 1) {
    const half = len >> 1, step = n / len;
    for (let i = 0; i < n; i += len) for (let k = 0; k < half; k++) {
      const wr = tw.c[k * step], wi = tw.s[k * step];
      const a = i + k, b = a + half;
      const xr = re[b] * wr - im[b] * wi, xi = re[b] * wi + im[b] * wr;
      re[b] = re[a] - xr; im[b] = im[a] - xi; re[a] += xr; im[a] += xi;
    }
  }
}

// One real signal -> series x0 + sum_k R cos(2π k t + PH), circles sorted by radius.
function analyse(sig) {
  const N = sig.length, H = N / 2;
  const re = Float64Array.from(sig), im = new Float64Array(N);
  fft(re, im);
  const amp = new Float64Array(H + 1), ph = new Float64Array(H + 1);
  for (let k = 0; k <= H; k++) {
    const m = Math.hypot(re[k], im[k]) / N;
    amp[k] = (k === 0 || k === H) ? m : 2 * m;      // a real signal splits its power between bin k and N-k
    ph[k] = Math.atan2(im[k], re[k]);
  }
  const order = [];
  for (let k = 1; k <= H; k++) order.push(k);
  order.sort((a, b) => amp[b] - amp[a]);
  const K = Float64Array.from(order), R = new Float64Array(H), PH = new Float64Array(H), rank = new Int32Array(H + 1);
  order.forEach((k, i) => { R[i] = amp[k]; PH[i] = ph[k]; rank[k] = i; });
  return { N, H, x0: re[0] / N, K, R, PH, amp, ph, rank, re, im };
}

// Signal made of only the M largest circles, by inverse FFT of the kept bins.
function partialSignal(a, M) {
  const N = a.N, re = new Float64Array(N), im = new Float64Array(N);
  re[0] = a.re[0]; im[0] = -a.im[0];
  for (let i = 0; i < M; i++) {
    const k = a.K[i];
    re[k] = a.re[k]; im[k] = -a.im[k];
    if (k < a.H) { re[N - k] = a.re[N - k]; im[N - k] = -a.im[N - k]; }
  }
  fft(re, im);                          // inverse = conj(fft(conj(X))) / N
  for (let n = 0; n < N; n++) re[n] /= N;
  return re;
}

// ---------- 6. one chain: returns the tip's coordinate A (the pen's x or y).
// pos (optional) receives each centre as (A, B) in the chain's own plane: A = sum R cos, B = sum R sin.
function chainTip(a, t, M, pos) {
  let A = a.x0, B = 0;
  if (pos) { pos[0] = A; pos[1] = B; }
  for (let i = 0; i < M; i++) {
    const ang = TAU * a.K[i] * t + a.PH[i];
    A += a.R[i] * Math.cos(ang); B += a.R[i] * Math.sin(ang);
    if (pos) { pos[2 * i + 2] = A; pos[2 * i + 3] = B; }
  }
  return A;
}

// ---------- pipeline ----------
function params() {
  return { blur: +$('p-blur').value, thr: +$('p-thr').value, minLen: +$('p-min').value, tol: +$('p-tol').value };
}

function runVectorize() {
  const t0 = performance.now(), p = params();
  const { w, h } = S;
  if (!S.blurCache || S.blurCache.sigma !== p.blur) S.blurCache = { sigma: p.blur, img: gaussBlur(S.gray, w, h, p.blur) };
  const blurred = S.blurCache.img;
  const raw = S.mode === 'edges' ? cannyEdges(blurred, w, h, p.thr / 100) : darkMask(blurred, w, h, p.thr / 100 * 255);
  const sk = dropCorners(thin(raw, w, h), w, h);
  let skPx = 0; for (let i = 0; i < sk.length; i++) skPx += sk[i];
  const chains = traceStrokes(sk, w, h);
  let strokes = chains.map(c => rdp(c, p.tol)).filter(c => c.length >= 2 && polyLen(c) >= Math.max(p.minLen, 1.5));
  const total = strokes.length;
  if (total > MAX_STROKES) {
    strokes = strokes.map(s => [polyLen(s), s]).sort((a, b) => b[0] - a[0]).slice(0, MAX_STROKES).map(a => a[1]);
  }
  S.skeleton = sk; S.strokes = strokes;
  let pts = 0; for (const s of strokes) pts += s.length;
  S.stats.vec = { ms: performance.now() - t0, skPx, chains: chains.length, strokes: strokes.length, dropped: total - strokes.length, pts };
}

function runConnect() {
  if (S.strokes.length === 0) { S.path = null; S.stats.conn = null; return; }
  const t0 = performance.now();
  const path = connect(S.strokes);
  path.ms = performance.now() - t0;
  S.path = path;
}

function runSignals() {
  const path = S.path;
  if (!path || path.P < 2) { S.sig = null; S.ax = S.ay = null; return; }
  S.sig = resample(path, S.N);
  if (!(S.sig.L > 1e-6)) { S.sig = null; S.ax = S.ay = null; return; }
  const t0 = performance.now();
  S.ax = analyse(S.sig.x);
  S.ay = analyse(S.sig.y);
  S.stats.fftMs = performance.now() - t0;
  // timing over a few repeats, less noisy than the single run above
  const t1 = performance.now(); const reps = 5;
  for (let r = 0; r < reps; r++) { analyse(S.sig.x); analyse(S.sig.y); }
  S.stats.fftRepMs = (performance.now() - t1) / reps;
  // with every circle the series must give back the resampled path
  let worst = 0; const step = Math.max(1, S.N / 256);
  for (let n = 0; n < S.N; n += step) {
    const t = n / S.N;
    worst = Math.max(worst, Math.hypot(chainTip(S.ax, t, S.ax.H) - S.sig.x[n], chainTip(S.ay, t, S.ay.H) - S.sig.y[n]));
  }
  S.stats.maxErr = worst;
  updateM();
}

function updateM() {
  if (!S.ax) return;
  const H = S.ax.H;
  const cnt = f => Math.max(1, Math.min(H, Math.round(Math.exp(f * Math.log(H)))));
  S.Mx = cnt(S.fx); S.My = cnt(S.fy);
  $('o-mx').textContent = S.Mx + ' / ' + H;
  $('o-my').textContent = S.My + ' / ' + H;
  S.partX = partialSignal(S.ax, S.Mx);
  S.partY = partialSignal(S.ay, S.My);
  buildCams();
  resetTrail();
  drawSignals(); drawSpectra();
}

// ---------- drawing: vectorize views ----------
function xform(cw, ch) {
  const s = Math.min(cw / S.w, ch / S.h);
  return { s, ox: (cw - S.w * s) / 2, oy: (ch - S.h * s) / 2 };
}
function faintImage(g, xf, alpha) {
  g.globalAlpha = alpha; g.imageSmoothingEnabled = true;
  g.drawImage(S.work, xf.ox, xf.oy, S.w * xf.s, S.h * xf.s);
  g.globalAlpha = 1;
}
function emptyMsg(g, w, h, msg) {
  g.fillStyle = S.theme.muted; g.font = '13px system-ui, sans-serif'; g.textAlign = 'center';
  g.fillText(msg, w / 2, h / 2); g.textAlign = 'left';
}

function layoutViews() {
  if (!S.w) return;
  for (const id of ['cv-src', 'cv-edge', 'cv-poly']) {
    const c = $(id), cw = c.clientWidth;
    c.style.height = Math.min(360, Math.round(cw * S.h / S.w)) + 'px';
  }
  const c = $('cv-path');
  c.style.height = Math.min(520, Math.round(c.clientWidth * S.h / S.w)) + 'px';
}

function drawVectorViews() {
  const th = S.theme;
  { // image
    const { g, w, h } = prep($('cv-src')), xf = xform(w, h);
    faintImage(g, xf, 1);
  }
  { // thinned detector output
    const { g, w, h } = prep($('cv-edge')), xf = xform(w, h);
    faintImage(g, xf, 0.18);
    const m = document.createElement('canvas'); m.width = S.w; m.height = S.h;
    const mg = m.getContext('2d'), id = mg.createImageData(S.w, S.h);
    const col = th.dark ? [230, 231, 234] : [29, 31, 35];
    for (let i = 0; i < S.skeleton.length; i++) if (S.skeleton[i]) {
      id.data[4 * i] = col[0]; id.data[4 * i + 1] = col[1]; id.data[4 * i + 2] = col[2]; id.data[4 * i + 3] = 255;
    }
    mg.putImageData(id, 0, 0);
    g.imageSmoothingEnabled = true;
    g.drawImage(m, xf.ox, xf.oy, S.w * xf.s, S.h * xf.s);
    if (!S.stats.vec.skPx) emptyMsg(g, w, h, 'nothing found: lower the threshold');
  }
  { // polylines
    const { g, w, h } = prep($('cv-poly')), xf = xform(w, h);
    faintImage(g, xf, 0.18);
    g.lineWidth = 1.4; g.lineJoin = 'round'; g.strokeStyle = th.pen;
    for (const s of S.strokes) {
      g.beginPath();
      s.forEach((p, i) => i ? g.lineTo(xf.ox + p[0] * xf.s, xf.oy + p[1] * xf.s) : g.moveTo(xf.ox + p[0] * xf.s, xf.oy + p[1] * xf.s));
      g.stroke();
    }
    g.fillStyle = th.ink;
    for (const s of S.strokes) for (const p of s) g.fillRect(xf.ox + p[0] * xf.s - 1, xf.oy + p[1] * xf.s - 1, 2, 2);
    if (!S.strokes.length) emptyMsg(g, w, h, 'no strokes: lower the threshold or min stroke');
  }
  const v = S.stats.vec;
  $('st-vec').textContent = v.strokes + ' strokes, ' + v.pts + ' vertices from ' + v.skPx + ' skeleton pixels' +
    (v.dropped ? ' (' + v.dropped + ' shortest dropped to keep ' + MAX_STROKES + ')' : '') + ', ' + fmt(v.ms, 0) + ' ms.';
}

function drawConnectView() {
  const th = S.theme, { g, w, h } = prep($('cv-path')), xf = xform(w, h);
  faintImage(g, xf, 0.12);
  const path = S.path;
  if (!path) { emptyMsg(g, w, h, 'nothing to connect'); $('st-conn').textContent = 'No strokes, so no path.'; return; }
  const X = v => xf.ox + v * xf.s, Y = v => xf.oy + v * xf.s, P = path.P;
  const showB = $('c-bridges').checked;
  // strokes, coloured by position along the path
  g.lineWidth = 2; g.lineJoin = 'round';
  let i = 0;
  while (i < P) {
    let j = i; g.beginPath(); g.moveTo(X(path.x[i]), Y(path.y[i]));
    while (!path.br[j]) { j++; g.lineTo(X(path.x[j]), Y(path.y[j])); }
    g.strokeStyle = 'hsl(' + Math.round(210 * i / P) + ',' + (th.dark ? '75%,65%' : '70%,40%') + ')';
    g.stroke();
    i = j + 1;
  }
  if (showB) {
    g.strokeStyle = th.bridge; g.lineWidth = 1.2; g.setLineDash([4, 3]);
    for (let k = 0; k < P; k++) if (path.br[k]) {
      const m = (k + 1) % P;
      g.beginPath(); g.moveTo(X(path.x[k]), Y(path.y[k])); g.lineTo(X(path.x[m]), Y(path.y[m])); g.stroke();
    }
    g.setLineDash([]);
  }
  g.fillStyle = th.ink; g.beginPath(); g.arc(X(path.x[0]), Y(path.y[0]), 4, 0, TAU); g.fill();
  const pct = 100 * path.bridgeLen / (path.bridgeLen + path.strokeLen || 1);
  $('st-conn').textContent = path.strokes + ' strokes joined by ' + path.strokes + ' bridges: ' + fmt(path.strokeLen, 0) + ' px of drawing + ' +
    fmt(path.bridgeLen, 0) + ' px of bridge (' + fmt(pct, 1) + '% of the path). Ordering shortened the bridges from ' +
    fmt(path.greedyCost, 0) + ' to ' + fmt(path.finalCost, 0) + ' px (nearest neighbour, then 2-opt + Or-opt), ' + fmt(path.ms, 0) + ' ms. Black dot: start.';
}

// ---------- drawing: signals and spectra ----------
function niceTicks(lo, hi, n) {
  const span = hi - lo, raw = span / n, p = Math.pow(10, Math.floor(Math.log10(raw)));
  const step = [1, 2, 5, 10].map(m => m * p).find(s => s >= raw) || raw;
  const out = []; for (let v = Math.ceil(lo / step) * step; v <= hi + 1e-9; v += step) out.push(v);
  return out;
}

function timePlot(id, sig, part, colour) {
  const th = S.theme, { g, w, h } = prep($(id));
  const L = 40, R = 6, T = 8, B = 18;
  let lo = Infinity, hi = -Infinity;
  for (let i = 0; i < sig.length; i++) { lo = Math.min(lo, sig[i]); hi = Math.max(hi, sig[i]); }
  const pad = (hi - lo) * 0.06 || 1; lo -= pad; hi += pad;
  const px = t => L + t * (w - L - R), py = v => T + (hi - v) / (hi - lo) * (h - T - B);
  g.font = '11px system-ui, sans-serif'; g.fillStyle = th.muted; g.strokeStyle = th.grid; g.lineWidth = 1;
  for (const v of niceTicks(lo, hi, 3)) { g.beginPath(); g.moveTo(L, py(v)); g.lineTo(w - R, py(v)); g.stroke(); g.fillText(String(Math.round(v)), 4, py(v) + 4); }
  [[0, '0'], [0.5, 'T/2'], [1, 'T']].forEach(([t, s]) => { g.textAlign = t === 1 ? 'right' : t === 0 ? 'left' : 'center'; g.fillText(s, px(t), h - 4); });
  g.textAlign = 'left';
  const draw = (arr, st, wd, dash) => {
    g.beginPath(); g.strokeStyle = st; g.lineWidth = wd; g.setLineDash(dash);
    for (let i = 0; i <= arr.length; i++) { const v = arr[i % arr.length]; i ? g.lineTo(px(i / arr.length), py(v)) : g.moveTo(px(0), py(v)); }
    g.stroke(); g.setLineDash([]);
  };
  draw(sig, colour, 1.6, []);
  if (part) draw(part, th.ref, 1.4, [4, 3]);
}

function drawSignals() {
  if (!S.sig) { for (const id of ['sig-x', 'sig-y']) { const { g, w, h } = prep($(id)); emptyMsg(g, w, h, 'no signal'); } return; }
  timePlot('sig-x', S.sig.x, S.partX, S.theme.x);
  timePlot('sig-y', S.sig.y, S.partY, S.theme.y);
  $('o-n').textContent = S.N;
  $('st-sig').textContent = S.N + ' samples, closed path length ' + fmt(S.sig.L, 0) + ' px, one sample every ' + fmt(S.sig.L / S.N, 2) + ' px of path.';
}

function specPlots(a, ida, idp, colour, M) {
  const th = S.theme;
  const H = a.H, lh = Math.log(H) || 1;
  const L = 44, R = 8, T = 8;
  { // amplitude, log-log
    const { g, w, h } = prep($(ida)), B = 18;
    const top = Math.pow(10, Math.ceil(Math.log10(a.R[0] * 1.001))), dec = 4;
    const px = k => L + Math.log(k) / lh * (w - L - R);
    const py = v => T + (1 - Math.max(0, (Math.log10(v / top) + dec)) / dec) * (h - T - B);
    g.font = '11px system-ui, sans-serif'; g.fillStyle = th.muted; g.strokeStyle = th.grid; g.lineWidth = 1;
    for (let d = 0; d <= dec; d++) {
      const v = top * Math.pow(10, -d), y = py(v);
      g.beginPath(); g.moveTo(L, y); g.lineTo(w - R, y); g.stroke();
      g.fillText(v >= 1 ? String(v) : String(+v.toFixed(4)), 4, y + 4);
    }
    const ticks = [1, 10, 100, 1000].filter(k => k <= H);
    if (ticks[ticks.length - 1] !== H) ticks.push(H);
    ticks.forEach(k => { g.textAlign = k === H ? 'right' : 'center'; g.fillText(String(k), px(k), h - 4); });
    g.textAlign = 'left';
    const y0 = py(top * Math.pow(10, -dec));
    for (const used of [false, true]) {
      g.strokeStyle = used ? colour : th.axis; g.globalAlpha = used ? 1 : 0.5; g.lineWidth = 1.2;
      g.beginPath();
      for (let k = 1; k <= H; k++) {
        if ((a.rank[k] < M) !== used) continue;
        const x = px(k); g.moveTo(x, y0); g.lineTo(x, py(a.amp[k]));
      }
      g.stroke();
    }
    g.globalAlpha = 1;
    g.fillStyle = th.muted; g.textAlign = 'right'; g.fillText('k (turns per lap) →', w - R, T + 8); g.textAlign = 'left';
  }
  { // phase, log k, linear degrees; only where the amplitude is not numerical noise
    const { g, w, h } = prep($(idp)), B = 16;
    const px = k => L + Math.log(k) / lh * (w - L - R);
    const py = v => T + (1 - (v + Math.PI) / TAU) * (h - T - B);
    g.font = '11px system-ui, sans-serif'; g.fillStyle = th.muted; g.strokeStyle = th.grid; g.lineWidth = 1;
    for (const d of [-180, 0, 180]) { const y = py(d * Math.PI / 180); g.beginPath(); g.moveTo(L, y); g.lineTo(w - R, y); g.stroke(); g.fillText(String(d), 4, y + 4); }
    const cut = a.R[0] * 1e-3;
    for (const used of [false, true]) {
      g.strokeStyle = used ? colour : th.axis; g.globalAlpha = used ? 1 : 0.5; g.lineWidth = 1.2;
      g.beginPath();
      for (let k = 1; k <= H; k++) {
        if ((a.rank[k] < M) !== used || a.amp[k] < cut) continue;
        const x = px(k); g.moveTo(x, py(0)); g.lineTo(x, py(a.ph[k]));
      }
      g.stroke();
    }
    g.globalAlpha = 1;
  }
}

function drawSpectra() {
  if (!S.ax) {
    for (const id of ['sp-xa', 'sp-xp', 'sp-ya', 'sp-yp']) { const { g, w, h } = prep($(id)); emptyMsg(g, w, h, 'no spectrum'); }
    $('st-fft').textContent = ''; return;
  }
  specPlots(S.ax, 'sp-xa', 'sp-xp', S.theme.x, S.Mx);
  specPlots(S.ay, 'sp-ya', 'sp-yp', S.theme.y, S.My);
  $('st-fft').textContent = 'Iterative radix-2 FFT, O(N log N): both transforms of N = ' + S.N + ' took ' + fmt(S.stats.fftRepMs, 2) +
    ' ms. Mean position (DC): x₀ = ' + fmt(S.ax.x0, 1) + ' px, y₀ = ' + fmt(S.ay.x0, 1) + ' px. Phase is shown only for circles above 0.1% of the largest. ' +
    'With all ' + S.ax.H + ' circles the pen is within ' + fmt(S.stats.maxErr, 6) + ' px of the resampled path.';
}

// ---------- 6. stage ----------
function layoutStage() {
  const c = $('stage'), W = c.parentElement.clientWidth;
  if (!S.w) return null;
  const M0 = Math.max(70, Math.min(190, Math.round(W * 0.27)));
  const dwMax = W - M0, dhMax = Math.min(dwMax, 640);
  const s = Math.min(dwMax / S.w, dhMax / S.h);
  const dw = S.w * s, dh = S.h * s, H = Math.round(M0 + dh);
  c.style.height = H + 'px'; $('cam').style.height = H + 'px';
  const dpr = window.devicePixelRatio || 1;
  const pw = Math.round(W * dpr), ph = Math.round(H * dpr);
  const lay = { W, H, M0, s, dw, dh, dpr };
  if (c.width !== pw || c.height !== ph || !S.trail || S.trail.width !== pw || S.trail.height !== ph) {
    c.width = pw; c.height = ph;
    S.trail = document.createElement('canvas'); S.trail.width = pw; S.trail.height = ph;
    $('cam').width = pw; $('cam').height = ph;
    S.ctrail = document.createElement('canvas'); S.ctrail.width = pw; S.ctrail.height = ph;
    S.prev = S.cprev = null;
  }
  S.layout = lay;
  return lay;
}

function resetTrail() {
  if (S.trail) S.trail.getContext('2d').clearRect(0, 0, S.trail.width, S.trail.height);
  if (S.ctrail) S.ctrail.getContext('2d').clearRect(0, 0, S.ctrail.width, S.ctrail.height);
  S.prev = S.cprev = null; S.t = 0;
  requestDraw();
}

// pen position in the picture, from the two independent chains
function penAt(t, Mx, My) { return [chainTip(S.ax, t, Mx), chainTip(S.ay, t, My === undefined ? Mx : My)]; }

function stepTrail(lay, pen, t, cam) {
  const g = (cam ? S.ctrail : S.trail).getContext('2d'), prev = cam ? S.cprev : S.prev;
  g.setTransform(lay.dpr, 0, 0, lay.dpr, 0, 0);
  const sx = lay.M0 + pen[0] * lay.s, sy = lay.M0 + pen[1] * lay.s;
  if (prev) {
    const bridge = S.sig.b[Math.floor(t * S.N + 0.5) % S.N] === 1;
    if (!(bridge && $('c-penup').checked)) {
      g.beginPath(); g.moveTo(prev[0], prev[1]); g.lineTo(sx, sy);
      g.strokeStyle = S.theme.pen; g.globalAlpha = bridge ? 0.45 : 1; g.lineWidth = bridge ? 1 : 2; g.lineCap = 'round'; g.stroke();
      g.globalAlpha = 1;
    }
  }
  if (cam) S.cprev = [sx, sy]; else S.prev = [sx, sy];
}

const posX = [], posY = [];
function drawStage() {
  const lay = S.layout || layoutStage();
  const c = $('stage');
  if (!lay) return;
  const g = c.getContext('2d'), th = S.theme;
  g.setTransform(lay.dpr, 0, 0, lay.dpr, 0, 0);
  g.clearRect(0, 0, lay.W, lay.H);
  const { M0, s, dw, dh, W, H } = lay;
  g.strokeStyle = th.line; g.lineWidth = 1; g.strokeRect(M0 + 0.5, M0 + 0.5, dw - 1, dh - 1);
  if (!S.sig) { emptyMsg(g, W, H, 'nothing to draw yet: no closed path'); return; }
  if (!S.layout.ready) S.layout.ready = true;

  if ($('c-ghost').checked) {         // the target: the resampled path, bridges fainter
    g.save(); g.translate(M0, M0); g.scale(s, s);
    g.lineWidth = 1 / s; g.strokeStyle = th.ref; g.globalAlpha = 0.45;
    const sg = S.sig, N = sg.N;
    g.beginPath();
    for (let n = 0; n <= N; n++) {
      const i = n % N, prevB = sg.b[(n - 1 + N) % N];
      if (n === 0 || prevB) g.moveTo(sg.x[i], sg.y[i]); else g.lineTo(sg.x[i], sg.y[i]);
    }
    g.stroke(); g.restore();
  }
  g.drawImage(S.trail, 0, 0, lay.W, lay.H);

  const a = S.ax, b = S.ay, Mx = S.Mx, My = S.My;
  const need = 2 * (Math.max(Mx, My) + 1);
  if (posX.length < need) { posX.length = need; posY.length = need; }
  const px = chainTip(a, S.t, Mx, posX), py = chainTip(b, S.t, My, posY);
  // screen mapping: X chain plane (A, B) -> (M0 + A s, M0/2 + B s);  Y chain rotated a quarter turn -> (M0/2 - B s, M0 + A s)
  const X = (A, B) => [M0 + A * s, M0 / 2 + B * s];
  const Y = (A, B) => [M0 / 2 - B * s, M0 + A * s];
  const pen = [M0 + px * s, M0 + py * s];

  const chain = (pos, M, map, clip, colour) => {
    g.save();
    g.beginPath(); g.rect(clip[0], clip[1], clip[2], clip[3]); g.clip();
    if ($('c-circles').checked) {
      g.strokeStyle = colour; g.globalAlpha = 0.28; g.lineWidth = 1; g.beginPath();
      for (let i = 0; i < M; i++) {
        const r = (pos === posX ? a : b).R[i] * s;
        if (r < 0.6) break;                    // sorted by radius, so the rest are smaller
        const p = map(pos[2 * i], pos[2 * i + 1]);
        g.moveTo(p[0] + r, p[1]); g.arc(p[0], p[1], r, 0, TAU);
      }
      g.stroke();
      g.globalAlpha = 0.9; g.lineWidth = 1.2; g.beginPath();
      for (let i = 0; i <= M; i++) { const p = map(pos[2 * i], pos[2 * i + 1]); i ? g.lineTo(p[0], p[1]) : g.moveTo(p[0], p[1]); }
      g.stroke();
    }
    g.restore();
  };
  chain(posX, Mx, X, [M0, 0, W - M0, M0], th.x);
  chain(posY, My, Y, [0, M0, M0, dh], th.y);

  // tips and lines from each tip to the pen
  const tx = X(px, posX[2 * Mx + 1]), ty = Y(py, posY[2 * My + 1]);
  g.save(); g.setLineDash([4, 3]); g.lineWidth = 1;
  g.strokeStyle = th.x; g.beginPath(); g.moveTo(tx[0], tx[1]); g.lineTo(pen[0], pen[1]); g.stroke();
  g.strokeStyle = th.y; g.beginPath(); g.moveTo(ty[0], ty[1]); g.lineTo(pen[0], pen[1]); g.stroke();
  g.restore();
  g.fillStyle = th.x; g.beginPath(); g.arc(tx[0], tx[1], 3, 0, TAU); g.fill();
  g.fillStyle = th.y; g.beginPath(); g.arc(ty[0], ty[1], 3, 0, TAU); g.fill();
  g.fillStyle = th.pen; g.strokeStyle = th.ink; g.lineWidth = 1.2;
  g.beginPath(); g.arc(pen[0], pen[1], 4.5, 0, TAU); g.fill(); g.stroke();
  g.fillStyle = th.muted; g.font = '11px system-ui, sans-serif';
  g.fillText('X chain → pen x', M0 + 6, 14); g.fillText('Y chain', 6, M0 + 14); g.fillText('→ pen y', 6, M0 + 27);

  const cx = (S.t * 100) + '%';
  $('cur-x').style.left = 'calc(40px + (100% - 46px) * ' + S.t + ')';
  $('cur-y').style.left = 'calc(40px + (100% - 46px) * ' + S.t + ')';
}

// ---------- 7. cams ----------
// One cam per axis, one turn per lap. Follower displacement s(θ) = x(θ) - min x, from the band-limited series (the
// same circles as above). Pitch curve r(θ) = Rb + s(θ). A knife-edge follower rides the pitch curve itself; a roller of
// radius ρ rides the inner offset of it, and the offset has a cusp (undercut) where the pitch curve is convex with
// radius of curvature below ρ. Pressure angle (radial translating follower): tan φ = r'(θ) / r(θ).
function camProfile(part, rho, roller) {
  const N = part.length, da = TAU / N;
  let mn = Infinity, mx = -Infinity;
  for (let i = 0; i < N; i++) { mn = Math.min(mn, part[i]); mx = Math.max(mx, part[i]); }
  const smax = mx - mn, Rb = 0.6 * smax + 10;
  const sd = new Float64Array(N), r = new Float64Array(N);
  for (let i = 0; i < N; i++) { sd[i] = part[i] - mn; r[i] = Rb + sd[i]; }
  const outline = new Path2D(), pitch = new Path2D(), bad = new Path2D();
  const px = new Float64Array(N), py = new Float64Array(N), uc = new Uint8Array(N);
  let maxPress = 0, nUc = 0;
  for (let i = 0; i < N; i++) {
    const a = i * da, c = Math.cos(a), sn = Math.sin(a), r0 = r[i];
    const rp = (r[(i + 1) % N] - r[(i - 1 + N) % N]) / (2 * da);
    const rpp = (r[(i + 1) % N] - 2 * r0 + r[(i - 1 + N) % N]) / (da * da);
    const Tx = rp * c - r0 * sn, Ty = rp * sn + r0 * c, T2 = Tx * Tx + Ty * Ty, T = Math.sqrt(T2);
    const kappa = (r0 * r0 + 2 * rp * rp - r0 * rpp) / (T2 * T);
    maxPress = Math.max(maxPress, Math.atan(Math.abs(rp) / r0));
    let x = r0 * c, y = r0 * sn;
    if (i) pitch.lineTo(x, y); else pitch.moveTo(x, y);
    if (roller) { x -= rho * Ty / T; y += rho * Tx / T; }        // offset inward along the outward normal (Ty, -Tx)/|T|
    px[i] = x; py[i] = y;
    uc[i] = roller && kappa > 0 && 1 - rho * kappa < 0 ? 1 : 0;
    nUc += uc[i];
    if (i) outline.lineTo(x, y); else outline.moveTo(x, y);
  }
  outline.closePath(); pitch.closePath();
  for (let i = 0; i < N; i++) if (uc[i]) {
    const j = (i + 1) % N;
    bad.moveTo(px[i], py[i]); bad.lineTo(px[j], py[j]);
  }
  return { N, min: mn, s: sd, r, Rb, smax, outline, pitch, bad, px, py, uc, nUc, maxPress: maxPress * 180 / Math.PI };
}

function buildCams() {
  if (!S.partX) { S.cam = null; return; }
  const roller = S.follower === 'roller';
  S.cam = { x: camProfile(S.partX, S.rho, roller), y: camProfile(S.partY, S.rho, roller), roller };
  const cx = S.cam.x, cy = S.cam.y;
  const bad = [];
  if (cx.nUc) bad.push('X');
  if (cy.nUc) bad.push('Y');
  const press = (n, c) => n + ' cam: max pressure angle ' + fmt(c.maxPress, 0) + '°';
  $('st-cam').innerHTML = press('X', cx) + ', ' + press('Y', cy) + (cx.maxPress > 30 || cy.maxPress > 30 ? ' (above about 30° a real follower jams)' : '') + '. ' +
    (roller ? (bad.length ? '<span class="bad">Undercut on the ' + bad.join(' and ') + ' cam (red): it cannot be made with a roller of radius ' + S.rho + ' px.</span>'
                          : 'No undercut with a roller of radius ' + S.rho + ' px.')
            : 'Knife-edge follower: exact for any profile.');
}

// pen position read back from the two cam profiles (linear interpolation between the N profile samples)
function camAt(cm, t) {
  const u = t * cm.N, i = Math.floor(u), f = u - i;
  return cm.min + (1 - f) * cm.s[i % cm.N] + f * cm.s[(i + 1) % cm.N];
}
function camPen(t) { return [camAt(S.cam.x, t), camAt(S.cam.y, t)]; }

function drawCam() {
  const lay = S.layout, c = $('cam');
  if (!lay) return;
  const g = c.getContext('2d'), th = S.theme, cm = S.cam;
  g.setTransform(lay.dpr, 0, 0, lay.dpr, 0, 0);
  g.clearRect(0, 0, lay.W, lay.H);
  const { M0, s: sc, dw, dh, W, H } = lay;
  g.strokeStyle = th.line; g.lineWidth = 1; g.strokeRect(M0 + 0.5, M0 + 0.5, dw - 1, dh - 1);
  if (!cm) { emptyMsg(g, W, H, 'nothing to draw yet: no closed path'); return; }
  g.drawImage(S.ctrail, 0, 0, W, H);
  const th2 = TAU * S.t, rodEnd = 0.47 * M0;
  const [px, py] = camPen(S.t);
  const pen = [M0 + px * sc, M0 + py * sc];
  const one = (c1, cx, cy, psi, colour) => {
    const cs = 0.42 * M0 / (c1.Rb + c1.smax);
    const rc = c1.Rb + camAt({ N: c1.N, min: 0, s: c1.s }, S.t);        // pitch radius under the follower
    g.save();
    g.translate(cx, cy); g.rotate(psi - th2); g.scale(cs, cs);
    g.lineJoin = 'round'; g.lineWidth = 1.6 / cs; g.strokeStyle = colour;
    g.globalAlpha = 0.1; g.fillStyle = colour; g.fill(c1.outline); g.globalAlpha = 1;
    g.stroke(c1.outline);
    if (cm.roller) {
      g.globalAlpha = 0.5; g.lineWidth = 1 / cs; g.setLineDash([4 / cs, 3 / cs]); g.stroke(c1.pitch); g.setLineDash([]); g.globalAlpha = 1;
      g.lineWidth = 2.6 / cs; g.strokeStyle = th.bad; g.stroke(c1.bad);
    }
    g.lineWidth = 1 / cs; g.strokeStyle = th.muted;                      // a mark on the cam so its turning shows
    g.beginPath(); g.moveTo(0, 0); g.lineTo(c1.Rb * 0.35, 0); g.stroke();
    g.restore();
    g.fillStyle = th.ink; g.beginPath(); g.arc(cx, cy, 2.5, 0, TAU); g.fill();
    const ux = Math.cos(psi), uy = Math.sin(psi);
    const tx = cx + ux * cs * rc, ty = cy + uy * cs * rc;                // follower tip (knife) or roller centre
    g.strokeStyle = th.ink; g.lineWidth = 2; g.beginPath(); g.moveTo(tx, ty); g.lineTo(cx + ux * rodEnd, cy + uy * rodEnd); g.stroke();
    if (cm.roller) { g.beginPath(); g.arc(tx, ty, S.rho * cs, 0, TAU); g.fillStyle = th.panel; g.fill(); g.stroke(); }
    else { g.fillStyle = th.ink; g.beginPath(); g.arc(tx, ty, 3, 0, TAU); g.fill(); }
    g.fillStyle = th.ink; g.fillRect(cx + ux * rodEnd - 3, cy + uy * rodEnd - 3, 6, 6);
    return [cx + ux * rodEnd, cy + uy * rodEnd, cs];
  };
  // X cam above the picture, follower sliding horizontally; Y cam left of it, follower sliding vertically
  const ex = one(cm.x, M0 + dw / 2, M0 / 2, 0, th.x);
  const ey = one(cm.y, M0 / 2, M0 + dh / 2, Math.PI / 2, th.y);
  // lever from follower to slide bar (the cams are drawn smaller than the picture, so the motion is magnified)
  g.save(); g.setLineDash([4, 3]); g.lineWidth = 1;
  g.strokeStyle = th.x; g.beginPath(); g.moveTo(ex[0], ex[1]); g.lineTo(pen[0], ex[1]); g.stroke();
  g.strokeStyle = th.y; g.beginPath(); g.moveTo(ey[0], ey[1]); g.lineTo(ey[0], pen[1]); g.stroke();
  g.restore();
  // cross-slide: the X bar is vertical and slides sideways, the Y bar is horizontal and slides up and down
  g.lineWidth = 3; g.lineCap = 'round';
  g.strokeStyle = th.x; g.beginPath(); g.moveTo(pen[0], ex[1]); g.lineTo(pen[0], pen[1]); g.stroke();
  g.strokeStyle = th.y; g.beginPath(); g.moveTo(ey[0], pen[1]); g.lineTo(pen[0], pen[1]); g.stroke();
  g.lineCap = 'butt';
  g.fillStyle = th.pen; g.strokeStyle = th.ink; g.lineWidth = 1.2;
  g.beginPath(); g.arc(pen[0], pen[1], 4.5, 0, TAU); g.fill(); g.stroke();
  g.fillStyle = th.muted; g.font = '11px system-ui, sans-serif';
  g.fillText('X cam → slide bar → pen x', M0 + 6, 14); g.fillText('Y cam', 6, M0 + 14); g.fillText('→ pen y', 6, M0 + 27);
  S.camScale = [ex[2], ey[2]];
}

let drawQueued = false;
function requestDraw() { if (!drawQueued && !S.playing) { drawQueued = true; requestAnimationFrame(() => { drawQueued = false; drawStage(); drawCam(); }); } }

// ---------- animation loop ----------
let last = 0, frames = 0, fpsT = 0, drawAcc = 0;
function frame(now) {
  requestAnimationFrame(frame);
  const dt = Math.min(0.1, (now - last) / 1000); last = now;
  if (S.dirty) flush();
  if (!S.sig) { return; }
  const t0 = performance.now();
  if (S.playing && S.layout) {
    const dT = dt * S.speed / BASE_PERIOD;                  // fraction of a lap this frame
    const sub = Math.max(1, Math.ceil(dT * S.N / 1.5));     // at most ~1.5 samples of the path per step
    for (let i = 1; i <= sub; i++) {
      let t = S.t + dT / sub; if (t >= 1) t -= 1; if (t < 0) t += 1;
      S.t = t;
      stepTrail(S.layout, penAt(t, S.Mx, S.My), t);
      if (S.cam) stepTrail(S.layout, camPen(t), t, true);
    }
  }
  if (S.playing || drawQueued) { drawStage(); drawCam(); }
  drawAcc += performance.now() - t0;
  frames++;
  if (now - fpsT > 1000) {
    S.fps = frames * 1000 / (now - fpsT); S.drawMs = drawAcc / frames; frames = 0; drawAcc = 0; fpsT = now;
    $('st-stage').textContent = 'N = ' + S.N + ', X chain ' + S.Mx + ' circles, Y chain ' + S.My + ' circles (' + (S.Mx + S.My) + ' in all), lap = ' + fmt(BASE_PERIOD / S.speed, 1) + ' s. ' +
      fmt(S.fps, 0) + ' fps, ' + fmt(S.drawMs, 2) + ' ms of work per frame.';
  }
}

// ---------- scheduling ----------
function invalidate(level) { S.dirty = Math.max(S.dirty, level); if (!S.playing) { requestAnimationFrame(() => { if (S.dirty) flush(); }); } }
function flush() {
  const lvl = S.dirty; S.dirty = 0;
  if (!S.gray) return;
  readTheme();
  if (lvl >= 3) { runVectorize(); drawVectorViews(); }
  if (lvl >= 2) { runConnect(); drawConnectView(); }
  if (lvl >= 1) { runSignals(); drawSignals(); drawSpectra(); layoutStage(); resetTrail(); }
}

// ---------- samples and input ----------
function makeSample() {
  const c = document.createElement('canvas'); c.width = 640; c.height = 480;
  const g = c.getContext('2d', { willReadFrequently: true });     // software canvas: reads back reliably
  g.fillStyle = '#fff'; g.fillRect(0, 0, 640, 480);
  g.strokeStyle = '#000'; g.lineWidth = 7; g.lineJoin = 'round'; g.lineCap = 'round';
  const poly = pts => { g.beginPath(); pts.forEach((p, i) => i ? g.lineTo(p[0], p[1]) : g.moveTo(p[0], p[1])); g.stroke(); };
  poly([[90, 240], [90, 380], [260, 380], [260, 240]]);                  // house
  poly([[70, 245], [175, 150], [280, 245]]);
  poly([[150, 380], [150, 305], [200, 305], [200, 380]]);               // door
  g.strokeRect(105, 262, 32, 32);                                       // window
  g.beginPath(); g.arc(470, 125, 48, 0, TAU); g.stroke();               // sun
  for (let i = 0; i < 10; i++) {
    const a = TAU * i / 10;
    poly([[470 + 68 * Math.cos(a), 125 + 68 * Math.sin(a)], [470 + 92 * Math.cos(a), 125 + 92 * Math.sin(a)]]);
  }
  const wave = []; for (let x = 330; x <= 600; x += 6) wave.push([x, 330 + 36 * Math.sin((x - 330) / 270 * TAU * 2)]);
  poly(wave);                                                           // a sine wave
  return c;
}

function useSample() {
  setMode('lines', true);
  for (let attempt = 0; attempt < 5; attempt++) {
    const c = makeSample();
    setSource(c, c.width, c.height, 'built-in sample');
    let mn = 255; for (let i = 0; i < S.gray.length; i++) if (S.gray[i] < mn) mn = S.gray[i];
    if (mn < 128) break;                    // the sample has black lines: a blank read-back means the canvas failed, draw it again
  }
}

async function loadFile(file) {
  if (!file || !/^image\//.test(file.type || 'image/')) return;
  try {
    const bmp = await createImageBitmap(file);
    setMode('edges', true);
    setSource(bmp, bmp.width, bmp.height, file.name);
  } catch (e) {
    const rd = new FileReader();
    rd.onload = () => {
      const img = new Image();
      img.onload = () => { setMode('edges', true); setSource(img, img.naturalWidth, img.naturalHeight, file.name); };
      img.onerror = () => { $('src-name').textContent = 'could not read ' + file.name; };
      img.src = rd.result;
    };
    rd.readAsDataURL(file);
  }
}

function setMode(mode, resetThr) {
  S.mode = mode;
  document.querySelectorAll('#mode button').forEach(b => b.setAttribute('aria-pressed', String(b.dataset.mode === mode)));
  if (resetThr) { $('p-thr').value = mode === 'edges' ? 25 : 50; showOutputs(); }
}

function showOutputs() {
  $('o-blur').textContent = (+$('p-blur').value).toFixed(1);
  $('o-thr').textContent = $('p-thr').value + '%';
  $('o-min').textContent = $('p-min').value + ' px';
  $('o-tol').textContent = (+$('p-tol').value).toFixed(1) + ' px';
  $('o-n').textContent = S.N;
  S.speed = Math.pow(2, +$('p-speed').value / 50);
  $('o-speed').textContent = '×' + S.speed.toFixed(2);
}

// ---------- wiring ----------
function init() {
  readTheme();
  matchMedia('(prefers-color-scheme: dark)').addEventListener('change', () => { readTheme(); S.dirty = Math.max(S.dirty, 0); redrawAll(); });

  $('file').addEventListener('change', e => { loadFile(e.target.files[0]); e.target.value = ''; });
  $('btn-sample').addEventListener('click', useSample);
  const drop = $('drop');
  ['dragenter', 'dragover'].forEach(ev => drop.addEventListener(ev, e => { e.preventDefault(); drop.classList.add('over'); }));
  ['dragleave', 'drop'].forEach(ev => drop.addEventListener(ev, e => { e.preventDefault(); drop.classList.remove('over'); }));
  drop.addEventListener('drop', e => loadFile(e.dataTransfer.files[0]));
  // dropping anywhere on the page must not navigate away
  window.addEventListener('dragover', e => e.preventDefault());
  window.addEventListener('drop', e => { e.preventDefault(); if (e.dataTransfer.files[0]) loadFile(e.dataTransfer.files[0]); });

  document.querySelectorAll('#mode button').forEach(b => b.addEventListener('click', () => { setMode(b.dataset.mode, true); invalidate(3); }));
  for (const id of ['p-blur', 'p-thr', 'p-min', 'p-tol']) $(id).addEventListener('input', () => { showOutputs(); invalidate(3); });
  $('p-n').addEventListener('input', () => { S.N = Math.pow(2, +$('p-n').value); showOutputs(); invalidate(1); });
  $('c-bridges').addEventListener('change', () => { readTheme(); drawConnectView(); });

  const togglePlay = () => {
    S.playing = !S.playing;
    for (const id of ['b-play', 'b-play2']) { $(id).textContent = S.playing ? 'Pause' : 'Play'; $(id).setAttribute('aria-pressed', String(S.playing)); }
    if (S.playing) { last = performance.now(); }
  };
  $('b-play').addEventListener('click', togglePlay); $('b-play2').addEventListener('click', togglePlay);
  $('b-reset').addEventListener('click', resetTrail); $('b-reset2').addEventListener('click', resetTrail);
  document.querySelectorAll('#follower button').forEach(b => b.addEventListener('click', () => {
    S.follower = b.dataset.f;
    document.querySelectorAll('#follower button').forEach(x => x.setAttribute('aria-pressed', String(x === b)));
    buildCams(); requestDraw();
  }));
  $('p-rho').addEventListener('input', () => { S.rho = +$('p-rho').value; $('o-rho').textContent = S.rho + ' px'; buildCams(); requestDraw(); });
  // X and Y circle counts are independent; "link" makes a change to one move the other
  const setFrac = (axis, f) => {
    f = Math.max(0, Math.min(1, f));
    const ids = axis === 'x' ? ['fx', 'p-mx'] : ['fy', 'p-my'];
    S[ids[0]] = f; $(ids[1]).value = Math.round(f * 1000);
  };
  for (const axis of ['x', 'y']) {
    const other = axis === 'x' ? 'y' : 'x';
    const apply = f => { setFrac(axis, f); if ($('c-link').checked) setFrac(other, f); updateM(); };
    $('p-m' + axis).addEventListener('input', () => apply(+$('p-m' + axis).value / 1000));
    $('mp-' + axis).addEventListener('click', e => {
      const v = e.target.dataset && e.target.dataset.m;
      if (!v || !S.ax) return;
      const m = v === 'all' ? S.ax.H : Math.min(+v, S.ax.H);
      apply(Math.log(m) / Math.log(S.ax.H));
    });
  }
  $('c-link').addEventListener('change', () => { if ($('c-link').checked) { setFrac('y', S.fx); updateM(); } });
  $('p-speed').addEventListener('input', showOutputs);
  for (const id of ['c-circles', 'c-ghost', 'c-penup']) $(id).addEventListener('change', () => { if (id === 'c-penup') resetTrail(); requestDraw(); });
  $('c-penup').addEventListener('change', () => { drawStage(); });

  let rz = 0;
  window.addEventListener('resize', () => { cancelAnimationFrame(rz); rz = requestAnimationFrame(() => { layoutViews(); layoutStage(); resetTrail(); redrawAll(); }); });
  showOutputs();
  useSample();
  requestAnimationFrame(frame);
}

function redrawAll() {
  if (!S.gray) return;
  readTheme();
  if (S.stats.vec) { drawVectorViews(); drawConnectView(); drawSignals(); drawSpectra(); }
  drawStage(); drawCam();
}

window.epicycles = { S, camPen, chainTip, penAt, analyse, fft, resample, connect, partialSignal, loadFile, flush, drawStage };
init();
})();
