/* Genoarch engine (JavaScript port of core/genoarch_core/engine.py).
 * Runs fully in the browser (or Node) on the compact index. Every reduction that
 * feeds a ranking runs in the same order as the Python engine, so disease scores
 * and null statistics are bit-identical; see scripts/parity_check.mjs.
 * Research use only: not a medical device, not for diagnosis. */
(function (root, factory) {
  const api = factory();
  if (typeof module === "object" && module.exports) module.exports = api;
  root.GenoarchEngine = api;
})(typeof globalThis !== "undefined" ? globalThis : this, function () {
  "use strict";
  const ENGINE_VERSION = "0.1.0";
  const MODES = ["asymmetric", "blended", "symmetric"];
  const UNIVERSES = { omim: ["OMIM", "DECIPHER"], all: ["OMIM", "DECIPHER", "ORPHA"] };
  const RESEARCH_USE_ONLY = "Research use only. Genoarch is not a medical device and must not be used for diagnosis, " +
    "treatment or any clinical decision.";
  const DEFAULTS = { mode: "blended", universe: "omim", T: 1.0, beta: 0.0, w: 0.1, alpha: 0.05, null_samples: 100,
    ic_bins: 10, eps: 0.05, q_max: 0.99, eig_top_k: 50, eig_min_ic: 1.0, hyp_top_diseases: 10, hyp_min_ic: 3.0,
    hyp_min_terms: 2 };
  const IC_SCALE = 1000000, FREQ_SCALE = 10000, UNKNOWN_FREQUENCY = 0.5;
  const PHENOTYPIC_ABNORMALITY = "HP:0000118";

  // ---------------------------------------------------------------- PRNG
  function fnv1a32(text) {
    let h = 0x811c9dc5 >>> 0;
    const bytes = new TextEncoder().encode(text);
    for (let i = 0; i < bytes.length; i++) h = Math.imul(h ^ bytes[i], 0x01000193) >>> 0;
    return h >>> 0;
  }
  function mulberry32(seed) {
    let a = seed | 0;
    return function () {
      a = (a + 0x6d2b79f5) | 0;
      let t = Math.imul(a ^ (a >>> 15), 1 | a);
      t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
      return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
  }

  // ------------------------------------------------------------- helpers
  const pad7 = (n) => "HP:" + String(n).padStart(7, "0");
  function undelta(arr) { const out = new Int32Array(arr.length); let acc = 0; for (let i = 0; i < arr.length; i++) { acc += arr[i]; out[i] = acc; } return out; }
  function bsearch(arr, v) { let lo = 0, hi = arr.length - 1; while (lo <= hi) { const m = (lo + hi) >> 1; const x = arr[m]; if (x === v) return true; if (x < v) lo = m + 1; else hi = m - 1; } return false; }
  function seqSum(arr) { let s = 0; for (let i = 0; i < arr.length; i++) s += arr[i]; return s; }
  function fixed4(x) { return Number(x).toFixed(4); }
  function h2(q) { return q > 0 && q < 1 ? -(q * Math.log2(q) + (1.0 - q) * Math.log2(1.0 - q)) : 0; }
  function normalizeHpoId(raw) {
    if (raw === null || raw === undefined) return null;
    let s = String(raw).trim();
    if (!s) return null;
    if (s.toLowerCase().startsWith("http")) { const m = s.match(/\/obo\/([A-Za-z]+)_(\d+)$/); if (m) s = m[1] + ":" + m[2]; }
    const m = s.toUpperCase().match(/^(?:HP)[:_]?(\d{1,7})$/);
    return m ? "HP:" + m[1].padStart(7, "0") : null;
  }
  class QueryError extends Error {}

  // ------------------------------------------------------- knowledge base
  function buildKB(doc) {
    if (!doc || doc.format !== "genoarch-index") throw new Error("not a genoarch-index document");
    const T = doc.terms, nT = T.hp.length;
    const kb = { meta: doc.meta || {}, nT };
    kb.termIds = T.hp.map(pad7);
    kb.termIndex = new Map(kb.termIds.map((t, i) => [t, i]));
    kb.names = T.name; kb.synonyms = T.syn;
    kb.parents = T.par.map((p) => Int32Array.from(p).sort());
    kb.ic = new Float64Array(nT); kb.cnt = new Int32Array(nT);
    for (let i = 0; i < nT; i++) { kb.ic[i] = T.ic6[i] / IC_SCALE; kb.cnt[i] = T.cnt[i]; }
    kb.children = Array.from({ length: nT }, () => []);
    for (let i = 0; i < nT; i++) for (const p of kb.parents[i]) kb.children[p].push(i);
    kb.children = kb.children.map((c) => Int32Array.from(c).sort());
    if (T.anc) kb.ancestors = T.anc.map(undelta);
    else {
      const memo = new Array(nT);
      const anc = (i) => { if (memo[i]) return memo[i]; const s = new Set([i]); for (const p of kb.parents[i]) for (const a of anc(p)) s.add(a); memo[i] = Int32Array.from(s).sort(); return memo[i]; };
      for (let i = 0; i < nT; i++) anc(i);
      kb.ancestors = memo;
    }
    kb.ancOff = new Int32Array(nT + 1);
    for (let i = 0; i < nT; i++) kb.ancOff[i + 1] = kb.ancOff[i] + kb.ancestors[i].length;
    kb.ancFlat = new Int32Array(kb.ancOff[nT]);
    for (let i = 0; i < nT; i++) kb.ancFlat.set(kb.ancestors[i], kb.ancOff[i]);
    kb.altMap = new Map(doc.alt.map(([a, i]) => [pad7(a), i]));
    kb.otherTerms = new Map(doc.other.map((n, k) => [pad7(n), (doc.other_names || [])[k] || ""]));
    kb.root = kb.termIndex.get(PHENOTYPIC_ABNORMALITY);

    const D = doc.diseases, nD = D.id.length;
    kb.nD = nD; kb.diseaseIds = D.id; kb.diseaseNames = D.name;
    kb.diseaseIndex = new Map(D.id.map((d, i) => [d, i]));
    kb.diseaseSource = D.id.map((d) => d.split(":")[0]);
    kb.annTerms = new Array(nD); kb.annFreqs = new Array(nD);
    kb.annLen = new Int32Array(nD); kb.annOff = new Int32Array(nD + 1);
    for (let d = 0; d < nD; d++) {
      const t = undelta(D.ann[d]), f = new Float64Array(t.length);
      for (let j = 0; j < t.length; j++) f[j] = D.f4[d][j] < 0 ? NaN : D.f4[d][j] / FREQ_SCALE;
      kb.annTerms[d] = t; kb.annFreqs[d] = f; kb.annLen[d] = t.length; kb.annOff[d + 1] = kb.annOff[d] + t.length;
    }
    const nA = kb.annOff[nD];
    kb.flatAnnDisease = new Int32Array(nA); kb.flatAnnTerm = new Int32Array(nA); kb.flatAnnFreqEff = new Float64Array(nA);
    for (let d = 0; d < nD; d++) {
      const o = kb.annOff[d];
      for (let j = 0; j < kb.annLen[d]; j++) {
        kb.flatAnnDisease[o + j] = d; kb.flatAnnTerm[o + j] = kb.annTerms[d][j];
        const f = kb.annFreqs[d][j]; kb.flatAnnFreqEff[o + j] = Number.isNaN(f) ? UNKNOWN_FREQUENCY : f;
      }
    }
    // propagated annotation sets and postings (diseases per term, ascending)
    kb.prop = new Array(nD);
    const mark = new Int32Array(nT).fill(-1);
    const postCount = new Int32Array(nT);
    for (let d = 0; d < nD; d++) {
      const acc = [];
      for (const t of kb.annTerms[d]) for (const a of kb.ancestors[t]) if (mark[a] !== d) { mark[a] = d; acc.push(a); }
      const arr = Int32Array.from(acc).sort(); kb.prop[d] = arr;
      for (const a of arr) postCount[a]++;
    }
    kb.postOff = new Int32Array(nT + 1);
    for (let i = 0; i < nT; i++) kb.postOff[i + 1] = kb.postOff[i] + postCount[i];
    kb.postFlat = new Int32Array(kb.postOff[nT]);
    const fill = kb.postOff.slice(0, nT);
    for (let d = 0; d < nD; d++) for (const a of kb.prop[d]) kb.postFlat[fill[a]++] = d;

    const G = doc.genes, nG = G.sym.length;
    kb.nG = nG; kb.geneSymbols = G.sym; kb.geneNcbi = G.ncbi;
    kb.geneDiseases = G.dis.map(undelta); kb.geneAssoc = G.assoc;
    kb.diseaseGenes = Array.from({ length: nD }, () => []);
    for (let g = 0; g < nG; g++) for (const d of kb.geneDiseases[g]) kb.diseaseGenes[d].push(g);
    kb.kGenes = new Int32Array(nD);
    for (let d = 0; d < nD; d++) kb.kGenes[d] = kb.diseaseGenes[d].length;
    kb._ascCache = new Map(); kb._descCache = new Map(); kb._descendants = new Map();
    return kb;
  }

  // -------------------------------------------------------------- engine
  class Engine {
    constructor(doc, calibration) {
      this.kb = buildKB(doc);
      this.calibration = calibration || {};
      this.version = ENGINE_VERSION;
      this._universes = new Map();
      this._search = null;
    }
    termLabel(i) { return { id: this.kb.termIds[i], name: this.kb.names[i] }; }
    diseaseLabel(d) { return { id: this.kb.diseaseIds[d], name: this.kb.diseaseNames[d] }; }
    isAncestor(a, t) { return bsearch(this.kb.ancestors[t], a); }
    ancByIc(i) {
      const kb = this.kb; let hit = kb._ascCache.get(i); if (hit) return hit;
      hit = Array.from(kb.ancestors[i]).sort((x, y) => (kb.ic[x] - kb.ic[y]) || (x - y));
      kb._ascCache.set(i, hit); return hit;
    }
    ancByIcDesc(i) {
      const kb = this.kb; let hit = kb._descCache.get(i); if (hit) return hit;
      hit = Array.from(kb.ancestors[i]).sort((x, y) => (kb.ic[y] - kb.ic[x]) || (x - y));
      kb._descCache.set(i, hit); return hit;
    }
    descendants(i) {
      const kb = this.kb; let hit = kb._descendants.get(i); if (hit) return hit;
      const seen = new Set([i]), stack = [i];
      while (stack.length) { const c = stack.pop(); for (const k of kb.children[c]) if (!seen.has(k)) { seen.add(k); stack.push(k); } }
      hit = Int32Array.from(seen).sort(); if (kb._descendants.size < 4096) kb._descendants.set(i, hit); return hit;
    }
    resolve(raw) {
      const kb = this.kb, tid = normalizeHpoId(raw);
      if (tid === null) return [null, "invalid"];
      if (kb.termIndex.has(tid)) return [kb.termIndex.get(tid), "ok"];
      if (kb.altMap.has(tid)) return [kb.altMap.get(tid), "alt"];
      if (kb.otherTerms.has(tid)) return [null, "not_phenotype"];
      return [null, "unknown"];
    }

    settings(o) {
      o = o || {};
      const cal = this.calibration;
      const mode = o.mode || cal.default_mode || DEFAULTS.mode;
      const universe = o.universe || cal.default_universe || DEFAULTS.universe;
      if (!MODES.includes(mode)) throw new QueryError("mode must be one of " + MODES.join(", "));
      if (!(universe in UNIVERSES)) throw new QueryError("universe must be one of omim, all");
      const fitted = ((cal.params || {})[mode] || {})[universe] || {};
      const s = Object.assign({}, DEFAULTS);
      for (const [k, v] of Object.entries(cal.constants || {})) if (k in DEFAULTS) s[k] = v;
      s.mode = mode; s.universe = universe;
      s.T = Number(o.T != null ? o.T : (fitted.T != null ? fitted.T : DEFAULTS.T));
      s.beta = Number(o.beta != null ? o.beta : (fitted.beta != null ? fitted.beta : DEFAULTS.beta));
      if (mode === "blended") s.w = Number(o.w != null ? o.w : (fitted.w != null ? fitted.w : DEFAULTS.w));
      else s.w = mode === "symmetric" ? 0.5 : 0.0;
      if (!(s.w >= 0 && s.w <= 1)) throw new QueryError("w must be in [0, 1]");
      if (!(s.beta >= 0)) throw new QueryError("beta must be >= 0");
      s.null_samples = Math.trunc(Number(o.null_samples != null ? o.null_samples : s.null_samples));
      s.alpha = Number(o.alpha != null ? o.alpha : s.alpha);
      s.calibrated = Object.keys(fitted).length > 0 && o.T == null && o.beta == null && o.w == null;
      if (!(s.T > 0)) throw new QueryError("T must be positive");
      if (!(s.null_samples >= 0 && s.null_samples <= 1000)) throw new QueryError("null_samples must be in [0, 1000]");
      return s;
    }

    universe(name, exclude) {
      const kb = this.kb;
      const excl = Array.from(new Set((exclude || []).map(Number))).sort((a, b) => a - b);
      if (!excl.length && this._universes.has(name)) return this._universes.get(name);
      const prefixes = new Set(UNIVERSES[name]);
      const mask = new Uint8Array(kb.nD);
      for (let d = 0; d < kb.nD; d++) mask[d] = prefixes.has(kb.diseaseSource[d]) ? 1 : 0;
      for (const d of excl) mask[d] = 0;
      const uidx = []; for (let d = 0; d < kb.nD; d++) if (mask[d]) uidx.push(d);
      const inPool = new Uint8Array(kb.nT);
      for (let i = 0; i < kb.flatAnnTerm.length; i++) if (mask[kb.flatAnnDisease[i]]) inPool[kb.flatAnnTerm[i]] = 1;
      const pool = []; for (let t = 0; t < kb.nT; t++) if (inPool[t]) pool.push(t);
      pool.sort((x, y) => (kb.ic[x] - kb.ic[y]) || (x - y));
      const nb = Math.max(1, Math.min(DEFAULTS.ic_bins, pool.length));
      const bins = [], edges = [];
      for (let b = 0; b < nb; b++) {
        const lo = Math.floor((b * pool.length) / nb), hi = Math.floor(((b + 1) * pool.length) / nb);
        const bin = pool.slice(lo, hi); bins.push(bin); edges.push(kb.ic[bin[bin.length - 1]]);
      }
      const gDis = [], gOff = [], gIds = [];
      for (let g = 0; g < kb.nG; g++) {
        const ds = []; for (const d of kb.geneDiseases[g]) if (mask[d]) ds.push(d);
        if (ds.length) { gOff.push(gDis.length); for (const d of ds) gDis.push(d); gIds.push(g); }
      }
      const dgDis = [], dgGene = [];
      for (const d of uidx) for (const g of kb.diseaseGenes[d]) { dgDis.push(d); dgGene.push(g); }
      const disWithGenes = Array.from(new Set(dgDis)).sort((a, b) => a - b);
      const u = { name, mask, uidx: Int32Array.from(uidx), bins, edges, gDis: Int32Array.from(gDis), gOff: Int32Array.from(gOff),
        gIds: Int32Array.from(gIds), dgDis: Int32Array.from(dgDis), dgGene: Int32Array.from(dgGene), disWithGenes, excluded: excl };
      if (!excl.length) this._universes.set(name, u);
      return u;
    }

    prepare(present, absent) {
      const kb = this.kb, warnings = [], unresolved = [];
      const resolveAll = (items, role) => {
        const out = [];
        for (const raw of items || []) {
          const [i, status] = this.resolve(String(raw));
          if (i === null) { unresolved.push({ input: String(raw), status, role }); continue; }
          if (status === "alt") warnings.push(`${raw} was mapped to its current id ${kb.termIds[i]}.`);
          out.push(i);
        }
        return Array.from(new Set(out)).sort((a, b) => a - b);
      };
      const pres = resolveAll(present, "present");
      let absn = resolveAll(absent, "absent");
      const implied = new Set();
      for (const p of pres) for (const a of kb.ancestors[p]) implied.add(a);
      for (const e of absn) if (implied.has(e)) warnings.push(`${kb.termIds[e]} was marked absent but is implied by a present term; ignored as absent.`);
      absn = absn.filter((e) => !implied.has(e));
      const presUsed = pres.filter((p) => !pres.some((q) => q !== p && this.isAncestor(p, q)));
      const absUsed = absn.filter((e) => !absn.some((f) => f !== e && this.isAncestor(f, e)));
      if (presUsed.length < pres.length) warnings.push("Present terms that are ancestors of other present terms were folded into the more specific term.");
      return { present: pres, absent: absn, present_used: presUsed, absent_used: absUsed, unresolved, warnings };
    }

    // ------------------------------------------------------------ scoring
    bmRow(p, row) {
      const kb = this.kb; row.fill(0);
      for (const a of this.ancByIc(p)) {
        const ica = kb.ic[a];
        if (ica > 0) for (let k = kb.postOff[a]; k < kb.postOff[a + 1]; k++) row[kb.postFlat[k]] = ica;
      }
      return row;
    }
    bmMatrix(terms) { return terms.map((p) => this.bmRow(p, new Float64Array(this.kb.nD))); }
    sDP(terms) {
      const kb = this.kb, inQ = new Uint8Array(kb.nT);
      for (const p of terms) for (const a of kb.ancestors[p]) inQ[a] = 1;
      const m = new Float64Array(kb.nT);
      for (let t = 0; t < kb.nT; t++) {
        let best = 0;
        for (let k = kb.ancOff[t]; k < kb.ancOff[t + 1]; k++) { const a = kb.ancFlat[k]; const v = inQ[a] ? kb.ic[a] : 0.0; if (v > best) best = v; }
        m[t] = best;
      }
      const sums = new Float64Array(kb.nD);
      for (let i = 0; i < kb.flatAnnTerm.length; i++) sums[kb.flatAnnDisease[i]] += m[kb.flatAnnTerm[i]];
      for (let d = 0; d < kb.nD; d++) sums[d] = sums[d] / kb.annLen[d];
      return sums;
    }
    absentFrequencies(e) {
      const kb = this.kb, inDesc = new Uint8Array(kb.nT);
      for (const t of this.descendants(e)) inDesc[t] = 1;
      const fe = new Float64Array(kb.nD);
      for (let d = 0; d < kb.nD; d++) {
        let best = 0;
        for (let k = kb.annOff[d]; k < kb.annOff[d + 1]; k++) { const v = inDesc[kb.flatAnnTerm[k]] ? kb.flatAnnFreqEff[k] : 0.0; if (v > best) best = v; }
        fe[d] = best;
      }
      return fe;
    }
    penalty(absent) {
      const kb = this.kb, pen = new Float64Array(kb.nD);
      for (const e of absent) { const fe = this.absentFrequencies(e); const ice = kb.ic[e]; for (let d = 0; d < kb.nD; d++) pen[d] += ice * fe[d]; }
      return pen;
    }
    components(present, absent, needDp) {
      const nD = this.kb.nD, bm = this.bmMatrix(present), acc = new Float64Array(nD);
      for (const row of bm) for (let d = 0; d < nD; d++) acc[d] += row[d];
      const n = Math.max(present.length, 1);
      for (let d = 0; d < nD; d++) acc[d] = acc[d] / n;
      return { bm, s_pd: acc, s_dp: needDp ? this.sDP(present) : null, pen: absent.length ? this.penalty(absent) : null, n };
    }
    static combine(comp, mode, beta, w) {
      const nD = comp.s_pd.length, out = new Float64Array(nD);
      for (let d = 0; d < nD; d++) {
        let base;
        if (mode === "symmetric") base = 0.5 * (comp.s_pd[d] + comp.s_dp[d]);
        else if (mode === "blended") base = (1.0 - w) * comp.s_pd[d] + w * comp.s_dp[d];
        else base = comp.s_pd[d];
        out[d] = beta && comp.pen ? base - beta * comp.pen[d] / comp.n : base;
      }
      return out;
    }
    order(score, uidx) { return Array.from(uidx).sort((a, b) => (score[b] - score[a]) || (a - b)); }
    posterior(score, uidx, T) {
      let smax = -Infinity; for (const d of uidx) if (score[d] > smax) smax = score[d];
      const w = new Float64Array(uidx.length);
      for (let i = 0; i < uidx.length; i++) w[i] = Math.exp((score[uidx[i]] - smax) / T);
      const Z = seqSum(w), post = new Float64Array(this.kb.nD);
      for (let i = 0; i < uidx.length; i++) post[uidx[i]] = w[i] / Z;
      return post;
    }

    // --------------------------------------------------------- null model
    nullTest(q, s, u, comp, score) {
      const kb = this.kb, R = s.null_samples, uidx = u.uidx, nU = uidx.length, nD = kb.nD;
      if (R < 20) return { samples: R, p_value: null, abstain: false, alpha: s.alpha, z: null };
      const seedText = `genoarch-null|${s.mode}|${s.universe}|${fixed4(s.beta)}|${fixed4(s.w)}|` +
        q.present_used.map((p) => kb.termIds[p]).join(",") + "|" + q.absent_used.map((e) => kb.termIds[e]).join(",");
      const rng = mulberry32(fnv1a32(seedText));
      const termBins = q.present_used.map((p) => {
        let icp = -Infinity; for (const a of kb.ancestors[p]) if (kb.ic[a] > icp) icp = kb.ic[a];
        const b = u.edges.findIndex((e) => icp <= e); return b >= 0 ? b : u.edges.length - 1;
      });
      const n = comp.n, S = new Float64Array(R * nU), row = new Float64Array(nD), acc = new Float64Array(nD);
      for (let r = 0; r < R; r++) {
        const terms = termBins.map((b) => { const pool = u.bins[b]; return pool[Math.floor(rng() * pool.length)]; });
        acc.fill(0);
        for (const t of terms) { this.bmRow(t, row); for (let d = 0; d < nD; d++) acc[d] += row[d]; }
        const spd = new Float64Array(nD); for (let d = 0; d < nD; d++) spd[d] = acc[d] / n;
        const nc = { s_pd: spd, s_dp: s.mode !== "asymmetric" ? this.sDP(terms) : null, pen: comp.pen, n };
        const sc = Engine.combine(nc, s.mode, s.beta, s.w);
        for (let j = 0; j < nU; j++) S[r * nU + j] = sc[uidx[j]];
      }
      const s1 = new Float64Array(nU), s2 = new Float64Array(nU);
      for (let r = 0; r < R; r++) for (let j = 0; j < nU; j++) { const x = S[r * nU + j]; s1[j] += x; s2[j] += x * x; }
      const mean = new Float64Array(nU), vr = new Float64Array(nU);
      for (let j = 0; j < nU; j++) { mean[j] = s1[j] / R; vr[j] = Math.max(s2[j] / R - mean[j] * mean[j], 0.0); }
      const srt = Float64Array.from(vr).sort(), k = srt.length;
      const v0 = k % 2 ? srt[(k - 1) / 2] : (srt[k / 2 - 1] + srt[k / 2]) / 2;
      const zObs = new Float64Array(nU); let jBest = 0;
      for (let j = 0; j < nU; j++) { zObs[j] = (score[uidx[j]] - mean[j]) / Math.sqrt(vr[j] + v0 + 1e-12); if (zObs[j] > zObs[jBest]) jBest = j; }
      const zMax = zObs[jBest], zNullMax = new Float64Array(R), rawMax = new Float64Array(R);
      for (let r = 0; r < R; r++) {
        let zm = -Infinity, rm = -Infinity;
        for (let j = 0; j < nU; j++) {
          const x = S[r * nU + j];
          const m1 = (s1[j] - x) / (R - 1);
          const v1 = Math.max((s2[j] - x * x) / (R - 1) - m1 * m1, 0.0);
          const z = (x - m1) / Math.sqrt(v1 + v0 + 1e-12);
          if (z > zm) zm = z; if (x > rm) rm = x;
        }
        zNullMax[r] = zm; rawMax[r] = rm;
      }
      let exceed = 0, exceedRaw = 0, obsRaw = -Infinity;
      for (let j = 0; j < nU; j++) if (score[uidx[j]] > obsRaw) obsRaw = score[uidx[j]];
      for (let r = 0; r < R; r++) { if (zNullMax[r] >= zMax) exceed++; if (rawMax[r] >= obsRaw) exceedRaw++; }
      const pValue = (1 + exceed) / (R + 1);
      const zFull = new Float64Array(nD).fill(NaN); for (let j = 0; j < nU; j++) zFull[uidx[j]] = zObs[j];
      let znm = -Infinity; for (const z of zNullMax) if (z > znm) znm = z;
      return { samples: R, alpha: s.alpha, p_value: pValue, abstain: pValue > s.alpha,
        statistic: "max standardized score (max-T, leave-one-out, pooled variance floor)", z_max: zMax, variance_floor: v0,
        z_max_disease: this.diseaseLabel(uidx[jBest]), null_z_max_mean: seqSum(zNullMax) / R, null_z_max_max: znm,
        p_value_raw_max: (1 + exceedRaw) / (R + 1), observed_max_score: obsRaw, null_max_score_mean: seqSum(rawMax) / R, _z: zFull };
    }

    // ----------------------------------------------------------- evidence
    evidence(d, q) {
      const kb = this.kb, prop = kb.prop[d], ann = kb.annTerms[d], freqs = kb.annFreqs[d];
      const rows = q.present_used.map((p) => {
        const mica = this.ancByIcDesc(p).find((a) => bsearch(prop, a));
        let best = -1;
        for (let j = 0; j < ann.length; j++) if (this.isAncestor(mica, ann[j]) && (best < 0 || kb.ic[ann[j]] > kb.ic[ann[best]])) best = j;
        const f = best < 0 || Number.isNaN(freqs[best]) ? null : freqs[best];
        return { term: this.termLabel(p), match: best >= 0 ? this.termLabel(ann[best]) : null, mica: this.termLabel(mica), ic: kb.ic[mica], frequency: f };
      });
      const absentRows = [];
      for (const e of q.absent_used) {
        const desc = this.descendants(e); let best = -1, bf = 0.0;
        for (let j = 0; j < ann.length; j++) if (bsearch(desc, ann[j])) { const fe = Number.isNaN(freqs[j]) ? 0.5 : freqs[j]; if (fe > bf) { best = j; bf = fe; } }
        if (best >= 0) absentRows.push({ term: this.termLabel(e), annotated: this.termLabel(ann[best]), frequency: bf, ic: kb.ic[e] });
      }
      return { disease: this.diseaseLabel(d), present: rows, absent: absentRows };
    }

    // ------------------------------------------------ next-best phenotypes
    nextPhenotypes(q, s, post, order, nOut) {
      const kb = this.kb, top = order.slice(0, s.eig_top_k);
      if (!top.length) return [];
      const piRaw = Float64Array.from(top, (d) => post[d]), Z = seqSum(piRaw), pi = piRaw.map((x) => x / Z);
      const qmaps = [], cand = new Set();
      for (const d of top) {
        const qm = new Map(), ann = kb.annTerms[d], fr = kb.annFreqs[d];
        for (let j = 0; j < ann.length; j++) {
          const fe = Number.isNaN(fr[j]) ? 0.5 : fr[j];
          for (const a of kb.ancestors[ann[j]]) { const cur = qm.has(a) ? qm.get(a) : -1.0; if (fe > cur) qm.set(a, fe); }
          cand.add(ann[j]);
        }
        qmaps.push(qm);
      }
      const blocked = new Set();
      for (const p of q.present) for (const a of kb.ancestors[p]) blocked.add(a);
      for (const e of q.absent) for (const t of this.descendants(e)) blocked.add(t);
      const terms = Array.from(cand).sort((a, b) => a - b).filter((t) => !blocked.has(t) && kb.ic[t] >= s.eig_min_ic);
      if (!terms.length) return [];
      const col = new Map(terms.map((t, j) => [t, j])), K = top.length, C = terms.length;
      const Q = new Float64Array(K * C).fill(s.eps);
      qmaps.forEach((qm, i) => { for (const [t, v] of qm) { const j = col.get(t); if (j !== undefined) Q[i * C + j] = Math.min(Math.max(v, s.eps), s.q_max); } });
      const pYes = new Float64Array(C), cond = new Float64Array(C);
      for (let i = 0; i < K; i++) for (let j = 0; j < C; j++) { const qq = Q[i * C + j]; pYes[j] += pi[i] * qq; cond[j] += pi[i] * h2(qq); }
      const eig = new Float64Array(C); for (let j = 0; j < C; j++) eig[j] = h2(pYes[j]) - cond[j];
      const idx = Array.from({ length: C }, (_, j) => j);
      const key = (j) => Math.floor(eig[j] * 1e9 + 0.5);
      idx.sort((a, b) => (key(b) - key(a)) || (terms[a] - terms[b]));
      return idx.slice(0, nOut).map((j) => {
        const sup = [], against = [];
        for (let i = 0; i < K; i++) { const v = Q[i * C + j]; if (v >= 0.5 && sup.length < 3) sup.push(top[i]); if (v === s.eps && against.length < 3) against.push(top[i]); }
        return Object.assign(this.termLabel(terms[j]), { eig_bits: eig[j], p_present: pYes[j], ic: kb.ic[terms[j]],
          supports: sup.map((d) => this.diseaseLabel(d)), against: against.map((d) => this.diseaseLabel(d)) });
      });
    }

    // ---------------------------------------------------------- hypotheses
    hypotheses(q, s, u, comp, post, order, nOut) {
      const kb = this.kb;
      if (!u.gIds.length || !q.present_used.length) return [];
      const top = order.slice(0, s.hyp_top_diseases), excluded = new Set();
      for (const d of top) for (const g of kb.diseaseGenes[d]) excluded.add(g);
      const nR = comp.bm.length, nC = u.gIds.length, gbm = new Float64Array(nR * nC);
      for (let r = 0; r < nR; r++) {
        const row = comp.bm[r];
        for (let c = 0; c < nC; c++) {
          const lo = u.gOff[c], hi = c + 1 < nC ? u.gOff[c + 1] : u.gDis.length;
          let best = row[u.gDis[lo]]; for (let k = lo + 1; k < hi; k++) if (row[u.gDis[k]] > best) best = row[u.gDis[k]];
          gbm[r * nC + c] = best;
        }
      }
      const acc = new Float64Array(nC), nEv = new Int32Array(nC);
      for (let r = 0; r < nR; r++) for (let c = 0; c < nC; c++) { const v = gbm[r * nC + c]; const hit = v >= s.hyp_min_ic; acc[c] += hit ? v : 0.0; if (hit) nEv[c]++; }
      const score = acc.map((x) => x / comp.n);
      const need = Math.min(s.hyp_min_terms, q.present_used.length);
      const cols = [];
      for (let c = 0; c < nC; c++) if (!excluded.has(u.gIds[c]) && nEv[c] >= need) cols.push(c);
      cols.sort((a, b) => (score[b] - score[a]) || (nEv[b] - nEv[a]) || (u.gIds[a] - u.gIds[b]));
      return cols.slice(0, nOut).map((c, k) => {
        const g = u.gIds[c], udis = Array.from(kb.geneDiseases[g]).filter((d) => u.mask[d]), ev = [];
        q.present_used.forEach((p, r) => {
          if (!(gbm[r * nC + c] >= s.hyp_min_ic)) return;
          let mica = null, via = null;
          for (const a of this.ancByIcDesc(p)) { const hitD = udis.find((d) => bsearch(kb.prop[d], a)); if (hitD !== undefined) { mica = a; via = hitD; break; } }
          ev.push({ term: this.termLabel(p), mica: this.termLabel(mica), ic: kb.ic[mica], via_disease: this.diseaseLabel(via) });
        });
        const known = udis.slice().sort((a, b) => (post[b] - post[a]) || (a - b)).slice(0, 3);
        return { rank: k + 1, symbol: kb.geneSymbols[g], ncbi: kb.geneNcbi[g], score: score[c], n_evidence: nEv[c], evidence: ev,
          known_diseases: known.map((d) => this.diseaseLabel(d)),
          label: "Hypothesis: shares specific phenotypes with this patient but is not linked to any of the " +
            `top ${top.length} ranked diseases. Not an established association; research follow-up only.` };
      });
    }

    resolveDiseases(ids) {
      const found = [], missing = [];
      for (const raw of ids || []) { const d = this.kb.diseaseIndex.get(String(raw).trim()); if (d === undefined) missing.push(String(raw)); else found.push(d); }
      return [Array.from(new Set(found)).sort((a, b) => a - b), missing];
    }

    // ------------------------------------------------------------- analyze
    analyze(present, absent, opts) {
      opts = opts || {};
      const kb = this.kb, s = this.settings(opts), q = this.prepare(present, absent || []);
      const parts = opts.parts || ["rank", "next", "hypotheses"];
      const topK = opts.top_k != null ? opts.top_k : 20, nNext = opts.n_next != null ? opts.n_next : 10;
      const nHyp = opts.n_hyp != null ? opts.n_hyp : 10, evidenceK = opts.evidence_k != null ? opts.evidence_k : 10;
      if (!q.present_used.length) throw new QueryError("At least one valid present HPO term (Phenotypic abnormality subtree) is required.");
      const [excl, missing] = this.resolveDiseases(opts.exclude_diseases);
      for (const m of missing) q.warnings.push(`Excluded disease ${m} is not in the knowledge base.`);
      const u = this.universe(s.universe, excl);
      if (u.uidx.length < 2) throw new QueryError("The disease universe is empty after exclusions.");
      const comp = this.components(q.present_used, q.absent_used, s.mode !== "asymmetric");
      const score = Engine.combine(comp, s.mode, s.beta, s.w);
      const order = this.order(score, u.uidx);
      const post = this.posterior(score, u.uidx, s.T);
      let H = 0; for (const d of u.uidx) { const p = post[d]; H += p > 0 ? p * Math.log2(p) : 0; }
      const entropy = -H;
      const gpost = new Float64Array(kb.nG);
      for (let i = 0; i < u.dgDis.length; i++) gpost[u.dgGene[i]] += post[u.dgDis[i]] / kb.kGenes[u.dgDis[i]];
      const geneCand = Array.from(new Set(u.dgGene)).sort((a, b) => a - b);
      const gkey = (g) => Math.floor(gpost[g] * 1e12 + 0.5);
      const gorder = geneCand.sort((a, b) => (gkey(b) - gkey(a)) || (a - b));
      let withGenes = 0; for (const d of u.disWithGenes) withGenes += post[d];
      const noGene = 1.0 - withGenes;
      const result = {
        engine: "genoarch-engine " + this.version,
        query: { present: q.present.map((p) => this.termLabel(p)), absent: q.absent.map((e) => this.termLabel(e)),
          present_used: q.present_used.map((p) => this.termLabel(p)), absent_used: q.absent_used.map((e) => this.termLabel(e)),
          unresolved: q.unresolved, warnings: q.warnings },
        settings: { mode: s.mode, universe: s.universe, T: s.T, beta: s.beta, w: s.w, alpha: s.alpha, null_samples: s.null_samples,
          calibrated: s.calibrated, n_universe: u.uidx.length, excluded_diseases: excl.map((d) => kb.diseaseIds[d]) },
        research_use_only: RESEARCH_USE_ONLY,
      };
      const top = order.slice(0, topK);
      if (parts.includes("rank")) {
        const diseases = top.map((d, r) => Object.assign({ rank: r + 1 }, this.diseaseLabel(d), { score: score[d], posterior: post[d],
          genes: kb.diseaseGenes[d].slice().sort((a, b) => a - b).map((g) => kb.geneSymbols[g]) }));
        const genes = gorder.slice(0, topK).map((g, r) => {
          const linked = Array.from(kb.geneDiseases[g]).filter((d) => u.mask[d]).sort((a, b) => (post[b] - post[a]) || (a - b));
          return { rank: r + 1, symbol: kb.geneSymbols[g], ncbi: kb.geneNcbi[g], posterior: gpost[g],
            diseases: linked.slice(0, 3).map((d) => Object.assign(this.diseaseLabel(d), { posterior: post[d] })), n_diseases: linked.length };
        });
        const p1 = post[top[0]], p2 = top.length > 1 ? post[top[1]] : 0.0;
        const nul = this.nullTest(q, s, u, comp, score);
        const zf = nul._z; delete nul._z;
        diseases.forEach((row, i) => { row.z_null = !zf || Number.isNaN(zf[top[i]]) ? null : zf[top[i]]; });
        result.diseases = diseases;
        result.genes = genes;
        result.uncertainty = { top1_posterior: p1, top1_margin: p1 - p2, entropy_bits: entropy,
          entropy_norm: u.uidx.length > 1 ? entropy / Math.log2(u.uidx.length) : 0.0,
          gene_top1_posterior: gorder.length ? gpost[gorder[0]] : 0.0, no_known_gene_mass: noGene, open_world: nul };
        result.evidence = top.slice(0, evidenceK).map((d) => this.evidence(d, q));
      }
      if (parts.includes("next")) result.next_phenotypes = this.nextPhenotypes(q, s, post, order, nNext);
      if (parts.includes("hypotheses")) result.hypotheses = this.hypotheses(q, s, u, comp, post, order, nHyp);
      return result;
    }

    // -------------------------------------------------------------- search
    search(query, limit) {
      limit = limit || 20;
      const kb = this.kb, qs = String(query || "").trim().toLowerCase();
      if (!qs) return [];
      const tok = (x) => x.match(/[a-z0-9]+/g) || [];
      if (!this._search) {
        this._search = kb.names.map((nm, i) => { const name = nm.toLowerCase(), syns = kb.synonyms[i].map((x) => x.toLowerCase());
          return [name, syns, tok(name), syns.map(tok)]; });
      }
      const out = [], seen = new Set();
      const m = qs.match(/^(?:hp[:_]?\s*)?(\d{1,7})$/);
      if (m && (qs.startsWith("hp") || m[1].length === 7)) { const [i] = this.resolve("HP:" + m[1].padStart(7, "0")); if (i !== null) { out.push([-1, 0, 0, i, null]); seen.add(i); } }
      const qtok = tok(qs);
      const allPrefix = (words) => qtok.length > 0 && qtok.every((t) => words.some((w) => w.startsWith(t)));
      this._search.forEach(([name, syns, ntok, stoks], i) => {
        if (seen.has(i)) return;
        let tier = null, syn = null, j;
        if (name === qs) tier = 0;
        else if ((j = syns.indexOf(qs)) >= 0) { tier = 1; syn = kb.synonyms[i][j]; }
        else if (name.startsWith(qs)) tier = 2;
        else if ((j = syns.findIndex((x) => x.startsWith(qs))) >= 0) { tier = 3; syn = kb.synonyms[i][j]; }
        else if (allPrefix(ntok)) tier = 4;
        else if ((j = stoks.findIndex((tk) => allPrefix(tk))) >= 0) { tier = 5; syn = kb.synonyms[i][j]; }
        else if (name.includes(qs)) tier = 6;
        else if ((j = syns.findIndex((x) => x.includes(qs))) >= 0) { tier = 7; syn = kb.synonyms[i][j]; }
        if (tier !== null) out.push([tier, -kb.cnt[i], kb.names[i].length, i, syn]);
      });
      out.sort((a, b) => (a[0] - b[0]) || (a[1] - b[1]) || (a[2] - b[2]) || (a[3] - b[3]));
      return out.slice(0, limit).map(([tier, , , i, syn]) => Object.assign(this.termLabel(i), { synonym: syn, ic: kb.ic[i], n_diseases: kb.cnt[i],
        match: tier < 0 ? "id" : (syn ? "synonym" : "name") }));
    }
  }

  async function decodeIndexBytes(bytes) {
    let u8 = bytes instanceof Uint8Array ? bytes : new Uint8Array(bytes);
    if (u8[0] === 0x1f && u8[1] === 0x8b) {
      if (typeof DecompressionStream === "undefined") throw new Error("This browser cannot decompress the index (DecompressionStream missing).");
      const ds = new DecompressionStream("gzip");
      const stream = new Blob([u8]).stream().pipeThrough(ds);
      u8 = new Uint8Array(await new Response(stream).arrayBuffer());
    }
    return JSON.parse(new TextDecoder().decode(u8));
  }
  function base64ToBytes(b64) { const bin = atob(b64.replace(/\s+/g, "")); const u8 = new Uint8Array(bin.length); for (let i = 0; i < bin.length; i++) u8[i] = bin.charCodeAt(i); return u8; }

  return { Engine, QueryError, ENGINE_VERSION, MODES, UNIVERSES, DEFAULTS, RESEARCH_USE_ONLY, fnv1a32, mulberry32, normalizeHpoId,
    decodeIndexBytes, base64ToBytes };
});
