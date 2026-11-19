/* Genoarch web app: routing, term entry, rendering. Analysis runs in this browser via engine.js. */
(function () {
  "use strict";
  const GE = window.GenoarchEngine;
  const ROUTES = ["home", "analyze", "method", "pricing", "about", "legal"];
  const EXAMPLES = [
    { key: "marfan", label: "Marfan-like", present: ["HP:0001166", "HP:0001083", "HP:0002616", "HP:0001519", "HP:0000545"], absent: [] },
    { key: "seizures", label: "Febrile and myoclonic seizures", present: ["HP:0002373", "HP:0002123", "HP:0001251", "HP:0001263"], absent: [] },
    { key: "nonspecific", label: "Nonspecific delay", present: ["HP:0001263", "HP:0001249"], absent: [] },
  ];
  const $ = (sel) => document.querySelector(sel);
  const esc = (s) => String(s == null ? "" : s).replace(/[&<>"']/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c]));
  const nf = (x) => Number(x).toLocaleString("en-US");
  const pct = (p) => (p == null ? "n/a" : p < 0.001 ? "<0.1%" : p < 0.1 ? (p * 100).toFixed(1) + "%" : Math.round(p * 100) + "%");
  const num = (x, d) => (x == null ? "n/a" : Number(x).toFixed(d));
  const pval = (p) => (p == null ? "n/a" : p < 0.001 ? p.toExponential(1) : p.toFixed(3));
  const hpoLink = (id) => `https://hpo.jax.org/browse/term/${encodeURIComponent(id)}`;
  const state = { engine: null, calibration: {}, method: null, terms: [], excluded: [], result: null, selected: 0,
    useServer: false, timer: null, token: 0, betaTouched: false };

  // ------------------------------------------------------------------ routing
  function route() {
    let r = (location.hash || "").replace(/^#/, "");
    if (!ROUTES.includes(r)) r = "home";
    document.querySelectorAll(".view").forEach((v) => { v.hidden = v.dataset.view !== r; });
    document.querySelectorAll(".nav a").forEach((a) => { if (a.dataset.route === r) a.setAttribute("aria-current", "page"); else a.removeAttribute("aria-current"); });
    if (r === "method") renderMethod();
  }

  // ------------------------------------------------------------------ data
  async function loadIndex() {
    const emb = document.getElementById("ga-index");
    if (emb && emb.textContent.trim()) return GE.decodeIndexBytes(GE.base64ToBytes(emb.textContent));
    const res = await fetch("data/index.json.gz");
    if (!res.ok) throw new Error(`The HPO index could not be loaded (HTTP ${res.status}).`);
    return GE.decodeIndexBytes(await res.arrayBuffer());
  }
  async function loadJSON(id, url) {
    const emb = document.getElementById(id);
    if (emb && emb.textContent.trim()) return JSON.parse(emb.textContent);
    try { const r = await fetch(url); return r.ok ? r.json() : null; } catch (e) { return null; }
  }

  // ------------------------------------------------------------------ arch
  function archSVG(res) {
    const ev = res.evidence && res.evidence[0];
    const rows = (ev ? ev.present : []).slice(0, 14);
    const ow = res.uncertainty.open_world, gene = res.genes[0];
    const n = rows.length, cx = 200, cy = 232, rIn = 112, gap = 1.6, key = 8.5;
    const left = Math.ceil(n / 2), right = n - left;
    const pol = (r, a) => [cx + r * Math.cos((a * Math.PI) / 180), cy - r * Math.sin((a * Math.PI) / 180)];
    const sector = (r0, r1, a0, a1) => {
      const [x1, y1] = pol(r1, a0), [x2, y2] = pol(r1, a1), [x3, y3] = pol(r0, a1), [x4, y4] = pol(r0, a0);
      return `M${x1.toFixed(2)} ${y1.toFixed(2)}A${r1} ${r1} 0 0 1 ${x2.toFixed(2)} ${y2.toFixed(2)}L${x3.toFixed(2)} ${y3.toFixed(2)}A${r0} ${r0} 0 0 0 ${x4.toFixed(2)} ${y4.toFixed(2)}Z`;
    };
    const blocks = [];
    const place = (count, from, to, startIdx, orderFromFoot) => {
      if (!count) return;
      const w = (from - to - gap * (count - 1)) / count;
      for (let k = 0; k < count; k++) {
        const a0 = from - k * (w + gap), a1 = a0 - w, row = rows[startIdx + k];
        const t = 16 + (34 * Math.min(row.ic, 9.5)) / 9.5;
        const [tx, ty] = pol(rIn + t / 2, (a0 + a1) / 2);
        blocks.push(`<path class="block" style="--i:${orderFromFoot(k)};opacity:${(0.38 + (0.62 * Math.min(row.ic, 9.5)) / 9.5).toFixed(3)}" d="${sector(rIn, rIn + t, a0, a1)}"><title>${esc(row.term.name)}: shares ${esc(row.mica.name)} (IC ${num(row.ic, 2)})</title></path>` +
          `<text class="num" x="${tx.toFixed(1)}" y="${ty.toFixed(1)}">${startIdx + k + 1}</text>`);
      }
    };
    place(left, 180, 90 + key + gap, 0, (k) => k);
    place(right, 90 - key - gap, 0, left, (k) => right - 1 - k);
    let ksClass = "keystone match", line1 = gene ? gene.symbol : "", line2 = gene ? `${pct(gene.posterior)} calibrated confidence` : "";
    if (ow.p_value == null) { ksClass = "keystone match"; line2 += ", open-world test off"; }
    else if (ow.abstain) { ksClass = "keystone empty"; line1 = "No confident match"; line2 = `abstain, p = ${pval(ow.p_value)}`; }
    const ks = `<path class="${ksClass}" style="--i:${Math.max(left, right) + 1}" d="${sector(rIn - 6, rIn + 52, 90 + key, 90 - key)}"><title>Keystone: top-ranked gene</title></path>`;
    const label = `<text class="gene" x="${cx}" y="${cy - 52}" ${ow.abstain ? 'style="font-size:17px"' : ""}>${esc(line1)}</text>` +
      `<text class="conf${ow.abstain ? " abstain" : ""}" x="${cx}" y="${cy - (ow.abstain ? 32 : 28)}">${esc(line2)}</text>`;
    const desc = ow.abstain ? "The arch is left without a keystone: no known disease fits these phenotypes better than chance." :
      `Arch of ${n} phenotype blocks; keystone is ${gene ? gene.symbol : "the top gene"} at ${pct(gene ? gene.posterior : 0)}.`;
    return `<svg class="arch" viewBox="0 0 400 250" role="img" aria-label="${esc(desc)}"><line class="ground" x1="20" y1="${cy}" x2="380" y2="${cy}"/>${blocks.join("")}${ks}${label}</svg>` +
      `<ol class="arch-legend">${rows.map((r, i) => `<li><b>${i + 1}</b>${esc(r.term.name)} <span>IC ${num(r.ic, 1)}</span></li>`).join("")}</ol>`;
  }

  // ------------------------------------------------------------------ analyze UI
  function presentIds() { return state.terms.filter((t) => t.status === "present").map((t) => t.id); }
  function absentIds() { return state.terms.filter((t) => t.status === "absent").map((t) => t.id); }

  function addTerm(id, status) {
    const eng = state.engine; if (!eng) return;
    const [i, st] = eng.resolve(id);
    if (i === null) return;
    const cid = eng.kb.termIds[i];
    const existing = state.terms.find((t) => t.id === cid);
    if (existing) existing.status = status || existing.status;
    else state.terms.push({ id: cid, name: eng.kb.names[i], status: status || "present", alt: st === "alt" ? id : null });
    renderTerms(); schedule();
  }
  function renderTerms() {
    const list = $("#term-list");
    list.innerHTML = state.terms.map((t, k) => `<li class="term ${t.status === "absent" ? "is-absent" : ""}">
      <span class="t-name">${esc(t.name)}</span>
      <span class="toggle" role="group" aria-label="${esc(t.name)} status">
        <button type="button" data-k="${k}" data-status="present" aria-pressed="${t.status === "present"}">Present</button><button type="button" data-k="${k}" data-status="absent" aria-pressed="${t.status === "absent"}">Absent</button>
      </span>
      <button type="button" class="remove" data-remove="${k}" aria-label="Remove ${esc(t.name)}">×</button>
      <a class="t-id" href="${hpoLink(t.id)}">${esc(t.id)}</a></li>`).join("");
    const np = presentIds().length, na = absentIds().length;
    $("#term-count").textContent = state.terms.length ? `${np} present, ${na} absent` : "";
    $("#term-empty").hidden = np > 0;
  }
  function renderExcluded() {
    $("#excluded-box").hidden = !state.excluded.length;
    $("#excluded-list").innerHTML = state.excluded.map((d, k) => `<li>${esc(d.name)} <button type="button" class="remove" data-unexclude="${k}" aria-label="Restore ${esc(d.name)}">×</button></li>`).join("");
  }

  function currentOptions() {
    const universe = document.querySelector('input[name="universe"]:checked').value;
    const o = { mode: $("#mode-select").value, universe, null_samples: Number($("#null-select").value),
      exclude_diseases: state.excluded.map((d) => d.id), top_k: 20 };
    if (state.betaTouched && $("#beta-input").value !== "") o.beta = Number($("#beta-input").value);
    return o;
  }
  function fittedBeta() {
    const o = { mode: $("#mode-select").value, universe: document.querySelector('input[name="universe"]:checked').value };
    return state.engine ? state.engine.settings(o).beta : 0;
  }
  function schedule() { clearTimeout(state.timer); state.timer = setTimeout(run, 180); }

  async function run() {
    if (!state.engine) return;
    const present = presentIds(), absent = absentIds();
    if (!present.length) { $("#results").hidden = true; setStatus("Add at least one present phenotype to see a ranking."); return; }
    const token = ++state.token;
    setStatus("Analyzing…");
    await new Promise((r) => setTimeout(r, 20));
    try {
      let res;
      const opts = currentOptions();
      if (state.useServer) {
        const r = await fetch("analyze", { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify(Object.assign({ present, absent }, opts)) });
        const j = await r.json();
        if (!r.ok) throw new Error(typeof j.detail === "string" ? j.detail : "The server rejected the request.");
        res = j;
      } else res = state.engine.analyze(present, absent, opts);
      if (token !== state.token) return;
      state.result = res; state.selected = 0;
      renderResults();
    } catch (e) { setStatus(e.message || String(e), true); $("#results").hidden = true; }
  }
  function setStatus(msg, isError) { const s = $("#status"); s.textContent = msg; s.classList.toggle("error", !!isError); s.hidden = !msg; }

  function renderResults() {
    const res = state.result, u = res.uncertainty, ow = u.open_world, s = res.settings;
    const warn = res.query.warnings.concat(res.query.unresolved.map((x) => `${x.input} was not used (${x.status.replace("_", " ")}).`));
    setStatus(warn.join(" "));
    $("#results").hidden = false;
    $("#result-arch").innerHTML = archSVG(res);
    const owBadge = ow.p_value == null ? `<span class="badge off">Open-world test off</span>` : ow.abstain ?
      `<span class="badge abstain">Abstain: no known disease fits better than chance</span> <span class="muted">p = ${pval(ow.p_value)}, ${ow.samples} random patients</span>` :
      `<span class="badge match">Known-disease match</span> <span class="muted">p = ${pval(ow.p_value)}, ${ow.samples} random patients</span>`;
    const g = res.genes[0], d = res.diseases[0];
    $("#verdicts").innerHTML = [
      ["Top gene", g ? `<span class="value">${esc(g.symbol)}</span> <span class="muted">${pct(g.posterior)} calibrated confidence</span>` : "none"],
      ["Top disease", d ? `<span class="value">${esc(d.name)}</span> <span class="muted">${pct(d.posterior)}, ${esc(d.id)}</span>` : "none"],
      ["Open-world test", owBadge],
      ["Spread", `<span class="value">${num(u.entropy_bits, 2)} bits</span> <span class="muted">entropy, ${pct(u.entropy_norm)} of the maximum; top-1 margin ${pct(u.top1_margin)}</span>`],
      ["Without a linked gene", `<span class="value">${pct(u.no_known_gene_mass)}</span> <span class="muted">of the disease probability</span>`],
      ["Settings", `<span class="muted">${esc(modeName(s.mode))}, ${nf(s.n_universe)} diseases, T ${num(s.T, 3)}, absent weight ${num(s.beta, 2)}${s.calibrated ? ", fitted values" : ", user override"}</span>`],
    ].map(([l, v]) => `<div class="verdict"><span class="label">${l}</span>${v}</div>`).join("");
    $("#gene-list").innerHTML = res.genes.slice(0, 15).map((x) => `<li><span class="r">${x.rank}</span><span class="n">${esc(x.symbol)}</span><span class="p">${pct(x.posterior)}</span>
      <span class="bar" aria-hidden="true"><span style="width:${Math.max(0, Math.min(100, x.posterior * 100)).toFixed(2)}%"></span></span>
      <span class="sub">${x.diseases.slice(0, 2).map((y) => esc(y.name)).join("; ")}${x.n_diseases > 2 ? `; ${x.n_diseases - 2} more` : ""}</span></li>`).join("");
    $("#disease-list").innerHTML = res.diseases.slice(0, 15).map((x, k) => `<li class="selectable ${k === state.selected ? "selected" : ""}" data-select="${k}" tabindex="0" aria-label="${esc(x.name)}, ${pct(x.posterior)}. Show evidence">
      <span class="r">${x.rank}</span><span class="n">${esc(x.name)}</span><span class="p">${pct(x.posterior)}</span>
      <span class="bar" aria-hidden="true"><span style="width:${Math.max(0, Math.min(100, x.posterior * 100)).toFixed(2)}%"></span></span>
      <span class="sub">${esc(x.id)}${x.genes.length ? ", " + esc(x.genes.slice(0, 4).join(", ")) : ", no linked gene"}${x.z_null != null ? `, z ${num(x.z_null, 1)}` : ""}</span>
      <span class="row-actions"><button type="button" class="link-button" data-exclude="${k}">Rule out</button></span></li>`).join("");
    $("#next-list").innerHTML = (res.next_phenotypes || []).map((x) => `<li><span class="n">${esc(x.name)}</span>
      <span class="act"><button type="button" class="small-button" data-add="${esc(x.id)}" data-status="present">Mark present</button><button type="button" class="small-button absent" data-add="${esc(x.id)}" data-status="absent">Mark absent</button></span>
      <span class="meta"><a href="${hpoLink(x.id)}">${esc(x.id)}</a>, ${num(x.eig_bits, 2)} bits expected gain, present in ${pct(x.p_present)} under the current ranking${x.supports.length ? `; supports ${esc(x.supports[0].name)}` : ""}</span></li>`).join("") ||
      `<li class="muted">No informative phenotype is left to ask about for the top diseases.</li>`;
    $("#hyp-list").innerHTML = (res.hypotheses || []).map((h) => `<li><div class="head"><span class="n">${esc(h.symbol)}</span><span class="tag">Hypothesis</span>
      <span class="muted">${h.n_evidence} matching terms, mean matched IC ${num(h.score, 2)}</span></div>
      <div class="ev">${h.evidence.map((e) => `${esc(e.term.name)} <span class="muted">(shares ${esc(e.mica.name)}, IC ${num(e.ic, 1)}, via ${esc(e.via_disease.name)})</span>`).join("; ")}</div>
      <div class="known">Known links: ${h.known_diseases.map((y) => esc(y.name)).join("; ")}</div></li>`).join("") ||
      `<li class="muted">No gene outside the top diseases matches two or more specific phenotypes of this patient.</li>`;
    renderEvidence();
  }
  function modeName(m) { return { blended: "Blended best-match average", asymmetric: "Patient-to-disease best-match average", symmetric: "Symmetric best-match average" }[m] || m; }
  function renderEvidence() {
    const res = state.result, ev = res.evidence[state.selected];
    if (!ev) { $("#ev-table").innerHTML = ""; return; }
    $("#ev-sub").textContent = `Why ${ev.disease.name} (${ev.disease.id}) matches. Select another disease in the list to compare.`;
    const pres = ev.present.map((r) => `<tr><td>${esc(r.term.name)}</td><td>${r.match ? esc(r.match.name) : "none"}</td><td>${esc(r.mica.name)}</td><td class="num">${num(r.ic, 2)}</td><td class="num">${r.frequency == null ? "unknown" : pct(r.frequency)}</td></tr>`).join("");
    const abs = ev.absent.length ? `<tr class="section"><td colspan="5">Absent phenotypes annotated to this disease</td></tr>` +
      ev.absent.map((r) => `<tr><td>${esc(r.term.name)}</td><td>${esc(r.annotated.name)}</td><td></td><td class="num">${num(r.ic, 2)}</td><td class="num">${pct(r.frequency)}</td></tr>`).join("") : "";
    $("#ev-table").innerHTML = `<thead><tr><th scope="col">Patient phenotype</th><th scope="col">Closest disease annotation</th><th scope="col">Shared ancestor</th><th scope="col" class="num">IC (nats)</th><th scope="col" class="num">Frequency in disease</th></tr></thead><tbody>${pres}${abs}</tbody>`;
  }

  // ------------------------------------------------------------------ search box
  function setupSearch() {
    const input = $("#term-search"), list = $("#term-results");
    let items = [], active = -1, timer = null;
    const close = () => { list.hidden = true; input.setAttribute("aria-expanded", "false"); active = -1; input.removeAttribute("aria-activedescendant"); };
    const paint = () => {
      list.innerHTML = items.map((h, k) => `<li id="opt-${k}" role="option" aria-selected="${k === active}" data-k="${k}"><span class="s-name">${esc(h.name)}</span>
        <span class="s-meta">${esc(h.id)}${h.synonym ? `, matched “${esc(h.synonym)}”` : ""}, ${nf(h.n_diseases)} diseases</span></li>`).join("");
      if (active >= 0) input.setAttribute("aria-activedescendant", `opt-${active}`);
    };
    const choose = (k) => { const h = items[k]; if (!h) return; addTerm(h.id, "present"); input.value = ""; close(); input.focus(); };
    input.addEventListener("input", () => {
      clearTimeout(timer);
      timer = setTimeout(() => {
        if (!state.engine) return;
        items = state.engine.search(input.value, 12); active = items.length ? 0 : -1;
        if (!items.length) { close(); return; }
        paint(); list.hidden = false; input.setAttribute("aria-expanded", "true");
      }, 90);
    });
    input.addEventListener("keydown", (e) => {
      if (list.hidden) return;
      if (e.key === "ArrowDown") { active = Math.min(items.length - 1, active + 1); paint(); e.preventDefault(); }
      else if (e.key === "ArrowUp") { active = Math.max(0, active - 1); paint(); e.preventDefault(); }
      else if (e.key === "Enter") { choose(active); e.preventDefault(); }
      else if (e.key === "Escape") close();
    });
    list.addEventListener("mousedown", (e) => { const li = e.target.closest("li[data-k]"); if (li) { e.preventDefault(); choose(Number(li.dataset.k)); } });
    input.addEventListener("blur", () => setTimeout(close, 120));
  }

  // ------------------------------------------------------------------ controls
  function setupControls() {
    $("#example-buttons").innerHTML = EXAMPLES.map((x) => `<button type="button" class="small-button" data-example="${x.key}">${esc(x.label)}</button>`).join("");
    document.addEventListener("click", (e) => {
      const t = e.target;
      if (t.closest("[data-example]")) { loadExample(t.closest("[data-example]").dataset.example); return; }
      if (t.dataset.status && t.dataset.k !== undefined) { state.terms[Number(t.dataset.k)].status = t.dataset.status; renderTerms(); schedule(); return; }
      if (t.dataset.remove !== undefined) { state.terms.splice(Number(t.dataset.remove), 1); renderTerms(); schedule(); return; }
      if (t.dataset.add) { addTerm(t.dataset.add, t.dataset.status); return; }
      if (t.dataset.exclude !== undefined) { const d = state.result.diseases[Number(t.dataset.exclude)]; state.excluded.push({ id: d.id, name: d.name }); renderExcluded(); schedule(); e.stopPropagation(); return; }
      if (t.dataset.unexclude !== undefined) { state.excluded.splice(Number(t.dataset.unexclude), 1); renderExcluded(); schedule(); return; }
      const row = t.closest("[data-select]");
      if (row) { state.selected = Number(row.dataset.select); document.querySelectorAll("#disease-list li").forEach((li, k) => li.classList.toggle("selected", k === state.selected)); renderEvidence(); }
    });
    document.addEventListener("keydown", (e) => {
      const row = e.target.closest && e.target.closest("[data-select]");
      if (row && (e.key === "Enter" || e.key === " ") && e.target === row) { e.preventDefault(); row.click(); }
    });
    $("#clear-terms").addEventListener("click", () => { state.terms = []; state.excluded = []; renderTerms(); renderExcluded(); schedule(); });
    document.querySelectorAll('input[name="universe"]').forEach((r) => r.addEventListener("change", () => { if (!state.betaTouched) $("#beta-input").value = fittedBeta(); schedule(); }));
    $("#mode-select").addEventListener("change", () => { if (!state.betaTouched) $("#beta-input").value = fittedBeta(); schedule(); });
    $("#beta-input").addEventListener("input", () => { state.betaTouched = true; schedule(); });
    $("#null-select").addEventListener("change", schedule);
    $("#engine-select").addEventListener("change", (e) => { state.useServer = e.target.value === "server"; schedule(); });
    $("#reset-settings").addEventListener("click", () => { applyDefaults(); schedule(); });
    $("#copy-json").addEventListener("click", async () => {
      const note = $("#copy-note");
      try { await navigator.clipboard.writeText(JSON.stringify(state.result, null, 1)); note.textContent = "Copied the results to the clipboard."; }
      catch (err) { note.textContent = "The browser blocked clipboard access. Use the API or the Python package to export results."; }
    });
  }
  function applyDefaults() {
    const s = state.engine.settings({});
    const radio = document.querySelector(`input[name="universe"][value="${s.universe}"]`); if (radio) radio.checked = true;
    $("#mode-select").value = s.mode; $("#beta-input").value = s.beta; $("#null-select").value = String(s.null_samples);
    state.betaTouched = false;
  }
  function loadExample(key, quiet) {
    const ex = EXAMPLES.find((x) => x.key === key); if (!ex || !state.engine) return;
    state.terms = []; state.excluded = [];
    ex.present.forEach((id) => addTermQuiet(id, "present")); ex.absent.forEach((id) => addTermQuiet(id, "absent"));
    renderTerms(); renderExcluded(); if (!quiet) schedule();
  }
  function addTermQuiet(id, status) {
    const [i] = state.engine.resolve(id); if (i === null) return;
    state.terms.push({ id: state.engine.kb.termIds[i], name: state.engine.kb.names[i], status });
  }
  async function detectServer() {
    if (location.protocol === "file:") return;
    try {
      const r = await fetch("health", { headers: { Accept: "application/json" } });
      const j = r.ok ? await r.json() : null;
      if (j && j.status === "ok") $("#engine-switch").hidden = false;
    } catch (e) { /* static hosting: no API */ }
  }

  // ------------------------------------------------------------------ method page
  function renderMethod() {
    const m = state.method, box = $("#method-body");
    if (!m) { box.innerHTML = state.engine ? `<p class="status">The method card is not available in this build.</p>` : `<p class="status">Loading the method card…</p>`; return; }
    if (box.dataset.done) return;
    const md = m.method, pm = m.parameters, data = m.data, ev = m.evaluation || {};
    const dc = ev.default_config || {};
    const defs = [["Information content", md.information_content], ["Similarity", md.similarity],
      ["Blended mode (default)", md.modes.blended], ["Asymmetric mode", md.modes.asymmetric], ["Symmetric mode", md.modes.symmetric],
      ["Absent phenotypes", md.absent_terms], ["Disease probability", md.posterior], ["Gene probability", md.gene_posterior],
      ["Uncertainty", md.uncertainty], ["Open-world test", md.open_world], ["Next phenotypes", md.next_phenotypes], ["Gene hypotheses", md.hypotheses]];
    const fitted = Object.entries(pm.fitted || {}).flatMap(([mode, us]) => Object.entries(us).map(([u, p]) =>
      `<tr><td>${esc(mode)}</td><td>${esc(u)}</td><td class="num">${num(p.T, 4)}</td><td class="num">${num(p.beta, 2)}</td><td class="num">${num(p.w, 2)}</td></tr>`)).join("");
    const metric = (label, v, fmt) => `<tr><td>${label}</td><td class="num">${v == null ? "n/a" : fmt(v)}</td></tr>`;
    const f3 = (v) => Number(v).toFixed(3);
    const ow = ev.open_world || {}, sim = (ev.next_phenotype_simulation || {});
    const rel = (ev.reliability_10 || []).filter((b) => b.n > 0);
    const relSVG = rel.length ? reliabilitySVG(rel) : "";
    const allCfg = Object.entries(ev.all_configs || {}).map(([k, v]) => `<tr><td>${esc(k)}</td><td class="num">${f3(v.gene_top1)}</td><td class="num">${f3(v.gene_top10)}</td><td class="num">${f3(v.disease_top1)}</td><td class="num">${f3(v.disease_top10)}</td><td class="num">${f3(v.ece_15)}</td></tr>`).join("");
    box.innerHTML = `<p class="ruo-strong">${esc(m.research_use_only)}</p>
      <p>${esc(m.intended_use)}</p>
      <h2>Not intended for</h2><ul>${m.not_intended_for.map((x) => `<li>${esc(x)}</li>`).join("")}</ul>
      <h2>Inputs</h2><p>${esc(m.inputs)}</p>
      <h2>Definitions</h2><dl class="defs">${defs.map(([k, v]) => `<dt>${esc(k)}</dt><dd>${esc(v)}</dd>`).join("")}</dl>
      <h2>Fitted parameters</h2>
      <p>Default: ${esc(pm.default_mode)} mode, ${esc(pm.default_universe)} disease sources. Fitted on ${esc((pm.fitted_on || {}).dataset || "n/a")}, ${nf((pm.fitted_on || {}).n_cases || 0)} cases from ${nf((pm.fitted_on || {}).n_diseases || 0)} diseases (${esc((pm.fitted_on || {}).split_rule || "")}).</p>
      <div class="table-scroll"><table class="metrics"><thead><tr><th>Mode</th><th>Sources</th><th class="num">T</th><th class="num">Absent weight</th><th class="num">w</th></tr></thead><tbody>${fitted}</tbody></table></div>
      <p class="muted">Constants: ${Object.entries(pm.constants || {}).map(([k, v]) => `${esc(k)} ${esc(v)}`).join(", ")}.</p>
      <h2>Evaluation</h2>
      <p>Held-out test slice of Phenopacket Store with no disease shared with the calibration slice: ${nf(ev.n_cases || 0)} cases, ${nf(ev.n_diseases || 0)} diseases. Default configuration. These are solved phenopacket cases used as a retrospective proxy for undiagnosed cases; no open cohort of truly undiagnosed patients is available to us, so the figures below describe solved, published cases. Log: <code>${esc(ev.log || "")}</code>.</p>
      <div class="table-scroll"><table class="metrics"><tbody>
        ${metric("Gene top-1 accuracy", dc.gene_top1, f3)}${metric("Gene top-3 accuracy", dc.gene_top3, f3)}${metric("Gene top-10 accuracy", dc.gene_top10, f3)}${metric("Gene mean reciprocal rank", dc.gene_mrr, f3)}
        ${metric("Disease top-1 accuracy", dc.disease_top1, f3)}${metric("Disease top-10 accuracy", dc.disease_top10, f3)}${metric("Disease mean reciprocal rank", dc.disease_mrr, f3)}
        ${metric("Negative log-likelihood of the true disease", dc.nll, f3)}${metric("Expected calibration error, 15 bins", dc.ece_15, f3)}${metric("Expected calibration error before temperature scaling", dc.ece_15_T1, f3)}
        ${metric("AUROC of confidence for a correct top-1", dc.auroc_conf_vs_correct, f3)}${metric("Share of cases with confidence at least 0.5", dc["coverage_conf_ge_0.5"], f3)}${metric("Top-1 accuracy within that share", dc["accuracy_conf_ge_0.5"], f3)}
      </tbody></table></div>
      ${relSVG}
      <h3>All configurations</h3>
      <div class="table-scroll"><table class="metrics"><thead><tr><th>Mode and sources</th><th class="num">Gene top-1</th><th class="num">Gene top-10</th><th class="num">Disease top-1</th><th class="num">Disease top-10</th><th class="num">ECE</th></tr></thead><tbody>${allCfg}</tbody></table></div>
      <h3>Open-world test</h3>
      <p>On ${nf(ow.n_cases || 0)} test cases with ${nf(ow.null_samples || 0)} random patients each (alpha ${esc(ow.alpha)}): abstained on ${pct(ow.abstain_rate_in_kb)} of cases whose disease is in the index, ${pct(ow.abstain_rate_gene_removed)} after removing every disease of the causal gene, and ${pct(ow.abstain_rate_random_patient)} of random phenotype sets.</p>
      <h3>Next-phenotype suggestions</h3>
      <p>Simulation starting from half of each patient’s terms (${nf((sim.eig_top1 || {}).n || 0)} cases): the top suggestion could be answered from the full record in ${nf(((sim.eig_top1 || {}).answer_present || 0) + ((sim.eig_top1 || {}).answer_absent || 0))} cases and improved the rank of the true disease in ${pct((sim.eig_top1 || {}).fraction_rank_improved_when_answered)} of those; a random candidate could be answered in ${nf(((sim.random_candidate || {}).answer_present || 0) + ((sim.random_candidate || {}).answer_absent || 0))} cases.</p>
      <h2>Data</h2>
      <p>HPO ontology release ${esc(data.hpo_ontology_version)}, annotation release ${esc(data.hpo_annotation_version)}: ${nf(data.counts.terms)} phenotype terms, ${nf(data.counts.diseases)} diseases, ${nf(data.counts.genes)} genes, ${nf(data.counts.annotations || 0)} annotations. ${esc(data.evaluation_dataset)}</p>
      <ul>${Object.entries(data.inputs || {}).map(([k, v]) => `<li><a href="${esc(v.url)}">${esc(k)}</a>, sha256 <code>${esc(v.sha256)}</code></li>`).join("")}</ul>
      <h2>Limitations</h2><ul>${m.limitations.map((x) => `<li>${esc(x)}</li>`).join("")}</ul>
      <h2>License</h2><p>${esc(m.license)}</p>
      <h2>AI-use disclosure</h2><p>${esc(m.ai_use_disclosure)}</p>`;
    box.dataset.done = "1";
  }
  function reliabilitySVG(bins) {
    const W = 360, H = 260, L = 44, B = 36, T = 12, R = 12, pw = W - L - R, ph = H - T - B;
    const x = (v) => L + v * pw, y = (v) => T + (1 - v) * ph;
    const ticks = [0, 0.25, 0.5, 0.75, 1];
    const grid = ticks.map((t) => `<line class="grid" x1="${x(0)}" y1="${y(t)}" x2="${x(1)}" y2="${y(t)}"/><text x="${L - 6}" y="${y(t) + 4}" text-anchor="end">${t}</text><text x="${x(t)}" y="${H - B + 16}" text-anchor="middle">${t}</text>`).join("");
    const pts = bins.map((b) => [x(b.mean_conf), y(b.accuracy), b.n]);
    return `<figure class="reliability"><svg viewBox="0 0 ${W} ${H}" role="img" aria-label="Reliability diagram of top-1 disease confidence on the test slice">
      ${grid}<line class="diag" x1="${x(0)}" y1="${y(0)}" x2="${x(1)}" y2="${y(1)}"/>
      <polyline class="ln" points="${pts.map((p) => p[0].toFixed(1) + "," + p[1].toFixed(1)).join(" ")}"/>
      ${pts.map((p) => `<circle class="pt" cx="${p[0].toFixed(1)}" cy="${p[1].toFixed(1)}" r="${(2.5 + Math.sqrt(p[2]) / 12).toFixed(1)}"><title>${p[2]} cases</title></circle>`).join("")}
      <text x="${L + pw / 2}" y="${H - 4}" text-anchor="middle">mean confidence in bin</text>
      <text x="12" y="${T + ph / 2}" text-anchor="middle" transform="rotate(-90 12 ${T + ph / 2})">observed accuracy</text></svg>
      <figcaption class="muted">Reliability of the top-1 disease confidence (10 bins; dot area grows with the number of cases; dashed line is perfect calibration).</figcaption></figure>`;
  }

  // ------------------------------------------------------------------ boot
  async function boot() {
    route();
    window.addEventListener("hashchange", () => { route(); window.scrollTo(0, 0); });
    setupSearch(); setupControls();
    try {
      const [doc, cal, mc] = await Promise.all([loadIndex(), loadJSON("ga-calibration", "data/calibration.json"), loadJSON("ga-method", "data/method_card.json")]);
      state.calibration = cal || {}; state.method = mc;
      state.engine = new GE.Engine(doc, state.calibration);
      const meta = state.engine.kb.meta, counts = meta.counts || {}, ann = (meta.annotation_header || {}).version;
      const fo = state.calibration.fitted_on || {};
      $("#home-data-line").textContent = `HPO release ${meta.ontology_version}, annotations ${ann}: ${nf(counts.terms)} phenotype terms, ${nf(counts.diseases)} diseases and ${nf(counts.genes)} genes.` +
        (fo.n_cases ? ` Confidence calibrated on ${nf(fo.n_cases)} Phenopacket Store cases.` : "");
      $("#footer-versions").textContent = `Genoarch engine ${GE.ENGINE_VERSION}, HPO ${meta.ontology_version}, annotations ${ann}`;
      applyDefaults();
      const ex = EXAMPLES[0];
      const hero = state.engine.analyze(ex.present, ex.absent, {});
      $("#hero-arch").innerHTML = archSVG(hero);
      loadExample("marfan", true);
      await run();
      if (!$("#method-body").dataset.done && !$('[data-view="method"]').hidden) renderMethod();
      detectServer();
    } catch (e) {
      setStatus(e.message || String(e), true);
      $("#home-data-line").textContent = e.message || String(e);
    }
  }
  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", boot); else boot();
})();
