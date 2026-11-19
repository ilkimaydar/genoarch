// Run the browser engine (web/js/engine.js) under Node on fixed patients.
// Usage: node scripts/parity_check.mjs <index.json.gz> <calibration.json> <patients.json> <out.json>
import { readFileSync, writeFileSync } from "node:fs";
import { gunzipSync } from "node:zlib";
import { createRequire } from "node:module";
import { dirname, join } from "node:path";
import { fileURLToPath } from "node:url";

const require = createRequire(import.meta.url);
const here = dirname(fileURLToPath(import.meta.url));
const { Engine } = require(join(here, "..", "web", "js", "engine.js"));
const [indexPath, calPath, patientsPath, outPath] = process.argv.slice(2);

let raw = readFileSync(indexPath);
if (raw[0] === 0x1f && raw[1] === 0x8b) raw = gunzipSync(raw);
const t0 = Date.now();
const engine = new Engine(JSON.parse(raw.toString("utf8")), calPath && calPath !== "-" ? JSON.parse(readFileSync(calPath, "utf8")) : {});
const loadMs = Date.now() - t0;
const spec = JSON.parse(readFileSync(patientsPath, "utf8"));
const patients = Array.isArray(spec) ? spec : spec.patients;
const queries = Array.isArray(spec) ? [] : spec.search || [];
const results = [], timings = [];
for (const p of patients) {
  const t = Date.now();
  let res;
  try { res = engine.analyze(p.present, p.absent, p.options || {}); } catch (e) { res = { error: String(e.message || e) }; }
  timings.push(Date.now() - t);
  results.push({ name: p.name, result: res });
}
const search = queries.map((q) => ({ query: q, result: engine.search(q, 10) }));
writeFileSync(outPath, JSON.stringify({ engine: "js", node: process.version, load_ms: loadMs, analyze_ms: timings, results, search }));
console.log(`js parity run: ${results.length} patients, load ${loadMs} ms, analyze median ${timings.slice().sort((a, b) => a - b)[Math.floor(timings.length / 2)]} ms`);
