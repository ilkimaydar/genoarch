# Deploying Genoarch

Three targets share one build. The static build needs no server; the Docker image adds the API.

| Target | What runs | Build |
| --- | --- | --- |
| Single file | `dist/genoarch.html` with the index embedded | `make static` |
| Static host (GitHub Pages, Hugging Face static Space, Cloudflare Pages) | `dist/site/` | `make static` |
| Container (Hugging Face Docker Space, any container host) | FastAPI + `web/` | `docker build .` |

## GitHub Pages

1. Repository Settings, Pages: set Source to "GitHub Actions".
2. Push to `main` (or run the workflow by hand). `.github/workflows/pages.yml` runs `scripts/build_static.py`
   on the committed index and calibration and deploys `dist/site`.
3. The site appears at `https://ilkimaydar.github.io/genoarch/`.

## Hugging Face Space `ilkimaydar/genoarch`

Option A, from GitHub Actions: add a repository secret `HF_TOKEN` (a Hugging Face token with write access),
then run the workflow "Publish static build to the Hugging Face Space". It creates the Space with
`sdk: static` if needed and uploads `dist/hf-space`.

Option B, by hand: create a Space named `genoarch` with the Static SDK, then upload the contents of
`dist/hf-space/` (its `README.md` carries the Space configuration). The Space serves `index.html`.

Option C, container: create a Docker Space and push this repository; the `Dockerfile` listens on port 7860,
which is the Spaces default.

## Docker

```sh
docker build -t genoarch:0.1.0 .
docker run --rm -p 8000:7860 genoarch:0.1.0     # http://127.0.0.1:8000, docs at /docs
```

Environment: `PORT` (default 7860), `GENOARCH_INDEX`, `GENOARCH_CALIBRATION`, `GENOARCH_WEB_DIR`.

## Cloudflare in front of genoarch.ai

Prerequisite: register `genoarch.ai` and add it to a Cloudflare account (free plan is enough), then change the
registrar's nameservers to the two Cloudflare nameservers shown in the dashboard.

### DNS

Static site on GitHub Pages:

| Type | Name | Content | Proxy |
| --- | --- | --- | --- |
| CNAME | `www` | `ilkimaydar.github.io` | Proxied |
| CNAME | `@` (flattened at the apex) | `ilkimaydar.github.io` | Proxied |

Then add `genoarch.ai` as the custom domain in the repository's Pages settings (this writes a `CNAME` file into
the Pages deployment) and enable "Enforce HTTPS" once the certificate is issued. GitHub also documents apex A
records; check its current list in the GitHub Pages documentation before using them (not verified here).

API on a container host (for example `api.genoarch.ai`): a proxied CNAME from `api` to the host name the
platform gives you. Hugging Face Space URLs (`https://ilkimaydar-genoarch.hf.space`) can be fronted the same way
only if the platform accepts the custom host name; check Hugging Face's current custom-domain support before
relying on it (not verified here).

### SSL / TLS

* Encryption mode: Full (strict). Both GitHub Pages and Hugging Face serve valid certificates.
* Edge Certificates: Always Use HTTPS on, Automatic HTTPS Rewrites on, minimum TLS 1.2, TLS 1.3 on.
* Enable HSTS only after HTTPS works on every subdomain you use.

### Caching rules

| Path | Rule | Reason |
| --- | --- | --- |
| `/data/index.json.gz` | Cache, edge TTL 1 day, browser TTL 1 hour | 1.45 MB, changes only with a new HPO release |
| `/js/*`, `/css/*` | Cache, edge TTL 1 day | change only on deploy; purge on deploy |
| `/` and `/index.html` | Cache with a short TTL (5 minutes) or respect origin headers | picks up new deploys quickly |
| `/analyze`, `/rank`, `/next-phenotypes`, `/hypotheses`, `/hpo/search`, `/health` | Bypass cache | per-request API responses |

Make sure Cloudflare does not re-compress or transform `index.json.gz`: the page reads the gzip bytes itself and
also accepts an already-decompressed body, so either behaviour works, but byte-for-byte caching is cheapest.
After each deploy, purge the cache for `/index.html`, `/js/*` and `/css/*` (or version the file names).

### Security

* WAF managed rules on; a rate-limiting rule on the API paths (for example 60 requests per minute per IP).
* No cookies or accounts exist in this version; the static app sends no phenotype data to any server.

## Release checklist

1. `make data index calibrate parity static test lint`
2. Check `dist/eval/parity/parity_report.json` (`all_pass: true`) and the evaluation summary.
3. Commit `web/data/{index.json.gz,calibration.json,method_card.json}` together with the code.
4. Push to `main`; Pages deploys automatically; run the Hugging Face workflow.
