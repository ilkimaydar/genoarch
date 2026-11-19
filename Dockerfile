# Genoarch: one image serving the FastAPI API and the static web app (web/).
# The compact index and calibration in web/data/ are committed, so the build needs no data download.
FROM python:3.12-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PORT=7860 \
    PYTHONPATH=/app/core:/app/api \
    GENOARCH_INDEX=/app/web/data/index.json.gz \
    GENOARCH_CALIBRATION=/app/web/data/calibration.json \
    GENOARCH_WEB_DIR=/app/web

WORKDIR /app
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt
COPY core ./core
COPY api ./api
COPY web ./web
RUN useradd --create-home --uid 1000 genoarch
USER genoarch

EXPOSE 7860
HEALTHCHECK --interval=30s --timeout=5s --start-period=20s \
  CMD python -c "import os, urllib.request; urllib.request.urlopen('http://127.0.0.1:%s/health' % os.environ.get('PORT', '7860'))" || exit 1
CMD ["sh", "-c", "exec uvicorn genoarch_api.main:app --host 0.0.0.0 --port ${PORT} --proxy-headers --forwarded-allow-ips='*'"]
