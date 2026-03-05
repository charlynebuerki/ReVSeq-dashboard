# ReVSeq Dashboard

Django dashboard for ReVSeq data visualization, with Auspice/Nextstrain trees.

## Recent Changes

- Refactored dashboard tree embedding to use a configurable base URL (`NEXTSTRAIN_BASE_URL`).
- Added deployment-safe Django host and iframe settings:
  - `DJANGO_ALLOWED_HOSTS`
  - `DJANGO_DEBUG`
  - `X_FRAME_OPTIONS = "SAMEORIGIN"`
- Added dedicated deploy compose file: `docker-compose.deploy.yml`.
- Standardized local vs hosted tree routing:
  - Local (no reverse proxy): `http://127.0.0.1:<port>`
  - Hosted (nginx reverse proxy): `/nextstrain`

## Local Run (No Docker)

Run Django:

```bash
export NEXTSTRAIN_BASE_URL=http://127.0.0.1:4001
python manage.py runserver
```

Run Auspice in a second terminal:

```bash
HOST=0.0.0.0 PORT=4001 npx auspice-revseq view --datasetDir dashboard/static/auspice
```

Open:
- `http://127.0.0.1:8000/home/`

## Local Run (Docker Compose)

`docker-compose.yml` is set up for local development.

```bash
docker compose up --build
```

Open:
- Dashboard: `http://localhost:8000/home/`
- Auspice: `http://localhost:4000/`

## Deployment (Docker + nginx)

### 1. Build and push versioned images

```bash
export REVSEQ_TAG=2026.03.02-6
docker build -f Dockerfile -t cbuerk/revseq-dashboard:$REVSEQ_TAG .
docker build -f Dockerfile_nextstrain -t cbuerk/revseq-nextstrain:$REVSEQ_TAG .
docker push cbuerk/revseq-dashboard:$REVSEQ_TAG
docker push cbuerk/revseq-nextstrain:$REVSEQ_TAG
```

### 2. Deploy on host

`docker-compose.deploy.yml` should include for `dashboard`:

```yaml
environment:
  - DJANGO_ALLOWED_HOSTS=revseq.charlynebuerki.com,localhost,127.0.0.1
  - DJANGO_DEBUG=0
  - NEXTSTRAIN_BASE_URL=/nextstrain
```

Then:

```bash
export REVSEQ_TAG=2026.03.02-6
docker compose -f docker-compose.deploy.yml pull
docker compose -f docker-compose.deploy.yml up -d --force-recreate
```

### 3. nginx reverse proxy (host)

Use a dedicated HTTPS server block for `revseq.charlynebuerki.com` with:

- `/` -> Django (`127.0.0.1:8000`)
- `/nextstrain/` -> Auspice (`127.0.0.1:4000`) with prefix rewrite
- `/dist/` -> Auspice static assets
- `/charon/` -> Auspice API
- explicit rewrite for dataset requests with `prefix=nextstrain/...`

Example (inside the `listen 443 ssl` block):

```nginx
location /nextstrain/ {
    rewrite ^/nextstrain/(.*)$ /$1 break;
    proxy_pass http://127.0.0.1:4000;
    proxy_set_header Host $host;
    proxy_set_header X-Forwarded-Proto $scheme;
    proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
}

location = /charon/getDataset {
    if ($arg_prefix ~ "^nextstrain/(.*)$") {
        proxy_pass http://127.0.0.1:4000/charon/getDataset?prefix=/$1;
    }
    proxy_pass http://127.0.0.1:4000/charon/getDataset?$query_string;
    proxy_set_header Host $host;
    proxy_set_header X-Forwarded-Proto $scheme;
    proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
}

location /dist/ {
    proxy_pass http://127.0.0.1:4000/dist/;
    proxy_set_header Host $host;
    proxy_set_header X-Forwarded-Proto $scheme;
    proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
}

location /charon/ {
    proxy_pass http://127.0.0.1:4000/charon/;
    proxy_set_header Host $host;
    proxy_set_header X-Forwarded-Proto $scheme;
    proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
}

location / {
    proxy_pass http://127.0.0.1:8000;
    proxy_set_header Host $host;
    proxy_set_header X-Forwarded-Proto $scheme;
    proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
}
```

After edits:

```bash
sudo nginx -t
sudo systemctl reload nginx
```

## Quick Health Checks

Hosted checks:

```bash
curl -I https://revseq.charlynebuerki.com/strain/Parainfluenza_3/
curl -I "https://revseq.charlynebuerki.com/nextstrain/HPIV-3?f_Node%20type=New"
curl -I "https://revseq.charlynebuerki.com/charon/getDataset?prefix=/HPIV-3"
curl -I "https://revseq.charlynebuerki.com/static/barplots/Parainfluenza_3_seq.html"
```

If charts appear blank, do a hard refresh (`Ctrl+Shift+R`) to clear cached iframe content.
