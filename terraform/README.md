# Terraform — qprimer-designer serving infra

Provisions the public web entry point for the qprimer-designer Streamlit app on Cloud Run:

```
Internet
   │ HTTPS (443) — Google-managed cert on the *.run.app URL
   ▼
Cloud Run service  (qprimer-designer, ingress=ALL, allUsers run.invoker)
   └─ CPU-only container image from Artifact Registry
      us-central1-docker.pkg.dev/sabeti-adapt/qprimer-designer/qprimer-designer
```

The site is **public** — no authentication. The team accepts a shared-tenant model;
results are not sensitive and (in this iteration) are **not retained durably** — each
Cloud Run instance keeps results on local disk only, which is lost on redeploy / scale
events. See "Deferred" below for adding GCS persistence later.

This stack creates:

- The **Artifact Registry** repo the CI `build-gar` job pushes the CPU image to (with
  cleanup policies: keep 20 recent versions; delete untagged after 7 days).
- A **runtime service account** (`qprimer-designer-run`) the Cloud Run services run as.
- The **production** service (`qprimer-designer`) and a **staging** service
  (`qprimer-designer-staging`), both public, with `cpu_idle=false`, session affinity,
  startup CPU boost, and scale-to-zero (`min_instances=0`).

It does **not** create the GPU image (that lives in GHCR), any GCS bucket, or a custom
domain mapping (see Deferred).

## Image split (why GAR is CPU-only)

Two images are built from one Dockerfile via the `TORCH_VARIANT` build arg:

- **GHCR** — multi-arch (amd64+arm64), **GPU**-enabled. Used by training and the
  Terra/batch CLI.
- **GAR** — amd64-only, **CPU**-only (slim, ~1.5–2 GB vs ~4.5 GB). This is what Cloud Run
  runs — Cloud Run has no GPU, so CUDA libraries would only bloat the image and slow cold
  starts.

## Prerequisites

1. **Deploy identity (already exists, shared).** CI authenticates via Workload Identity
   Federation as `gha-deployer@sabeti-adapt.iam.gserviceaccount.com`, which already holds
   `roles/run.admin`, `roles/artifactregistry.writer`, and `roles/iam.serviceAccountUser`
   at the project level — so it can deploy these services and act as the runtime SA with
   **no extra IAM**.
2. **GitHub repo secrets** on `broadinstitute/qprimer_designer` (same values carmen uses):
   - `GCP_WIF_PROVIDER` — the Workload Identity provider resource name.
   - `GCP_DEPLOY_SA` — `gha-deployer@sabeti-adapt.iam.gserviceaccount.com`.
3. **gcloud auth** for running terraform locally: `gcloud auth application-default login`.

## Apply

`terraform.tfvars` can be empty — all variables have sensible defaults.

```bash
cd terraform
terraform init
terraform plan
terraform apply
```

State is local and gitignored (see `.gitignore`); consider a GCS backend if more than one
person manages this. **Bootstrap order:** apply terraform first (creates the GAR repo +
service shells), *then* push a branch so CI can publish the CPU image and deploy revisions.

## Deferred (future changes)

- **GCS result persistence:** create a dedicated bucket with a 1-day lifecycle rule, mount
  it as a gen2 GCS volume on both services, grant the runtime SA `roles/storage.objectAdmin`,
  and set `QPRIMER_DATA_DIR` to the mount. The app already roots its data dirs at
  `QPRIMER_DATA_DIR`, so this is config-only on the app side. Makes "Past Results" survive
  redeploys and be shared across users.
- **Custom domain** (`qprimer-designer.sabeti.broadinstitute.org`): add a
  `google_cloud_run_domain_mapping`, verify the domain in `sabeti-adapt`, and add a CNAME to
  `ghs.googlehosted.com` in the `sabeti.broadinstitute.org` Cloud DNS zone.
