variable "project_id" {
  description = "GCP project that owns the Cloud Run services and Artifact Registry repo."
  type        = string
  default     = "sabeti-adapt"
}

variable "region" {
  description = "Region for the Cloud Run services and the (regional) Artifact Registry repo."
  type        = string
  default     = "us-central1"
}

variable "service_name" {
  description = "Cloud Run service name and the user-visible app slug."
  type        = string
  default     = "qprimer-designer"
}

variable "image" {
  description = "Fully-qualified container image to deploy. The CI workflow pushes a :<sha> tag; for steady-state deploys we point at :latest. This is the CPU-only GAR image (Cloud Run has no GPU)."
  type        = string
  default     = "us-central1-docker.pkg.dev/sabeti-adapt/qprimer-designer/qprimer-designer:latest"
}

variable "gar_repository_id" {
  description = "Artifact Registry repository ID for the qprimer-designer Docker image."
  type        = string
  default     = "qprimer-designer"
}

variable "cloud_run_min_instances" {
  description = "Minimum Cloud Run instances. 0 = scale to zero (accept cold starts for this low-frequency app)."
  type        = number
  default     = 0
}

variable "cloud_run_max_instances" {
  description = "Maximum Cloud Run instances. Caps blast radius of misuse."
  type        = number
  default     = 5
}

variable "cloud_run_cpu" {
  description = "CPU per Cloud Run instance. Heavier than carmen because the primer-design pipeline runs bowtie2 + MAFFT + torch in-process."
  type        = string
  default     = "4"
}

variable "cloud_run_memory" {
  description = "Memory per Cloud Run instance. Validate against a real design run for OOM headroom."
  type        = string
  default     = "4Gi"
}

variable "cloud_run_timeout_seconds" {
  description = "Per-request timeout. Typical design runs ~2 min after recent optimizations; 600s gives headroom."
  type        = number
  default     = 600
}

variable "cloud_run_concurrency" {
  description = "Max concurrent requests per container. Must be >1: Streamlit's persistent websocket plus file-upload XHRs are separate Cloud Run requests, and routing the upload to a fresh sessionless instance returns 400. Kept modest so a single instance runs only a few heavy pipeline jobs; max_instances absorbs additional load. Session affinity pins a browser to one instance."
  type        = number
  default     = 8
}
