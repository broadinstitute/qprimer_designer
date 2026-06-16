output "cloud_run_service_name" {
  value = google_cloud_run_v2_service.app.name
}

output "cloud_run_url" {
  description = "Direct *.run.app URL of the production service (reachable since ingress=ALL)."
  value       = google_cloud_run_v2_service.app.uri
}

output "staging_service_name" {
  description = "Cloud Run staging service. CI deploys per-branch tagged revisions to it."
  value       = google_cloud_run_v2_service.staging.name
}

output "staging_base_url" {
  description = "Base *.run.app URL of the staging service. Per-branch URLs prepend `<branch>---`."
  value       = google_cloud_run_v2_service.staging.uri
}

output "artifact_registry_repo" {
  description = "Artifact Registry repo the CI `build-gar` job pushes the CPU image to."
  value       = "${var.region}-docker.pkg.dev/${var.project_id}/${var.gar_repository_id}"
}

output "runtime_service_account" {
  value = google_service_account.runtime.email
}
