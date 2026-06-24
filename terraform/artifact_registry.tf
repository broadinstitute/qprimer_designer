###############################################################################
# Artifact Registry repository that holds the qprimer-designer CPU image.      #
# The CI workflow (.github/workflows/docker.yml, job `build-gar`) pushes here. #
# The GPU multi-arch image lives in GHCR, not here.                            #
###############################################################################

resource "google_artifact_registry_repository" "qprimer_designer" {
  project       = var.project_id
  location      = var.region
  repository_id = var.gar_repository_id
  description   = "CPU-only container image for the qprimer-designer Streamlit web app (Cloud Run)."
  format        = "DOCKER"

  cleanup_policy_dry_run = false

  cleanup_policies {
    id     = "keep-recent-versions"
    action = "KEEP"
    most_recent_versions {
      keep_count = 20
    }
  }

  cleanup_policies {
    id     = "delete-old-untagged"
    action = "DELETE"
    condition {
      tag_state  = "UNTAGGED"
      older_than = "604800s" # 7 days
    }
  }

  depends_on = [google_project_service.services]
}
