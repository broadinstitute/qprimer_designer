###############################################################################
# Staging Cloud Run service                                                   #
#                                                                             #
# A second, deliberately permissive service used by CI to deploy a tagged     #
# revision per branch push, so engineers can test changes against the real    #
# container before merging.                                                   #
#                                                                             #
# Per-branch revisions are created by .github/workflows/docker.yml            #
# (`staging-deploy`) using `gcloud run deploy --tag <branch> --no-traffic`,   #
# which yields URLs of the form                                               #
#   https://<branch>---<service>-<hash>.<region>.run.app                      #
# distinct per branch, with no traffic shifting on the base URL.              #
#                                                                             #
# CI is the source of truth for staging revisions AND their resource config   #
# (the gcloud flags in staging-deploy). terraform only owns the service shell #
# + IAM, hence lifecycle.ignore_changes on the whole template.                #
###############################################################################

resource "google_cloud_run_v2_service" "staging" {
  name                = "${var.service_name}-staging"
  project             = var.project_id
  location            = var.region
  ingress             = "INGRESS_TRAFFIC_ALL"
  deletion_protection = false

  template {
    service_account                  = google_service_account.runtime.email
    timeout                          = "${var.cloud_run_timeout_seconds}s"
    max_instance_request_concurrency = var.cloud_run_concurrency
    session_affinity                 = true

    scaling {
      min_instance_count = 0
      max_instance_count = var.cloud_run_max_instances
    }

    containers {
      image = var.image

      ports {
        container_port = 8080
      }

      resources {
        limits = {
          cpu    = var.cloud_run_cpu
          memory = var.cloud_run_memory
        }
        cpu_idle          = false
        startup_cpu_boost = true
      }
    }
  }

  # CI is the source of truth for staging revisions; terraform only bootstraps
  # the service shell. New revisions arrive via `gcloud run deploy --tag ...`.
  lifecycle {
    ignore_changes = [
      template,
      client,
      client_version,
    ]
  }

  depends_on = [google_project_service.services]
}

resource "google_cloud_run_v2_service_iam_member" "staging_public" {
  project  = google_cloud_run_v2_service.staging.project
  location = google_cloud_run_v2_service.staging.location
  name     = google_cloud_run_v2_service.staging.name
  role     = "roles/run.invoker"
  member   = "allUsers"
}
