###############################################################################
# Cloud Run service hosting the Streamlit web app (production)                  #
###############################################################################

resource "google_cloud_run_v2_service" "app" {
  name                = var.service_name
  project             = var.project_id
  location            = var.region
  ingress             = "INGRESS_TRAFFIC_ALL"
  deletion_protection = false

  template {
    service_account                  = google_service_account.runtime.email
    timeout                          = "${var.cloud_run_timeout_seconds}s"
    max_instance_request_concurrency = var.cloud_run_concurrency
    # Streamlit's session_state lives in the Python process; pin a browser to
    # one instance so file-upload XHRs and the websocket land together.
    session_affinity = true

    scaling {
      min_instance_count = var.cloud_run_min_instances
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
        # CPU always allocated (not throttled between requests): the pipeline
        # runs snakemake as a subprocess whose work must keep progressing while
        # the Streamlit request that launched it is in flight. (carmen can idle-
        # throttle because it barely computes; this app cannot.)
        cpu_idle          = false
        startup_cpu_boost = true
      }
    }
  }

  # CI deploys new images via `gcloud run deploy --image ...`; terraform owns the
  # service shell + non-image config. The image tag drifts from this default
  # between releases — that's expected and harmless.
  lifecycle {
    ignore_changes = [
      client,
      client_version,
      template[0].containers[0].image,
    ]
  }

  depends_on = [google_project_service.services]
}

# Public site: ingress=ALL above lets anyone reach the service; this binding
# lets anyone invoke. There is no auth tier — the team accepts a public,
# shared-tenant model (results are not sensitive and are not retained durably).
resource "google_cloud_run_v2_service_iam_member" "public_invoker" {
  project  = google_cloud_run_v2_service.app.project
  location = google_cloud_run_v2_service.app.location
  name     = google_cloud_run_v2_service.app.name
  role     = "roles/run.invoker"
  member   = "allUsers"
}

###############################################################################
# Custom domain mapping (Cloud Run v1 domain-mapping API, works with v2       #
# services). GCP auto-provisions and renews the TLS cert once the DNS records  #
# below are in place.                                                          #
#                                                                              #
# Prerequisite: the GCP identity running terraform must be a verified owner    #
# of sabeti.broadinstitute.org in Google Search Console. One-time step per     #
# domain, not per service.                                                      #
###############################################################################

resource "google_cloud_run_domain_mapping" "custom_domain" {
  count    = var.custom_domain != "" ? 1 : 0
  name     = var.custom_domain
  location = var.region
  project  = var.project_id

  metadata {
    namespace = var.project_id
  }

  spec {
    route_name = google_cloud_run_v2_service.app.name
  }

  depends_on = [google_cloud_run_v2_service.app]
}
