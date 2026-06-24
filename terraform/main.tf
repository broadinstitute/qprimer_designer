###############################################################################
# Required APIs                                                                #
###############################################################################

resource "google_project_service" "services" {
  for_each = toset([
    "run.googleapis.com",
    "artifactregistry.googleapis.com",
    "iamcredentials.googleapis.com",
  ])
  project            = var.project_id
  service            = each.value
  disable_on_destroy = false
}

###############################################################################
# Service account that the Cloud Run services run as                          #
#                                                                             #
# Distinct from the CI deploy identity (gha-deployer@sabeti-adapt), which is  #
# shared, pre-existing project infra and already holds run.admin +           #
# artifactregistry.writer + iam.serviceAccountUser at the project level — so  #
# it can deploy these services and act as this runtime SA with no extra IAM.  #
###############################################################################

resource "google_service_account" "runtime" {
  account_id   = "${var.service_name}-run"
  display_name = "Runtime SA for ${var.service_name} Cloud Run service"
  project      = var.project_id
}
