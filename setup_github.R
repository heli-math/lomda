#!/usr/bin/env Rscript
# =============================================================================
# setup_github.R
# Run this script once to initialise git and push lomda to GitHub.
# =============================================================================
# Prerequisites:
#   1. Create an empty repo at https://github.com/yourusername/lomda
#   2. Install gh CLI or have your PAT configured:
#      usethis::gh_token_help()
#   3. Run: Rscript setup_github.R  (from inside the lomda/ directory)
# =============================================================================

# ---- Step 1: Install devtools / usethis if needed --------------------------
if (!requireNamespace("devtools",  quietly = TRUE)) install.packages("devtools")
if (!requireNamespace("usethis",   quietly = TRUE)) install.packages("usethis")
if (!requireNamespace("roxygen2",  quietly = TRUE)) install.packages("roxygen2")

# ---- Step 2: Generate documentation ----------------------------------------
message("Generating roxygen documentation...")
roxygen2::roxygenise()

# ---- Step 3: Check the package ---------------------------------------------
message("Running R CMD check (this may take a minute)...")
devtools::check(quiet = FALSE)

# ---- Step 4: Initialise git ------------------------------------------------
if (!file.exists(".git")) {
  system("git init")
  system("git add -A")
  system('git commit -m "Initial commit: lomda v0.1.0"')
}

system('git config --global user.name "He Li"')
system('git config --global user.email "lihestat@outlook.com"')
system("git add -A")
system('git commit -m "Initial commit: lomda v0.1.0"')
system("git log --oneline")
system("git branch -M main")
system("git push -u origin main")


# ---- Step 5: Push to GitHub ------------------------------------------------
# Replace 'yourusername' with your actual GitHub username.
gh_user <- Sys.getenv("GITHUB_USER",
                       unset = "heli-math")  # <-- change this
remote  <- paste0("https://github.com/", gh_user, "/lomda.git")

system(paste("git remote add origin", remote))
system("git branch -M main")
system("git push -u origin main")

message("
Done! Your lomda package is now on GitHub.

Next steps:
  * Visit https://github.com/", gh_user, "/lomda
  * Enable GitHub Actions for automated R-CMD-check
  * Add a pkgdown site with: usethis::use_pkgdown_github_pages()
  * Submit to CRAN when ready: devtools::release()
")
