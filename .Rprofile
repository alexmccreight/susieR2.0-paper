# Auto-load shared path + aesthetic helpers when R starts at the repository root.
# Figure scripts also `source("R/paths.R")` defensively, so this is a convenience
# rather than a requirement.
local({
  if (file.exists("R/paths.R")) {
    try(source("R/paths.R"), silent = TRUE)
    try(source("R/aesthetics.R"), silent = TRUE)
  }
})
