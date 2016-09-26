# 'file' must already be defined
deps1 <- gsub("\\)","",gsub("library\\(","",readLines(file)[grep("library\\(", readLines(file))]))
deps2 <- gsub("\\)","",gsub("require\\(","",readLines(file)[grep("require\\(", readLines(file))]))
deps <- unique(c(deps1, deps2))

deps <- deps[!(deps %in% installed.packages())]
if (length(deps) > 0) {
  command <- paste0('install.packages(c("',
                    paste(deps, collapse='","'),'"))')
  eval(parse(text=command))
} else {
  cat("All R dependencies already installed!\n")
}
