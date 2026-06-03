indir <- "/Users/zhuchendi/Desktop/TB234/pNS/metadata/"

files <- list.files(
  indir,
  pattern = "\\.csv$",
  full.names = TRUE
)

for (f in files) {
  message("Processing: ", basename(f))
  
  df <- read.csv(f, check.names = FALSE)
  
  required <- c(
    "pNS",
    "OBSERVED_SYN",
    "OBSERVED_NSY",
    "EXPECTED_SYN",
    "EXPECTED_NSY"
  )
  
  if (!all(required %in% colnames(df))) {
    warning("Skip file due to missing columns: ", basename(f))
    next
  }
  
  df$expected_nonsyn <- df$OBSERVED_SYN *
    (df$EXPECTED_NSY / df$EXPECTED_SYN)
  
  df$P <- mapply(
    function(obs, exp) {
      if (is.na(obs) || is.na(exp) || exp <= 0) {
        return(NA)
      }
      ppois(
        obs - 1,
        lambda = exp,
        lower.tail = FALSE
      )
    },
    df$OBSERVED_NSY,
    df$expected_nonsyn
  )
  
  df$FDR <- p.adjust(df$P, method = "BH")
  
  df$Positive_Selection <- ifelse(
    df$pNS > 1 & df$FDR < 0.05,
    TRUE,
    FALSE
  )
  
  df$expected_nonsyn <- NULL
  
  outfile <- file.path(
    indir,
    paste0(
      tools::file_path_sans_ext(basename(f)),
      "_with_stats.csv"
    )
  )
  
  write.csv(df, outfile, row.names = FALSE)
}

