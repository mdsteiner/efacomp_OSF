library(EFAtools)

# create latex tables
names_package <- gsub(pattern = "_", replacement = "\\_",
                      names(population_models$loadings)[1:27], fixed = TRUE)
names_package_lab <- names(population_models$loadings)[1:27]
names_paper <- c("18|3|6", "6|3|6", "9|3|6", "12|3|6", "15|3|6", "18|3|3",
                 "18|3|9", "18|3|369b", "18|3|369w", "18|3|46|1c", "18|3|46|3c",
                 "12|3m|6", "18|3|6n", "6|3|369wb", "9|3|369wb", "12|3|369wb",
                 "15|3|369wb", "12|6|6", "18|6|6", "24|6|6", "30|6|6", "36|6|6",
                 "12|6|369wb", "18|6|369wb", "24|6|369wb", "30|6|369wb", "36|6|369w")
names_cors <- stringr::str_to_title(names(population_models$phis_3)) 



sink("output/tables_pattern_matrices.txt")

# tables pattern matrices
for (i in seq_len(27)) {
  temp_patt <- population_models$loadings[[i]]
  colwidth <- 12 / ncol(temp_patt)
  
  # create subsection
  cat("\\subsection{Case", names_paper[i] ,"}\n\n")
  
  # table header
  cat("\\begin{table}[htb!]\n\\caption{Population Pattern Matrix of Case",
      names_paper[i], "(", names_package[i], 
      "in the EFAtools R Package).}\n\\label{", paste0("tab:", names_package_lab[i]),
      "}\n\\begin{tabular}{", rep(paste0("C{", colwidth, "cm}"), ncol(temp_patt)),
      "} \\hline\n", paste(paste0("F", 1:ncol(temp_patt)), collapse = " & "),
      "\\\\ \\hline\n")
  
  # values
  cat(paste(apply(temp_patt, 1, function(x) paste(round(x, 3), collapse = " & ")),
            collapse = "\\\\\n"))

  # table footer
  cat(" \\\\ \\hline\n\\multicolumn{", ncol(temp_patt), "}{p{14.5cm}}{\\emph{Note:} Population pattern matrix used in the simulation analyses. To create the population correlation matrices, this pattern matrix was combined with the different population factor-intercorrelation matrices as specified in equation 1 in the main article.}\n\\end{tabular}\n\\end{table}")
  
  # add float barrier
  cat("\n\n\\FloatBarrier\n\n")
}

cat("\\section{Factor Intercorrelations}\n\n")

# tables factor intercorrelations (three factors)
for (i in seq_len(4)) {
  temp_patt <- population_models$phis_3[[i]]
  colwidth <- 12 / (ncol(temp_patt) + 1)
  
  # create subsection
  cat("\\subsection{", names_cors[i] ," Correlations (Three Factors)}\n\n")
  
  # table header
  cat("\\begin{table}[htb!]\n\\caption{Population Factor-Intercorrelation Matrix With",
      names_cors[i], "Correlations.}\n\\label{", paste0("tab:", names_cors[i],
      "3fac"), "}\n\\begin{tabular}{", rep(paste0("C{", colwidth, "cm}"), ncol(temp_patt) + 1),
      "} \\hline\n", paste(c(" ", paste0("F", 1:ncol(temp_patt))), collapse = " & "),
      "\\\\ \\hline\n")
  
  # values
  cat(paste(paste(rownames(temp_patt), "&"),
            apply(temp_patt, 1, function(x) paste(x, collapse = " & ")),
            collapse = "\\\\\n"))
  
  # table footer
  cat(" \\\\ \\hline\n\\multicolumn{", ncol(temp_patt)+1, "}{p{14.5cm}}{\\emph{Note:} Population factor-intercorrelation matrix used in the simulation analyses. To create the population correlation matrices, this factor-intercorrelation matrix was combined with the different population pattern matrices as specified in equation 1 in the main article.}\n\\end{tabular}\n\\end{table}")
  
  # add float barrier
  cat("\n\n\\FloatBarrier\n\n")
}

# tables factor intercorrelations (three factors)
for (i in seq_len(4)) {
  temp_patt <- population_models$phis_6[[i]]
  colwidth <- 12 / (ncol(temp_patt) + 1)
  
  # create subsection
  cat("\\subsection{", names_cors[i] ," Correlations (Six Factors)}\n\n")
  
  # table header
  cat("\\begin{table}[htb!]\n\\caption{Population Factor-Intercorrelation Matrix With",
      names_cors[i], "Correlations.}\n\\label{", paste0("tab:", names_cors[i],
                                                        "6fac"), "}\n\\begin{tabular}{", rep(paste0("C{", colwidth, "cm}"), ncol(temp_patt) + 1),
      "} \\hline\n", paste(c(" ", paste0("F", 1:ncol(temp_patt))), collapse = " & "),
      "\\\\ \\hline\n")
  
  # values
  cat(paste(paste(rownames(temp_patt), "&"),
            apply(temp_patt, 1, function(x) paste(round(x, 3), collapse = " & ")),
            collapse = "\\\\\n"))
  
  # table footer
  cat(" \\\\ \\hline\n\\multicolumn{", ncol(temp_patt)+1, "}{p{14.5cm}}{\\emph{Note:} Population factor-intercorrelation matrix used in the simulation analyses. To create the population correlation matrices, this factor-intercorrelation matrix was combined with the different population pattern matrices as specified in equation 1 in the main article.}\n\\end{tabular}\n\\end{table}")
  
  # add float barrier
  cat("\n\n\\FloatBarrier\n\n")
}

sink()

