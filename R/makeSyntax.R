makeSyntax = function(dataDir,
                      Hclass,
                      syntax,
                      maxClassSplit1,
                      nKeepVariables,
                      namesKeepVariables,
                      CC = 0) {

syntax[grep("infile", syntax)] = utils::capture.output(cat(paste0("infile \'", dataDir, "'")))
  syntax[grep("write", syntax)] =  paste0(
    "write = 'H", Hclass, "c", CC, "_sol", 1, ".csv'
    writeestimatedvalues='ev", CC, "_sol", 1, ".txt'
    outfile 'H", Hclass, "c", CC, "_sol", 1, ".txt' classification ",
    ifelse(nKeepVariables>0,
           paste0("keep ", paste0(namesKeepVariables, collapse = ", ")),
           ""),";"
  )

  # lm is the number of lines of one model.
  lm = which(syntax == "end model")[1] - which(syntax == "model")[1]
  syntaxModel = syntax[which(syntax == "model") :which(syntax == "end model")]
  newSyntax = syntax

  for(nClassSplit in 2:maxClassSplit1){
    newSyntax = c(newSyntax, syntaxModel)
    newSyntax[grep("write", newSyntax)][nClassSplit] =  paste0(
      "write = 'H", Hclass, "c", CC, "_sol", nClassSplit, ".csv'
      writeestimatedvalues='ev", CC, "_sol", nClassSplit, ".txt'
      outfile 'H", Hclass, "c", CC, "_sol", nClassSplit, ".txt' classification ",
      ifelse(nKeepVariables>0,
             paste0("keep ", paste0(namesKeepVariables, collapse = ", ")),
             ""),";"
    )
    newSyntax[grep("Cluster nominal", newSyntax)][nClassSplit] = paste0("      Cluster nominal ", nClassSplit, ";")
  }

  utils::write.table(newSyntax, paste0("LCT", CC, ".lgs"), row.names = FALSE, quote = FALSE, col.names = FALSE)
}
