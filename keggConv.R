#!/usr/bin/env Rscript
library(KEGGREST)
library(RCurl)
args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error

get_kegg <- function(uniprotID){
  ko.id <- head(keggConv("genes", paste0("uniprot:",uniprotID)))
  if (length(ko.id) == 0) {
    ko="unidentified"
    newline <- as.data.frame(unname(cbind(uniprotID,ko,"NA","unidentified","unidentified")))
  }else{
    ko <- unname(ko.id)
    pwy.id <- keggLink("pathway",ko)
    if (length(pwy.id) == 0) {
      newline <- as.data.frame(unname(cbind(uniprotID,ko,"NA","unidentified","unidentified")))
    }else{
      pwy.name <- getURL(paste0("http://togows.dbcls.jp/entry/pathway/",pwy.id,"/name"))
      if (length(pwy.name) == 0) {
        newline <- as.data.frame(unname(cbind(uniprotID,ko,pwy.id,"unidentified","unidentified")))
      }else{
        pwy.class <- getURL(paste0("http://togows.dbcls.jp/entry/pathway/",pwy.id,"/classes"))
        if (length(pwy.class) == 0) {
          newline <- as.data.frame(unname(cbind(uniprotID,ko,pwy.id,pwy.name,"unidentified")))
        }else{
          newline <- as.data.frame(unname(cbind(uniprotID,ko,pwy.id,pwy.name,pwy.class)))
        }
      }
    }
  }
  colnames(newline) <- c('uptID','ko','pwyID','pwyName','pwyCls')
  return(newline)
}

result=get_kegg(args[1])
trans <- t(cbind(result$ko,result$pwyName,result$pwyCls))
# cat(result$ko,result$pwyName,result$pwyCls)
cat(trans,sep = '|')
