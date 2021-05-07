> if (!requireNamespace("BiocManager", quietly = TRUE))
+     install.packages("BiocManager")
> BiocManager::install("biomaRt")

> library(biomaRt)
> listMarts()

> ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl")

#first table
>getBM(attributes=c("entrezgene_id","hgnc_symbol","ensembl_gene_id"), filters="mim_morbid_accession", values=c("603218", "604802", "143100", "606438", "607136"), mart=ensembl)

#second table
>getBM(attributes=c("hgnc_symbol","ensembl_gene_id", "ensembl_transcript_id"), filters="mim_morbid_accession", values=c("603218", "604802", "143100", "606438", "607136"), mart=ensembl)
