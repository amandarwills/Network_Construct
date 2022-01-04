# Network_Construct
Scripts for constructing gene and metabolite networks.

metabolite_filter.pl: Perl script to filter aligned peak lists from El Maven based on "good peak count" and peak quality scores.

DAM.R (differentially accumulated metabolites): R script to identify metabolites that are differentially accumulated at a specific timepoint between two conditions.

getGeneCoNet.R: R script to construct gene co-expression network

getMetCoNet.R: R script to construct metabolite co-expression network

metabolite_gene_KEGG_based_network Scripts for downloading and connecting the RHEA numbers to KEGG Reation IDs and to build a KEGG reaction network that we can annotate with info about the gene or metabolite/compound expression/accumulation profiles. This should hopefully be useful for connecting timeseries gene and metabolite data. KEGG gives structure to the data and allows us to see which gene expression patterns might be feeding into a given pool or metabolites.
