This repo contains scripts used for the publication "[Dynamic and complex transcription factor binding during an inducible response in yeast](https://www.ncbi.nlm.nih.gov/pubmed/19487574)" by Ni L et al. (2009) Genes Devel. 23(11)1351-63. Seven transcription factors involved in yeast's salt response were characterized for their binding sites at four time points after switching to a high salt medium. Expression time course data in yeast strains lacking either of two of these transcription factors were also analyzed to generate insigths on the transcriptional dynamics of a physiological response program.   

* chipchip_mednorm.R: 2D median-normalization of counts on the Nimblegen microarray. 
* chipchip_reformat.pl: Reads in two data files (for Red and Green channels) from a
     chipchip experiment. Outputs an output file with 6 cols: chr, pos, x, y, Red, Green.
* chipchip.R : Reads three replicates (2d-normalized) from a Chip-chip experiment. 
       A 7-probe moving window (375 nucleotides on the genome)
       is used to estimate the signal value of the middle probe and the probability 
       that it is significantly different from 0, using the Wilcoxon test.
* scoreruns.pl: Reads four column output from R script that analyzes chipchip data;
       Identifies probes having ratio and log p value above threshold.
       defines "runs" that are within maxgap1, then eliminates those runs
       that are less than minrun, them joins runs that are within maxgap2
       of each other. Calculates peaks, including peaks on peaks (shoulder peaks).
       Report flanking genes on the W and C strands (geneW1, geneW2, geneC1, geneC2).
* locate.pl: takes chromosomal coordinates of a ChIP chip "hit" and returns gene names, 
       common names(if available), and description, for both strands, where the hit is 
       within 1500 nucleotides upstream of the ORF start. Depending on whether the input 
       coordinates are completely within an ORF (on either strand), or not, it also returns a 
       a 1 or 0, respectively.
* load_hits.pl: insert into table ALLHITS expt name, and other data from hitlist.
* load_randomized_hitlist.pl: reads a chipchip hitlist, reads hit lengths in order,
        then generates randomized hitlist in the yeast genome having
        the same hit lengths.
* loadTFsites.pl: insert into TFSITES table chr#, start pos, end pos
        (peaks that are within 500 nt of each other are placed in the same row)
* loadTSsites.pl: insert into TSSITES table chr#, start pos, end pos, normalized expr,
        geneW1, geneW2, geneC1, geneC2
* loadyap4chipexpr.pl insert into YAP4TIMES table chr#, start pos, end pos, normalized expr.
* intersect.pl: Takes two hitlists from a chipchip expt.  Writes a
	a file with the number of overlaps for the first N hits in
	each list having the highest signals.
* intersect_clust.pl: Takes a chipchip hitlists and a clustered hitlist and 
        returns overlapping hits.
* isect_chipchip_expt.pl: Compare a sorted gene listed from a microarray experiment
        with a coordinate list from a chipchip experiment.  Converts
        the gene names to a coordinate list, then looks for overlaps
        in the chipchip hit list.
* intersect_graph.pl: Takes two hitlists from a chipchip expt.  Writes a
	a file with the number of overlaps for the first N hits in
	each list having the highest signals.
*yeast_randomized_hits: reads a chipchip hitlist, reads hit lengths in order,
       then generates random hits in the yeast genome having 
       the same hit lengths. Produces a hitlist with hit locations.

