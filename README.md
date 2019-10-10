
# GWAS of glycosuria

Date of publication - 2019-10-07

## Title

Common variation on chromosome 16 is associated with glycosuria in
pregnancy: Findings from a genome-wide association study in European
women.

## Citation

……

## Summary statistics

All summary statistcs are available on the [data.bris
website]().

## Summary

![](./manuscript/figures/Supplementary%20Figure%202.%20Manhattan%20plot%20of%20a%20GWAS%20of%20self-reported%20glycosuria%20in%20the%20third%20trimester%20of%20pregnancy%20in%20ALSPAC.png)

*Manhattan plot of a GWAS of self-reported glycosuria in the third
trimester of pregnancy in ALSPAC*

We conducted a genome-wide association study (GWAS) of glycosuria (sugar
in urine) in pregnant mothers from the Avon Longitudinal Study of
Parents and Children (ALSPAC). Due to a lack of available external data
sources replication was not possible, instead we performed a GWAS in the
Northern Finland Birth Cohort 1986 (NFBC1986) where we used mothers
phenotype and the mothers’ offsprings genotype. To estimate the maternal
effects from offspring genotypes we doubled the effect estimates and
standard errors of the GWAS results (see:
[PMID 27029810](https://www.ncbi.nlm.nih.gov/pubmed/?term=27029810),
[PMID 29030599](https://www.ncbi.nlm.nih.gov/pubmed/?term=29030599),
[PMID 9778168](https://www.ncbi.nlm.nih.gov/pubmed/?term=9778168))

In order to reproduce this analysis you will need access to
[ALSPAC](http://www.bristol.ac.uk/alspac/researchers/access/) and
[NFBC1986](https://www.oulu.fi/nfbc/node/19668).

 

All scripts used for this work are in
[`scripts`](https://github.com/mattlee821/001_glycosuria_GWAS/tree/master/scripts):

  - All of the scripts use relative file paths (`./my/file/path`) where
    `./` is the directory `001_glycosuria_GWAS`
  - This work was performed using Terminal on a Mac, and scripts were
    run on the University of Bristol High Performance Computer,
    [BlueCrystal 3](https://www.acrc.bris.ac.uk/acrc/phase3.htm)

 

The manuscript and all figures, tables and supplementary information are
in
[`manuscript`](https://github.com/mattlee821/001_glycosuria_GWAS/tree/master/manuscript/).

## Review

All information from the review stage of publication are available in
[`manuscript/review`](https://github.com/mattlee821/001_glycosuria_GWAS/tree/master/manuscript/review).
This includes reviewer comments and responses to reviewer comments. The
script used to perform reviewer requested additional analysis is in
[`scripts`](https://github.com/mattlee821/001_glycosuria_GWAS/tree/master/scripts)
and labelled
[`reviewer_analysis.R`](https://github.com/mattlee821/001_glycosuria_GWAS/tree/master/scripts/reviewer_analysis.R).
