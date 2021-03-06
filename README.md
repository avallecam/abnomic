
<!-- README.md is generated from README.Rmd. Please edit that file -->

# `abnomic`: antibody profile immunomics project

<!-- badges: start -->
<!-- badges: end -->

The goal of `abnomic` is to provide a reproducible workflow for the
undergraduate thesis project: “Comparación de la respuesta de
anticuerpos ante *Plasmodium vivax* en pacientes de la Amazonía Peruana
según su severidad y episodios previos mediante microarreglos de
proteínas”.

## Outputs

-   take a look at the approved
    [project](https://github.com/avallecam/abnomic/blob/main/01-project.pdf)

## To-Do

-   ( ) create a `renv` environment

<!-- - ( ) [MUST] Reproducibility array - request to ADi -->
<!-- - ( ) [MUST] Peptides list - request to ADi -->

## Reproducibility

Dataset required:

-   for `03-sevrcov.Rmd`
    -   `data-raw/malaria severa 28Nov2014_Ed.dta`
    -   `data-raw/severe_malaria_vivax SEP2016.dta`
    -   `data-raw/samples.csv`
-   for `04-abnomic.Rmd`
    -   `data-raw/ADi-NAMRU6_Data-samples.csv`
    -   `data-raw/samples.csv`
    -   `data-raw/PfPv500_NonpresentReactiveTargets.xlsx`
-   for `04-abnomic.Rmd` generated after PlasmoDB strategy
    -   `data/04-listgen.tsv`
    -   `data/04-listall.tsv`
    -   `data/04-listgit.tsv`
    -   `data/00-protn6c.xlsx`
    -   `data/00-protn6c.xlsx`

### plasmodb search strategy

**whole genome**

-   link:
    <https://plasmodb.org/plasmo/app/workspace/strategies/import/1247ea53c4a464a3>  
-   local route: data/20211013-02\_whole\_31ago2017-pvx\_genome.txt
-   characteristics:
    -   5665 genes
-   variables:
    -   Gene ID
    -   Transcript ID
    -   Product description
    -   Gene name or symbol

**all genes by ID in microarray chip**

-   link:
    <https://plasmodb.org/plasmo/app/workspace/strategies/import/b32a1aefa2a27bbd>  
-   local route:
    data/20211013-01\_allgen\_20ago2017-microarray\_chip.txt
-   characteristics:
    -   873 genes
    -   880 transcripts
-   variables:
    -   Same as previous, plus
    -   Ortholog count
    -   Gene variability (Total SNPs, NonSyn/Syn SNP ratio)
    -   Protein extracellular localization and secretory pathways (\# TM
        domains, SignalP scores)
    -   Gene ontology (components, functions, processes)
    -   Gene characteristics (\# exons and transcripts, length of
        transcript, CDS and protein)

<!-- ## Related projects -->
<!-- - paper: epidemiological analysis of original case-control study -->
<!-- - paper: severe pv malaria ab response, only pv | combined pv + pf analysis -->
<!-- - ~~paper: previous exposure ab response (self-report with lot of biases)~~ -->
<!-- - paper: drug treatment related ab response -->
