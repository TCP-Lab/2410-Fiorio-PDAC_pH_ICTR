# 2410-Fiorio-PDAC_pH_ICTR
##### FGFR2 isoforms in PDAC

## Background
The __fibroblast growth factor receptor 2__--___FGFR2___--exists in multiple
isoforms from alternative splicing events. In particular, it is known from
literature that
- under physiological conditions, the _isoform b_ (___FGFR2IIIb___) is typically
	expressed at epithelium level, while the _isoform c_ (___FGFR2IIIc___) is
	usually found at stromal/mesenchymal level;
- expression of the _c_ isoform at the epithelial level promotes carcinogenesis,
	 with special reference to pancreatic ductal adenocarcinoma (PDAC).

> __References__
> - [PMID 20094046](https://pubmed.ncbi.nlm.nih.gov/20094046/)
> - [PMID 22440254](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3349828/)
> - [PMID 23444225](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3644028/)

## Aim
It would be interesting to re-analyze at isoform-level the RNA-Seq data set from
our 2022 study about epithelial–to-mesenchymal transition (EMT) effects induced
by acute or prolonged exposure of PANC-1 cell line to acidic pH environment:
 
> Audero, M.M. _et al_. __Acidic Growth Conditions Promote Epithelial-to-Mesenchymal Transition to Select More Aggressive PDAC Cell Phenotypes In Vitro__. _Cancers_ 2023, 15, 2572. https://doi.org/10.3390/cancers15092572
> ([PMID 37174038](https://pmc.ncbi.nlm.nih.gov/articles/PMC10177299/))

The aim now is to investigate the expression patterns of FGFR2 isoforms in
PANC-1 cells as an _in vitro_ model for PDAC when a more aggressive phenotype
is triggered by acidosis.

## Methods
### From Reads to Counts
The $12 \times 2$ PE FASTQ files from
[PMID 37174038](https://doi.org/10.3390/cancers15092572)
went through a standard
[___x.FASTQ___](https://github.com/TCP-Lab/x.FASTQ)
pipeline for quality control, adapter and quality trimming, read alignment (by
STAR), and transcript abundance quantification (by RSEM). For the last two steps,
Genome assembly GRCh38 (_hg38_) was used, together with the
`Homo_sapiens.GRCh38.110.gtf` GTF file for annotation. _TPM_ and _expected count_
expression matrices were neventually assembled for both genes and isoforms,
using Ensembl ENSG and ENST, respectively, as primary IDs.

### Identifying the isoforms of interest
Searching for human _FGFR2_ gene through the _Ensembl genome browser_
([ENSG00000066468](https://www.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000066468;r=10:121478332-121598458) > `Summary` > `Show transcript table`)
returned 41 possible splice variants, among which 25 protein coding transcripts.
Of these, only 8 are _golden_ genes.

> [!TIP]
> __The Ensembl/Havana merge:__ Species which have both HAVANA and Ensembl gene
> annotation (i.e., Human, Mouse, Rat, and Zebrafish) undergo a merge of the two
> sets of gene models. A merged (or ___golden___) gene indicates that annotation
> was provided by both Ensembl and HAVANA. Where a transcript model is annotated
> only by Ensembl or HAVANA, it is displayed as an unmerged (or ___red___) model.
>
>__References__
> - [Gene annotation in Ensembl](https://www.ensembl.org/info/genome/genebuild/index.html)
> - [The Ensembl and Havana merge](https://www.ensembl.org/info/genome/genebuild/annotation_merge.html)
> - [The HAVANA team](https://www.sanger.ac.uk/project/manual-annotation/)

Nevertheless, it was not obvious to find which of those ENST IDs corresponded to
the two isoforms of interest identified by the two "common" names
___FGFR2IIIb___ and ___FGFR2IIIc___.

Thankfully, the __P21802 · FGFR2_HUMAN__ entry from the __UniProt__ database
features a specific section about protein
[Sequence & Isoforms](https://www.uniprot.org/uniprotkb/P21802/entry#sequences),
also reporting all their names and aliases (synonyms). Based on this information,
the correspondence between the common names of the isoforms and their official
UniProtKB/Swiss-Prot IDs was found to be as follows:
- __P21802-1__ (canonical sequence)
	Synonyms: BEK, __FGFR2IIIc__
- __P21802-3__
	Synonyms: BFR-1, __FGFR2IIIb__, KGFR

This in turn allowed the transcripts of interest to be identified in __Ensembl__
by the `UniProt Match` column (both of them being _gold_ protein coding biotypes).

| Common Name    | FGFR2IIIb         | FGFR2IIIc          |
| -------------- | ----------------- | ------------------ |
| Transcript ID  | [ENST00000457416.7](https://www.ensembl.org/Homo_sapiens/Transcript/Summary?db=core;g=ENSG00000066468;r=10:121478332-121598458;t=ENST00000457416) | [ENST00000358487.10](https://www.ensembl.org/Homo_sapiens/Transcript/Summary?db=core;g=ENSG00000066468;r=10:121478332-121598458;t=ENST00000358487)|
| Name           | FGFR2-215         | FGFR2-206          |
| Translation ID | ENSP00000410294.2 | ENSP00000351276.6  |
| CCDS           | CCDS7620          | CCDS31298          |
| UniProt Match  | P21802-3          | P21802-1           |
| RefSeq Match   | NM_022970.4       | NM_000141.5        |

> [!TIP]
> Get a schematic representations of exon-intron structure for the isoforms of
> interest following this path:
> [Ensembl `Gene` Tab](https://www.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000066468;r=10:121478332-121598458) >
> `Summary` > `Splice variants` > `Basic Gene Annotations from GENCODE 47` >
> _left click on a transcript_ > `Zoom on feature`. Alternatively, you can directly
> select the Ensembl `Location` Tab > `Region in detail` menu item. In the _Region
> Image_ pane, transcripts are drawn as boxes (exons) and lines connecting the boxes
> (introns). Filled boxes represent coding sequence and unfilled boxes (or portions
> of boxes) represent UnTranslated Regions (UTR). Here, tracks above the blue bar
> (_Contigs_) are on the forward strand of the chromosome, and tracks under the blue
> bar are on the reverse strand.

![Full transcripts](/figs/Human_10121478334_121598458_full_trimmed.png "Full transcripts")
![Zoom on swapped exons](/figs/Human_10121513723_121521615_zoom_trimmed.png "Zoom on swapped exons")


`Summary` > `Transcript comparison` > `Select transcripts` > `Download sequence` >
choose what to download among: the sequences of the mature (spliced) transcripts (as cDNA),
just the coding sequences (CDS), 5' UTRs, 3' UTRs, the list of exons, the list of introns,
the entire genomic sequences


### Kerblam! Workflows
Provided that the input data and metadata matrices are correctly named and
present in the `./data/in` folder, this
[Kerblam!](https://www.kerblam.dev/)
project runs 2 independent workflows.

```bash
kerblam run dea
```
can be used to rerun the same DEA published in
[PMID 37174038](https://doi.org/10.3390/cancers15092572),
but with the most up-to-date versions of the same software and R packages used
there (in times before Kerblam!). It may be interesting to compare the resulting
lists of DEGs and reflect on reproducibility...

```bash
kerblam run iso
```
is the workflow to be used for the isoform-level analysis of _FGFR2_ expression.

## Results
### RSEM Expected Counts
```
                          gene_level isoform_sum
4-days-pH-6_6-12-21_09_21         48       48.00
4-days-pH-6_6-16                  31       31.01
4-days-pH-6_6-7_1                 11       11.00
control-12                        36       36.01
control-13-21_09_21               27       27.00
control-7_4                       23       23.00
control-7_6                       22       22.00
pH-selected-p7_5                  94       94.00
pH-selected-p7_6                  59       59.00
pH-selected-p7_7                  53       53.00

>>> Ctrl
                    mean      SD
ENST00000358487    1.005    2.01
ENST00000457416    0.000    0.00

>>> Acute
                    mean      SD
ENST00000358487        0       0
ENST00000457416        0       0

>>> Selected
                    mean      SD
ENST00000358487     9.22  8.7795
ENST00000457416     0.00  0.0000
```
### TPMs
```
                          gene_level isoform_sum
4-days-pH-6_6-12-21_09_21       0.59        0.59
4-days-pH-6_6-16                0.60        0.60
4-days-pH-6_6-7_1               0.15        0.15
control-12                      0.30        0.30
control-13-21_09_21             1.02        1.00
control-7_4                     0.56        0.56
control-7_6                     0.31        0.31
pH-selected-p7_5                0.73        0.73
pH-selected-p7_6                0.43        0.43
pH-selected-p7_7                0.53        0.52

>>> Ctrl
                    mean      SD
ENST00000358487    0.005    0.01
ENST00000457416    0.000    0.00

>>> Acute
                    mean      SD
ENST00000358487        0       0
ENST00000457416        0       0

>>> Selected
                    mean      SD
ENST00000358487   0.0533  0.0551
ENST00000457416   0.0000  0.0000
```

## Discussion
___FGFR2IIIb___ (ENST00000457416) isoform does not appear to be expressed in any
experimental group (and this may be consistent with the fact that PANC-1 cells
are a model of epithelioid carcinoma).

___FGFR2IIIc___ (ENST00000358487) isoform would appear to be more expressed in
the acidic-pH _Selected_ group, consistent with the fact that this treatment
(pH-selection by long-term acidic pressure followed by recovery to pH 7.4)
favored EMT and correlated with a more aggressive tumor phenotype. However,
considering the standard deviations within groups, and applying even the most
"liberal" thresholds of 0.5 TPMs (or 10 counts) for lowly expressed genes, it is
hard to support the claim that the FGFR2IIIc isoform is actually expressed (or
more expressed) in pH-selected PANC-1 cells.

__Ultimately, we have just a (very) weak evidence that prolonged exposure of
PANC-1 cells to acidic pH environment induces the expression of isoform IIIc.__

The reasons for this could be that (i) the receptor is poorly expressed by
PANC-1 in general or (ii) the sequencing depth is not sufficient. However,
considering the average depth of ~70 Mreads/sample and the fact that FGFR2 gene
is detected at very low--but still significant--levels (10-100 counts; ~0.5
TPMs) in every sample, I lean toward the first hypothesis and do not think that
further increasing the sequencing depth would bring any benefit.

> __References__
> - [DESeq2 vignette](https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#pre-filtering)
> - [reddit: TPM counts cut off?](https://www.reddit.com/r/bioinformatics/comments/1dn19gn/tpm_counts_cut_off/?rdt=64484)

