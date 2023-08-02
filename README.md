# Data Harmonisation Procedures
Guidelines for merging genome-wide genotype data obtained from multi-way admixed populations.

## Steps
1. Processing individual datasets
2. Merging of datasets
3. Ancestry inference
4. Batch effect coffection

**Software required for procedures**
- PLINK
- BCFTOOLS
- KING
- SHAPEIT2
- ADMIXTURE
- PONG
- RFMix

## 01 - Processing individual datasets 
**Initial quality control**

Individual datasets must undergo independent quality control procedures prior to merging. Remove SNPs and individuals with missing information, monomorphic sites and SNPs deviating from Hardy-Weinberg equilibrium (HWE). Remove sex chromosomes and chromosome 23. Check for any phenotypic information of population under investigation, in order to remove any indivual who have any missing phenotypic information, such as age or sex or disease status.

``` 
plink --bfile filename --mind 0.1 --geno 0.05 --hwe 0.00001 -chr1-22 --make-bed --out outputfile
```

**Imputation**

Imputing separate datasets is recommended prior to merging to increase the number of overlapping SNPs between datasets. This step is crucial when merging separate datasets genotyped on disparate genotyping platforms. There are several publicly available imputation servers which offer different imputation algorithms and reference panels: 
- Sanger imputation server (https://imputation.sanger.ac.uk)
- TOPMed imputation server (https://imputation.biodatacatalyst.nhlbi.nih.gov)
- Michigan imputation server (https://imputationserver.sph.umich.edu/index.html)

The accuracy and quality of imputation performed by different imputation algorithms and reference panels on different populations is an area of ongoing research. The population structure of the tagret population and similarity to the imputation reference panel must be carefully considered before imputation. We recommend investigating the performance of different imputation algorithm and reference panel combinations to determine the best strategy for your unique target population. 

For our five-way admixed South African population, we employed the PBWT imputation algorithm and the African Genome Resource panels, made available through the Sanger imputation server. The server provides detailed instructions for preparing files for imputation (https://imputation.sanger.ac.uk/?instructions=1). 

**Post-imputation quality control** 

We recommend removing poorly imputed SNPs from individual datasets to minimise spurious associations in downstream analyses. To assess imputation quality we consider the INFO score quality metric. The INFO score ranges from 0 to 1, where values near 1 indicate that a SNP has been imputed with a high degree of certainty. BCFTOOLS can be used to calculate the INFO scores. 

```
for i in {1..22}; do bcftools +impute-info  ${i}.pbwt_reference_impute.vcf.gz  -Oz -o chr${i}_imputed_info.vcf.gz; done
```

We can choose an INFO score threshold to filter out poorly imputed sites. We recommend filtering out all SNPs with INFO scores less than 0.8. This can be done with BCFTOOLS. 

```
bcftools filter -e 'INFO/INFO<0.8' input.vcf.gz -Oz -o output.vcf.gz
```

We need to check and remove related indivduals within individual datasets before merging. The software KING and kinship coefficient is sufficient to identify close relationships between participants who have extensive population strcuture. Usually up to 2nd Degree relatedness, indicated by a kinship coefficient of 0.0884-0.177. 

```
king -b inputfile.bed --kinship --degree2 --prefix outputfile
```

Make a file with related individuals and remove. 

```
plink --bfile inputfile --remove related_individuals_file --make-bed --out outputfile
```

## 02 - Merging individual datasets

The merge function in PLINK appears to exclude variants that are common in both data sets, which is unexpected.The following steps extract the variants common in both file sets and merge the data:

1. Extract SNP IDs from BIM files and sort them numerically (it will be useful to convert SNP IDs to chr:bp format before doing this).

```
plink --bfile inputfile --set-all-var-ids @:# --make-bed --out outputfile
awk '{print $2}' inputfile.bim | sort > inputfile_SNPs_sorted.txt
```

2. Create a list of SNPs common to the datasets you wish to merge.

```
comm -12 inputfile1_sorted.txt inputfile2_sorted.txt > intersecting_snps.txt
```

3. Extract this list of SNPs from the separate datasets.

```
plink --bfile inputfile1 --extract intersecting_snps.txt --make-bed --out inputfile1_intersect_snps

plink --bfile inputfile2 --extract intersecting_snps.txt --make-bed --out inputfile2_intersect_snps
```

4. Merge the datasets.

```
plink --bfile inputfile1_intersect_snps --bmerge inputfile2_intersect_snps.bed inputfile2_intersect_snps.bim inputfile2_intersect_snps.fam --make-bed --out outputfile
```

Usual errors with merge are:
- Flipstrand - usually the case if not, remove the variants.
- Remove variants to to having 3+ alleles.

5. After merging, we need to check and remove related indivduals again.

```
king -b inputfile.bed --kinship --degree2 --prefix outputfile
```

Make a file with related individuals and remove. 

```
plink --bfile inputfile --remove related_individuals_file --make-bed --out outputfile
```

Additional QC measures can be applied: 

```
plink --bfile inputfile --ind 0.1 --geno 0.05 --hwe 0.00001 --make-bed --out outputfile
```

## 03 - Ancestry inference

Ancestry inference involves two steps: global ancestry inference and local ancestry inference. 

**Global ancestry inference**

1. Prepare reference file 

Download Reference data from the 1000GP website: 

```wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr*.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz```

*Populations used in 1000GP:* 

- GBR (40)
- CHB (40)
- YRI (40)

*Populations used in external databases:*

- Malay (40)
- Nama (40)

Before merging reference datasets, make sure the bim file format are the same for all datasets and all used the hg19 build for chromosome positions

Convert vcf files to plink binary format: 

``` plink --vcf inputfile --biallelic-only strict --recode --make-bed --out outputfile```


2. Merge reference and target datasets using steps outlined in 02 - Merging individual datasets.


3. Conduct LD pruning on the data before inferring global ancestry (we only want the haplotypes in tight LD blocks).

```
plink --bfile inputfile --indep-pairwise 50 10 0.1 
```

Remove pruned.out file fom data. 

```
plink --bfile inputfile --extract plink.prune.in --make-bed --out outputfile 
```


4. Run ADMIXTURE software on filtered binary files (see https://dalexander.github.io/admixture/). 

```
for k in {3..10}; do admixture --cv Target_Reference_merged_unrelated_LD_pruned.bed ${k} -j4 | tee log${k}; done
```


5. Visulaise results with PONG (see https://github.com/ramachandran-lab/pong/tree/master). 

The inputfiles required: 

- pong_filemap
- pop_order
- ind2pop

Run pong software: 

```
pong -m pong_filemap -n pop_order -i pop_ids
```

