# Data Harmonisation Procedures
Guidelines for merging genome-wide genotype data obtained from multi-way admixed populations.

## Steps
1. Processing individual datasets.
2. Merging of datasets.
3. Ancestry inference.
4. Batch effect correction.

**Software required for procedures**
- PLINK (https://www.cog-genomics.org/plink/2.0/)
- VCFTOOLS v0.1.17 (https://vcftools.github.io/index.html)
- BCFTOOLS (https://samtools.github.io/bcftools/bcftools.html)
- Eagle v2.4.1 (phases genotypes prior to imputation - https://alkesgroup.broadinstitute.org/Eagle/#x1-30002)
- Minimac4 (https://github.com/statgen/Minimac4)
- KING (https://www.kingrelatedness.com)
- SHAPEIT2 (https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html)
- ADMIXTURE (https://dalexander.github.io/admixture/)
- PONG (https://brown.edu/Research/Ramachandran_Lab/files/pong/pong-manual.pdf)
- RFMix (https://github.com/slowkoni/rfmix)

## 01 - Processing individual datasets 
**Initial quality control**

Individual datasets must undergo independent quality control procedures prior to merging. Remove SNPs and individuals with missing information, monomorphic sites (for which the MAF < 1/2n, where n is the number of samples in the dataset) and SNPs deviating from Hardy-Weinberg equilibrium (HWE). Remove sex chromosomes and chromosome 23. Check for any phenotypic information of population under investigation, in order to remove any indivual who have any missing phenotypic information, such as age or sex or disease status. In our study, we filtered out variants with a genotype call rate less than 95%, variants that deviated from the HWE p-value threshold of 0.00001. We also removed individuals with average genotype call rates less than 90%. These QC thresholds may be adjusted based on the specific requirements of your study. 

``` 
plink --bfile filename --mind 0.1 --geno 0.05 --hwe 0.00001 --maf -chr1-22 --make-bed --out outputfile
```

VCFTOOLS can be used to perform genotype QC on files in VCF file format. 

```
vcftools --vcf input_filename --max-missing 0.95 --hwe 0.00001 --maf --max-alleles 2 --recode --recode-INFO-all --out output_filename
```

Convert the PLINK binary files to VCF file format for imputation. 

```
plink --bfile binary_fileset --recode vcf-iid --out vcf_file 
```

**Special considerations for merging with WGS data**

Merging WGS data with array genotype data requires special considerations. Genotype data obtained from different cohorts should only be combined using intersecting or common SNPs across all groups to avoid high rates of missing data in the final merged dataset. The high coverage genotype data obtained from WGS has the ability to completely overlap with markers on any array. However, merging WGS and array data can be complicated by differences in genotyping technologies with a high potential for introducing batch effects. The following describe steps for preparing WGS data for merging with genotype data:

1.	Generate the list of SNPs to be extracted from the WGS data. The list should be in a textt format with the rsIDs or position of the SNP in the chromosome:basepair format (one SNP per line, no headers). Regardless of the method, ensure that the corresponding reference genome used for genotype calling across datasets is consistent.

2.	Extract the SNP list from WGS data using PLINK and convert the output file to binary PLINK file format.

```
plink --vcf input_filename.vcf --extract SNP_list.txt --make-bed --out output_filename
```

3.	Follow QC, imputation and post-imputation QC procedures outlined above.

**Imputation**

Imputing separate datasets is recommended prior to merging to increase the number of overlapping SNPs between datasets. This step is crucial when merging separate datasets genotyped on disparate genotyping platforms. There are several publicly available imputation servers which offer different imputation algorithms and reference panels: 
- Sanger imputation server (https://imputation.sanger.ac.uk)
- TOPMed imputation server (https://imputation.biodatacatalyst.nhlbi.nih.gov)
- Michigan imputation server (https://imputationserver.sph.umich.edu/index.html)

In the event that researchers have a local/personalized reference panel for imputation, or do not wish to use a remote imputation server, phasing and imputation can be performed locally. Here, we provide steps for performing phasing and imputation using Eagle and MiniMac4:

1.	Split target and reference VCFs into individual chromosome files using PLINK. 

```
for i in {1..22}; do plink --vcf TargetSamples.vcf --chr ${i} --recode vcf --out TargetSamples _chr${i}; done

for i in {1..22}; do plink --vcf ReferenceSamples.vcf --chr ${i} --recode vcf --out ReferenceSamples _chr${i}; done
```

2.	Phase target per-chromosome VCFs using Eagle and the phased reference VCF.

```
for i in {1..22}; do Eagle_v2.4.1/eagle --noImpMissing --vcfOutFormat=z --vcfTarget TargetSamples _chr${i}.vcf --geneticMapFile Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz --outPrefix TargetSamples _chr${i}.phased --vcfRef ReferenceSamples_chr${i}.vcf; done
```

3.	Impute missing genotypes in phased target per-chromosome VCFs using Minimac4 and the phased reference VCF.

```
for i in {1..22}; do minimac4 --refHaps ReferenceSamples_chr${i}.vcf --haps TargetSamples _chr${i}.phased.vcf.gz --allTypedSites --cpus 10 --prefix TargetSamples_imputed_chr${i}; done
```

The per-chromosome imputed VCF files can then be concatenated using BCFTools to produce one VCF (containing all chromosomes and individuals). 

```
bcftools concat TargetSamples_imputed_chr*.vcf -Oz -o allChromosomes_impute.vcf.gz
```

The accuracy and quality of imputation performed by different imputation algorithms and reference panels on different populations is an area of ongoing research. The population structure of the tagret population and similarity to the imputation reference panel must be carefully considered before imputation. We recommend investigating the performance of different imputation algorithm and reference panel combinations to determine the best strategy for your unique target population. 

For our five-way admixed South African population, we employed the PBWT imputation algorithm and the African Genome Resource panels, made available through the Sanger imputation server. The server provides detailed instructions for preparing files for imputation (https://imputation.sanger.ac.uk/?instructions=1). 

**Post-imputation quality control** 

We recommend removing poorly imputed SNPs from individual datasets to minimise spurious associations in downstream analyses. To assess imputation quality we consider the INFO score quality metric. The INFO score ranges from 0 to 1, where values near 1 indicate that a SNP has been imputed with a high degree of certainty. BCFTOOLS can be used to calculate the INFO scores. 

```
for i in {1..22}; do bcftools +impute-info  ${i}.imputed.vcf.gz  -Oz -o chr${i}_imputed_info.vcf.gz; done
```

We can choose an INFO score threshold to filter out poorly imputed sites. We recommend filtering out all SNPs with INFO scores less than 0.8. This can be done with BCFTOOLS. Again, the INFO score threshold may be adjusted based on the specific requirements of your study.

```
bcftools filter -e 'INFO/INFO<0.8' input.vcf.gz -Oz -o output.vcf.gz
```

We need to check and remove related indivduals within individual datasets before merging. The software KING and kinship coefficient is sufficient to identify close relationships between participants who have extensive population structure. Usually up to 2nd Degree relatedness, indicated by a kinship coefficient of 0.0884-0.177. Identical by descent (IBD) segments are rapidly inferred and individuals in relationship paris (such as parent/offpring or sibling pairs) are identified using KING:

```
king -b input_filename.bed --ibdseg
```

One individual from each relationship-pair identified should be chosen to be removed from the dataset. The remove_relatives.py script (developed as part of the PONDEROSA algorithm - https://github.com/williamscole/PONDEROSA/tree/master) compiles a list of unrelated individuals, given the .seg file generated by the above KING command and the .fam file from PLINK. Pairs of individuals who are less than 2nd degree relatives will be identified in the following command:

```
python remove_relatives.py None king.seg input_filename.fam 2nd
```

Use PLINK and the output list of unrelated individuals to filter related individuals from the dataset:

```
plink --bfile input_filename --keep unrelated_individuals.txt --make-bed --out output_filename
```

Repeat this procedure for each separate dataset.

## 02 - Merging individual datasets

If required, convert genomic positions to the same coordinate system (hg19, hg38, etc.). The LiftoverVcf tool from GATK can be used to “lift over” genetic positions in a VCF from one genome build to another. The LiftoverVCF tool forms part of the Picard command-line tools, which are provided in a single executable picard.jar file (which can be downloaded from https://github.com/broadinstitute/picard/releases/tag/3.1.1). Follow the instructions at https://broadinstitute.github.io/picard/ to download and install the Picard toolset. 

Once the toolset has been installed, run the following commands to convert VCF genomic coordinates from one genome build to another (e.g., from hg19 to hg38). 

```
java -jar picard.jar LiftoverVcf \\
     I=input.vcf \\
     O=lifted_over.vcf \\
     CHAIN=b37tohg38.chain \\
     REJECT=rejected_variants.vcf \\
     R=reference_sequence.fasta
```

The merge function in PLINK appears to exclude variants that are common in both data sets, which is unexpected. The following steps extract the variants common in both file sets and merge the data:

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

*Usual errors with merge are:*

*- Flipstrand - usually the case if not, remove the variants.*

*- Remove variants to to having 3+ alleles.*

5. After merging, we need to check and remove related indivduals again.

```
king -b input_filename.bed --ibdseg

python remove_relatives.py None king.seg input_filename.fam 2nd
```

Use PLINK and the output list of unrelated individuals to filter related individuals from the dataset:

```
plink --bfile input_filename --keep unrelated_individuals.txt --make-bed --out output_filename
```

Additional QC measures can be applied based on the requirements of your study: 

```
plink --bfile inputfile --ind 0.1 --geno 0.05 --hwe 0.00001 --make-bed --out outputfile
```

## 03 - Ancestry inference

Ancestry inference involves two steps: global ancestry inference and local ancestry inference. 

**Global ancestry inference**

1. Prepare reference file. 

Download Reference data from the 1000GP website: 

```wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr*.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz```

*Populations used in 1000GP:* 

- GBR (40)
- CHB (40)
- YRI (40)

*Populations used in external databases:*

- Malay (40)
- Nama (40)

Before merging reference datasets, make sure the bim file format are the same for all datasets and all used the hg19/GRCh37 build for chromosome positions.

Convert vcf files to plink binary format: 

```
plink --vcf inputfile --biallelic-only strict --recode --make-bed --out outputfile
```


2. Merge reference and target datasets using steps outlined in *02 - Merging individual datasets*.


3. Conduct LD pruning on the data before inferring global ancestry (we only want the haplotypes in tight LD blocks).

```
plink --bfile inputfile --indep-pairwise 50 10 0.1 
```

Remove ```pruned.out``` file fom data. 

```
plink --bfile inputfile --extract plink.prune.in --make-bed --out outputfile 
```


4. Run ADMIXTURE software on filtered binary files (see https://dalexander.github.io/admixture/). 

```
for k in {3..10}; do admixture --cv Target_Reference_merged_unrelated_LD_pruned.bed ${k} -j4 | tee log${k}; done
```


5. Visualise results with PONG (see https://github.com/ramachandran-lab/pong/tree/master). 

*The inputfiles required:* 

- pong_filemap
- pop_order
- ind2pop

Run pong software: 

```
pong -m pong_filemap -n pop_order -i pop_ids
```

*The global ancestry results obtained from this step is crucial to examine if the correct reference/ancestral populations are used for ancestry inference. RFMix will be used to refine global ancestry proportions.* 


**Local ancestry inference**

Per chromosome VCF files must be phased for ancestry inference using RFMix. Recombination maps are required to phase data, therefore inferring where each allele came from which parent. We used the African American genetic map build 37. This map was based on data from the HapMap consortium and includes genetic data from individuals with 80% West African ancestry and 20% European ancestry (see https://ftp.ncbi.nlm.nih.gov/hapmap/recombination/2011-01_phaseII_B37/).

1. Split data into individual chromosomes before phasing.

```
for i in {1..22}; do plink --bfile inputfile --chr ${i} --make-bed --out outputfile.chr${i}; done
```

2. Phase each chromosome separately using SHAPEIT2.

```
shapeit --input-bed inputfile.bed inputfile.bim inputfile.fam --input-map recombinationmap ouput-max outputfile.haps
```

3. Convert phased binary files back to VCF format.

```
shapeit -convert --input-haps inputfile.haps --output-vcf outputfile.vcf
```

4. Prepare reference and query VCF files required for RFMix.

```
vcftools --vcf inputfile.vcf --keep query_individuals.txt --recode --out outputfile_query.vcf
vcftools --vcf inputfile.vcf --keep reference_individuals.txt --recode --out outputfile_query.vcf
```

5. Run RFMix using default parameters and the ```---reanalyze-reference``` flag (see https://github.com/slowkoni/rfmix/blob/master/MANUAL.md).

```
rfmix -f inputfile_query.vcf -r inputfile_reference.vcf -m sample_map_SAC -g geneticmap -n 4 -G 15 -e 3 -o outputfile --chromosome --reanalyze-reference
```

*Outputfiles generated by RFMix*

- rfmix.Q. *The global ancestry proportions in this file can be used as covariates in the linear regression models in step 04 - Batch effect correction, or in other downstream analyses.*
- sis.tsv
- msp.tsv
- fb.tsv

## 04 - Batch effect correction 

A pseudo-case-control comparison method to detect batch effects involves coding case/control status by batch followed by running an association analysis testing each batch against all other batches. For example, the status of all samples from one dataset will be coded as a case, while the status of every other sample is to be coded as a control. A genome-wide association test will then be performed. This process will be repeated for each batch. If any single dataset has more positive signals compared to the other datasets then batch effects may be responsible for producing spurious results. If batch effects are present, the genomic inflation factor (λ) for these “pseudo-GWASs” will be greater than one. Batch effects can be resolved by removing those SNPs which pass the threshold for suggestive significance (P-value < 1x10-4) from the merged dataset, as these SNPs are affected by batch effects. 

1. Concatenate all per chromosome files together using BCFTOOLS.

```
bcftools concat inputfile_chr*.vcf -Ov -o output.vcf
```

2. Convert VCF files to plink binary format: 

```
plink --vcf inputfile --biallelic-only strict --recode --make-bed --out outputfile
```

*Steps for pseudo-case-control comparison.*

1. Code one batch/dataset as case in the FAM file.

2. Code all other batches/datasets as controls in the FAM file.

3. Run logistic regression model (use sex, age and global ancestry proportion covariates). You will require a file containing all the phenotype information of your cohort. Example Covariates_GA_PC.txt

```
plink --bfile inputfile --glm sex --covar Covariates_GA_PC.txt --covar-name AGE,AFR,EUR,SAN --covar-variance-standardize --adjust --ci 0.95 --out GWAS
```

4. Create a list of SNPs surpassing suggestive threshold (P-value < 10-4).

5. Repeat for each batch/dataset.

6. Remove SNPs affected by batch from merged dataset.

```
plink --bfile inputfile --extract BatchEffectedSnps.txt --make-bed --out outputfile
```

**Your dataset is now ready for downstream analyses!**
