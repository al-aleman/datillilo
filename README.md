![enter image description here](https://st3.depositphotos.com/9223738/19443/v/600/depositphotos_194439354-stock-illustration-desert-seamless-pattern-with-joshua.jpg)

Scripts and data for **"Divergent lineages in a young species: the case of Datilillo (*Yucca valida*), a broadly distributed plant from the Baja California Peninsula"**, [preprint here](https://www.biorxiv.org/content/10.1101/2023.05.22.541794v1). It is assumed that any required software, functions, libraries, packages (etc) are previously installed. Some arguments (e.g., the number of threads to use) should be adjusted to the system's specifications.

*Project description:* We examined the phylogeographic patterns of *Y. valida* throughout the species' geographic range. We hypothesised that past climatic fluctuations caused populations isolation, and since the species has a short-distance dispersal, we expected to find divergent lineages across its distribution. We genotyped 160 plants from 20 locations by NextRAD sequencing and aimed to i) identify the number of genetic lineages of *Y. valida* across its range, ii) reconstruct its populations' demographic history and iii) estimate the species' age.

### DATA ANALYSIS
**Raw data processing and SNP–calling**

[This](https://github.com/al-aleman/datillilo/blob/main/raw_reads_list.txt) is the list of the raw demultiplexed raw sequences, we're working on making them available via Dryad ASAP. The SNP dataset obtained after running the Stacks pipeline is already available [here](https://github.com/al-aleman/datillilo/blob/main/README.md#:~:text=structure/%0Acd%20structure-,Here,-is%20a%20compressed)

    # Quality check (pre-trimming) - should be done wherever the raw reads live
    mkdir quality
    fastqc -t 32 ./*.fastq.gz -o ./quality
    cd quality
    multiqc --interactive .
    # The multiqc_report.html file can be checked in a local browser

    # Quality-trimming + adapters removal + length forcing for Stacks
    cd ..
    mkdir cleandata    
    for i in $(ls *.fastq.gz | sed 's/.fastq.gz//')
    do
    bbduk.sh in=$i\.fastq.gz \
    out=./cleandata/$i\_cleaned.fastq.gz ktrim=r k=17 hdist=1 mink=8 \
    ref=/home/aaleman/DATA/software/bbmap/resources/nextera.fa.gz \
    trimq=10 minlen=110 ftr=119 ftl=10 qtrim=rl  trimq=10
    done

    # Quality check (post-trimming)
    cd cleandata
    mkdir quality
    fastqc -t 32 ./*.fastq.gz -o ./quality
    cd quality
    multiqc --interactive .
    # The multiqc_report.html file can be checked in a local browser

> The (**13**) samples below had less than one million raw reads and were
removed as a quality control filter before assemblying loci:
2458_CTAGTACG-ACTGCATA_S118_L002_R1_001_cleaned.fastq.gz
2458_TAAGGCGA-TAGATCGC_S2_L002_R1_001_cleaned.fastq.gz
2458_CTAGTACG-AGAGTAGA_S116_L002_R1_001_cleaned.fastq.gz
2458_CTAGTACG-GTAAGGAG_S117_L002_R1_001_cleaned.fastq.gz
2458_CTAGTACG-TAGATCGC_S113_L002_R1_001_cleaned.fastq.gz
2458_CTAGTACG-TATCCTCT_S115_L002_R1_001_cleaned.fastq.gz
2458_GTAGAGGA-GTAAGGAG_S61_L002_R1_001_cleaned.fastq.gz
1295_GTTGTGGC-AAGGAGTA_S625_L008_R1_001_cleaned.fastq.gz
1295_GTTGTGGC-CTAAGCCT_S628_L008_R1_001_cleaned.fastq.gz
2458_CGAGGCTG-CTAAGCCT_S95_L002_R1_001_cleaned.fastq.gz
2458_CTAGTACG-AAGGAGTA_S119_L002_R1_001_cleaned.fastq.gz
2458_CTAGTACG-CTCTCTAT_S114_L002_R1_001_cleaned.fastq.gz
2458_GCTACGCT-AAGGAGTA_S82_L002_R1_001_cleaned.fastq.gz

[valida.txt](https://github.com/al-aleman/datillilo/blob/main/valida.txt) is the popmap for Stacks (please verify that is tab- and not space-separated).

    # De-novo Stacks' assembly
    cd ..
    mkdir stacks
    denovo_map.pl --samples ./ --popmap ./valida.txt -o ./stacks -T 32 -M 3 -n 4 -m 8
    cd stacks
    # SNP-calling
    mkdir populations
    populations -P . -M ./../valida.txt -O ./populations/ -t 32 --vcf --min-maf 0.05 -r 0.5
    cd populations
    # populations.snps.vcf is our output of interest, we'll make a copy to modify freely knowing that this one is safe
    mkdir structure
    cp populations.snps.vcf ./structure/
    cd structure

Here is a compressed copy of [populations' output](https://github.com/al-aleman/datillilo/blob/main/populations.snps.vcf.gz).
We will change the sequencing identifiers for the sampling identifiers [with this script](https://github.com/al-aleman/datillilo/blob/main/replace.sh); additionaly  [popmap.txt](https://github.com/al-aleman/datillilo/blob/main/popmap.txt) and [GBS_SNP_filter.txt](https://github.com/al-aleman/datillilo/blob/main/GBS_SNP_filter.txt) will be required to use [GBS_SNP_filter.sh](https://github.com/laninsky/GBS_SNP_filter/tree/master)

    chmod +x replace.sh
    ./replace.sh

    # Chromosome-renamer (for Plink compatibility) 
    awk '/^#/ {print} !/^#/ {sub(/^([0-9]+)/, "chrom_&"); print}' populations.snps.vcf > valida.vcf
    
    # Ensure that GBS_SNP_filter.txt has 2 extra empty lines after #CHROM
    # GBS_SNP_filter.sh will be run simply to filter for one SNP/locus (prioritizing the SNP found in the most individuals -if a tie, then the SNP with the highest average coverage, and if a tie, a random SNP)
    # Our interest is the file *.oneSNP.vcf, therefore we do not care if things crash after this.
    GBS_SNP_filter.sh
    # I suggest taking valida.oneSNP.vcf to an independent folder as this is the dataset on which we can start making biological questions
    
 - [x] Raw data processing and SNP–calling is done!
 Here is a copy of [valida.oneSNP.vcf](https://github.com/al-aleman/datillilo/blob/main/valida.oneSNP.vcf)

**Genetic (nuclear) structure**
In a new folder, where (ideally) only *valida.oneSNP.vcf* lives

    # I prefer to use Plink v.1
    plink --vcf valida.oneSNP.vcf --allow-extra-chr --distance square --pca --make-bed
    # valida.oneSNP.dist and valida.oneSNP.dist.id are the pairwise genetic distance matrix results
    # valida.oneSNP.eigenval, and valida.oneSNP.eigenvec are the Principal Component Analysis (PCA) results
    
    # ADMIXTURE does not accept chromosome names that are not human chromosomes. We will just exchange the first column of the .bim file by 0
    awk '{$1="0";print $0}' valida.oneSNP.bim > valida.oneSNP.bim.tmp
    mv valida.oneSNP.bim.tmp valida.oneSNP.bim
        
    # ADMIXTURE RUN    
    for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
    do
    admixture --cv *.bed $K -j32 | tee log${K}.out
    done
    # The CV error can be seen with grep -h CV log*.out

*Results can be easily plotted using R:*

Two new popmaps will be needed, [pops.txt](https://github.com/al-aleman/datillilo/blob/main/pops.txt) and [PCA.txt](https://github.com/al-aleman/datillilo/blob/main/PCA.txt).

**First**, for the ADMIXTURE results, I like to use [Joana Meier's script](https://github.com/speciationgenomics/scripts/blob/master/plotADMIXTURE.r), it requires four arguments, the prefix for the ADMIXTURE output files (-p ), the file with the species information (-i ), the maximum number of K to be plotted (-k ), and a list with the populations or species separated by commas (-l <pop1,pop2...>). The list of populations provided with -l gives the order in which the populations or species shall be plotted. Example with K = 5 below:

    Rscript plotADMIXTURE.r -p valida.oneSNP -i pops.txt -k 5 -l Northern,Central, Southern

**Second**, this is how I plot the PCA:

    library(tidyverse)
    pca <- read.table("./valida.oneSNP.eigenvec")
    eigenval <- scan("./valida.oneSNP.eigenval")
    # set names
    pca <- pca[,-1]
    names(pca)[1] <- "ind"
    names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))
    spp_loc<-read.table("pops.txt", sep="\t", header=T)
    pca <- as_tibble(data.frame(pca, spp_loc))
    pca$Lineage <- factor(pca$Lineage, levels = c("Northern", "Central", "Southern"))
    pve <- data.frame(PC = 1:10, pve = eigenval/sum(eigenval)*100)
    a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
    a + ylab("Percentage variance explained") + theme_light() 
    cumsum(pve$pve)
    b <- ggplot(pca, aes(PC1, PC2, col = Lineage, shape=Lineage))
    b <- b + scale_color_manual(values = c("#006363","#FFB409","#D81B63")) #your colors here
    b <- b + scale_shape_manual(values=c(15, 16, 17))
    b <- b + theme_bw() + theme(panel.border =  element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +geom_hline(yintercept = 0, color="lightgray") +  geom_vline(xintercept = 0, color="lightgray")
    b + geom_point(size = 7)+ xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))

**Third**, for the NJ tree from the samples' pairwise genetic distance matrix, I join the IDs with the matrix by hand (*suggestions are welcome, if anyone knows a better way*), then transform it as below, and visualize the final tree in [iTOL 6](https://itol.embl.de/upload.cgi).

    library(ape)
    datamatrix = as.matrix(read.table("./valida.oneSNP.dist", header=T, fill=T))
    hc = hclust(as.dist(datamatrix), method="complete")
    phylotree = as.phylo(hc)
    plot(phylotree, type="fan")
    write.tree(phy=phylotree, file="NJ.nwk")

**Note: Genetic diversity** descriptors (HE, H.O., and FIS) for each lineage and the mean pairwise Weir and Cockerham's weighted values of genetic differentiation (standardised FST) were calculated with GENODIVE 3.0. To assess the isolation by distance null hypothesis, i.e., the correlation between the geographical and genetic distance between samples, a Mantel test was carried out correlating the individuals' Euclidean genetic distance and the pairwise distance between sampling locations (in kilometres, produced with the [Geographic Distance Matrix Generator](https://biodiversityinformatics.amnh.org/open_source/gdmg/)), after introducing a random jitter (+/- 0.05°) to samples belonging to shared sampling locations ([file here](https://github.com/al-aleman/datillilo/blob/main/coordinates_jitter.txt)). The significance of the Mantel test was evaluated with 99,999 permutations. Any required genetic data conversion was completed in PGDSpider 2.1.1.5. Alternatively (if not having acces to GENODIVE), Plink, VCFTools, Genepop, and/or Arlequin can be used to get these answers.

**Outlier SNPs and signatures of differentiation by local adaptation** analyses are performed in R, and the scripts can be found [here](https://github.com/al-aleman/datillilo/blob/main/PCAdapt_LEA.R). Here is the [matrix containing the environmental variables](https://github.com/al-aleman/datillilo/blob/main/valida.oneSNP.env) that was made with from WorldClim's Bioclimatic variables. From these results, **83 SNPs** were removed from the dataset to produce a [vcf with neutral SNPs](https://github.com/al-aleman/datillilo/blob/main/valida.neutral.vcf) for the demographic history analyses. To do this, I usually get the SNP list from the  [valida.oneSNP.vcf](https://github.com/al-aleman/datillilo/blob/main/valida.oneSNP.vcf) file (with VCFTools and the first two columns of the --missing-site output, for example), and then print the lines where the outliers are to make a flagged list of the positions to exclude.
