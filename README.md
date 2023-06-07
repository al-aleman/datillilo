![enter image description here](https://st3.depositphotos.com/9223738/19443/v/600/depositphotos_194439354-stock-illustration-desert-seamless-pattern-with-joshua.jpg)

Data and scripts for **"Divergent lineages in a young species: the case of Datilillo (*Yucca valida*), a broadly distributed plant from the Baja California Peninsula"**, [preprint here](https://www.biorxiv.org/content/10.1101/2023.05.22.541794v1). It is assumed that required software, functions, libraries, packages, modules (etc) were previously installed. Some arguments (e.g., the number of threads to use) should be adjusted to the system's specifications. 

> Note: Please read the software manuals if following these scripts for another research or dataset and adjust parameters accordingly.

**Project description:** We examined the phylogeographic patterns of *Y. valida* throughout the species' geographic range. We hypothesised that past climatic fluctuations caused populations isolation, and [since the species has a short-distance dispersal](https://academic.oup.com/biolinnean/article/136/2/364/6565171), we expected to find divergent lineages across its distribution. We genotyped 160 plants from 20 locations by NextRAD sequencing and aimed to i) identify the number of genetic lineages of *Y. valida* across its range, ii) reconstruct its populations' demographic history and iii) estimate the species' age using whole-chloroplast-genome data.

**Metadata:** Sampling and sequencing IDs plus sampling sites' coordinates can be accesed [here](https://github.com/al-aleman/datillilo/blob/main/data/sampling_metadata_yuccavalida.tsv).

---
### DATA ANALYSIS
### Raw data processing and SNP–calling

Note: We're currently uploading the second half of the demultiplexed raw sequences to Dryad, they'll be up shortly. The workflow and SNP data obtained after running the Stacks pipeline is already available [below](https://github.com/al-aleman/datillilo/blob/main/README.md#:~:text=structure/%0Acd%20structure-,Here,-is%20a%20compressed). 

Requirements: [FastQC](https://github.com/s-andrews/FastQC), [MultiQC](https://github.com/ewels/MultiQC) and [BBTools](https://github.com/kbaseapps/BBTools).

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
    ref=bbmap/resources/nextera.fa.gz \
    minlen=110 ftr=119 ftl=10 qtrim=rl  trimq=10
    done

    # Quality check (post-trimming)
    cd cleandata
    mkdir quality
    fastqc -t 32 ./*.fastq.gz -o ./quality
    cd quality
    multiqc --interactive .
    # The multiqc_report.html file can be checked in a local browser

These (**13**) [samples](https://github.com/al-aleman/datillilo/blob/main/data/flagged_raw_sequences.txt) had less than one million raw reads and were removed as a quality control filter before assemblying loci.

The file [valida.txt](https://github.com/al-aleman/datillilo/blob/main/data/valida.txt) is the popmap for running [Stacks](https://catchenlab.life.illinois.edu/stacks/) (please verify that is tab- and not space-separated). Note that Stacks' parameters were optimized according to [Paris et al. (2017)](https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.12775) and [Mastretta-Yanes et al. (2014)](https://onlinelibrary.wiley.com/doi/full/10.1111/1755-0998.12291). That involved iterating different ranges of values for the parameters *m* (the minimum number of identical reads required to create a stack), *M* (the number of mismatches allowed between loci on a single individual), and *n* (the number of mismatches allowed between loci when building the catalog). The Stacks' run below is based on the optimal set of parameters that maximized the amount of reliable information.

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

Here is a compressed copy of the [populations' output](https://github.com/al-aleman/datillilo/blob/main/data/populations.snps.vcf.gz).
We will change the sequencing identifiers for the sampling identifiers [with this script](https://github.com/al-aleman/datillilo/blob/main/scripts/replace.sh); additionaly [popmap.txt](https://github.com/al-aleman/datillilo/blob/main/data/popmap.txt) and [GBS_SNP_filter.txt](https://github.com/al-aleman/datillilo/blob/main/data/GBS_SNP_filter.txt) will be required to use [GBS_SNP_filter.sh](https://github.com/laninsky/GBS_SNP_filter/tree/master)

    chmod +x replace.sh
    ./replace.sh

    # Chromosome-renamer (for Plink compatibility) 
    awk '/^#/ {print} !/^#/ {sub(/^([0-9]+)/, "chrom_&"); print}' populations.snps.vcf > valida.vcf
    
    # Ensure that GBS_SNP_filter.txt has 2 extra empty lines after #CHROM
    # GBS_SNP_filter.sh will be run simply to filter for one SNP/locus
    # (prioritizing the SNP found in the most individuals -if a tie, then the SNP with the highest average coverage, and if a tie, a random SNP)
    # Our interest is the file *.oneSNP.vcf, therefore we do not care if things crash after this.
    ./GBS_SNP_filter.sh
    
    # I suggest taking valida.oneSNP.vcf to an independent folder as this is the dataset on which we can start making biological questions
    mkdir analyses
    mv valida.oneSNP.vcf analyses
    cd analyses
    
 - [x] Raw data processing and SNP–calling is done!
 Here is a copy of [valida.oneSNP.vcf](https://github.com/al-aleman/datillilo/blob/main/data/valida.oneSNP.vcf)

---
### Genetic (nuclear) structure
Requirements: [Plink 1.9](https://www.cog-genomics.org/plink/1.9/basic_stats), [ADMIXTURE](https://dalexander.github.io/admixture/), and [R](https://www.r-project.org/).

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

Results can be easily plotted using R (Two new popmaps will be needed, [pops.txt](https://github.com/al-aleman/datillilo/blob/main/data/pops.txt) and [PCA.txt](https://github.com/al-aleman/datillilo/blob/main/data/PCA.txt)):

For the ADMIXTURE results, I like to use [Joana Meier's script](https://github.com/speciationgenomics/scripts/blob/master/plotADMIXTURE.r), it requires four arguments: the prefix for the ADMIXTURE output files (-p ), the file with the species information (-i ), the maximum number of *K* to be plotted (-k ), and a list with the populations or species separated by commas (-l <pop1,pop2...>). The list of populations provided with -l gives the order in which the populations or species shall be plotted. Example with *K* = 5 below:

    Rscript plotADMIXTURE.r -p valida.oneSNP -i pops.txt -k 5 -l Northern,Central,Southern

This is how I plot the PCA in R:

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

For the NJ tree from the samples' pairwise genetic distance matrix, I join the IDs with the *.dist matrix by hand (suggestions are welcome if anyone knows a better way), then I transform it as below in R, and visualize the final tree in [iTOL 6](https://itol.embl.de/upload.cgi).

    library(ape)
    datamatrix = as.matrix(read.table("./valida.oneSNP.dist", header=T, fill=T))
    hc = hclust(as.dist(datamatrix), method="complete")
    phylotree = as.phylo(hc)
    plot(phylotree, type="fan")
    write.tree(phy=phylotree, file="NJ.nwk")

---
### Genetic diversity

Descriptors (HE, HO, and FIS) for each lineage and the mean pairwise Weir and Cockerham's weighted values of genetic differentiation (standardised FST) were calculated with [GENODIVE](https://www.bentleydrummer.nl/software/software/GenoDive.html). To assess the isolation by distance null hypothesis, i.e., the correlation between the geographical and genetic distance between samples, a Mantel test was carried out correlating the individuals' Euclidean genetic distance and the pairwise distance between sampling locations (in kilometres, produced with the [Geographic Distance Matrix Generator](https://biodiversityinformatics.amnh.org/open_source/gdmg/)), after introducing a random jitter (+/- 0.05°) to samples belonging to shared sampling locations ([file here](https://github.com/al-aleman/datillilo/blob/main/data/coordinates_jitter.txt)). The significance of the Mantel test was evaluated with 99,999 permutations. Genetic data conversion from *.vcf* to *.gen* was completed in [PGDSpider](http://www.cmpg.unibe.ch/software/PGDSpider/). Alternatively (if not having acces to GENODIVE), Plink, [VCFtools](https://vcftools.github.io/index.html), [Genepop](https://genepop.curtin.edu.au/), and/or [Arlequin](http://cmpg.unibe.ch/software/arlequin35/) can be used to run these analyses.

---
### Outlier SNPs and signatures of differentiation by local adaptation

Analyses are performed in R, and the scripts can be found [here](https://github.com/al-aleman/datillilo/blob/main/scripts/PCAdapt_LEA.R). Here is the [matrix containing the environmental variables](https://github.com/al-aleman/datillilo/blob/main/data/valida.oneSNP.env) that was made with from WorldClim's Bioclimatic data. From these results, **83 SNPs** were removed of the dataset (63 unique SNPs potentially associated with adaptation to multiple environmental variables identified by LFMM and 20 putative outlier SNPs detected by PCAdapt) to produce a [vcf with neutral SNPs](https://github.com/al-aleman/datillilo/blob/main/data/valida.neutral.vcf) for the demographic history analyses. To do this, I usually get the SNP list from the  [valida.oneSNP.vcf](https://github.com/al-aleman/datillilo/blob/main/data/valida.oneSNP.vcf) file (with VCFtools and the first two columns of the -*-missing-site* output, for example), and then print the lines where the outliers are to make a flagged list of the positions to exclude.

---
### Demographic history

Requirements: [δaδi](https://bitbucket.org/gutenkunstlab/dadi/src/master/), Daniel Portik's [δaδi_pipeline](https://github.com/dportik/dadi_pipeline) and [Stacks_pipeline](https://github.com/dportik/Stacks_pipeline), and [easySFS](https://github.com/isaacovercast/easySFS). With the neutral *vcf* as input (ideally, in a new folder), we're back to populations (Stacks). We can use the neutral dataset to get the *haplotypes.tsv file, which we will input to [Convert_tsv_to_dadi.py](https://github.com/dportik/Stacks_pipeline/blob/master/stacks-pipeline-scripts/Convert_tsv_to_dadi.py) and produce [the SNP file](https://github.com/al-aleman/datillilo/blob/main/scripts/%CE%B4a%CE%B4i_pipeline/SNPs_file_Northern_Central_Southern.txt) to run δaδi. 

> This is a friendly reminder to test which python version works best for you with each software.

    populations -V valida.neutral.vcf -O . -M pops.txt -t 32 --vcf
    
    # The script below will make the file SNPs_file_Northern_Central_Southern.txt, which will be needed to run Models_3D.py
    python Convert_tsv_to_dadi.py -i valida.neutral.haplotypes.tsv -o . -p popmap.txt
    
    # easySFS is a super-fast way to preview which projection is best, I do not use it for data conversion
    # I found 18,18,18 as a good trade-off to maximizing an even number of segregating sites (second number) and a balanced sample size (first number)
    python easySFS.py -i valida.neutral.vcf -p pops.txt --preview -a

[Here](https://github.com/al-aleman/datillilo/tree/main/scripts/%CE%B4a%CE%B4i_pipeline) is everything needed to test the four possible three-populations simple models of i) simultaneous, ii) admixed, and consecutive divergence from iii) North to South and iv) South to North, without gene flow or changes in populations' sizes. They can be easily run as below. WARNING: This is  slow. I'm sharing the [outputs of these analyses](https://github.com/al-aleman/datillilo/tree/main/scripts/%CE%B4a%CE%B4i_pipeline/results) as a consideration for time constrains. According to [the results' summary](https://github.com/al-aleman/datillilo/blob/main/scripts/%CE%B4a%CE%B4i_pipeline/results/Results_Summary_Extended.txt), the optimal demographic model was the simultaneous divergence of the three nuclear genetic lineages from one common ancestral population.

    python dadi_Run_3D_Set_North-South.py
    python dadi_Run_3D_Set_South-North.py
    python dadi_Run_3D_Set_Simultaneous.py
    python dadi_Run_3D_Set_Admixed.py
    
    # Once finished, results can be summarized by running
    Summarize_Outputs.py .

 - [x] Objectives i) (identifying the number of genetic lineages of *Y. valida* across its range) and ii) (reconstructing its populations' demographic history) are done!

---
### Whole chloroplast genome sequences recovery

This is a reference-guided workflow to reconstruct whole-chloroplast-genome sequences (as a by-product of the enzymatic fragmentation for the high-throughput sequencing libraries preparation without isolating cpDNA), as shown by [Aleman et al. (2023)](https://www.biorxiv.org/content/10.1101/2023.04.21.537876v1). We included raw-sequencing-data of 40 samples from [Arteaga et al. (2020)](https://www.frontiersin.org/articles/10.3389/fpls.2020.00685/full) (sequences list [here](https://github.com/al-aleman/datillilo/blob/main/data/sequences_capensis.txt), see publication for metadata) and used the chloroplast genome of *Y. schidigera* (GenBank: [NC_032714.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_032714.1)) as reference.

Requirements: [BBTools](https://github.com/kbaseapps/BBTools), [Bowtie](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml), [Samtools](http://www.htslib.org/), [ANGSD](http://www.popgen.dk/angsd/index.php/ANGSD),  this [stats compiler script](https://github.com/al-aleman/datillilo/blob/main/scripts/mapping_stats.sh),  [another ID replacer](https://github.com/al-aleman/datillilo/blob/main/scripts/another_replacer.sh), [a file renamer](https://github.com/al-aleman/datillilo/blob/main/scripts/fasta_renamer.sh), and optionally, [seqtk](https://github.com/lh3/seqtk). 

    # Quality-trimming + adapters removal
    cd ..
    mkdir cleandata    
    for i in $(ls *.fastq.gz | sed 's/.fastq.gz//')
    do
    bbduk.sh in=$i\.fastq.gz \
    out=./cleandata/$i\_cleaned.fastq.gz ktrim=r k=17 hdist=1 mink=8 \
    ref=bbmap/resources/nextera.fa.gz \
    trimq=10 minlen=100 qtrim=rl  trimq=10
    done
    cd cleandata
    
    # Indexing the reference genome (the fasta file should be inside the folder)
    bowtie2-build ./schidigera.fasta schidigera
    mkdir plastome
    
    # Read-mapping
    for i in $(ls *.fastq.gz | sed 's/.fastq.gz//')
    do 
    bowtie2 -x ./schidigera -U $i\.fastq.gz -S ./plastome/$i.sam --very-sensitive-local -p 32
    done
    cd plastome
    mkdir bam
    
    # Transform the .sam files into .bam
    for i in $(ls *.sam | sed 's/.sam//')
    do 
    samtools view -S -b $i\.sam -@ 32 | samtools sort -@ 32 > ./bam/$i.bam
    done
    cd bam
    mkdir chloroplast_sequences coverage flagstats
    
    # Mapping-statistics (merely informative)
    ./mapping_stats.sh
    # Let's remember to change the sequencing identifiers for the sampling identifiers
    chmod +x another_replacer.sh
    ./another_replacer.sh
    # Outputs are pretty intuitive; however, depth calculation is estimated considering sites with breath of coverage = 0 to be depth = 0 
    # This causes an underestimation. Solving it uniquely requires a simple calculation (real depth = (estimated depth * 100) / breath of coverage)
    
    # Nucleotide-calling 
    for i in $(ls *.bam | sed 's/.bam//')
    do
    angsd -dofasta 2 -doCounts 1 -i $i.bam -minQ 20 -minMapQ 20 -nThreads 32 -r NC_013823.1 -out ./chloroplast_sequences/$i.plast
    done
    cd chloroplast_sequences
    
    # Outputs of ANSGD will be called according to the sequencing identifiers, but the headers of every file will be based on the name of the reference genome. Let's change this.
    # First, change the sequencing identifiers for the sampling identifiers in the files' names
    chmod +x fasta_renamer.sh
    ./fasta_renamer.sh
    
    # Now the fasta headers
    for sequence in *.fa;
    do
    awk '/^>/{print ">" substr(FILENAME,1,length(FILENAME)-3); next} 1' $sequence > $sequence.fasta
    done
    
    # Let's bring all sequences (by species) on one multi-FASTA
    cat V*fasta > valida.fasta
    cat C*fasta > capensis.fasta
    
    # I prefer FASTA files to be one-liners
    seqtk seq valida.fasta > valida.fa
    seqtk seq capensis.fasta > capensis.fa

 - [x] Our whole-chloroplast-genome-sequences by species are ready!
 [Here](https://github.com/al-aleman/datillilo/blob/main/data/capensis_valida_fasta.zip) is a folder with the chloroplast sequences for each species, one species per file.

Then, based on the mapping statistics I subset the [sequences from *Y. valida* with a breadth of coverage >0.5](https://github.com/al-aleman/datillilo/blob/main/data/subset_fastaIDs.txt) and concatenated them with the chloroplast genome reference from *Y. schidigera* for the intraspecific phylogenetic analyses.

    seqtk subseq valida.fa subset_fastaIDs.txt > subset.fasta
    cat subset.fasta ./../../../schidigera.fasta > subset_root.fasta
    mkdir raxml
    seqtk seq subset_root.fasta > ./raxml/pre_beast.fa

Additionally,  using [this python script](https://github.com/al-aleman/datillilo/blob/main/scripts/consensus.py), I produced chloroplast consensus sequences for *Y. valida* and *Y. capensis*, using all the available sequences for each species.

    mkdir beast
    python consensus.py valida.fa valida_cons.fa
    python consensus.py capensis.fa capensis_cons.fa
    cat *_cons.fa > ./beast/consensus.fa

---
### Phylogenetic relationships and molecular clock analyses
Requirements: [RAxML-ng](https://github.com/amkozlov/raxml-ng), [BEAST](https://github.com/beast-dev/beast-mcmc/releases/tag/v1.8.4) (I found version 1.8.4 to be highly straightforward for replicating [Smith et al. (2021)](https://bsapubs.onlinelibrary.wiley.com/doi/full/10.1002/ajb2.1633)), and [MAFFT](https://mafft.cbrc.jp/alignment/server/).

    cd raxml
    raxml-ng --all --msa *fa --model GTR+I+G --prefix tree --seed 12345 \
    --outgroup Yucca_schidigera --bs-metric fbp,tbe --tree rand{1000} \
    --bs-trees autoMRE --threads 32
    # RAxML run (the outputs of interest are *.raxml.supportFBP and *.raxml.supportTBE)

RAxML results can be visualized on [iTOL 6](https://itol.embl.de/upload.cgi).

The consensus sequences were aligned to the chloroplast reference genomes of 15 species of the Agavoideae (*[Agave attenuata](https://www.ncbi.nlm.nih.gov/nuccore/NC_032696.1/)*, *[Beschorneria septentrionalis](https://www.ncbi.nlm.nih.gov/nuccore/NC_032699.1/)*, *[Camassia scilloides](https://www.ncbi.nlm.nih.gov/nuccore/NC_032700.1/)*, *[Chlorogalum pomeridianum](https://www.ncbi.nlm.nih.gov/nuccore/NC_032701.1/)*, *[Hesperaloe campanulata](https://www.ncbi.nlm.nih.gov/nuccore/NC_032702.1/)*, *[Hesperaloe parviflora](https://www.ncbi.nlm.nih.gov/nuccore/NC_032703.1/)*, *[Hesperocallis undulata](https://www.ncbi.nlm.nih.gov/nuccore/NC_032704.1/)*, *[Hesperoyucca whipplei](https://www.ncbi.nlm.nih.gov/nuccore/NC_032705.1/)*, *[Hosta ventricosa](https://www.ncbi.nlm.nih.gov/nuccore//NC_032706.1)*, *[Manfreda virginica](https://www.ncbi.nlm.nih.gov/nuccore/NC_032707.1/)*, *[Schoenolirion croceum](https://www.ncbi.nlm.nih.gov/nuccore/NC_032710.1/)*, *[Yucca brevifolia](https://www.ncbi.nlm.nih.gov/nuccore/NC_032711.1/)*, *[Yucca filamentosa](https://www.ncbi.nlm.nih.gov/nuccore/NC_032712.1/)*, *[Yucca queretaroensis](https://www.ncbi.nlm.nih.gov/nuccore/NC_032713.1/)*, and *[Yucca schidigera](https://www.ncbi.nlm.nih.gov/nuccore/NC_032714.1/)*)  using MAFFT on the CBRC server, applying the default settings plus the flag *--nzero*.

PGDSpider 2.1.1.5 was used to convert data from *.fasta* to *.nex*. The aligned [FASTA](https://github.com/al-aleman/datillilo/blob/main/data/smith_aleman.fa) and [NEXUS](https://github.com/al-aleman/datillilo/blob/main/data/smith_aleman.nex) files are available.

BEAUti was used to generate the BEAST XML file, with the following settings: A general time reversible model estimating the proportion of invariant sites (GTR + I) was employed. Gamma-distributed priors for the five substitution types and uniform priors (0-1) for the base frequencies and the proportion of invariant sites were established. An uncorrelated relaxed clock model with log-normal distributed rates was set up (mean = 1, standard deviation = 0.33), assuming a Yule process for the tree prior, starting from a random tree, establishing *Hosta ventricosa* as the outgroup, and forcing the ingroup's monophyly. A log-normal prior distribution on the age of the ingroup was specified to 14.5 million years (offset, mean = 1 standard deviation = 1.4), [based on the minimum age of the fossil *Protoyucca shadesii*](https://onlinelibrary.wiley.com/doi/full/10.1111/boj.12233), completing two independent runs of 60 million steps long with a sample every 6,000,000 generations. The runs were joined in LogCombiner 1.8.4, considering a 10% burn-in for each input. A maximum credibility tree was produced in TreeAnnotator 1.8.4, and each node's means, medians, and 95% confidence intervals were outputted. The final tree was visualized in FigTree. Here is the [XML-input for BEAST](https://github.com/al-aleman/datillilo/blob/main/data/smith_aleman.xml); and it can be run from the command-line as:

    beast -seed 12345 -threads 32 *xml

 - [x] Objective iii) estimating the species' age using whole-chloroplast-genome data is done! The only thing left is intepreting these results.

### END OF DATA ANALYSIS

*Thinking-chair time*
![enter image description here](https://images-wixmp-ed30a86b8c4ca887773594c2.wixmp.com/f/1032b378-9e4c-4e5b-9344-03d0fafe3154/ddzq14w-74e76891-57a8-4e2c-804a-32410994090d.png/v1/fill/w_1192,h_670,q_70,strp/blue_s_clues_and_you_living_room_by_jack1set2_ddzq14w-pre.jpg?token=eyJ0eXAiOiJKV1QiLCJhbGciOiJIUzI1NiJ9.eyJzdWIiOiJ1cm46YXBwOjdlMGQxODg5ODIyNjQzNzNhNWYwZDQxNWVhMGQyNmUwIiwiaXNzIjoidXJuOmFwcDo3ZTBkMTg4OTgyMjY0MzczYTVmMGQ0MTVlYTBkMjZlMCIsIm9iaiI6W1t7ImhlaWdodCI6Ijw9NzAyIiwicGF0aCI6IlwvZlwvMTAzMmIzNzgtOWU0Yy00ZTViLTkzNDQtMDNkMGZhZmUzMTU0XC9kZHpxMTR3LTc0ZTc2ODkxLTU3YTgtNGUyYy04MDRhLTMyNDEwOTk0MDkwZC5wbmciLCJ3aWR0aCI6Ijw9MTI0OCJ9XV0sImF1ZCI6WyJ1cm46c2VydmljZTppbWFnZS5vcGVyYXRpb25zIl19.QvIJO1CX4j-Tg0ebPf7OAI0cCL4a3LVNwP_noSFkP28)
