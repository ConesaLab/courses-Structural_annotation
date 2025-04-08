# Sructural annotation course

This tutorial will guide you through the porcess of creating an structural annotation using a reference genome and long reads RNA data from a PacBio sequencing experiment. 


# 1. Gene model training

AUGUSTUS has multiple gene models that can be used to produce an initial annotation on a given genome. However, this models are tailored to specific species, being adapted to their genomes and gene content. Thus, in the case of having a species that is not present within Augustus model, we need to create our own gene model. In order to achieve this, in this tutorial we will create two main resources: 

1. A gene model from a set of high-confidence genes
2. Hints derived from an lr-RNASeq experiment

## 1.1 Preprocessing and gene set creation

The first step in the workflow, and where our tutorial begins, is to create a set of high confidence genes, with known their coding sequence (CDS).These genes can originate from various sources, such as a previous version of the organism's genome, a closely related species, or predictions based on the genome under study. For the purposes of this tutorial, we will tackle the most complex scenario: deriving high-confidence genes directly from our own genome. While this may appear redundant, this step is essential to build an initial gene model that will serve as a foundation for evidence-driven gene prediction.

To achieve this, we will utilize **BUSCO** (Benchmarking Universal Single-Copy Orthologs), a widely used tool designed to assess the completeness of genome assemblies and annotated gene sets. BUSCO identifies genes based on their evolutionary conservation as single-copy orthologs across specific phylogenetic branches. These conserved genes are expected to be present in all species within a lineage and typically exist as single copies. We can take advantage of this well conserved genes to produce an initial training set for AUGUSTUS based on the single-copy core genes that we find in our genome.

BUSCO offers multiple lineage-specific datasets tailored to different lineages of organisms. It is highly important to select the dataset that is the closest to our target species, to have as many genes as possible in the initial search. In order to check which datasets are available in BUSCO, and the one that is the most suitable for a given experiment, you can run the following command:

```bash
busco --list-dataset
```

This command will produce a list of the different datasets that BUSCO has available in its database. You will have to select the closest to you sample's taxonomy. 

❓**Trivia: Which dataset would you select for the tutorial sample?**  
<details><summary>Solution</summary>
Here is some more text that was hidden before.
</details><br>

Once the dataset is selected, BUSCO has multiple modes to be run in, and different software that will do the gene search. First of all, BUSCO can be run in `genome`, `proteome` or `transcriptome` mode, which is dependant on the type of input given. In our case, we will use the `genome` mode, as we want to do a full genome search, but you will see that later we have to change for the quality assessment of the final annotaiton. 

When it comes to the gene prediction tool, for eukaryotes there are three options: `miniprot`. `augustus` and `metaeuk`. Each of this programs is recommended for a specific situation, being miniprot the default in eukaryote genomes for its speed and accuracy. In our case, we will use miniprot, for the sake of time, as it is a gene mapper, rather than predictior. Since we only want to find the core genes from the list, miniprot will be more than enough. After running BUSCO, ignore the quality results of the genome, as we are only interested in the genes that were found by the tool. The sequences of those genes can be directly found within BUSCO intermediate files. To gather all the found sequences and change the genes. 

```bash
busco command
busco_gather_aa.py <busco_dir>
```

# 1.2 Gene set creation

For this following part, once we have all the protein sequences gathered in one file, we will have to cluster them to eliminate redundancy. Augustus will . To achieve this, we will use  `cd-hit` to create gene clusters with 80% of similarity, and selecting the cluster representative for the final gene set. These gene set then will be turned to GeneBank format, so Augustus can work with it propperly. One key aspect in this transformation is adding a flanking region to the genes. This is done so Augustus has more context about the region that surrounds the genes, both upstream and downstream. You want to capture information about how the intergenic region looks like, but without accidentally including coding regions in there. This might be hard to estimate some times and dependant on the genome (due to the gene content). To play it safe, we will use a flanking region of 1000bp, as we saw in [Paniagua et al. 2025](https://genome.cshlp.org/content/early/2025/03/04/gr.279864.124.long)

```bash
# Clustering
cd-hit -o $dir/complete_buscos.cdhit -c 0.8 -i {input} -p 1 -d 0 -T 4 -M 48000 
grep ">" $dir/complete_buscos.cdhit | cut -f2 -d">" | cut -f1 > {output}

# Concatenating the genes and converting to GeneBank format
concatenate_GFF.py <input>
flanking_region=1000
gff2gbSmallDNA.pl {input.gff} {input.genome} $flanking_region {output} &> {log}

```

It has been observed that AUGUSTUS may suffer from having too many genes in the training set. In these cases, the more you feed the model is not always the merrier. As it is usually said "trash in, trahs out". That is why, for instance, we only search for highly conserved genes to do the training. Also, there is a limit of genes where the model reaches a plateau in its prediction, which is dependant on the number of genes used for training. In our most recent assessment ([Paniagua et al., 2025](https://genome.cshlp.org/content/early/2025/03/04/gr.279864.124.long)), it was determined that a flanking size of 1000bp and more than 2000 genes yielded great results, with 5000 genes being the best option (anything above that did not improve the results). That parameter can be changed in the following script:

```bash
subset_genes.py <gene_file> <gene_number>
```

## 1.3 Augustus training

Once the gene set is prepared, lets dive into AUGUSTUS. Augustus has some models already precomputed using their curated set of genes, and this models are ready to use after installing Augusutus in your machine. You can check those models available by running `augustus --species=help`. **❓Trivia: Which species has the most models available?**

<details><summary>Solution</summary>
_Coprinus cinereus_ has 4 models 
</details>

In our case, we are going to build our own model. For that, first we will have to create the species profile in Augustus. That can be done automatically by running:

```bash
new_species.pl --species=<species_name>
```

With that species created, we can go ahead and run the initial training of the model:

```bash
etraining --species=<species_name> <input_genes> > <training_results>
```
This initial trainign will update the parameters on the profile of our species under `$AUGUSTUS_CONFIG_PATH/species/<species_name>`. List the diretory and take a look at the files, specifically at the metaparameters. <Insert explanation about augustus metaparemeters, and how they are not that important righ now>. We will have to modify some parameters now to account for the "bad genes". These are those genes that might have a premature stop codon in their sequence. We will get rid of them and update the model so Augustus is more accurate. The first step is to find out which genes had a premature stop codon and then, eliminate them from the input that was given to Augustus training step and retrain the model. If we run augusuts `etraining` again, the old model and its parameter will be updated and overwritten.

```bash
grep 'in sequence' {input} | cut -f7 -d' ' | sed s/://g | sort -u > {output}
fitlerGenes.pl <bad_genes> <good_genes> > <filtered_genes>
etraining --species={params.name} {input} > {output}
```
Finally, we will modify the stop codon frequency

```bash
tail -6 <second_model> | head -3 > <freq.out>
modify_SC_freq.py
```


# 2. _Ab initio_ prediction

Now that the gene model for our target genome has been created and curated, we will perform the _ab initio_ prediction of the annotaiton. This means that we will only use the gene model to predict the genes and create an structural annotation. Augutus _ab inito_ prediction works by employing a **generalized Hidden Markov Model**, statistically modelling the structure of genes within the genomic sequence, using only the provided gene model as reference.

<details>
<summary>Breakdown of the AUGUSTUS process</summary>

1. **States represent genomic features**: The GHMM defines different states that correspond to various parts of a gene and the surrounding genomic regions, such as exons, introns, start codons, stop codons, splice sites, and intergenic regions.

2. **Emission probabilities**: Each state in the model is associated with emission probabilities, which define the likelihood of that state producing a particular DNA sequence. For example, coding exon states will have probabilities that reflect the typical nucleotide composition of coding sequences.

3. **Transition probabilities**: The model also includes transition probabilities, which determine the likelihood of moving from one state to another in a biologically meaningful way. For instance, an intron state is likely to be followed by a splice site state.

4. **Finding the optimal parse**: Given an input DNA sequence, AUGUSTUS uses algorithms (similar to the Viterbi algorithm) to find the most probable sequence of states, also known as a parse, that could have generated that sequence. This optimal parse represents the predicted gene structure, identifying the locations and boundaries of genes, exons, and introns.

5. **Statistical modeling of gene components**: AUGUSTUS probabilistically models various characteristics of genes, including the sequence around splice sites, the branch point region, the start and stop codons, the length distribution of exons and introns, and even the distribution of the number of exons per gene. This allows the model to evaluate the likelihood of different potential gene structures based on learned statistical patterns from training data.

6. **No external evidence in ab initio prediction**: In its pure ab initio mode, AUGUSTUS relies solely on the information present within the input DNA sequence and the statistical models it has been trained on to predict gene structures.

Essentially, AUGUSTUS uses a sophisticated statistical framework to "decode" the genomic sequence and identify the most likely arrangement of gene components based on inherent sequence features and learned patterns.
</details><br>

Now, lets dive into the prediction itself:

```bash
augustus --species={params.name} {input.genome} --protein=on --codingseq=on > {output}
```

We specify the options of `protein` and `codingseq`, so Augustus will include them as commentaries in the text and can be later extracted to generate the proteome for the quality control. 

❓**Trivia: How many genes have been predicted?**  
<details><summary>Solution</summary> 
Include number
</details>

## 2.1 Quality Control

We will use three tools for the quality control of the annotation created by Augustus. These tools will be BUSCO, Omark and AGAT. 

1. [AGAT](https://github.com/NBISweden/AGAT) is a suite of scripts and modules that  has multiple functinons, which can be obtained via `agat --tools`. In our case, it will help us extract the basic information about the annotation, such as the number of genes. This can be run directly on the gtf COMPLETE AGAT ONCE IT IS INCLUDED IN THE PIPELIN

2. [BUSCO](https://busco.ezlab.org/busco_userguide.html#protein-mode) run in protein mode will give a fast approximation of the number of core genes that were predicted, offering a measure of completness to the annotation. 

3. [OMARk](https://github.com/DessimozLab/OMArk) is a software similar to BUSCO, since it produces a quality assessment of the proteome based on the completeness and consitency. However, it has a twist, as it is able to detect contamination from closely related species.


However, before running OMARk and BUSCO, we need to prepare the proteome from Augustus output. In order to do so, Augustus has some predefined perl script that extracts the protein sequences from the main file (thus the `protein=on` flag in the prediction).

```bash
getAnnoFasta.pl {input.annot} --seqfile={input.reference}
```
### Running OMARk

Running OMARk requires two steps (three if the database has not been downloaded). Similar to Augustus, you need to select the database that suits your genome the best. 

**:question: Trivia: Which database should we select?**
<details><summary>Solution</summary>
To be filled
</details><br>

OMARk relies on knowing where each protein in the query proteome maps within the precomputed gene families (HOGs) of the OMA database. To achieve this, [OMAmer](https://github.com/DessimozLab/omamer) first assigns proteins to HOGs using a fast, alignment-free k-mer-based method that compares the k-mer content of query proteins to those in the OMA database. This preprocessing step is essential for OMArk to analyze homologous relationships and taxonomic origins, enabling it to assess completeness, consistency, and detect contamination efficiently. For that reason, first we will run OMAmer before proceeding with OMARk.

```bash
# Download the database (optional)
wget <database_url> -O {output}
# OMAmer
omamer --db {input.omark_db} --query {input.proteome} --out {output}
omark -f {input.omamer} -d {input.omark_db} -o $(dirname {output})
```

### Running BUSCO

Busco is much more simple to run. It needs the same lineage and paramters as we used with the whole genome, but this time we will use the proteome as input. 

**Try to do it on your own :wink:**

`
<details>
<summary>BUSCO command</summary>
```bash
busco -i {input.proteome} -o {output} -l {params.lineage} \
            -m proteins -c {threads} --force --download_path {dir.busco_dir}
```
</details><br>

### Running AGAT

TO BE FILLED

# 3. _Evidence driven_ annotation

In this final section, lnc-RNA sequencing data will be integrated as part of the prediction step to improve the quality of the annotation. We will incorporate this data as hints, which will be used to guide AUGUSTUS to favor gene models consistent with the experimental data, resulting in more reliable structural and functional annotation, particularly for complex eukaryotic genomes and non-model organisms with limited existing genomic information. Pre-processing of lr-RNA seq data, such as transcript assembly and filtering, is often necessary to generate high-confidence evidence for AUGUSTUS.

## 3.1 Preprocessing of lr-RNA seq data

We will assume that we have a final file with the processed RNA long reads. The first step will be to map our reads to the reference genome, and produce the initial transcriptome. For that, we will use [IsoQuant](https://github.com/ablab/IsoQuant)

```bash

```