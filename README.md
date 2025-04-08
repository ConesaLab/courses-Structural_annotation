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

 Question -> Which dataset would you select for the tutorial sample?
<details><summary>Solution</summary>
Here is some more text that was hidden before.
</details>

Once the dataset is selected, BUSCO has multiple modes to be run in, with the main difference being the software that will do the gene search. There are three options: `miniprot`. `augustus` and `metaeuk`. Each of this programs is recommended for a specific situation, being miniprot the default in eukaryote genomes for its speed and accuracy. In our case, we will use miniprot, for the sake of time. Since we want to find asdjflkas. Finally, we will run BUSCO in genome mode, as we want to find gene sin the whole genome. 

```bash
busco command
```

After running BUSCO, ignore the quality results of the genome, as we are only interested in the genes that were found by the tool. The sequences of those genes can be directly obtained from BUSCO intermediate files. TO gather tyhem all and change the genes to a format that later AUGUSTUS can use, you will have to run the following script:

```busco
busco_gather_aa.py <busco_dir>
```

This script basically gathers all the found genes protein sequences that were saved in individual files into one big file with all the genes that were found. Then, using `cd-hit` we will create clusters of genes, to reduce the redundancy of genes that might be hihgly similar, thus avoiding overfitting the model with a specific gene type. To the cluster representatives, some processig is applied to obtain a gene bank file. This is similar to a gff in nature, but with small style modifications, so AUGUSTUS is able to work with it. It essentially indicates the start and end of each gene, as well as the CDS (CHECK THIS PART). Here comes one of the parts that can be modified based on the experiment, or tweaked in case you are not happy with the prediction. (INCLUDE THE FLANKING REGION)

It has been observed that AUGUSTUS may suffer from having too many genes in the training set. In these cases, the more you feed the model is not always the merrier. As it is usually said "trash in, trahs out". That is why, for instance, we only search for highly conserved genes to do the training. Also, there is a limit of genes where the model reaches a plateau in its prediction, which is dependant on the number of genes used for training. In our most recent assessment (Paniagua et al., 2025)[https://genome.cshlp.org/content/early/2025/03/04/gr.279864.124.long], it was determined that a flanking size of 1000bp and more than 2000 genes yielded great results, with 5000 genes being the best option (anything above that did not improve the results). That parameter can be changed in the following script:

```bash
subset_genes.py <gene_file> <gene_number>
```

## 1.2 Augustus training

Once the gene set is prepared, lets dive into AUGUSTUS. Augustus has some models already precomputed using their curated set of genes, and this models are ready to use after installing Augusutus in your machine. You can check those models available by running `augustus --species=help`. **Trivia: Which species has the most models available?**

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
This initial trainign will update the parameters on the profile of our species under `$AUGUSTUS_CONFIG_PATH/species/<species_name>`. List the diretory and take a look at the files, specifically at the metaparameters. <INsert explanation about augustus metaparemeters, and how they are not that important righ now>. We will have to modify some parameters now to account for the "bad genes". These are those genes that might have a premature stop codon in their sequence. We will get rid of them and update the model so Augustus is more accurate. The first step is to find out which genes had a premature stop codon and then, eliminate them from the input that was given to Augustus training step and retrain the model. If we run augusuts `etraining` again, the old model and its parameter will be updated and overwritten.

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


