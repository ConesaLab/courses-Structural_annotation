# Sructural annotation course

This tutorial will guide you through the porcess of creating an structural annotation using a reference genome and long reads RNA data from a PacBio sequencing experiment. 


# 1. Gene model training

## 1.1 Preprocessing and gene set creation

The first step in the workflow, where our tutorial begins, is to create a set of high confidence genes, with their CDS known. This might come from many sources, such as another version of our organism, a closely related species or a prediction on our own genome. For the sake of this tutorial, and choosing the most complex case, we will use a set of genes that comes from our own genome. Eventhough this migth seem redundandt, we will first predict a set of high confidence genes in order to create a gene model to do the _evidence-driven_ prediction. To fulfill this goal, we will use **BUSCO.** BUSCO is a tool that was designed to produce an assessment of the quality of an assembled genome or proteome, by determning how many genes are present from a set of core genes. BUSCO has multiple datasets that contain genes that should be present in all the species that are under the evolutionary point of a certain bracnh of a phylogeitc tree. Thus, it is important to select among the datasets that busco has available the one that suits each species the most. That can be done through the followign command:

```bash
busco --list-dataset
```

This command will produce a list of the different datasets that BUSCO has available in its database. You will have to select the closest to you sample's taxonomy. 

## Question -> Which dataset would you select for the tutorial sample?
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


