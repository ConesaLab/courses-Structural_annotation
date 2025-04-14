# Structural annotation course

This tutorial will guide you through the process of creating an structural annotation using a reference genome and long reads RNA data from a PacBio sequencing experiment. 

<details>
<summary> 📖 Theoretical background</summary>

The first step in understanding a genome typically involves **structural annotation**—the process of identifying protein-coding genes and their associated features. One core method used in this phase is **_ab initio_ gene prediction**, which relies solely on the genomic sequence itself. This approach uses statistical models, such as **Hidden Markov Models (HMMs)**, trained to detect signal sensors—including **splice sites**, **start codons**, and **stop codons**-as well as content sensors, like **codon usage patterns** that are characteristic of coding regions.

While **_ab initio_ prediction** is valuable because it does not require prior experimental data, it often struggles with **lower accuracy**, especially in capturing **complete gene structures**, **untranslated regions (UTRs)**, and **alternative isoforms**. Moreover, its effectiveness is heavily dependent on the availability of **species-specific training models**.

To address these limitations, **evidence-based genome annotation** combines experimental data with _ab initio_ predictions to enhance accuracy. Among the most powerful sources of evidence is **long-read RNA sequencing (lr-RNA seq)**. Unlike traditional short-read RNA-seq, lr-RNA seq can capture **full-length transcripts**, offering direct insights into exon-intron boundaries, alternative splicing, transcription start and end sites, and UTRs.

When used as extrinsic evidence or **"hints"**, lr-RNA seq data can dramatically improve the performance of _ab initio_ gene finders like AUGUSTUS. This integration results in **more accurate and complete gene models**, enabling the identification of **novel isoforms** and providing a deeper understanding of the transcriptome. Such an approach is especially valuable for the **annotation of non-model organisms**, where genomic resources are often limited.

</details><br>

We will use a range of tools to obtain both an _ab initio_ and an evidence-driven annotation. The tools that we will use are:

1. **BUSCO** to produce a set of high-confidence genes and assess the completeness of the genome and the annotation
2. **AUGUSTUS** to generate the gene model and both annotations
3. **AGAT** to extract the proteome from the annotation
4. **cd-hit** will clusters the high-confidence genes to eliminate redundancy
5. **IsoQuant** to produce the lr-RNA seq transcriptome from the raw lr-RNA seq data
6. **OMARk** to assess the completeness and consistency of the proteome
7. **SQANTI3** to filter the raw transcriptome and eliminate low quality isoforms

## Table of contents

- [0. Prerequisites](#0-prerequisites)
  - [Software installation](#software-installation)
  - [Data download](#data-download)
- [1. Gene model training](#1-gene-model-training)
  - [1.1 Preprocessing](#11-preprocessing)  
  - [1.2 Gene model creation](#12-gene-model-creation)
  - [1.3 Augustus training](#13-augustus-training)
- [2. _Ab initio_ prediction](#2-ab-initio-prediction)
  - [2.1 Quality Control](#21-quality-control)
- [3. _Evidence driven_ annotation](#3-evidence-driven-annotation)
  - [3.1 Preprocessing of lr-RNA seq data](#31-preprocessing-of-lr-rna-seq-data)
  - [3.2 Hint creation](#32-hint-creation)
  - [3.3 Final evidence-driven annotation](#33-final-evidence-driven-annotation)



# 0. Prerequisites 

Before starting the tutorial, it is key to have a clean and organized working environment. The initial step, even before processing any data is to prepare the working environment. In bioinformatics, an organized workspace is vital, so when you come after some time to your project, you can find and understand what you were doing, rather than spend hours searching through weirdly named directories. It is important to always create three directories:


- scripts: all the scripts will be stored here, with meaningful names
- data: Raw data will go in here and, if you want and need, databases
- results: Create a sub directory for every different process you do. If you run a process multiple times with different parameters, include them in the directory name, so you will differentiate them in the future.

```bash
mkdir scripts
mkdir data
```

### Software installation

Most of the tools can be directly installed using conda (or mamba), with the exception of SQANTI3, which needs to be downloaded from GitHub and the scripts added to the path. This is simple and can be done by following [this tutorial](https://github.com/ConesaLab/SQANTI3/wiki/Dependencies-and-installation).

<details><summary>🛠️ SQANTI3 Installation </summary>

```bash
wget https://github.com/ConesaLab/SQANTI3/releases/download/v5.3.6/SQANTI3_v5.3.6.zip
mkdir -p tools/sqanti3
unzip SQANTI3_v5.3.6.zip -d tools/sqanti3
# Install the SQANTI3 conda environment
conda env create -f tools/sqanti3/SQANTI3.conda_env.yml
```
</details><br>

For the rest of the tools, you can install their environments, which you will find in the directory `tools/conda_envs`. You can install them by running the following command:

```bash
mamba create env -n <tool_name> -f tools/conda_envs/<tool_name>.yml
```

### Data download

For the tutorial, we will use a mouse 
> TODO: Fill this once the dataset is decided upon


# 1. Gene model training

In any annotation pipeline, the first and most important task is to decide on the **gene model** that will be used by the prediction tool. AUGUSTUS has multiple gene models that can be used to produce an initial annotation on a given genome. However, this models are tailored to specific species, being adapted to their genomes and gene content. Thus, in the case of having a species that is not present within Augustus model, we need to create our own gene model. In order to achieve this, in this tutorial we will create two main resources: 

1. A gene model from a set of high-confidence genes
2. Hints derived from an lr-RNASeq experiment

The final structural annotation will be run through different tools, like the previously mentioned BUSCO, AGAT and OMARk. These tools will be used to assess the quality of the annotation, both in terms of completeness, consistency and number of genes.

## 1.1 Preprocessing

The first step in the workflow, and where our tutorial begins, is to create a set of high confidence genes, with known their coding sequence (CDS).These genes can originate from various sources, such as a previous version of the organism's genome, a closely related species, or predictions based on the genome under study. For the purposes of this tutorial, we will tackle the most complex scenario: deriving high-confidence genes directly from our own genome. While this may appear redundant, this step is essential to build an initial gene model that will serve as a foundation for evidence-driven gene prediction.

To achieve this, we will utilize **BUSCO** (Benchmarking Universal Single-Copy Orthologs), a widely used tool designed to assess the completeness of genome assemblies and annotated gene sets. BUSCO identifies genes based on their evolutionary conservation as single-copy orthologs across specific phylogenetic branches. These conserved genes are expected to be present in all species within a lineage and typically exist as single copies. We can take advantage of this well conserved genes to produce an initial training set for AUGUSTUS based on the single-copy core genes that we find in our genome.

BUSCO offers multiple lineage-specific datasets tailored to different lineages of organisms. It is highly important to select the dataset that is the closest to our target species, to have as many genes as possible in the initial search. In order to check which datasets are available in BUSCO, and the one that is the most suitable for a given experiment, you can run the following command:

```bash
busco --list-dataset
```

This command will produce a list of the different datasets that BUSCO has available in its database. You will have to select the closest to your sample's taxonomy. 

❓**Trivia: Which dataset would you select for the tutorial sample?**  
<details><summary>Solution</summary>
Here is some more text that was hidden before.S
</details><br>

Once the dataset is selected, BUSCO has multiple modes to be run in, and different software that will do the gene search. First of all, BUSCO can be run in `genome`, `proteome` or `transcriptome` mode, which is dependant on the type of input given. In our case, we will use the `genome` mode, as we want to do a full genome search, but you will see that later we have to change for the quality assessment of the final annotation. 

When it comes to the gene prediction tool, for eukaryotes there are three options: `miniprot`. `augustus` and `metaeuk`. Each of this programs is recommended for a specific situation, being miniprot the default in eukaryote genomes for its speed and accuracy. In our case, we will use miniprot, for the sake of time and resources, as it is a gene mapper, rather than a predictor. Since we only want to find the core genes from the list, miniprot will be more than enough. After running BUSCO, ignore the quality results of the genome, as we are only interested in the genes that were found by the tool. The sequences of those genes can be directly found within BUSCO intermediate files. To gather all the found sequences and change the genes. 

```bash
busco command
busco_gather_aa.py <busco_dir>
```

# 1.2 Gene model creation

For this following part, once we have all the protein sequences gathered in one file, we will have to cluster them to eliminate redundancy. Augustus will . To achieve this, we will use  `cd-hit` to create gene clusters with 80% of similarity, and selecting the cluster representative for the final gene set. These gene set then will be turned to GeneBank format, so Augustus can work with it properly. One key aspect in this transformation is adding a flanking region to the genes. This is done so Augustus has more context about the region that surrounds the genes, both upstream and downstream. You want to capture information about how the intergenic region looks like, but without accidentally including coding regions in there. This might be hard to estimate some times and dependant on the genome (due to the gene content). To play it safe, we will use a flanking region of 1000bp, as we saw in [Paniagua et al. 2025](https://genome.cshlp.org/content/early/2025/03/04/gr.279864.124.long)

```bash
# Clustering
cd-hit -o $dir/complete_buscos.cdhit -c 0.8 -i {input} -p 1 -d 0 -T 4 -M 48000 
grep ">" $dir/complete_buscos.cdhit | cut -f2 -d">" | cut -f1 > {output}

# Concatenating the genes and converting to GeneBank format
concatenate_GFF.py <input>
flanking_region=1000
gff2gbSmallDNA.pl {input.gff} {input.genome} $flanking_region {output} &> {log}

```

It turns out that when it comes to training AUGUSTUS, more is not always better. In fact, having too many genes in the training set can actually hurt performance — classic case of "trash in, trash out." That’s why we focus on selecting only highly conserved genes for training, to make sure we’re feeding the model with the highest quality data. On top of that, there’s a point where adding more genes doesn’t really improve predictions — the model sort of hits a plateau. In our latest evaluation ([Paniagua et al., 2025](https://genome.cshlp.org/content/early/2025/03/04/gr.279864.124.long)), we found that using a flanking size of 1000 bp and over 2000 genes worked really well, with 5000 genes being the sweet spot. Beyond that, adding more genes didn’t make much of a difference. You can adjust this parameter in the following script:

```bash
subset_genes.py <gene_file> <gene_number>
```

## 1.3 Augustus training

Once the final gene set is ready, lets dive into AUGUSTUS. Augustus has some models that already are precomputed using their curated set of genes, and this models are ready to use after installing Augusutus in your machine. You can check those models available by running `augustus --species=help`. **❓Trivia: Which species has the most models available?**

<details><summary>Solution</summary>
_Coprinus cinereus_ has 4 models 
</details>

In our case, we are going to build our own model. For that, first we will have to create the species profile in Augustus. That can be done automatically by running:

```bash
new_species.pl --species=<species_name>
```

With the profile created, we can go ahead and run the initial training of the model:

```bash
etraining --species=<species_name> <input_genes> > <training_results>
```
This initial training will update the parameters on the profile of our species under `$AUGUSTUS_CONFIG_PATH/species/<species_name>`. List the directory and take a look at the files, specifically at the metaparameters. <Insert explanation about augustus metaparemeters, and how they are not that important right now>. We will have to modify some parameters now to account for the "bad genes". These are those genes that might have a premature stop codon in their sequence. We will get rid of them and update the model so Augustus is more accurate. The first step is to find out which genes had a premature stop codon and then, eliminate them from the input that was given to Augustus training step and retrain the model. If we run Augusuts `etraining` again, the old model and its parameter will be updated and overwritten.

```bash
grep 'in sequence' {input} | cut -f7 -d' ' | sed s/://g | sort -u > {output}
filterGenes.pl <bad_genes> <good_genes> > <filtered_genes>
etraining --species={params.name} {input} > {output}
```
Finally, we will modify the stop codon frequency

```bash
tail -6 <second_model> | head -3 > <freq.out>
modify_SC_freq.py
```


# 2. _Ab initio_ prediction

Now that the gene model for our target genome has been created and curated, we will perform the _ab initio_ prediction of the annotation. This means that we will only use the gene model to predict the genes and create an structural annotation. Augustus _ab initio_ prediction works by employing a **generalized Hidden Markov Model**, statistically modelling the structure of genes within the genomic sequence, using only the provided gene model as reference.

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

1. [AGAT](https://github.com/NBISweden/AGAT) is a suite of scripts and modules that  has multiple functions, which can be obtained via `agat --tools`. In our case, it will help us extract the basic information about the annotation, such as the number of genes. 

2. [BUSCO](https://busco.ezlab.org/busco_userguide.html#protein-mode) run in protein mode will give a fast approximation of the number of core genes that were predicted, offering a measure of completeness to the annotation. 

3. [OMARk](https://github.com/DessimozLab/OMArk) is a software similar to BUSCO, since it produces a quality assessment of the proteome based on the completeness and consitency. However, it has a twist, as it is able to detect contamination from closely related species.

### Running AGAT

We will use two scripts from AGAT: `agat_convert_sp_gxf2gxf.pl` and `agat_sp_statistics.pl`. The first one will convert the Augustus output to GFF3 format, which is the standard format for an annotation, eliminating all the extra information that Augustus adds to the output. The second one will give us a summary of the annotation, including the number of genes, exons, introns and other features. 

```bash
agat_convert_sp_gxf2gxf.pl -g {input.annot} -o {output.gff3}
agat_sp_statistics.pl --gff {input.gff3} -o {output}
```

Take a look at AGAT's statistics, as it will give you a summary of the annotation. **:question: Trivia: How many genes have been predicted?**
>TODO: Include the number of genes predicted
<details><summary>Solution</summary>
Include number
</details><br>


Before running OMARk and BUSCO, we need to prepare the proteome from Augustus output. In order to do so, Augustus has some predefined perl script that extracts the protein sequences from the main file (thus the `protein=on` flag in the prediction).

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

Busco is much more simple to run. It needs the same lineage and parameters as we used with the whole genome, but this time we will use the proteome as input. 

**Try to do it on your own :wink:**

<details>
<summary>BUSCO command</summary>
```bash
busco -i {input.proteome} -o {output} -l {params.lineage} \
            -m proteins -c {threads} --force --download_path {dir.busco_dir}
```
</details><br>

>TODO: Include busco questions, about the percentage of the completeness and that

# 3. _Evidence driven_ annotation

In this final section, lnc-RNA sequencing data will be integrated as part of the prediction step to improve the quality of the annotation. We will incorporate this data as hints, which will be used to guide AUGUSTUS to favor gene models consistent with the experimental data, resulting in more reliable structural and functional annotation, particularly for complex eukaryotic genomes and non-model organisms with limited existing genomic information. Pre-processing of lr-RNA seq data, such as transcript assembly and filtering, is often necessary to generate high-confidence evidence for AUGUSTUS.

## 3.1 Preprocessing of lr-RNA seq data

We will assume that we have a final file with the processed RNA long reads. There are many ways in which we can process them and obtain a final transcriptome. In this course, we will use [IsoQuant](https://github.com/ablab/IsoQuant), as it usually yields the most consistent results, with the cost of misperforming in novel transcripts discovery. However, for our purposes, it will be better to have a more conservative approach. Another pipeline for this task could be [IsoSeq3](https://isoseq.how/), which is more complex, as processes the reads from their initial subread state, shows more diversity and less consistency among biological replicates. 

With this in mind, lets run IsoQuant 😄

```bash
isoquant.py --reference /PATH/TO/reference_genome.fasta \
  --fastq /PATH/TO/sample1.fastq.gz  \
  --data_type pacbio -o OUTPUT_FOLDER
```

With the raw transcriptome in ready, it is time to curate it before use. For that, we will use [SQANTI3](https://github.com/ConesaLab/SQANTI3). SQANTI3 is a tool designed for the quality control, curation and annotation of lon-read transcriptomes, specifically designed for lr-RNA data. The main module that we will use is the Quality Control module (SQANTI3 QC). This module operates by classifying all the isoforms within a transcriptome into one of the possible structural categories. As well, it has the ability to integrate a variety of orthogonal data into the results and the classification process, such as short-read data or CAGE-seq data. If you want to know more about it, you can check out the [wiki](https://github.com/ConesaLab/SQANTI3/wiki)

<details>
<summary><strong>SQANTI3 structural categories</strong></summary>

1. <strong>Full-Splice-Match (FSM):</strong> Transcript models where all splice junctions perfectly match a known reference transcript.  
2. <strong>Incomplete-Splice-Match (ISM):</strong> Transcript models that match a consecutive subset of splice junctions of a known reference transcript.  
3. <strong>Novel-In-Catalog (NIC):</strong> Transcript models with at least one splice junction not present in the reference annotation, but formed by known splice sites.  
4. <strong>Novel-Not-In-Catalog (NNC):</strong> Transcript models with at least one splice junction that utilizes a novel splice site not present in the reference annotation.  
5. <strong>Antisense:</strong> Transcript models that align to a gene locus but on the opposite strand to all annotated transcripts of that gene.  
6. <strong>Fusion:</strong> Transcript models whose exons align to two or more distinct gene loci.  
7. <strong>Genic genomic:</strong> Transcript models composed of exons that align within a gene locus but do not reconstruct any known or novel splicing pattern.  
8. <strong>Intergenic:</strong> Transcript models whose exons align to genomic regions outside of any annotated gene.  

</details><br>
In order to run SQANTI3 we need three mandatory inputs: 

1. Raw transcriptome, either in gtf or fasta format
2. Reference genome, in fasta format
3. Reference annotation, in gtf format

In our case, we already have the transcriptome and the reference genome, and we will use the _ab initio_ annotation as the reference annotation. This is because we do not care right now about the structural categories, but we want the other information that SQANTI3 provides, such as if any isoform has Retrotrasncriptase switching (RT-Switching) or non-canonical junctions. We will use these parameters as the initial filters of our transcriptome.

```bash 
conda activate sqanti3
sqanti3_qc.py results/isoquant/OUT/OUT.transcript_models.gtf \
    $reference_gtf $reference_genome \
    -d results/sqanti -o isoforms \
    -t 5 --report skip 

# SQANTI filtering
sqanti3_filter.py rules results/sqanti/isoforms_classification.txt --gtf results/sqanti/isoforms_corrected.gtf \
    --skip_report -j data/isoform_filter.json \
    -d results/sqanti  -o isoforms 
```

<!---
TODO: Explain a bit SQANTI filter and add the explanation of the json file
-->

In the second part, we use the filter module of SQANTI. This module offers a way to remove potential artifacts from your transcriptome data using a **rules-based approach**. This method relies on you, the user, to define specific criteria based on the attributes of your transcripts as determined by the SQANTI3 QC step.

In order to curate this transcriptome, since we cannot rely on the structural categories, because the refernce comes from an _ab initio_ prediction, we will only filter using external information. In the case of having short read data, it would be used here to select isoforms in which the junctions are supported.

<details><summary><strong>Filter rules</strong></summary>
The filtering rules are defined in a JSON file, which contains the following:

```json
{
  "rest":[{
      "perc_A_downstream_TTS":[0,59],
      "all_canonical":"canonical",
      "RTS_stage":"FALSE",
      "exons": 2
    }
  ]
}
```
### Explicación de las reglas del filtro SQANTI3

Las reglas de filtro que has proporcionado en formato JSON se aplicarán a **todas las categorías estructurales de isoformas para las cuales no se hayan definido reglas específicas**. Esto se debe a que están bajo la clave `"rest"` [1].

Esta configuración define **una única regla** para la categoría "rest". Para que una isoforma sea considerada un **"Isoform"** (y por lo tanto, sea conservada por el filtro), **debe cumplir con todos los siguientes requisitos simultáneamente** [2]:

*   **`"perc_A_downstream_TTS"` debe estar dentro del intervalo de 0 a 59 (inclusive)** [3]. Esto significa que el porcentaje de adeninas en la secuencia genómica inmediatamente posterior al sitio de terminación de la transcripción (TTS) debe ser menor del 60% [4]. Esta condición se utiliza para **filtrar posibles artefactos de intrapriming** [4, 5].
*   **`"all_canonical"` debe ser igual a `"canonical"`** [3]. Esto indica que **todas las uniones de empalme (splice junctions) de la isoforma deben ser canónicas** [6, 7]. Las uniones canónicas son aquellas que siguen las secuencias GT-AG, GC-AG o AT-AC [8].
*   **`"RTS_stage"` debe ser igual a `"FALSE"`** [3]. Esto significa que la isoforma **no debe haber sido marcada como sospechosa de ser un artefacto de cambio de transcriptasa reversa (RT-Switching)** [6, 9].
*   **`"exons"` debe ser igual a `2`** [3]. Esto especifica que la isoforma **debe tener exactamente dos exones**.

En resumen, cualquier isoforma que no pertenezca a una categoría estructural con reglas definidas explícitamente en el archivo JSON **será considerada un artefacto y descartada a menos que cumpla simultáneamente con todas estas cuatro condiciones**: tener un bajo porcentaje de adeninas aguas abajo del TTS, tener todas sus uniones de empalme canónicas, no ser sospechosa de ser un artefacto de RT-Switching y tener exactamente dos exones [2, 10].

Es importante recordar que **si una isoforma tiene un valor faltante (`NA` o similar) en cualquiera de estas columnas, se considerará que no cumple con la condición y será tratada como un artefacto** [10].

</details><br>

**:question: Trivia: How many different genes have been left after filtering?**

<details><summary>Solution</summary>
To be filled
</details><br>

# 3.2 Hint creation

The final step of the transcriptome preprocessing is to generate a hints file, which is the format that Augustus will use to merge the lr-RNA evidence with the previous gene model. 

```bash
# Extract the filtered isoforms
grep -f results/sqanti/isoforms_inclusion-list.txt results/sqanti/isoforms_corrected.gtf.cds.gff > results/sqanti/isoforms_corrected_filtered.gtf.cds.gff

# Hint creation
tmp_dir=results/hints/tmp
mkdir -p $tmp_dir
grep -P  "\t(CDS|exon)\t" results/sqanti/isoforms_corrected_filtered.gtf.cds.gff | gtf2gff.pl --printIntron --out=$tmp_dir/tmp.gff
# Remove gene_id and change transcript id for grp_id
sed -i 's/gene_id[^;]*;//g' $tmp_dir/tmp2.gff
sed -i 's/transcript_id \\"/grp=/g' $tmp_dir/tmp2.gff
# Add the source
cat $tmp_dir/tmp2.gff | sed "s/\\";/;pri=1;src=PB/g" > {output}
rm -r $tmp_dir

```

The hints file creation follows a complex process, which is described in the following steps:

<details>
<summary>Hint creation process</summary>

1. <strong>Filtering the input GTF file:</strong> The script filters the input GTF file to extract only the relevant features (CDS, exon, or intron) based on the specified UTR option. It uses <code>grep</code> to search for lines containing the specified feature types and then converts the GTF format to GFF format using <code>gtf2gff.pl</code>. The output is stored in a temporary directory.  

2. <strong>Removing gene_id and changing transcript_id:</strong> The script removes the <code>gene_id</code> field from the GFF file and replaces the <code>transcript_id</code> field with <code>grp_id</code>. This is done using <code>sed</code> to perform in-place text replacements.  

3. <strong>Adding source information:</strong> The script adds source information to the GFF file by replacing the <code>"</code> character with a custom string that includes the source (<code>PB</code>) and priority (<code>pri=1</code>). This is done using <code>sed</code> again.  

4. <strong>Outputting the final hints file:</strong> The final hints file is created by redirecting the modified GFF content to the specified output file. The temporary directory is removed afterward.  
</details> <br>


# 3.3 Final evidence-driven annotation

Now that we have the hints file, we can run Augustus again, but this time with the `--hints` flag. This will allow Augustus to use the hints file as a guide for the prediction. 

```bash
augustus --species={params.name} {input.genome} --protein=on --codingseq=on \
  --hintsfile={input.hints} --extrinsicCfgFile={params.cfg} > {output}
```
The `--extrinsicCfgFile` parameter is used to specify the configuration file that contains the parameters for the hints. This file can be copied directly from the Augustus configuration directory and you should only add your source under the `[SOURCES]` section. The hints file will be used to guide the prediction, and the configuration file will specify how to use the hints. In our case, just copy the code below.

<details>
<summary>Configuration file</summary>

```ini
# extrinsic information configuration file for AUGUSTUS
# include with --extrinsicCfgFile=filename
# date: 2025-04-08
# Pablo Atienza (pablo.atienza@csic.es)


# source of extrinsic information:
# M manual anchor (required)
# P protein database hit
# E est database hit
# C combined est/protein database hit
# D Dialign
# R retroposed genes
# T transMapped refSeqs
# PB PacBio (long reads, circular consensus)

[SOURCES]
M RM PB

#
# individual_liability: Only unsatisfiable hints are disregarded. By default this flag is not set
# and the whole hint group is disregarded when one hint in it is unsatisfiable.
# 1group1gene: Try to predict a single gene that covers all hints of a given group. This is relevant for
# hint groups with gaps, e.g. when two ESTs, say 5' and 3', from the same clone align nearby.
#
[SOURCE-PARAMETERS]
PB individual_liability
#   feature        bonus         malus   gradelevelcolumns
#		r+/r-
#
# the gradelevel colums have the following format for each source
# sourcecharacter numscoreclasses boundary    ...  boundary    gradequot  ...  gradequot
# 

[GENERAL]
      start     1          1  M    1  1e+100  RM  1     1    PB    1       1
       stop     1          1  M    1  1e+100  RM  1     1    PB    1       1
        tss     1          1  M    1  1e+100  RM  1     1    PB    1       1
        tts     1          1  M    1  1e+100  RM  1     1    PB    1       1
        ass     1          1  M    1  1e+100  RM  1     1    PB    1       1
        dss     1          1  M    1  1e+100  RM  1     1    PB    1       1
   exonpart     1       0.98  M    1  1e+100  RM  1     1    PB    1       1e5
       exon     1          1  M    1  1e+100  RM  1     1    PB    1       1e10
 intronpart     1          1  M    1  1e+100  RM  1     1    PB    1       1e5  
     intron     1         .1  M    1  1e+100  RM  1     1    PB    1       1e10  
    CDSpart     1          1  M    1  1e+100  RM  1     1    PB    1       1e5
        CDS     1          1  M    1  1e+100  RM  1     1    PB    1       1e15
    UTRpart     1          1  M    1  1e+100  RM  1     1    PB    1       1
        UTR     1          1  M    1  1e+100  RM  1     1    PB    1       1
     irpart     1          1  M    1  1e+100  RM  1     1    PB    1       1
nonexonpart     1          1  M    1  1e+100  RM  1  1.15    PB    1       1
  genicpart     1          1  M    1  1e+100  RM  1     1    PB    1       1


#
# Explanation: see original extrinsic.cfg file
#
```

</details><br>

# 3.4 Quality control

The quality control of the final annotation is the same as before, but this time we will use the new annotation file as input for BUSCO and OMARk. 

Try to do it in your own 😉

<details>
<summary>Click here if desperate (or lazy)</summary>

```bash
# Proteome conversion
getAnnoFasta.pl {input.annot} 
# OMARk
omamer --db {input.omark_db} --query {input.proteome} --out {output}
omark -f {input.omamer} -d {input.omark_db} -o $(dirname {output})
# BUSCO
busco -i {input.proteome} -o {output} -l {params.lineage} \
            -m proteins -c {threads} --force --download_path {dir.busco_dir}
# AGAT
```
</details><br>

# 4. Conclusion

With this, you have reached the end of this tutorial. You have learned how to create a gene model from scratch, how to use it to predict genes in a genome and how to use lr-RNA seq data to improve the prediction. You have also learned how to assess the quality of the annotation using different tools.

Now, compare the results of the first and second annotation. What are the main differences? Do you think that the lr-RNA seq data improved the prediction? Why?

