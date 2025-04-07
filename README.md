# Sructural annotation course

This tutorial will guide you through the porcess of creating an structural annotation using a reference genome and long reads RNA data from a PacBio sequencing experiment. 


# 1. Gene model training

The first step in the workflow, where our tutorial begins, is to create a set of high confidence genes, with their CDS known. This might come from many sources, such as another version of our organism, a closely related species or a prediction on our own genome. For the sake of this tutorial, and choosing the most complex case, we will use a set of genes that comes from our own genome. Eventhough this migth seem redundandt, we will first predict a set of high confidence genes in order to create a gene model to do the _evidence-driven_ prediction. To fulfill this goal, we will use **BUSCO.** BUSCO is a tool that was designed to produce an assessment of the quality of an assembled genome or proteome, by determning how many genes are present from a set of core genes. BUSCO has multiple datasets that contain genes that should be present in all the species that are under the evolutionary point of a certain bracnh of a phylogeitc tree. Thus, it is important to select among the datasets that busco has available the one that suits each species the most. That can be done through the followign command:

```bash
busco --list-dataset
```

This command will produce a list of the different datasets that BUSCO has available in its database. You will have to select the closest to you sample's taxonomy. 

## Question -> Which dataset would you select for the tutorial sample?
<details><summary>Click this!</summary>
Here is some more text that was hidden before.
</details>
