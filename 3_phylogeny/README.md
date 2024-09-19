# Phylogenetic inference

## Plastid phylogeny

1. DNA sequences of individual plastid genes were extracted from the genBank annotation using `get_mtgenes_from_gb.py`.

2. Alignments were generated using MAFFT-linsi (available in the `Raw` folder) and manually inspected in Geneious (available in the `Round3_pasta_aln` folder).

3. A subset of 25 plastid genes with conserved sequences in holoparasites were used to produce a concatenation-based phylogeny in IQTREE.
```
iqtree -s ptG25_sp44_pasta.fasta -p ptG25_sp44_pasta.nex -m MFP+MERGE
```

## Mitochondrial phylogeny

Similar to the plastid phylogeny pipeline. See the folder `mt_aln_genetrees` for individual gene alignments and gene trees.
