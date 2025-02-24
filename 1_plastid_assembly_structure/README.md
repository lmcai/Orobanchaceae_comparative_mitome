# Plastid genome assembly

1. Plastid genomes were assembled using GetOrganelle using the following commands:
```
get_organelle_from_reads.py -1 forward.fq -2 reverse.fq -o prefix -R 15 -k 21,45,65,85,105 -F embplant_pt
```

2. Annotation was completed using Geneious

3. To summarize the length of IR regions and single copy regions, the script `IR_len_from_GenBank.py` was applied to the genbank annotation from Geseq to get the length of IR regions. For holoparasite that IR cannot be easily annotated, reference sequences from `IR_SSR_LSR_reference.fas` was used to identify IR, LSR, and SSR.
