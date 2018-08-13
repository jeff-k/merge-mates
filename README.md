# merge-mates

Merge mated read pairs from next-generation sequencing data

## Synopsis

```bash
merge-mates -f sample-r1.fastq -r sample-r2.fastq -o output_prefix
```

### Mating

Mating determines the position within the read that optimally represents
overlap, if it exists. There are three possible conditions:

* Overlap
* Read-through
* Non-overlap

These conditions are represented with a signed Option type. Negative
overlap represents read-through.

### Merging

Merging two overlapping reads requires a strategy to mend their
disagreements. For next-gen data, there will also be quality strings to
merge.

## Related tools

Merging read pairs is an established practice for paired-end next
generation sequencing data. A summary of existing tools is provided here:

* PEAR (2013)
* COPE (2012)
* PANDAseq (2012)
* FLASH (2011)
* SHE-RA (2010)

Others: BIPES, iTAG.
