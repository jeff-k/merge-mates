# merge-mates

![overlap with read-through](doc/readthru.png)

Merge mated read pairs from next-generation sequencing data using
a _semiglobal_ overlap alignment without gap penalties.

## Synopsis

```bash
merge-mates sample-r1.fastq sample-r2.fastq -o output_prefix
```

### Mating

Mating determines the position within the read that optimally represents
overlap, if it exists. There are three possible conditions:

* **Overlap**, where the fragment size is less than twice the read length.
  These reads can be mated into a single fragment.
* **Read-through**, when the fragment size is less than the read length.
  In these cases overhanging bases on the 3' ends will not be meaningful.
* **Non-overlap**, where the fragment is size is greater than the read
  length. These reads cannot be mated. 

These conditions are represented with a signed Option type (negative
overlap represents read-through).

### Merging

Merging two overlapping reads requires a strategy to mend their
disagreements. For next-gen data, there will also be quality strings to
merge.

## Application Notes

### Reporting and QC

merge-mates may be run with the `-s` switch to additionally output summary
statistics to stderr.

TODO: output summary of adapter sequences left over from read through
reads.

### Platform dependent options

The utility of merge-mates increases with read length. 150bp paired end
reads will not produce interesting results.

TODO: allow manual configuration of read direction.

## Related tools

Merging read pairs is an established practice for paired-end next
generation sequencing data. A summary of existing tools is provided here:

* PEAR (2013)
* COPE (2012)
* PANDAseq (2012)
* FLASH (2011)
* SHE-RA (2010)

Others: BIPES, iTAG.
