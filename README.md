# NXF - TGS ONT PIPELINE

### About
Nextflow pipeline for processing raw Oxford Nanopore run data - merge/rename files, generate report, run ONT workflows 
(wf-clone-validation, wf-bacterial-genomes or wf-amplicon).

The required inputs for the pipeline are:
- path to samplesheet - path to a `csv` file with columns `sample`, `barcode` and `user`. Other columns can also be there.
- path to raw ONT data (`path/to/fastq_pass`)

### Running
Only `nextflow` and `docker` are required. The pipelines runs the following steps:
1. `merge_reads` - merge all reads belonging to a barcode and rename according to the provided samplesheet
2. `report` - generate a `csv` and `html` reports about read counts, quality etc. One report per user is generated.
3. `assembly` - perform assembly using one of the epi2me pipelines wf-clone-validation, wf-bacterial-genomes or wf-amplicon (default wf-clone-validation)

```bash
nextflow run angelovangel/nxf-tgs --samplesheet path/to/samplesheet.csv --fastq path/to/fastq_pass
```

It is possible to run only `merge_reads` or `merge_reads` + `report`. For this, use the `-entry` parameter (note the single dash) in nextflow like so:

```bash
nextflow run angelovangel/nxf-tgs --samplesheet path/to/samplesheet.csv --fastq path/to/fastq_pass -entry report
```
This will run `merge_reads` + `report` and no asembly.