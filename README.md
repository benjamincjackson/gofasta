    
# gofasta




gofasta is a command-line utility developed to handle SARS-CoV-2 alignments, but should be useful more generally for handling consensus genomes of any microbial pathogen. It has been used on datasets of millions of sequences, including by the United Kingdom's [daily SARS-CoV-2 genome processing pipelines](https://github.com/COG-UK), [Pangolin](https://github.com/cov-lineages/pangolin) and [Civet](https://github.com/artic-network/civet).

If you use gofasta in your work, please cite [https://doi.org/10.1093/bioinformatics/btac424](https://doi.org/10.1093/bioinformatics/btac424).

## Installation

You can install gofasta from Bioconda:

`conda install bioconda::gofasta`

If you have [Go](https://go.dev/) installed, you can run `go install github.com/virus-evolution/gofasta@latest` to build a binary of the latest release locally (maybe in `~/go/bin/`). Branches, commits or tags other than `latest` can be built in the same way.

You can also build the contents of this repository:

```
git clone https://github.com/virus-evolution/gofasta.git
cd gofasta
go build -o gofasta
```

## Commands

For a full list of commands and options, you can run `gofasta` with the `-h` flag at the command line: `gofasta -h`,  `gofasta sam -h`, `gofasta sam variants -h`, etc. 

e.g., for `gofasta updown list -h`:

<details><summary><b>Click to see the help message</b></summary>

```
❯ gofasta updown list -h
Generate input CSV files for gofasta updown topranking

Example usage:

	gofasta updown list -r reference.fasta -q alignment.fasta -o mutationlist.csv

Non-ATGC nucleotides are not recommended in the --reference, and --reference and --query must be aligned to
the same thing.

--outfile is a CSV-format file with the columns: query,SNPs,ambiguities,SNPcount,ambcount. There is one row
for each sequence in --query. SNPs is a "|"-delimited list of SNPs relative to --reference. ambiguities is a
"|"-delimited list of ranges (1-based, inclusive) of tracts of ambiguities (anything that isn't ATGC).

Usage:
  gofasta updown list [flags]

Flags:
  -q, --query string     Alignment of sequences to parse, in fasta format (default "stdin")
  -o, --outfile string   Output to write (default "stdout")
  -t, --threads int      Number of threads (and CPUs) to use: 0, the default, for all of them
      --strict           Stop with an error at the first record that can't be processed, instead of warning and skipping it
  -h, --help             help for list

Global Flags:
  -r, --reference string   Reference sequence, in fasta format, which is treated as the root of the imaginary tree
```
</details>

## Usage

Any input file (or stdin) can be gzipped: gofasta recognises gzip by its content, so the file can be called anything, and decompresses it as it reads it. Output is not compressed.

### Sam to fasta format conversion

[minimap2](https://github.com/lh3/minimap2) provides fast and reliable pairwise alignments between SARS-CoV-2 consensus genomes. gofasta provides utility functions to convert the sam format alignments it outputs into fasta format.

Pairwise alignments are less computationally demanding than multiple sequence alignment-proper and they offer an efficient approximation to it when sequences are closely related.

<details><summary><b>Making a single alignment with multiple sequences</b></summary>
</br>

This is the pipeline currently used by [Pangolin](https://github.com/cov-lineages/pangolin):

```
minimap2 -a -x asm20 --sam-hit-only --secondary=no --score-N=0 MN908947.fa unaligned.consensus.fasta -o aligned.sam

gofasta sam toMultiAlign -s aligned.sam --start 266 --end 29674 --pad -o aligned.fasta
```

`MN908947.fa` is the reference genome for SARS-CoV-2 (here is its [Genbank accession](https://www.ncbi.nlm.nih.gov/nuccore/MN908947.3)), and `unaligned.consensus.fasta` contains all your consensus genomes (in one file). In this case we've clipped the alignment to just the coding sequence (`--start 266 --end 29674`) but replaced the trimmed regions with `N`s to retain the reference length in the output (`--pad`). <b>Insertions relative to the reference are discarded from the output file so everything is the same (== reference) length.</b>


If minimap2 splits one of your genomes into several alignments (a primary alignment and supplementary alignments, which share a name in the sam file), gofasta combines them into one sequence. Where two of them overlap on the reference, the alignment of the first one in the file is used for the overlapping stretch (a message on stderr says if the alignments disagree there), so that, for example, an insertion that is present in both is only counted once. Secondary alignments and unmapped records are ignored.

minimap2 and gofasta can both make use of multiple threads. Additionally, minimap2 writes to stdout by default and gofasta reads from stdin if we don't specify a file for `-s`, so we can avoid writing the intermediate sam to disk entirely. Both of these should speed things up:

```
minimap2 -t8 -a -x asm20 --score-N=0 MN908947.fa unaligned.consensus.fasta | gofasta sam toma -t2 > aligned.fasta
```

We give minimap2 more threads here because it's doing more work. `toma` is just an alias for `toMultiAlign`.

<i>But I don't want to have to write all this code every time I want to align something</i>. That's understandable. In which case you could define a shell function in your `~/.zshrc` or `~/.bashrc` file, something like:

```
function align() {
  minimap2 -a -x asm20 --score-N=0 ~/path/to/MN908947.fa $1 | gofasta sam toma
}
```

and then you can run `align unaligned.consensus.fasta > aligned.fasta` in your terminal.

</details>

<details><summary><b>Making pairwise alignments (including insertions)</b></summary>
</br>
gofasta will also convert the minimap2 alignment into pairwise alignments, writing one file per consensus genome, including insertions relative to the reference and writing the reference itself to the output files:

```
minimap2 -a -x asm20 --score-N=0 MN908947.fa unaligned.consensus.fasta |\
  gofasta sam toPairAlign -r MN908947.fa -o outputdir
```

`outputdir` is a directory which will be created if it does not already exist. Individual fasta files will be written to `outputdir/*fasta`. Filenames are derived from the fasta descriptions (with `/` replaced by `_`).

Just as with `toMultiAlign`, you can clip the alignments to coordinates of your choosing. For example, if you are only interested in the Spike region:

```
minimap2 -a -x asm20 --score-N=0 MN908947.fa unaligned.consensus.fasta |\
  gofasta sam topa -r MN908947.fa --start 21563 --end 25384 -o outputdir
```

`--start` and `--end` are 1-based inclusive, and they are <i>in reference coordinates</i>. So you will get the region corresponding to the reference Spike regardless of insertions or deletions anywhere in your query genome.

There is no `--pad` option for `toPairAlign` because we don't expect all the sequences to be reference length.

</details>

### Searching for neighbours by genetic distance

gofasta provides two utilities to search for the closest genetic neighbours of a number of query sequences among a set of target sequences. In both cases the queries are loaded into memory and the targets are streamed from disk, so the target file can be arbitrarily large.

<details><summary><b>Traditional distance measures</b></summary>
</br>

Use `gofasta closest` to find the closest neighbours by traditional genetic distance measures:

```
❯ gofasta closest --help | sed -n '/Flags/,$p'
Flags:
  -t, --threads int       Number of threads (and CPUs) to use: 0, the default, for all of them
      --query string      Alignment of sequences to find neighbours for, in fasta format
      --target string     Alignment of sequences to search for neighbours in, in fasta format
  -m, --measure string    Which distance measure to use (raw, snp or tn93) (default "raw")
  -n, --number int        (Optional) the closest n sequences to each query will be returned
  -d, --max-dist string   (Optional) return all sequences less than or equal to this distance away (inf for no limit)
  -o, --outfile string    The output file to write (default "stdout")
      --table             Write a long-form table of the output
      --strict            Stop with an error at the first record that can't be processed, instead of warning and skipping it
  -h, --help              help for closest
```
 The available distance measures are `raw` (the default) which is the number of nucleotide differences <i>per site</i> ; `snp`, which is the total number of nucleotide differences, and `tn93`, which is [Tamura and Nei's (1993) evolutionary distance](https://academic.oup.com/mbe/article/10/3/512/1016366).

raw distance is calculated as: `l_d / (l_d + l_s)`, where `l_d` is the number of sites which are certainly different between two sequences and `l_s` is the number of sites which are certainly the same.

snp distance = `l_d`. For both the raw and snp distances, ambiguous IUPAC codes are treated as the set of bases that they represent.

tn93 distance is calculated according to equation (7) in [the paper](https://academic.oup.com/mbe/article/10/3/512/1016366). Only `ATGC` bases are considered when calculating this measure.

The routine is parallelised across queries, so there is no point setting `-t` greater than the number of sequences in `--query`.

</details>

<details><summary><b>Directional snp-distance</b></summary>
</br>

Use `gofasta updown topranking` to find the closest neighbours by (directional) snp distance. This is the routine used by [Civet3](https://github.com/artic-network/civet):

```
❯ gofasta updown topranking -h | sed -n '/Flags/,$p'
Flags:
  -q, --query string             File with sequences to find neighbours for. Either the CSV output of gofasta updown list, or an alignment in fasta format
  -t, --target string            File of sequences to look for neighbours in. Either the CSV output of gofasta updown list, or an alignment in fasta format
  -o, --outfile string           CSV-format file of closest neighbours to write (default "stdout")
      --table                    write a long-form table of the output
      --ignore string            Optional plain text file of IDs to ignore in the target file when searching for neighbours
      --dist-all int             Maximum allowed SNP-distance between target and query sequence in any direction. Overrides the settings below
      --dist-up int              Maximum allowed SNP-distance from query for sequences in the parent bin
      --dist-down int            Maximum allowed SNP-distance from query for sequences in the child bin
      --dist-side int            Maximum allowed SNP-distance from query for sequences in the sibling bin
      --size-total int           Max number of neighbours to find (attempts to split equally between same/up/down/side). A hard limit
      --size-up int              Max number of closest parent sequences to find, if size-total not specified. A soft limit unless --no-fill
      --size-down int            Max number of closest child sequences to find, if size-total not specified. A soft limit unless --no-fill
      --size-side int            Max number of closest sibling sequences to find, if size-total not specified. A soft limit unless --no-fill
      --size-same int            Max number of identical sequences to find, if size-total not specified. A soft limit unless --no-fill
      --threshold-pair float32   Up to this proportion of consequential sites is allowed to be ambiguous in either sequence for each pairwise comparison (default 0.1)
      --threshold-target int     Target can have at most this number of ambiguities to be considered (default 10000)
      --dist-push int            Push the --dist boundaries outwards so that bins have at least these many closest SNP-distances for which there are neighbours, where possible
      --no-fill                  Don't make up for a shortfall in any of --size-up, -down, -side or -same by increasing the count for other bins
  -h, --help                     help for topranking

Global Flags:
  -r, --reference string   Reference sequence, in fasta format, which is treated as the root of the imaginary tree
```

This routine has the potential to be faster than traditional genetic distance measures. It uses patterns of derived mutations to make inferences about the likely phylogenetic relationships between closely related sequences. Briefly, sequences are compressed so that they are represented by only their ATGC nucleotide differences from a common reference sequence, which is treated like the root of an imaginary tree. SNPs can either be unique to the query sequence, unique to the target sequence, or present in the intersection of the two. SNPs present in the intersection are interpreted as representing shared ancestry between query and target, whereas an excess of SNPs in either the query or target set can be interpreted to give directionality relative to a root sequence. See [the Civet preprint](https://www.medrxiv.org/content/10.1101/2021.12.13.21267267v1) for a detailed explanation of the scheme (Figure 2 and Figure S1). 

The input `--query` and `--target` files can either be alignments in fasta format, or they can be csv-format files produced by `gofasta updown list` (or one of each). Using the csv-format files should be faster to the extent that they are quicker to read from disk compared to alignments, which initially contain the information for every site.

An example of command-line use and more explanation is available by running `gofasta updown topranking --help`.

</details>

### Annotating mutations

Use `gofasta variants` to extract nucleotide changes and indels relative to a reference sequence from a multiple sequence alignment, and, given an annotation of the reference, amino acid changes as well; `gofasta sam variants` does the same for alignments in sam format (and always needs an annotation).

<details><summary><b>Nucleotide changes</b></summary>
</br>

Without an annotation, `gofasta variants` lists the nucleotide changes, deletions and insertions in a multiple sequence alignment (`--query`) compared to a reference sequence. The reference can be a sequence in a separate fasta file (`--reference`) of the same width as the alignment (i.e. aligned to the same thing), or a record of the alignment (`--reference-id`, its ID).

The basic usage creates a csv file with a header and a line for each sequence in `--query`. The first column is the sequence name, and the second column is a '|'-delimited list of changes: `nuc:` followed by the reference allele, the 1-based position in reference coordinates and the query allele; `del:POS:LEN` for a deletion starting at reference position POS; `ins:POS:LEN` for an insertion after reference position POS. Deletions at the start and end of the alignment (uncovered ends) are not reported.

IUPAC ambiguity codes are treated as the set of bases that they represent, and only certainly-different changes are reported. For example an output of `nuc:A101S` is possible, but `nuc:A101W` is not.
```
❯ gofasta variants -r MN908947.fa -q aligned.fasta -o variants.csv
❯
❯ head -n2 variants.csv
query,mutations
query1,nuc:T670G|nuc:G4184A|nuc:C4321T|nuc:C9344T|nuc:A9424G|nuc:C9534T|nuc:C10198T|nuc:G10447A|nuc:C10449A|nuc:G12160A|nuc:C12880T|nuc:C14408T|nuc:C15714T|nuc:C17410T|nuc:C19955T|nuc:A20055G|nuc:T21570G|nuc:C21618T|nuc:G21987A|nuc:T22200G|nuc:G22578A|nuc:T22679C|nuc:C22686T|nuc:A22688G|nuc:A23403G|nuc:C23525T|nuc:T23599G|nuc:C23604A|nuc:C23854A|nuc:G23948T|nuc:T24469A|nuc:C25000T|nuc:C26060T|nuc:C26270T|nuc:G27382C|nuc:A27383T|nuc:T27384C|nuc:G27788T|nuc:C27807T|nuc:A28271T|nuc:C28311T|nuc:C28724T|nuc:G28881A|nuc:G28882A|nuc:G28883C|nuc:A29510C
```

If you invoke `--aggregate`, the proportion of each change in the whole alignment is written.
```
❯ gofasta variants -r MN908947.fa -q aligned.fasta --aggregate
mutation,frequency
nuc:C44T,0.250000000
nuc:C241T,0.750000000
nuc:T670G,0.916666667
nuc:C1314T,0.083333333
nuc:C1613A,0.083333333
nuc:C1684T,0.083333333
nuc:C2790T,0.833333333
nuc:C3037T,0.833333333
...
```

</details>

<details><summary><b>Amino acid, indel and neutral nucleotide changes</b></summary>
</br>

If you provide an annotation, gofasta can also annotate amino acid changes relative to a reference sequence. The annotation can be provided (with the suffix `.gb`, `.gbk` or `.gbff`) in [genbank flat file format](https://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html) or (with the suffix `.gff` or `.gff3`) [gff version 3 format](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md). Examples of both for SARS-CoV-2 are available under `resources/` in this repository. Give it with `--annotation`, which takes the format from the suffix, or with `--genbank` or `--gff` to say the format whatever the file is called (`--gff stdin` reads it from stdin, if nothing else is).

The two relevant routines are `gofasta variants` (for annotating mutations in alignments in fasta format) and `gofasta sam variants` (for annotating mutations in alignments in sam format). They should give the same output for the same alignment and the same annotation. Multiple sequence alignments in fasta format don't need to be in reference coordinates for `gofasta variants`, but if they aren't, a sequence in the same space as the annotation must be present in the alignment. If the alignment is being read from stdin, this sequence must be the first sequence in the alignment, but doesn't have to be if the file is being read from disk. The reference sequence in fasta format needs to be provided to `gofasta sam variants` unless it is present in your annotation. As usual, run either command with the `-h` flag for example command lines and detailed help.

For a genbank format annotation, the annotation will be parsed such that the genome is split into protein-coding regions based on CDS features, and intergenic regions (everything that isn't in CDS). Mutations are then annotated with `ins` (insertion), `del` (deletion), `aa` (amino acid change) or `nuc` (a nucleotide change in intergenic sequence, or a synonymous change in a CDS). Where CDSs overlap, a change is reported as an `aa` for each CDS in which it changes the amino acid, and also as a `nuc` if it is synonymous in any of them: for example, a change that alters an amino acid in gene A but is synonymous in an overlapping gene B is reported both as `aa:A:...` and as `nuc:...`.

gff format annotation gives you more flexibility for naming amino acid changes. Currently, the annotation will be parsed such that the genome is split into protein-coding regions based on feature lines whose `type` (in column 3) is either `CDS` or `mature_protein_region_of_CDS`, and intergenic regions (everything else). For the purposes of annotating amino acids, `CDS` or `mature_protein_region_of_CDS` feature lines that have a `Name=something` tag,value pair in the attributes column (column 9) will be represented in the output, under that name. Thus you can define regions as protein-coding using a `CDS` feature line (for example orf1a in SARS-CoV-2) but annotate amino acid changes in its constituent protein products using `mature_protein_region_of_CDS` feature lines with `Name=` attributes. [See the example](https://github.com/virus-evolution/gofasta/blob/master/resources/sarscov2-reduced.gff). Other attributes such as `gene=` are not used, so with GFF3 files from NCBI, whose `Name=` is the protein accession, amino acid changes are reported as e.g. `aa:YP_009724390.1:D614G`; edit the `Name=` values (or use the GenBank file) to report gene names instead.

The reading frame of a CDS is given by the phase (column 8) of its 5'-most line in a gff file (the last line of a CDS on the minus strand), or by its `/codon_start` in a genbank file: translation starts that many bases (phase), or one fewer (`/codon_start`), from its 5' end. A CDS that starts with a partial codon in this way may also end with one, which is ignored. A CDS that starts with a complete codon must be a whole number of codons long. Changes in the bases of these partial codons are reported as `nuc`.

Examples of the output formats:

	ins:2028:3 - a 3-base insertion immediately after (1-based) position 2028 in reference coordinates
	del:11288:9 - a 9-base deletion whose first missing nucleotide is at (1-based) position 11288 in reference coordinates
	aa:S:D614G - the amino acid at (1-based) residue 614 in the S gene is a D in the reference and a G in this sequence
	aa:nsp12:P323L - the amino acid at (1-based) residue 323 in the rdrp gene is a P in the reference and an L in this sequence
	nuc:C3037T - the nucleotide at (1-based) position 3037 in reference coordinates is a C in the reference and a T in this sequence
	aa:S:H69- - the amino acid at (1-based) residue 69 in the S gene is deleted in this sequence
	aa:S:214EPE - the amino acids EPE are inserted after (1-based) residue 214 in the S gene in this sequence

In-frame insertions and deletions in a CDS are also described at the amino acid level, residue by residue, in the coordinates of the reference protein: a deleted residue is written with `-`, inserted residues after the number of the residue they follow, and a residue changed next to an indel as a substitution. For example, the 6-base deletion at 22029 in the S gene of the SARS-CoV-2 Delta variant is `aa:S:E156G|del:22029:6|aa:S:F157-|aa:S:R158-`. The amino acid changes are worked out from the codons that the indel touches, as aligned, so they follow the alignment; an indel in a run of repeated residues is placed as far towards the C-terminus as possible. A change in the bases of those codons is part of the amino acid change rather than a separate `nuc` entry, and `--append-snps` lists the nucleotide indel (and any such changes) in the parentheses. Nothing after a stop codon is reported (an insertion that brings in a stop gives, e.g., `aa:N:Y87*`). Indels that shift the reading frame, come after one that does, run past the end of the CDS, touch its start or stop codon, or touch missing data (`N`, or the uncovered ends of an alignment) are reported only as `ins` or `del`.

As without an annotation, the default mode writes a csv with one line per query sequence, and each sequence's mutations in the second column. Use `--aggregate` to get the overall frequencies of mutations in the alignment(s).

So, for example, you can find the frequencies of all the amino acid changes at residue 681 in the Spike gene, and the nucleotide changes underlying them, from the sample of SARS-CoV-2 sequences in `aligned.fasta` like:

```
❯ gofasta variants -q aligned.fasta --annotation MN908947.gb --aggregate --append-snps | grep "^aa:S:P681"
aa:S:P681H(nuc:C23604A),0.004000000
aa:S:P681R(nuc:C23604G),0.983000000
```

or find which sequences have `P681H`:
```
❯ gofasta variants -q aligned.fasta --annotation MN908947.gb | grep "S:P681H" | cut -d, -f1
COGUK/PHEC-XXXX107/PHEC
COGUK/PHEC-XXXX003/PHEC
COGUK/PHEC-XXXX996/PHEC
COGUK/PHEC-XXXX544/PHEC
```

</details>

## Context, limitations and alternatives

Alternatives to minimap2 for pairwise viral genome alignment exist. Notably, [Nextalign](https://github.com/nextstrain/nextclade) [(Aksamentov et al. 2021)](https://joss.theoj.org/papers/10.21105/joss.03773.pdf) can use a genome annotation to apply a reading-frame-aware gap penalty, and will perform translation and amino acid alignment to call amino acid mutations. 

gofasta’s functions for finding the closest genetic neighbours to a set of query sequences were designed for the specific use-case of querying a large dataset of target genomes by genetic distance to place a smaller set of focal sequences in context, without having to infer all the phylogenetic relationships between hundreds of thousands or millions of sequences. An alternative approach is currently provided by an up-to-date a [maximum parsimony phylogenetic tree consisting of publicly available SARS-CoV-2 genomes](https://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/) [(McBroome et al. 2021)](https://academic.oup.com/mbe/article/38/12/5819/6361626) [and associated software](https://github.com/yatisht/usher) which can be used for placing new samples and querying the dataset [(Turakhia et al. 2021)](https://www.nature.com/articles/s41588-021-00862-7).

Other software exists for calculating pairwise genetic distances, including [tn93](https://github.com/veg/tn93) and [pairsnp](https://github.com/gtonkinhill/pairsnp).

## Acknowledgements and third party licences

gofasta incorporates [bíogo hts](https://github.com/biogo/hts), [cobra](https://github.com/spf13/cobra) and [pflag](https://github.com/spf13/pflag), each distributed under its own licence. Their licence notices are reproduced in `NOTICE`; `gofasta licences` prints a summary and `gofasta licences -l` prints them in full.

gofasta uses a slightly modified version of the bit-level coding scheme for nucleotides by Emmanuel Paradis, which is used in [ape](https://doi.org/10.1093/bioinformatics/btg412). For a description of it, see [here](https://emmanuelparadis.github.io/misc/BitLevelCodingScheme.html).

The functions with SAM files as input were written to handle pairwise alignments between assembled SARS-CoV-2 genomes from [minimap2](https://github.com/lh3/minimap2).

