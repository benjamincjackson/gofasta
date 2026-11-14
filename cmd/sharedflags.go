package cmd

import (
	"github.com/spf13/cobra"
)

// variantsOutputOptions holds the flags that say which mutations variants and sam variants report, and how
type variantsOutputOptions struct {
	start, end int
	aggregate  bool
	threshold  float64
	appendSNP  bool
	groupAA    bool
}

// addVariantsOutputFlags adds the flags of variantsOutputOptions to c (variants and sam variants), bound to o
func addVariantsOutputFlags(c *cobra.Command, o *variantsOutputOptions) {
	c.Flags().IntVar(&o.start, "start", -1, "Only report variants with a reference base at or after this position")
	c.Flags().IntVar(&o.end, "end", -1, "Only report variants with a reference base at or before this position")
	c.Flags().BoolVar(&o.aggregate, "aggregate", false, "Report the proportions of each change")
	c.Flags().Float64Var(&o.threshold, "threshold", 0.0, "If --aggregate, only report changes with a freq greater than or equal to this value")
	c.Flags().BoolVar(&o.appendSNP, "append-snps", false, "Report the codon's SNPs in parenthesis after each amino acid mutation")
	c.Flags().BoolVar(&o.groupAA, "group-aa", false, "Report the amino acid changes from one indel as one entry (aa:S:HV69-- instead of aa:S:H69-|aa:S:V70-)")
	hideDefault(c, "start", "end")
}

// checkVariantsOutputFlags checks the flags addVariantsOutputFlags adds
func checkVariantsOutputFlags(cmd *cobra.Command) error {
	return firstError(
		floatAtLeast(cmd, 0, "threshold"),
		startNotAfterEnd(cmd),
	)
}

// mutationFormatsHelp describes the mutations variants and sam variants report, for their long help
const mutationFormatsHelp = `Mutations are annotated with ins (insertion), del (deletion), aa (amino acid change) or nuc (a nucleotide
change that isn't in a codon that is represented by an amino acid change). The formats are:

	ins:2028:3 - a 3-base insertion immediately after (1-based) position 2028 in reference coordinates
	del:11288:9 - a 9-base deletion whose first missing nucleotide is at (1-based) position 11288 in reference coordinates
	aa:S:D614G - the amino acid at (1-based) residue 614 in the S gene is a D in the reference and a G in this sequence
	aa:S:H69- - the amino acid at (1-based) residue 69 in the S gene is an H in the reference and is deleted in this sequence
	aa:S:214EPE - the amino acids EPE are inserted immediately after (1-based) residue 214 in the S gene in this sequence
	nuc:C3037T - the nucleotide at (1-based) position 3037 in reference coordinates is a C in the reference and a T in this sequence

With --group-aa, the amino acid changes from one indel are one entry: aa:S:HV69-- for aa:S:H69-|aa:S:V70-.

An indel that shifts the reading frame of a coding sequence is reported only as an indel, as are the indels
after it in that coding sequence, unless indels within 3 codons of it restore the frame. Changes of single
bases are still reported codon by codon.
`

// trimOptions holds the flags with which sam toMultiAlign and sam toPairAlign trim and wrap their output
type trimOptions struct {
	start, end int
	wrap       int
}

// addTrimFlags adds --start, --end and --wrap to c (sam toMultiAlign, with --pad, and sam toPairAlign),
// bound to o
func addTrimFlags(c *cobra.Command, o *trimOptions, withPad bool) {
	pad := ""
	if withPad {
		pad = ", or replaced with N if --pad"
	}
	c.Flags().IntVar(&o.start, "start", -1, "1-based first position (in reference coordinates) to keep in the output. Bases before it are omitted"+pad)
	c.Flags().IntVar(&o.end, "end", -1, "1-based last position (in reference coordinates) to keep in the output. Bases after it are omitted"+pad)
	c.Flags().IntVarP(&o.wrap, "wrap", "w", -1, "Wrap the output alignment to this number of nucleotides wide. Omit this option not to wrap the output.")
	hideDefault(c, "start", "end", "wrap")
}

// hideDefault stops the help from showing the named flags' defaults, which only mean "not given"
func hideDefault(c *cobra.Command, names ...string) {
	for _, n := range names {
		c.Flags().Lookup(n).DefValue = "0"
	}
}

// addOutfileFlag adds --outfile (-o, default stdout) to c, for writing what (e.g. "the variants")
func addOutfileFlag(c *cobra.Command, what string) {
	c.Flags().StringP("outfile", "o", "stdout", "Where to write "+what)
}
