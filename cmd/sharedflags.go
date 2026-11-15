package cmd

import (
	"runtime"

	"github.com/spf13/cobra"
	"github.com/spf13/pflag"
)

// threadsUsage is the description of every command's -t, --threads flag, whose default, 0, means all CPUs
const threadsUsage = "Number of threads (and CPUs) to use (with no flag, all CPUs are used)"

// addThreadsFlag adds -t, --threads to flags (a command's flags, or a parent's persistent flags), setting
// *threads
func addThreadsFlag(flags *pflag.FlagSet, threads *int) {
	flags.IntVarP(threads, "threads", "t", 0, threadsUsage)
}

// resolveThreads returns the number of threads to use for -t threads (0 for all CPUs), and limits the
// CPUs gofasta uses to it, so that -t also limits the CPUs used (e.g. -t 1 on a shared machine)
func resolveThreads(threads int) int {
	n := runtime.NumCPU()
	if threads > 0 && threads < n {
		runtime.GOMAXPROCS(threads)
		return threads
	}
	return n
}

// addStrictFlag adds the --strict flag, which each command that skips records it can't process has, to c,
// setting *strict
func addStrictFlag(c *cobra.Command, strict *bool) {
	c.Flags().BoolVar(strict, "strict", false, "Stop with an error at the first record that can't be processed, instead of warning and skipping it")
}

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
	c.Flags().IntVar(&o.start, "start", -1, "Only report mutations with a reference base at or after this position")
	c.Flags().IntVar(&o.end, "end", -1, "Only report mutations with a reference base at or before this position")
	c.Flags().BoolVar(&o.aggregate, "aggregate", false, "Report the frequency of each mutation")
	c.Flags().Float64Var(&o.threshold, "threshold", 0.0, "If --aggregate, only report mutations with a freq greater than or equal to this value")
	c.Flags().BoolVar(&o.appendSNP, "append-snps", false, "Report the codon's SNPs in parenthesis after each amino acid change")
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

With --group-aa, the amino acid changes from one indel are one entry: aa:S:HV69-- instead of aa:S:H69-|aa:S:V70-.

An indel that shifts the reading frame of a coding sequence is reported only as an indel, as are the indels
after it in that coding sequence, unless indels within 3 codons of it restore the frame. Changes of single
bases are still reported codon by codon.
`

// hideDefault stops the help from showing the named flags' defaults, which only mean "not given"
func hideDefault(c *cobra.Command, names ...string) {
	for _, n := range names {
		c.Flags().Lookup(n).DefValue = "0"
	}
}

// requireFlags marks the named flags of c required, and says so in their help
func requireFlags(c *cobra.Command, names ...string) {
	for _, n := range names {
		c.MarkFlagRequired(n)
		f := c.Flags().Lookup(n)
		f.Usage += " (required)"
	}
}

// addOutfileFlag adds --outfile (-o, default stdout) to c, for writing what (e.g. "the variants")
func addOutfileFlag(c *cobra.Command, what string) {
	c.Flags().StringP("outfile", "o", "stdout", "Where to write "+what)
}
