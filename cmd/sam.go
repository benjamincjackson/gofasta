package cmd

import (
	"fmt"

	"github.com/spf13/cobra"
)

func newSamCmd() *cobra.Command {
	c := &cobra.Command{
		Use:   "sam",
		Short: "Convert SAM-format alignments to fasta format, and extract variants from SAM-format alignments",
		Long:  `Convert SAM-format alignments to fasta format, and extract variants from SAM-format alignments`,
		Args:  cobra.MinimumNArgs(1),
		RunE: func(cmd *cobra.Command, args []string) error {
			// Args is MinimumNArgs(1), so there is always an argument here, and it is not a subcommand
			// (cobra would have run that instead)
			return fmt.Errorf("unknown command %q for %q", args[0], cmd.CommandPath())
		},
	}

	c.AddCommand(newToMultiAlignCmd(), newToPairAlignCmd(), newSamVariantsCmd())
	return c
}

// addSamfileFlag adds -s, --samfile (default stdin) to c, a sam command
func addSamfileFlag(c *cobra.Command) {
	c.Flags().StringP("samfile", "s", "stdin", "SAM file to read")
}

// samReferenceUsage describes the --reference of the sam commands that read one
const samReferenceUsage = "Reference fasta file used to make the --samfile alignment"

// trimOptions holds the flags with which sam toMultiAlign and sam toPairAlign trim and wrap their output
type trimOptions struct {
	start, end int
	wrap       int
}

// addTrimFlags adds --start and --end to c (sam toMultiAlign, with --pad, and sam toPairAlign), bound to o
func addTrimFlags(c *cobra.Command, o *trimOptions, withPad bool) {
	pad := ""
	if withPad {
		pad = ", or replaced with N if --pad"
	}
	c.Flags().IntVar(&o.start, "start", -1, "1-based first position (in reference coordinates) to keep in the output. Bases before it are omitted"+pad)
	c.Flags().IntVar(&o.end, "end", -1, "1-based last position (in reference coordinates) to keep in the output. Bases after it are omitted"+pad)
	hideDefault(c, "start", "end")
}

// addWrapFlag adds -w, --wrap to c (sam toMultiAlign and sam toPairAlign), bound to o
func addWrapFlag(c *cobra.Command, o *trimOptions) {
	c.Flags().IntVarP(&o.wrap, "wrap", "w", -1, "Wrap the output alignment to this number of nucleotides wide. Omit this option not to wrap the output.")
	hideDefault(c, "wrap")
}
