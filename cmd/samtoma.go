package cmd

import (
	"github.com/spf13/cobra"

	"github.com/virus-evolution/gofasta/pkg/gfio"
	"github.com/virus-evolution/gofasta/pkg/sam"
)

// toMultiAlignOptions holds the values of the flags of the command newToMultiAlignCmd returns
type toMultiAlignOptions struct {
	start  int
	end    int
	pad    bool
	wrap   int
	strict bool
}

// newToMultiAlignCmd builds the sam toMultiAlign command; s holds the sam command's flags, which its subcommands share
func newToMultiAlignCmd(s *samOptions) *cobra.Command {
	var o toMultiAlignOptions
	c := &cobra.Command{
		Use:     "toMultiAlign",
		Aliases: []string{"tomultialign", "toma"},
		Short:   "Convert a SAM file to a multiple alignment in fasta format",
		Long: `Convert a SAM file to a multiple alignment in fasta format

Insertions relative to the reference are omitted, so all sequences in the output are the same ( = reference)
length.

Example usage:
	gofasta sam toMultiAlign -s aligned.sam -o aligned.fasta

If you want, you can trim (and optionally pad) the output alignment to coordinates of your choosing:
	gofasta sam toMultiAlign -s aligned.sam --start 266 --end 29674 --pad -o aligned.fasta

If input and output files are not specified, the behaviour is to read the sam file from stdin and write the
fasta file to stdout, e.g.:
	minimap2 -a -x asm20 --score-N=0 reference.fasta unaligned.fasta | gofasta sam toMultiAlign > aligned.fasta`,

		PreRunE: func(cmd *cobra.Command, args []string) error {
			return firstError(
				intAtLeast(cmd, 0, "threads"),
				outputNotInput(cmd, "outfile", "samfile", "reference"),
			)
		},
		RunE: func(cmd *cobra.Command, args []string) (err error) {

			samIn, err := openIn(cmd, "samfile")
			if err != nil {
				return err
			}
			defer samIn.Close()

			out, err := openOut(cmd, "outfile")
			if err != nil {
				return err
			}
			defer gfio.CloseOut(out, &err)

			return sam.ToMultiAlign(samIn, out, o.wrap, o.start, o.end, o.pad, resolveThreads(s.threads), o.strict)
		},
	}

	c.Flags().IntVar(&o.start, "start", -1, "1-based first nucleotide position to retain in the output. Bases before this position are omitted, or are replaced with N if --pad")
	c.Flags().IntVar(&o.end, "end", -1, "1-based last nucleotide position to retain in the output. Bases after this position are omitted, or are replaced with N if --pad")
	c.Flags().BoolVar(&o.pad, "pad", false, "If --start and/or --end, replace the trimmed-out regions with Ns, else replace external deletions with Ns")
	c.Flags().StringP("outfile", "o", "stdout", "Where to write the alignment")
	c.Flags().IntVarP(&o.wrap, "wrap", "w", -1, "Wrap the output alignment to this number of nucleotides wide. Omit this option not to wrap the output.")

	addStrictFlag(c, &o.strict)
	c.Flags().SortFlags = false
	return c
}
