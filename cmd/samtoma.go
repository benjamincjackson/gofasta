package cmd

import (
	"errors"

	"github.com/spf13/cobra"

	"github.com/virus-evolution/gofasta/pkg/gfio"
	"github.com/virus-evolution/gofasta/pkg/sam"
)

// toMultiAlignOptions holds the values of the flags of the command newToMultiAlignCmd returns
type toMultiAlignOptions struct {
	outfile   string
	start     int
	end       int
	pad       bool
	wrap      int
	trim      bool
	trimStart int
	trimEnd   int // --trim, --trimstart and --trimend are hidden, kept for backwards compatibility
	strict    bool
}

// s holds the values of the sam command's flags, which its subcommands share
func newToMultiAlignCmd(s *samOptions) *cobra.Command {
	var o toMultiAlignOptions
	c := &cobra.Command{
		Use:     "toMultiAlign",
		Aliases: []string{"tomultialign", "toma"},
		Short:   "Convert a SAM file to a multiple alignment in fasta format",
		Long: `Convert a SAM file to a multiple alignment in fasta format

Insertions relative to the reference are omitted, so all sequences in the output are the same ( = reference) length.

Example usage:
	gofasta sam toMultiAlign -s aligned.sam -o aligned.fasta

If you want, you can trim (and optionally pad) the output alignment to coordinates of your choosing:
	gofasta sam toMultiAlign -s aligned.sam --start 266 --end 29674 --pad -o aligned.fasta

If input and output files are not specified, the behaviour is to read the sam file from stdin and write
the fasta file to stdout, e.g.:
	minimap2 -a -x asm20 --score-N=0 reference.fasta unaligned.fasta | gofasta sam toMultiAlign > aligned.fasta`,

		PreRunE: func(cmd *cobra.Command, args []string) error {
			return firstError(
				intAtLeast(cmd, 1, "threads"),
				outputNotInput(cmd, "fasta-out", "samfile", "reference"),
				trimFlagsErr(&o),
			)
		},
		RunE: func(cmd *cobra.Command, args []string) (err error) {

			/*
				Start of trimming argument reconciliation to maintain backwards compatibility
				(trimFlagsErr, from PreRunE, has checked the two sets of flags are not combined)
			*/
			if o.trimStart != -1 {
				o.start = o.trimStart + 1
			}

			if o.trimEnd != -1 {
				o.end = o.trimEnd
			}
			/*
				End of trimming argument reconciliation to maintain backwards compatibility
			*/

			samIn, err := openIn(cmd, "samfile")
			if err != nil {
				return err
			}
			defer samIn.Close()

			out, err := openOut(cmd, "fasta-out")
			if err != nil {
				return err
			}
			defer gfio.CloseOut(out, &err)

			err = sam.ToMultiAlign(samIn, out, o.wrap, o.start, o.end, o.pad, s.threads, o.strict)

			return
		},
	}

	c.Flags().IntVarP(&o.start, "start", "", -1, "1-based first nucleotide position to retain in the output. Bases before this position are omitted, or are replaced with N if --pad")
	c.Flags().IntVarP(&o.end, "end", "", -1, "1-based last nucleotide position to retain in the output. Bases after this position are omitted, or are replaced with N if --pad")
	c.Flags().BoolVarP(&o.pad, "pad", "", false, "If --start and/or --end, replace the trimmed-out regions with Ns, else replace external deletions with Ns")
	c.Flags().StringVarP(&o.outfile, "fasta-out", "o", "stdout", "Where to write the alignment")
	c.Flags().IntVarP(&o.wrap, "wrap", "w", -1, "Wrap the output alignment to this number of nucleotides wide. Omit this option not to wrap the output.")

	c.Flags().BoolVarP(&o.trim, "trim", "", false, "Trim the alignment")
	c.Flags().IntVarP(&o.trimStart, "trimstart", "", -1, "Start coordinate for trimming (0-based, half open)")
	c.Flags().IntVarP(&o.trimEnd, "trimend", "", -1, "End coordinate for trimming (0-based, half open)")
	c.Flags().MarkHidden("trim")
	c.Flags().MarkHidden("trimstart")
	c.Flags().MarkHidden("trimend")

	addStrictFlag(c, &o.strict)
	c.Flags().SortFlags = false
	return c
}

// junk:

// trimFlagsErr checks that the old trimming flags (--trim, --trimstart, --trimend) are not combined
// with --start and --end, which replace them
func trimFlagsErr(o *toMultiAlignOptions) error {
	if (o.trim || o.trimStart != -1 || o.trimEnd != -1) && (o.start != -1 || o.end != -1) {
		return errors.New(`As of version 1.0.0 --start and --end replace --trim, --trimstart and --trimend 
The old flags and their behaviour are retained for backwards compatibility, but you shouldn't combine the two sets of flags.
Note the change of coordinate system if moving from old to new flags.

`)
	}
	return nil
}
