package cmd

import (
	"github.com/spf13/cobra"

	"github.com/virus-evolution/gofasta/pkg/gfio"
	"github.com/virus-evolution/gofasta/pkg/snps"
)

// snpsOptions holds the values of the flags of the command newSNPsCmd returns
type snpsOptions struct {
	reference string
	query     string
	outfile   string
	hardGaps  bool
	aggregate bool
	threshold float64
	strict    bool
}

func newSNPsCmd() *cobra.Command {
	var o snpsOptions
	c := &cobra.Command{
		Use:   "snps",
		Short: "Find snps relative to a reference",
		Long: `Find snps relative to a reference.

Example usage:
	gofasta snps -r reference.fasta -q alignment.fasta -o snps.csv

reference.fasta and alignment.fasta must be the same width.

With the default settings the output is a csv-format file with one line per query sequence, and two columns:
'query' and 'SNPs', the second of which is a "|"-delimited list of snps in that query.

If you set --aggregate (and optionally a --threshold) it will return the SNPs present in the entire sample
(whose frequency is equal to/above --threshold) and their frequencies.

Setting --hard-gaps treats alignment gaps as different from {ATGC}.

If query and outfile are not specified, the behaviour is to read the query alignment
from stdin and write the snps file to stdout, e.g. you could do this:
	cat alignment.fasta | gofasta snps -r reference.fasta > snps.csv`,

		PreRunE: func(cmd *cobra.Command, args []string) error {
			return firstError(
				floatAtLeast(cmd, 0, "threshold"),
				outputNotInput(cmd, "outfile", "query", "reference"),
			)
		},
		RunE: func(cmd *cobra.Command, args []string) (err error) {

			query, err := openIn(cmd, "query")
			if err != nil {
				return err
			}
			defer query.Close()

			ref, err := openIn(cmd, "reference")
			if err != nil {
				return err
			}
			defer ref.Close()

			out, err := openOut(cmd, "outfile")
			if err != nil {
				return err
			}
			defer gfio.CloseOut(out, &err)

			err = snps.SNPs(ref, query, o.hardGaps, o.aggregate, o.threshold, out, o.strict)

			return
		},
	}

	c.Flags().StringVarP(&o.reference, "reference", "r", "", "Reference sequence, in fasta format")
	c.Flags().StringVarP(&o.query, "query", "q", "stdin", "Alignment of sequences to find snps in, in fasta format")
	c.Flags().StringVarP(&o.outfile, "outfile", "o", "stdout", "Output to write")
	c.Flags().BoolVarP(&o.hardGaps, "hard-gaps", "", false, "Don't treat alignment gaps as missing data")
	c.Flags().BoolVarP(&o.aggregate, "aggregate", "", false, "Report the proportions of each change")
	c.Flags().Float64VarP(&o.threshold, "threshold", "", 0.0, "If --aggregate, only report snps with a freq greater than or equal to this value")

	c.Flags().Lookup("hard-gaps").NoOptDefVal = "true"
	c.Flags().Lookup("aggregate").NoOptDefVal = "true"

	addStrictFlag(c, &o.strict)
	c.Flags().SortFlags = false
	return c
}
