package cmd

import (
	"github.com/spf13/cobra"

	"github.com/virus-evolution/gofasta/pkg/gfio"
	"github.com/virus-evolution/gofasta/pkg/updown"
)

// updownListOptions holds the values of the flags of the command newUpdownListCmd returns
type updownListOptions struct {
	threads int
	query   string
	outfile string
	strict  bool
}

func newUpdownListCmd() *cobra.Command {
	var o updownListOptions
	c := &cobra.Command{
		Use:   "list",
		Short: "Generate input CSV files for gofasta updown topranking",
		Long: `Generate input CSV files for gofasta updown topranking

Example usage:

	gofasta updown list -r reference.fasta -q alignment.fasta -o mutationlist.csv

Non-ATGC nucleotides are not recommended in the --reference, and --reference and --query must
be aligned to the same thing.

--outfile is a CSV-format file with the columns: query,SNPs,ambiguities,SNPcount,ambcount. There is one row
for each sequence in --query. SNPs is a "|"-delimited list of SNPs relative to --reference. ambiguities is
a "|"-delimited list of ranges (1-based, inclusive) of tracts of ambiguities (anything that isn't ATGC).
`,

		PreRunE: func(cmd *cobra.Command, args []string) error {
			return firstError(
				inheritedRequired(cmd, "reference"),
				intAtLeast(cmd, 0, "threads"),
				outputNotInput(cmd, "outfile", "query", "reference"),
			)
		},
		RunE: func(cmd *cobra.Command, args []string) (err error) {

			ref, err := openIn(cmd, "reference")
			if err != nil {
				return err
			}
			defer ref.Close()

			query, err := openIn(cmd, "query")
			if err != nil {
				return err
			}
			defer query.Close()

			out, err := openOut(cmd, "outfile")
			if err != nil {
				return err
			}
			defer gfio.CloseOut(out, &err)

			err = updown.List(ref, query, out, resolveThreads(o.threads), o.strict)

			return
		},
	}

	c.Flags().StringVarP(&o.query, "query", "q", "stdin", "Alignment of sequences to parse, in fasta format")
	c.Flags().StringVarP(&o.outfile, "outfile", "o", "stdout", "Output to write")

	addThreadsFlag(c.Flags(), &o.threads)
	addStrictFlag(c, &o.strict)
	c.Flags().SortFlags = false
	return c
}
