package cmd

import (
	"github.com/spf13/cobra"
	"io"

	"github.com/virus-evolution/gofasta/pkg/gfio"
	"github.com/virus-evolution/gofasta/pkg/variants"
)

// variantsOptions holds the values of the flags of the command newVariantsCmd returns
type variantsOptions struct {
	query       string
	referenceID string // a record of the query alignment
	reference   string // a fasta file
	threads     int
	variantsOutputOptions
	annotationOptions
	strict bool
}

func newVariantsCmd() *cobra.Command {
	var o variantsOptions
	c := &cobra.Command{
		Use:   "variants",
		Short: "Annotate mutations relative to a reference from a multiple sequence alignment in fasta format",
		Long: `Annotate mutations relative to a reference from a multiple sequence alignment in fasta format

Example usage:

	gofasta variants -q alignment.fasta -a MN908947.gb -r MN908947.fasta > variants.csv
	gofasta variants -q alignment.fasta -a MN908947.gb --reference-id MN908947.3 > variants.csv
	gofasta variants -q alignment.fasta -r MN908947.fasta > variants.csv

The reference is a fasta file of its own, the same width as the alignment (--reference), or a record of the
alignment (--reference-id; if the alignment is read from stdin, it must be the first record). Without either,
the sequence in the annotation file is used, in which case the alignment must be in the same coordinates.
Without an annotation, only nucleotide changes, deletions and insertions are reported, and a reference must be
given.

gff-format annotations must be valid version 3 files. See github.com/virus-evolution/gofasta for more details
of the format.

If --query and --outfile are not given, the alignment is read from stdin and the variants are written to
stdout.

You can use --aggregate to report the overall proportions of each mutation in the alignment, and --threshold
to filter on frequency.

` + mutationFormatsHelp,
		PreRunE: func(cmd *cobra.Command, args []string) error {
			return firstError(
				intAtLeast(cmd, 0, "threads"),
				checkVariantsOutputFlags(cmd),
				checkAnnotationFlags(cmd, false),
				referenceFlags(cmd, o.annotationOptions.given(cmd)),
				outputNotInput(cmd, "outfile", "query", "reference", "annotation", "genbank", "gff"),
			)
		},
		RunE: func(cmd *cobra.Command, args []string) (err error) {

			msa, err := openIn(cmd, "query")
			if err != nil {
				return err
			}
			defer msa.Close()

			stdin := gfio.IsStdin(o.query)

			var ref io.ReadCloser // nil unless --reference
			if o.reference != "" {
				ref, err = openIn(cmd, "reference")
				if err != nil {
					return err
				}
				defer ref.Close()
			}

			var anno io.ReadCloser // nil, and annoSuffix "", without an annotation
			var annoSuffix string
			if o.annotationOptions.given(cmd) {
				anno, annoSuffix, err = o.annotationOptions.open(cmd)
				if err != nil {
					return err
				}
				defer anno.Close()
			}

			out, err := openOut(cmd, "outfile")
			if err != nil {
				return err
			}
			defer gfio.CloseOut(out, &err)

			return variants.Variants(msa, stdin, o.referenceID, ref, anno, annoSuffix, out, o.start, o.end, o.aggregate, o.threshold, o.appendSNP, o.groupAA, resolveThreads(o.threads), o.strict)
		},
	}

	c.Flags().StringVarP(&o.query, "query", "q", "stdin", "Alignment of the sequences to annotate, in fasta format")
	c.Flags().StringVarP(&o.reference, "reference", "r", "", "Reference sequence, in a fasta file of its own, the same width as the alignment")
	c.Flags().StringVar(&o.referenceID, "reference-id", "", "ID of the reference record in the alignment (instead of --reference)")
	addAnnotationFlags(c, &o.annotationOptions)
	addOutfileFlag(c, "the variants")
	addVariantsOutputFlags(c, &o.variantsOutputOptions)
	addThreadsFlag(c.Flags(), &o.threads)

	addStrictFlag(c, &o.strict)
	c.Flags().SortFlags = false
	return c
}
