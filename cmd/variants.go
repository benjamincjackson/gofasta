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
	outfile     string
	threads     int
	aggregate   bool
	threshold   float64
	appendSNP   bool
	start       int
	end         int
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

Mutations are annotated with ins (insertion), del (deletion), aa (amino acid change) or nuc (a nucleotide
change that isn't in a codon that is represented by an amino acid change). The formats are:

	ins:2028:3 - a 3-base insertion immediately after (1-based) position 2028 in reference coordinates
	del:11288:9 - a 9-base deletion whose first missing nucleotide is at (1-based) position 11288 in reference coordinates
	aa:S:D614G - the amino acid at (1-based) residue 614 in the S gene is a D in the reference and a G in this sequence
	nuc:C3037T - the nucleotide at (1-based) position 3037 in reference coordinates is a C in the reference and a T in this sequence

Frame-shifting mutations in coding sequence are reported as indels but are ignored for subsequent amino-acids
in the alignment.
`,
		PreRunE: func(cmd *cobra.Command, args []string) error {
			return firstError(
				intAtLeast(cmd, 0, "threads"),
				floatAtLeast(cmd, 0, "threshold"),
				startNotAfterEnd(cmd),
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

			var refIn io.Reader // (an untyped nil if there is no --reference)
			if ref != nil {
				refIn = ref
			}
			err = variants.Variants(msa, stdin, o.referenceID, refIn, anno, annoSuffix, out, o.start, o.end, o.aggregate, o.threshold, o.appendSNP, resolveThreads(o.threads), o.strict)

			return
		},
	}

	c.Flags().StringVarP(&o.query, "query", "q", "stdin", "Alignment of the sequences to annotate, in fasta format")
	c.Flags().StringVarP(&o.reference, "reference", "r", "", "Reference sequence, in a fasta file of its own, the same width as the alignment")
	c.Flags().StringVar(&o.referenceID, "reference-id", "", "ID of the reference record in the alignment (instead of --reference)")
	addAnnotationFlags(c, &o.annotationOptions)
	c.Flags().StringVarP(&o.outfile, "outfile", "o", "stdout", "Name of the file of variants to write")
	c.Flags().IntVarP(&o.start, "start", "", -1, "Only report variants with a reference base at or after this position")
	c.Flags().IntVarP(&o.end, "end", "", -1, "Only report variants with a reference base at or before this position")
	c.Flags().BoolVarP(&o.aggregate, "aggregate", "", false, "Report the proportions of each change")
	c.Flags().Float64VarP(&o.threshold, "threshold", "", 0.0, "If --aggregate, only report changes with a freq greater than or equal to this value")
	c.Flags().BoolVarP(&o.appendSNP, "append-snps", "", false, "Report the codon's SNPs in parenthesis after each amino acid mutation")
	addThreadsFlag(c.Flags(), &o.threads)

	c.Flags().Lookup("aggregate").NoOptDefVal = "true"
	c.Flags().Lookup("append-snps").NoOptDefVal = "true"

	addStrictFlag(c, &o.strict)
	c.Flags().SortFlags = false
	return c
}
