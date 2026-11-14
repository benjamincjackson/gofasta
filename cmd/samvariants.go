package cmd

import (
	"io"

	"github.com/spf13/cobra"

	"github.com/virus-evolution/gofasta/pkg/gfio"
	"github.com/virus-evolution/gofasta/pkg/sam"
)

// samVariantsOptions holds the values of the flags of the command newSamVariantsCmd returns
type samVariantsOptions struct {
	variantsOutputOptions
	annotationOptions
	strict bool
}

// newSamVariantsCmd builds the sam variants command; s holds the sam command's flags, which its subcommands share
func newSamVariantsCmd(s *samOptions) *cobra.Command {
	var o samVariantsOptions
	c := &cobra.Command{
		Use:   "variants",
		Short: "Annotate mutations relative to a reference from an alignment in sam format",
		Long: `Annotate mutations relative to a reference from an alignment in sam format

Example usage:
	gofasta sam variants -s aligned.sam -r reference.fasta -a annotation.gb -o variants.csv
	gofasta sam variants -s aligned.sam -r reference.fasta -a annotation.gff -o variants.csv

If input sam and output csv files are not specified, the behaviour is to read the sam from stdin and write the
variants to stdout.

--reference should be the same sequence that was used to generate the sam file, and should be in the same
coordinates as the --annotation. You don't have to provide a file to --reference if your annotation has the
fasta record in it.

gff-format annotations must be valid version 3 files. See github.com/virus-evolution/gofasta for more details
of the format.

` + mutationFormatsHelp,

		PreRunE: func(cmd *cobra.Command, args []string) error {
			return firstError(
				intAtLeast(cmd, 0, "threads"),
				checkVariantsOutputFlags(cmd),
				checkAnnotationFlags(cmd, true),
				outputNotInput(cmd, "outfile", "samfile", "reference", "annotation", "genbank", "gff"),
			)
		},
		RunE: func(cmd *cobra.Command, args []string) (err error) {

			samIn, err := openIn(cmd, "samfile")
			if err != nil {
				return err
			}
			defer samIn.Close()

			refFromFile := false
			var ref io.ReadCloser // without --reference, the annotation's sequence is the reference
			if s.reference != "" {
				ref, err = openIn(cmd, "reference")
				if err != nil {
					return err
				}
				defer ref.Close()
				refFromFile = true
			}

			anno, annoSuffix, err := o.annotationOptions.open(cmd)
			if err != nil {
				return err
			}
			defer anno.Close()

			out, err := openOut(cmd, "outfile")
			if err != nil {
				return err
			}
			defer gfio.CloseOut(out, &err)

			return sam.Variants(samIn, ref, refFromFile, anno, annoSuffix, out, o.start, o.end, o.aggregate, o.threshold, o.appendSNP, o.groupAA, resolveThreads(s.threads), o.strict)
		},
	}

	addAnnotationFlags(c, &o.annotationOptions)
	addOutfileFlag(c, "the variants")
	addVariantsOutputFlags(c, &o.variantsOutputOptions)

	addStrictFlag(c, &o.strict)
	c.Flags().SortFlags = false

	return c
}
