package cmd

import (
	"io"

	"github.com/spf13/cobra"

	"github.com/virus-evolution/gofasta/pkg/gfio"
	"github.com/virus-evolution/gofasta/pkg/sam"
)

// samVariantsOptions holds the values of the flags of the command newSamVariantsCmd returns
type samVariantsOptions struct {
	outfile   string
	aggregate bool
	threshold float64
	appendSNP bool
	start     int
	end       int
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

Mutations are annotated with ins (insertion), del (deletion), aa (amino acid change) or nuc (a nucleotide
change that isn't in a codon that is represented by an amino acid change). The formats are:

	ins:2028:3 - a 3-base insertion immediately after (1-based) position 2028 in reference coordinates
	del:11288:9 - a 9-base deletion whose first missing nucleotide is at (1-based) position 11288 in reference coordinates
	aa:S:D614G - the amino acid at (1-based) residue 614 in the S gene is a D in the reference and a G in this sequence
	aa:S:H69- - the amino acid at (1-based) residue 69 in the S gene is an H in the reference and is deleted in this sequence
	aa:S:214EPE - the amino acids EPE are inserted immediately after (1-based) residue 214 in the S gene in this sequence
	nuc:C3037T - the nucleotide at (1-based) position 3037 in reference coordinates is a C in the reference and a T in this sequence

Frame-shifting mutations in coding sequence are reported as indels but are ignored for subsequent amino-acids
in the alignment.
`,

		PreRunE: func(cmd *cobra.Command, args []string) error {
			return firstError(
				intAtLeast(cmd, 0, "threads"),
				floatAtLeast(cmd, 0, "threshold"),
				startNotAfterEnd(cmd),
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

			err = sam.Variants(samIn, ref, refFromFile, anno, annoSuffix, out, o.start, o.end, o.aggregate, o.threshold, o.appendSNP, resolveThreads(s.threads), o.strict)

			return err
		},
	}

	addAnnotationFlags(c, &o.annotationOptions)
	c.Flags().StringVarP(&o.outfile, "outfile", "o", "stdout", "Where to write the variants")
	c.Flags().IntVarP(&o.start, "start", "", -1, "Only report variants with a reference base at or after this position")
	c.Flags().IntVarP(&o.end, "end", "", -1, "Only report variants with a reference base at or before this position")
	c.Flags().BoolVarP(&o.aggregate, "aggregate", "", false, "Report the proportions of each change")
	c.Flags().Float64VarP(&o.threshold, "threshold", "", 0.0, "If --aggregate, only report changes with a freq greater than or equal to this value")
	c.Flags().BoolVarP(&o.appendSNP, "append-snps", "", false, "Report the codon's SNPs in parenthesis after each amino acid mutation")

	c.Flags().Lookup("aggregate").NoOptDefVal = "true"
	c.Flags().Lookup("append-snps").NoOptDefVal = "true"

	addStrictFlag(c, &o.strict)
	c.Flags().SortFlags = false

	return c
}
