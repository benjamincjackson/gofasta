package cmd

import (
	"errors"
	"os"

	"github.com/spf13/cobra"

	"github.com/virus-evolution/gofasta/pkg/gfio"
	"github.com/virus-evolution/gofasta/pkg/variants"
)

func init() {
	rootCmd.AddCommand(newVariantsCmd())
}

// variantsOptions holds the values of the flags of the command newVariantsCmd returns
type variantsOptions struct {
	msa        string
	reference  string
	annotation string
	outfile    string
	threads    int
	aggregate  bool
	threshold  float64
	appendSNP  bool
	start      int
	end        int
	genbank    string // the hidden --genbank, kept for backwards compatibility
	strict     bool
}

func newVariantsCmd() *cobra.Command {
	var o variantsOptions
	c := &cobra.Command{
		Use:   "variants",
		Short: "Annotate mutations relative to a reference from a multiple sequence alignment in fasta format",
		Long: `Annotate mutations relative to a reference from a multiple sequence alignment in fasta format

Example usage:

	./gofasta variants --msa alignment.fasta --annotation MN908947.gb --reference MN908947.3 > variants.csv
	./gofasta variants --msa another.fasta --annotation MN908947.gff --reference Wuhan-Hu-1 > variants.csv

--reference is the name of the reference record in --msa. If you are reading the --msa from stdin, it must
be the first sequence. If you don't provide a --reference the program will try to use the fasta record in the
annotation file, in which case the --msa must be in the same coordinates.

gff-format annotations must be valid version 3 files. See github.com/virus-evolution/gofasta for more details
of the format.

If input --msa and output csv files are not specified, the behaviour is to read the alignment from stdin and write
the variants to stdout.

You can use --aggregate to report the overall proportions of each mutation in the --msa, and --threshold to filter on 
frequency.

Mutations are annotated with ins (insertion), del (deletion), aa (amino acid change) or nuc (a nucleotide change that
isn't in a codon that is represented by an amino acid change). The formats are:

	ins:2028:3 - a 3-base insertion immediately after (1-based) position 2028 in reference coordinates
	del:11288:9 - a 9-base deletion whose first missing nucleotide is at (1-based) position 11288 in reference coordinates
	aa:S:D614G - the amino acid at (1-based) residue 614 in the S gene is a D in the reference and a G in this sequence
	nuc:C3037T - the nucleotide at (1-based) position 3037 in reference coordinates is a C in the reference and a T in this sequence

Frame-shifting mutations in coding sequence are reported as indels but are ignored for subsequent amino-acids in the alignment.	
`,
		PreRunE: func(cmd *cobra.Command, args []string) error {
			// --reference is the ID of a record in the msa, not a file
			return firstError(
				intAtLeast(cmd, 1, "threads"),
				floatAtLeast(cmd, 0, "threshold"),
				startNotAfterEnd(cmd),
				outputNotInput(cmd, "outfile", "msa", "annotation", "genbank"),
			)
		},
		RunE: func(cmd *cobra.Command, args []string) (err error) {

			msa, err := openIn(cmd, "msa")
			if err != nil {
				return err
			}
			defer msa.Close()

			stdin := false
			if o.msa == "stdin" {
				stdin = true
			}

			// some backwards compatibility wrangling of --genbank vs --annotation
			var anno *os.File
			var annoSuffix string
			if o.genbank != "" {
				if o.annotation != "" {
					return errors.New("--annotation replaces --genbank")
				}
				anno, err = openIn(cmd, "genbank")
				if err != nil {
					return err
				}
				annoSuffix = "gb"
			} else {
				annoSuffix, err = annotationFormat(o.annotation)
				if err != nil {
					return err
				}
				anno, err = openIn(cmd, "annotation")
				if err != nil {
					return err
				}
			}
			defer anno.Close()

			out, err := openOut(cmd, "outfile")
			if err != nil {
				return err
			}
			defer gfio.CloseOut(out, &err)

			err = variants.Variants(msa, stdin, o.reference, anno, annoSuffix, out, o.start, o.end, o.aggregate, o.threshold, o.appendSNP, o.threads, o.strict)

			return
		},
	}

	c.Flags().StringVarP(&o.msa, "msa", "", "stdin", "Multiple sequence alignment in fasta format")
	c.Flags().StringVarP(&o.reference, "reference", "r", "", "The ID of the reference record in the msa")
	c.Flags().StringVarP(&o.annotation, "annotation", "a", "", "GenBank or GFF3 format annotation file. Must have suffix .gb, .gbk, .gbff, .gff or .gff3")
	c.Flags().StringVarP(&o.outfile, "outfile", "o", "stdout", "Name of the file of variants to write")
	c.Flags().IntVarP(&o.start, "start", "", -1, "Only report variants with a reference base at or after this position")
	c.Flags().IntVarP(&o.end, "end", "", -1, "Only report variants with a reference base at or before this position")
	c.Flags().BoolVarP(&o.aggregate, "aggregate", "", false, "Report the proportions of each change")
	c.Flags().Float64VarP(&o.threshold, "threshold", "", 0.0, "If --aggregate, only report changes with a freq greater than or equal to this value")
	c.Flags().BoolVarP(&o.appendSNP, "append-snps", "", false, "Report the codon's SNPs in parenthesis after each amino acid mutation")
	c.Flags().IntVarP(&o.threads, "threads", "t", 1, "Number of threads to use")

	c.Flags().Lookup("aggregate").NoOptDefVal = "true"
	c.Flags().Lookup("append-snps").NoOptDefVal = "true"

	c.Flags().StringVarP(&o.genbank, "genbank", "", "", "Genbank format annotation")
	c.Flags().MarkHidden("genbank")

	addStrictFlag(c, &o.strict)
	c.Flags().SortFlags = false
	return c
}
