package cmd

import (
	"errors"
	"io"

	"github.com/spf13/cobra"
)

// annotationOptions holds the annotation flags of variants and sam variants: exactly one of them is given
type annotationOptions struct {
	annotation string // format from the file's suffix
	genbank    string // GenBank, whatever the file is called
	gff        string // GFF3, whatever the file is called
}

// annotationFlagNames are the names of the annotation flags, in the order they are described
var annotationFlagNames = []string{"annotation", "genbank", "gff"}

// addAnnotationFlags adds --annotation, --genbank and --gff to c, bound to a
func addAnnotationFlags(c *cobra.Command, a *annotationOptions) {
	c.Flags().StringVarP(&a.annotation, "annotation", "a", "", "GenBank or GFF3 annotation file, its format taken from its suffix: .gb, .gbk, .gbff, .gff or .gff3 (optionally followed by .gz)")
	c.Flags().StringVar(&a.genbank, "genbank", "", "GenBank annotation file, whatever it is called (instead of --annotation)")
	c.Flags().StringVar(&a.gff, "gff", "", "GFF3 annotation file, whatever it is called (instead of --annotation)")
}

// checkAnnotationFlags checks that at most one annotation flag was given, and if required, that one was
func checkAnnotationFlags(cmd *cobra.Command, required bool) error {
	given := 0
	for _, name := range annotationFlagNames {
		if cmd.Flags().Changed(name) {
			given++
		}
	}
	switch {
	case given == 0 && required:
		return errors.New("give an annotation with one of --annotation, --genbank or --gff")
	case given > 1:
		return errors.New("give only one of --annotation, --genbank and --gff")
	}
	return nil
}

// given reports whether an annotation flag was given
func (a annotationOptions) given(cmd *cobra.Command) bool {
	for _, name := range annotationFlagNames {
		if cmd.Flags().Changed(name) {
			return true
		}
	}
	return false
}

// open opens the annotation and returns it with its format ("gb" or "gff"); checkAnnotationFlags has
// checked that at most one of the flags was given, and given has said that one was
func (a annotationOptions) open(cmd *cobra.Command) (io.ReadCloser, string, error) {
	name, format := "annotation", ""
	switch {
	case cmd.Flags().Changed("genbank"):
		name, format = "genbank", "gb"
	case cmd.Flags().Changed("gff"):
		name, format = "gff", "gff"
	default:
		var err error
		if format, err = annotationFormat(a.annotation); err != nil {
			return nil, "", err
		}
	}
	f, err := openIn(cmd, name)
	return f, format, err
}
