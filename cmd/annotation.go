package cmd

import (
	"errors"
	"io"
	"strings"

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

// checkAnnotationFlags checks that exactly one annotation flag was given
func checkAnnotationFlags(cmd *cobra.Command) error {
	given := 0
	for _, name := range annotationFlagNames {
		if cmd.Flags().Changed(name) {
			given++
		}
	}
	switch {
	case given == 0:
		return errors.New("give an annotation with one of --annotation, --genbank or --gff")
	case given > 1:
		return errors.New("give only one of --annotation, --genbank and --gff")
	}
	return nil
}

// open opens the annotation and returns it with its format ("gb" or "gff"); checkAnnotationFlags has
// checked that exactly one of the flags was given
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

// oneStdin checks that at most one of the named input flags is stdin (its default may be)
func oneStdin(cmd *cobra.Command, names ...string) error {
	stdin := make([]string, 0)
	for _, name := range names {
		if f := cmd.Flag(name); f != nil && f.Value.String() == "stdin" {
			stdin = append(stdin, "--"+name)
		}
	}
	if len(stdin) > 1 {
		return errors.New("only one input can be read from stdin, but " + strings.Join(stdin, " and ") + " are both stdin")
	}
	return nil
}
