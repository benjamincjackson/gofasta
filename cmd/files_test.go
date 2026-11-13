package cmd

import (
	"testing"

	"github.com/spf13/cobra"
)

// An input file that can't be opened is reported with the flag that named it
func TestOpenInNamesTheFlag(t *testing.T) {
	c := &cobra.Command{Use: "test"}
	var reference, msa string
	c.Flags().StringVarP(&reference, "reference", "r", "", "")
	c.Flags().StringVar(&msa, "msa", "", "")
	c.Flags().Set("reference", "not/a/file.whatever")
	c.Flags().Set("msa", "not/a/file.whatever")

	if _, err := openIn(c, "reference"); err == nil || err.Error() != "open -r / --reference not/a/file.whatever: no such file or directory" {
		t.Error(err)
	}
	if _, err := openIn(c, "msa"); err == nil || err.Error() != "open --msa not/a/file.whatever: no such file or directory" {
		t.Error(err)
	}
}
