package cmd

import (
	"errors"

	"github.com/spf13/cobra"
)

// newSNPsCmd builds the hidden snps command, which says what replaces it: it accepts any arguments and
// flags, so that an old command line gets this message rather than a flag error
func newSNPsCmd() *cobra.Command {
	return &cobra.Command{
		Use:                "snps",
		Hidden:             true,
		DisableFlagParsing: true,
		RunE: func(cmd *cobra.Command, args []string) error {
			cmd.SilenceUsage = true
			return errors.New(`gofasta snps was removed in version 2: gofasta variants does the same without an annotation.
Instead of
	gofasta snps -r reference.fasta -q alignment.fasta
use
	gofasta variants --reference-fasta reference.fasta --msa alignment.fasta
which reports changes as nuc:A123T (and deletions and insertions as del: and ins:), in reference coordinates`)
		},
	}
}
