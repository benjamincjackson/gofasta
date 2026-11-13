package cmd

import (
	"fmt"

	"github.com/spf13/cobra"
)

// updownOptions holds the values of the flags of the command newUpdownCmd returns, which its subcommands share
type updownOptions struct {
	reference string
}

func newUpdownCmd() *cobra.Command {
	var o updownOptions
	c := &cobra.Command{
		Use:   "updown",
		Short: "Get pseudo-tree-aware catchments for query sequences from alignments",
		Long:  `Get pseudo-tree-aware catchments for query sequences from alignments`,
		Args:  cobra.MinimumNArgs(1),
		RunE: func(cmd *cobra.Command, args []string) error {
			// Args is MinimumNArgs(1), so there is always an argument here, and it is not a subcommand
			// (cobra would have run that instead)
			return fmt.Errorf("unknown command %q for %q", args[0], cmd.CommandPath())
		},
	}

	c.PersistentFlags().StringVarP(&o.reference, "reference", "r", "", "Reference sequence, in fasta format, which is treated as the root of the imaginary tree")
	c.AddCommand(newUpdownListCmd(), newToprankingCmd(&o))
	return c
}
