package cmd

import (
	"fmt"

	"github.com/spf13/cobra"
)

func newUpdownCmd() *cobra.Command {
	c := &cobra.Command{
		Use:   "updown",
		Short: updownHelpShort,
		Long:  updownHelpShort,
		Args:  cobra.MinimumNArgs(1),
		RunE: func(cmd *cobra.Command, args []string) error {
			// Args is MinimumNArgs(1), so there is always an argument here, and it is not a subcommand
			// (cobra would have run that instead)
			return fmt.Errorf("unknown command %q for %q", args[0], cmd.CommandPath())
		},
	}

	c.AddCommand(newUpdownListCmd(), newToprankingCmd())
	return c
}

const updownHelpShort = "Get pseudo-phylogenetic-tree-aware neighbours by genetic distance from alignments"
