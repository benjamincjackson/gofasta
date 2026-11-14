package cmd

import (
	"fmt"

	"github.com/spf13/cobra"
)

func newSamCmd() *cobra.Command {
	c := &cobra.Command{
		Use:   "sam",
		Short: "Do things with sam files",
		Long:  `Do things with sam files`,
		Args:  cobra.MinimumNArgs(1),
		RunE: func(cmd *cobra.Command, args []string) error {
			// Args is MinimumNArgs(1), so there is always an argument here, and it is not a subcommand
			// (cobra would have run that instead)
			return fmt.Errorf("unknown command %q for %q", args[0], cmd.CommandPath())
		},
	}

	c.AddCommand(newToMultiAlignCmd(), newToPairAlignCmd(), newSamVariantsCmd())
	return c
}
