package cmd

import (
	"fmt"

	"github.com/spf13/cobra"
)

// samOptions holds the values of the flags of the command newSamCmd returns, which its subcommands share
type samOptions struct {
	threads   int
	file      string
	reference string
}

func newSamCmd() *cobra.Command {
	var o samOptions
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

	addThreadsFlag(c.PersistentFlags(), &o.threads)
	c.PersistentFlags().StringVarP(&o.file, "samfile", "s", "stdin", "SAM file to read")
	c.PersistentFlags().StringVarP(&o.reference, "reference", "r", "", "Reference fasta file used to generate the sam file")
	c.AddCommand(newToMultiAlignCmd(&o), newToPairAlignCmd(&o), newSamVariantsCmd(&o))
	return c
}
