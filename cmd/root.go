package cmd

import (
	"os"

	"github.com/spf13/cobra"
)

// newRootCmd builds the gofasta command and all its subcommands, each with its flags at their defaults
func newRootCmd() *cobra.Command {
	root := &cobra.Command{
		Use:   "gofasta",
		Short: "Command-line utilities for genomic epidemiology research",
		Long: `Command-line utilities for genomic epidemiology research
		
If you use gofasta in your work, please cite:

Jackson B (2022). gofasta: command-line utilities for genomic epidemiology research. Bioinformatics 38 (16), 4033-4035
https://doi.org/10.1093/bioinformatics/btac424`,
		Version: "1.2.3",
	}
	root.AddCommand(newClosestCmd(), newSNPsCmd(), newVariantsCmd(), newSamCmd(), newUpdownCmd(), newLicencesCmd())
	silenceUsageOnceRunning(root)
	return root
}

// silenceUsageOnceRunning makes each command that does the work (one without subcommands) print only the
// error, not the usage, if it fails once it is running: the flags and arguments were valid (they were
// parsed, and checked by PreRunE), so the problem is with the input, and the usage would bury the error.
// Errors in the flags or arguments, and unknown subcommands, still print the usage.
func silenceUsageOnceRunning(c *cobra.Command) {
	for _, sub := range c.Commands() {
		silenceUsageOnceRunning(sub)
	}
	if c.HasSubCommands() || c.RunE == nil {
		return
	}
	run := c.RunE
	c.RunE = func(cmd *cobra.Command, args []string) error {
		cmd.SilenceUsage = true
		return run(cmd, args)
	}
}

// addStrictFlag adds the --strict flag, which each command that skips records it can't process has, to c,
// setting *strict
func addStrictFlag(c *cobra.Command, strict *bool) {
	c.Flags().BoolVar(strict, "strict", false, "Stop with an error at the first record that can't be processed, instead of warning and skipping it")
}

// Execute builds the gofasta command and runs it with the program's arguments.
func Execute() {
	// cobra prints the error (and, for errors in the flags or arguments, the usage) to stderr itself;
	// stdout is reserved for data
	if err := newRootCmd().Execute(); err != nil {
		os.Exit(1)
	}
}
