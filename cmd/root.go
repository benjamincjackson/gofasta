package cmd

import (
	"errors"
	"strings"

	"github.com/spf13/pflag"

	"os"

	"github.com/virus-evolution/gofasta/pkg/gfio"

	"github.com/spf13/cobra"
)

// newRootCmd builds the gofasta command and all its subcommands, each with its flags at their defaults
func newRootCmd() *cobra.Command {
	root := &cobra.Command{
		Use:   "gofasta",
		Short: "Command-line utilities for genomic epidemiology research",
		Long: `Command-line utilities for genomic epidemiology research

Files are named by their paths; "stdin" or "-" reads an input from stdin, and "stdout" or "-" writes an output
to stdout. Any input can be gzipped. (--annotation can't be stdin, as its format comes from its suffix: use
--genbank or --gff.)

If you use gofasta in your work, please cite:

Jackson B (2022). gofasta: command-line utilities for genomic epidemiology research. Bioinformatics 38 (16), 4033-4035
https://doi.org/10.1093/bioinformatics/btac424`,
		Version: "2.0.0pre",
	}
	root.AddCommand(newClosestCmd(), newSNPsCmd(), newVariantsCmd(), newSamCmd(), newUpdownCmd(), newLicencesCmd())
	root.CompletionOptions.DisableDefaultCmd = true
	requiredFlagsFirst(root)
	silenceUsageOnceRunning(root)
	return root
}

// requiredFlagsFirst makes each command with a PreRunE check its required flags (MarkFlagRequired), and
// that at most one input is stdin (oneStdin), before the checks in PreRunE: cobra checks them only after PreRunE, so a check there that reads a missing flag
// would report something less clear than the flag being missing.
func requiredFlagsFirst(c *cobra.Command) {
	for _, sub := range c.Commands() {
		requiredFlagsFirst(sub)
	}
	if c.PreRunE == nil {
		return
	}
	pre := c.PreRunE
	c.PreRunE = func(cmd *cobra.Command, args []string) error {
		if err := cmd.ValidateRequiredFlags(); err != nil {
			return err
		}
		if err := oneStdin(cmd); err != nil {
			return err
		}
		return pre(cmd, args)
	}
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

// notInputs are the string flags that don't name an input file: outputs, and values that are not file
// names; every other string flag names an input
var notInputs = map[string]bool{"outfile": true, "outpath": true, "reference-id": true, "measure": true, "max-dist": true}

// oneStdin checks that at most one of the command's inputs is stdin ("stdin" or "-", given or by
// default): every string flag except those in notInputs counts
func oneStdin(cmd *cobra.Command) error {
	stdin := make([]string, 0)
	cmd.InheritedFlags() // merges the parent's persistent flags into cmd.Flags()
	cmd.Flags().VisitAll(func(f *pflag.Flag) {
		if f.Value.Type() == "string" && !notInputs[f.Name] && gfio.IsStdin(f.Value.String()) {
			stdin = append(stdin, "--"+f.Name)
		}
	})
	switch {
	case len(stdin) == 2:
		return errors.New("only one input can be read from stdin, but " + stdin[0] + " and " + stdin[1] + " are both stdin")
	case len(stdin) > 2:
		return errors.New("only one input can be read from stdin, but " + strings.Join(stdin[:len(stdin)-1], ", ") + " and " + stdin[len(stdin)-1] + " are all stdin")
	}
	return nil
}

// Execute builds the gofasta command and runs it with the program's arguments.
func Execute() {
	// cobra prints the error (and, for errors in the flags or arguments, the usage) to stderr itself;
	// stdout is reserved for data
	if err := newRootCmd().Execute(); err != nil {
		os.Exit(1)
	}
}
