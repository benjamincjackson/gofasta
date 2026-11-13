package cmd

import (
	"io"
	"os"

	"github.com/spf13/cobra"

	"github.com/virus-evolution/gofasta/pkg/gfio"
)

// flagLabel names a flag as errors about its value do: "-q / --query", or "--query" with no shorthand
func flagLabel(cmd *cobra.Command, name string) string {
	f := cmd.Flag(name)
	if f.Shorthand == "" {
		return "--" + f.Name
	}
	return "-" + f.Shorthand + " / --" + f.Name
}

// openIn opens the file given by the named flag for reading ("stdin" for stdin), decompressing it if it
// is gzipped
func openIn(cmd *cobra.Command, name string) (io.ReadCloser, error) {
	return gfio.OpenIn(cmd.Flag(name).Value.String(), flagLabel(cmd, name))
}

// openOut creates the file given by the named flag ("stdout" for stdout)
func openOut(cmd *cobra.Command, name string) (*os.File, error) {
	return gfio.OpenOut(cmd.Flag(name).Value.String())
}
