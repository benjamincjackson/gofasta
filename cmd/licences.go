package cmd

import (
	"fmt"
	"strings"

	"github.com/spf13/cobra"
)

// Notice is the content of the NOTICE file (the licence notices of the third-party software compiled
// into gofasta), set by package main, which embeds it
var Notice string

func newLicencesCmd() *cobra.Command {
	var long bool
	c := &cobra.Command{
		Use:     "licences",
		Aliases: []string{"licenses"},
		Short:   "Print licence information",
		Long:    "Print a summary of licence information, or with -l, the full text, including all third-party notices",
		Run: func(cmd *cobra.Command, args []string) {
			out := cmd.OutOrStdout()
			if long {
				fmt.Fprint(out, Notice)
				return
			}
			fmt.Fprintln(out, noticeSummary(Notice))
			fmt.Fprintln(out, "\nUse \"gofasta licences -l\" to print the full text, including all third-party notices.")
		},
	}
	c.Flags().BoolVarP(&long, "long", "l", false, "Print the full licence information, including all third-party notices")
	return c
}

// noticeSummary is the opening statement of the notice: everything before the first third-party notice
func noticeSummary(notice string) string {
	summary, _, _ := strings.Cut(notice, "\n---")
	return strings.TrimRight(summary, "\n")
}
