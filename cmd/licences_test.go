package cmd

import (
	"bytes"
	"testing"
)

const testNotice = "x is licensed under the MIT License; see LICENSE.\n\nIt also uses y.\n\n" +
	"--------------------------------------------------------------------------------\n" +
	"BSD 3-Clause License, used by:\n  y v1.0.0\n----\n\nfull text\n"

// licences prints the notice up to its first third-party notice; with -l, all of it
func TestLicences(t *testing.T) {
	saved := Notice
	Notice = testNotice
	defer func() { Notice = saved }()

	summary := "x is licensed under the MIT License; see LICENSE.\n\nIt also uses y.\n\nUse \"gofasta licences -l\" to print the full text, including all third-party notices.\n"
	for _, tc := range []struct {
		args []string
		want string
	}{
		{[]string{"licences"}, summary},
		{[]string{"licenses"}, summary},
		{[]string{"licences", "-l"}, testNotice},
		{[]string{"licences", "--long"}, testNotice},
	} {
		root := newRootCmd()
		out := new(bytes.Buffer)
		root.SetOut(out)
		root.SetArgs(tc.args)
		if err := root.Execute(); err != nil {
			t.Fatalf("%v: %v", tc.args, err)
		}
		if out.String() != tc.want {
			t.Errorf("%v: got\n%q\nwant\n%q", tc.args, out.String(), tc.want)
		}
	}
}
