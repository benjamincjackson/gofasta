/*
Package gfio provides io functionality, including to/from stdin/stdout,
and helpful error messages for bad file paths (naming, for example, the
command-line option that gave the path)
*/
package gfio

import (
	"errors"
	"io/fs"
	"os"
)

// parseInErr returns a modified *fs.PathError-style error, augmented with the name of what was being
// opened (e.g. the command-line flag that gave the path)
func parseInErr(err error, label string) error {
	switch x := err.(type) {
	case *fs.PathError:
		return errors.New(x.Op + " " + label + " " + x.Path + ": " + x.Err.Error())
	default:
		return err
	}
}

// OpenIn opens path for reading, or returns stdin if path is "stdin". An error names label (e.g. the
// command-line flag that gave the path, such as "-q / --query").
func OpenIn(path, label string) (*os.File, error) {
	if path == "stdin" {
		return os.Stdin, nil
	}
	f, err := os.Open(path)
	if err != nil {
		return f, parseInErr(err, label)
	}
	return f, nil
}

// OpenOut creates path for writing, or returns stdout if path is "stdout".
func OpenOut(path string) (*os.File, error) {
	if path == "stdout" {
		return os.Stdout, nil
	}
	return os.Create(path)
}

// CloseOut closes an output file and, if the close fails (the final flush of buffered data to a full disk,
// for example, is only reported by Close), stores the error in *err unless an error is already there.
// Use it as: defer gfio.CloseOut(out, &err), with err a named return value.
func CloseOut(f *os.File, err *error) {
	if cerr := f.Close(); cerr != nil && *err == nil {
		*err = cerr
	}
}
