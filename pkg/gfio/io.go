/*
Package gfio provides io functionality, including to/from stdin/stdout,
and helpful error messages for bad file paths (naming, for example, the
command-line option that gave the path)
*/
package gfio

import (
	"bufio"
	"compress/gzip"
	"errors"
	"io"
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

// IsStdin reports whether path means stdin: "stdin" or "-"
func IsStdin(path string) bool { return path == "stdin" || path == "-" }

// IsStdout reports whether path means stdout: "stdout" or "-"
func IsStdout(path string) bool { return path == "stdout" || path == "-" }

// OpenIn opens path for reading, or stdin if path is "stdin" or "-". Gzipped input (recognised by its first two
// bytes, so whatever the file is called) is decompressed as it is read; concatenated gzip members, as
// bgzip writes, are read as one stream. An uncompressed file is returned as its *os.File; a gzipped file
// can be moved back to its start with Seek(0, io.SeekStart). An error names label (e.g. the command-line
// flag that gave the path, such as "-q / --query").
func OpenIn(path, label string) (io.ReadCloser, error) {
	if IsStdin(path) {
		br := bufio.NewReader(os.Stdin)
		if isGzip(br) {
			zr, err := gzip.NewReader(br)
			if err != nil {
				return nil, errors.New("reading gzipped stdin: " + err.Error())
			}
			return readCloser{zr, os.Stdin}, nil
		}
		return readCloser{br, os.Stdin}, nil
	}
	f, err := os.Open(path)
	if err != nil {
		return nil, parseInErr(err, label)
	}
	br := bufio.NewReader(f)
	if !isGzip(br) {
		// nothing has been consumed from f that a caller needs: go back to its start and return it as it is
		if _, err := f.Seek(0, io.SeekStart); err != nil {
			f.Close()
			return nil, err
		}
		return f, nil
	}
	zr, err := gzip.NewReader(br)
	if err != nil {
		f.Close()
		return nil, errors.New("reading gzipped " + label + " " + path + ": " + err.Error())
	}
	return &gzipFile{f: f, br: br, zr: zr}, nil
}

// isGzip reports whether br starts with the gzip magic number
func isGzip(br *bufio.Reader) bool {
	b, err := br.Peek(2)
	return err == nil && b[0] == 0x1f && b[1] == 0x8b
}

// readCloser reads from one reader and closes another
type readCloser struct {
	io.Reader
	io.Closer
}

// gzipFile is a gzipped file being decompressed as it is read
type gzipFile struct {
	f  *os.File
	br *bufio.Reader
	zr *gzip.Reader
}

func (g *gzipFile) Read(p []byte) (int, error) { return g.zr.Read(p) }

func (g *gzipFile) Close() error { return g.f.Close() }

// Seek supports only moving back to the start of the decompressed data (offset 0, io.SeekStart)
func (g *gzipFile) Seek(offset int64, whence int) (int64, error) {
	if offset != 0 || whence != io.SeekStart {
		return 0, errors.New("a gzipped file can only be read again from its start")
	}
	if _, err := g.f.Seek(0, io.SeekStart); err != nil {
		return 0, err
	}
	g.br.Reset(g.f)
	if err := g.zr.Reset(g.br); err != nil {
		return 0, err
	}
	return 0, nil
}

// OpenOut creates path for writing, or returns stdout if path is "stdout" or "-".
func OpenOut(path string) (*os.File, error) {
	if IsStdout(path) {
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
