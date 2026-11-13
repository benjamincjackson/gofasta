package cmd

import (
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/spf13/cobra"
	"github.com/spf13/pflag"
)

// findCmd returns the command at path (e.g. "sam variants") in a newly built gofasta command, so its
// flags are at their defaults
func findCmd(t *testing.T, path string) *cobra.Command {
	t.Helper()
	c, _, err := newRootCmd().Find(strings.Fields(path))
	if err != nil || c.CommandPath() != "gofasta "+path {
		t.Fatalf("no command %q", path)
	}
	return c
}

// preRun sets the given flags of a new command at path and runs its PreRunE. A flag marked required that
// is not in set is given a placeholder value, so each case tests only what it sets.
func preRun(t *testing.T, path string, set map[string]string) error {
	t.Helper()
	c := findCmd(t, path)
	c.InheritedFlags() // merges the parent's persistent flags into c.Flags(), so they can be set
	c.Flags().VisitAll(func(f *pflag.Flag) {
		if _, ok := set[f.Name]; !ok && len(f.Annotations[cobra.BashCompOneRequiredFlag]) > 0 {
			if err := c.Flags().Set(f.Name, "required-placeholder"); err != nil {
				t.Fatalf("setting --%s: %v", f.Name, err)
			}
		}
	})
	for k, v := range set {
		if err := c.Flags().Set(k, v); err != nil {
			t.Fatalf("setting --%s: %v", k, err)
		}
	}
	return c.PreRunE(c, nil)
}

func TestFlagValidation(t *testing.T) {
	dir := t.TempDir()
	f := filepath.Join(dir, "in.fa")
	if err := os.WriteFile(f, []byte(">a\nA\n"), 0o644); err != nil {
		t.Fatal(err)
	}
	other := filepath.Join(dir, "out.csv")

	type tc struct {
		path string
		set  map[string]string
	}

	bad := []tc{
		{"closest", map[string]string{"threads": "-1"}},
		{"closest", map[string]string{"number": "-3"}},
		{"closest", map[string]string{"query": f, "outfile": f}},
		{"closest", map[string]string{"max-dist": "-1"}},
		{"closest", map[string]string{"max-dist": "-0.5"}},
		{"closest", map[string]string{"max-dist": "nan"}},
		{"closest", map[string]string{"max-dist": "x"}},
		{"snps", map[string]string{"threshold": "-0.5"}},
		{"snps", map[string]string{"reference": f, "outfile": f}},
		{"variants", map[string]string{"threads": "0"}},
		{"variants", map[string]string{"threshold": "-0.2"}},
		{"variants", map[string]string{"msa": f, "outfile": f}},
		{"variants", map[string]string{"start": "9", "end": "5"}},
		{"sam variants", map[string]string{"start": "9", "end": "5"}},
		{"sam variants", map[string]string{"threads": "0"}},
		{"sam variants", map[string]string{"threshold": "-0.5"}},
		{"sam variants", map[string]string{"reference": f, "outfile": f}},
		{"sam toMultiAlign", map[string]string{"samfile": f, "fasta-out": f}},
		{"sam toPairAlign", map[string]string{"threads": "-1"}},
		{"updown topranking", map[string]string{"size-up": "-2"}},
		{"updown topranking", map[string]string{"dist-push": "-1"}},
		{"updown topranking", map[string]string{"threshold-pair": "-0.1"}},
		{"updown topranking", map[string]string{"ignore": f, "outfile": f}},
		{"updown topranking", map[string]string{"query": "q.csv", "target": "t.csv"}}, // no --size-* or --dist-*
		{"updown topranking", map[string]string{"size-up": "2", "query": "q.txt", "target": "t.csv"}},
		{"updown topranking", map[string]string{"size-up": "2", "query": "q.fasta", "target": "t.csv"}}, // no --reference
		{"closest", map[string]string{"measure": "bogus"}},
		{"updown list", map[string]string{"query": f, "outfile": f}},
		{"updown list", map[string]string{}},     // no --reference
		{"sam toPairAlign", map[string]string{}}, // no --reference
	}
	for _, b := range bad {
		if err := preRun(t, b.path, b.set); err == nil {
			t.Errorf("%s %v: expected an error", b.path, b.set)
		}
	}

	good := []tc{
		{"closest", map[string]string{}},
		{"closest", map[string]string{"threads": "0"}},
		{"closest", map[string]string{"number": "5", "table": "true"}},
		{"closest", map[string]string{"table": "true"}}, // ignored without -n/-d
		{"closest", map[string]string{"max-dist": "0"}},
		{"closest", map[string]string{"max-dist": "inf"}},
		{"closest", map[string]string{"query": f, "outfile": other}},
		{"snps", map[string]string{"threshold": "0.5"}},
		{"variants", map[string]string{"start": "5", "end": "5", "threshold": "0.1"}},
		{"variants", map[string]string{"start": "5"}},
		{"variants", map[string]string{"end": "5"}},
		{"sam variants", map[string]string{"start": "5"}},
		{"variants", map[string]string{"reference": f, "outfile": f}}, // --reference is a record ID here
		{"sam variants", map[string]string{"threads": "4", "reference": f, "outfile": other}},
		{"sam toMultiAlign", map[string]string{}},
		{"updown topranking", map[string]string{"size-up": "2", "threshold-pair": "0.2", "query": "q.csv", "target": "t.csv"}},
		{"updown topranking", map[string]string{"size-total": "2", "query": "q.fasta", "target": "t.csv", "reference": "r.fasta"}},
		{"closest", map[string]string{"measure": "TN93"}},
		{"updown list", map[string]string{"reference": f}},
		{"sam toPairAlign", map[string]string{"reference": f}},
	}
	for _, g := range good {
		if err := preRun(t, g.path, g.set); err != nil {
			t.Errorf("%s %v: unexpected error %v", g.path, g.set, err)
		}
	}
}

func TestAnnotationFormat(t *testing.T) {
	for path, want := range map[string]string{
		"a.gb": "gb", "a.gbk": "gb", "a.GBFF": "gb", "dir.x/a.gff": "gff", "a.gff3": "gff", "A.GFF3": "gff",
		"a.gff.gz": "gff", "a.gb.GZ": "gb",
	} {
		if got, err := annotationFormat(path); err != nil || got != want {
			t.Errorf("annotationFormat(%q) = %q, %v; want %q", path, got, err, want)
		}
	}
	for _, path := range []string{"a.gtf", "a", "a.gtf.gz", "a.gz", ""} {
		if _, err := annotationFormat(path); err == nil {
			t.Errorf("annotationFormat(%q): want an error", path)
		}
	}
}
