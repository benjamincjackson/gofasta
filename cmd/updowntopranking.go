package cmd

import (
	"bufio"
	"errors"
	"io"
	"path/filepath"
	"strings"

	"github.com/spf13/cobra"

	"github.com/virus-evolution/gofasta/pkg/gfio"
	"github.com/virus-evolution/gofasta/pkg/updown"
)

// toprankingOptions holds the values of the flags of the command newToprankingCmd returns
type toprankingOptions struct {
	threads         int
	query           string
	target          string
	ignore          string
	outfile         string
	table           bool
	sizeTotal       int
	sizeUp          int
	sizeDown        int
	sizeSide        int
	sizeSame        int
	noFill          bool
	distAll         int
	distUp          int
	distDown        int
	distSide        int
	distPush        int
	thresholdPair   float32
	thresholdTarget int
	strict          bool
}

// newToprankingCmd builds the updown topranking command; u holds the updown command's flags, which its subcommands share
func newToprankingCmd(u *updownOptions) *cobra.Command {
	var o toprankingOptions
	c := &cobra.Command{
		Use:   "topranking",
		Short: "Get pseudo-tree-aware catchments for query sequences from alignments",
		Long: `Get pseudo-tree-aware catchments for query sequences from alignments

Example usage:
	gofasta updown topranking --query smallquery.fasta -r WH04.fasta --target mutationlist.csv --size-total 1000 -o catchment.csv

For each sequence in --query, this routine finds the closest sequences by SNP-distance in --target, binned
according to whether they are likely children, parents, or siblings of, or on a polytomy with, the query
sequence. It does this by comparing SNPs relative to a common reference sequence which is imagined to be the
root of the tree.

--query and --target can either be alignments in fasta format, or the CSV output of gofasta updown list, or
one of each. They must have file extensions .csv, .fasta, .fa, .fna or .fas (in any case, optionally followed
by .gz). If either is an alignment, you must provide --reference, and this should be the same sequence that
was used by gofasta updown list.

Use the --dist flags to filter on SNP-distances in each direction. As long as you haven't also used any --size
flags, the program will return all the targets that are equal or less than the specified SNP-distance(s) away.
If --dist-push is invoked with an integer (i) argument, the program will push the SNP-distance boundaries to
cover the all sequences that are the closest i SNP-distances away, and then it will return all the neighbours
at those distances in that bin.

If you provide a number to --size-total, the output will try to include this many closest sequences to the
query, split evenly between the four bins. If one or more bins has a shortfall, the sizes of the other bins
will increase until --size-total is met, if possible, unless you use --no-fill.

You can provide --size-up, -down, -side and -same, instead of --size-total, if you want different proportions
for each bin. The program will aim to provide the sum of these numbers in total in the output, and will make
up for a shortfall in one bin by increasing the count of the other bins where possible, unless --no-fill.

You can combine the two types of flag (size and dist), to return only the closest n sequences under a set
distance (as long as you haven't also invoked --dist-push).
`,
		PreRunE: func(cmd *cobra.Command, args []string) error {
			return firstError(
				intAtLeast(cmd, 0, "threads", "size-total", "size-up", "size-down", "size-side", "size-same",
					"dist-all", "dist-up", "dist-down", "dist-side", "dist-push", "threshold-target"),
				floatAtLeast(cmd, 0, "threshold-pair"),
				outputNotInput(cmd, "outfile", "query", "target", "reference", "ignore"),
				updown.CheckTopRankingArgs(o.sizeTotal, o.sizeUp, o.sizeDown, o.sizeSide, o.sizeSame,
					o.distAll, o.distUp, o.distDown, o.distSide, o.distPush),
				toprankingInputTypesErr(&o, u.reference),
			)
		},
		RunE: func(cmd *cobra.Command, args []string) (err error) {

			qtype, ttype, err := toprankingInputTypes(&o, u.reference)
			if err != nil {
				return err
			}

			// targets to ignore (potentially):
			ignoreArray := make([]string, 0)
			if len(o.ignore) != 0 {
				f, err := openIn(cmd, "ignore")
				if err != nil {
					return err
				}
				defer f.Close()
				s := bufio.NewScanner(f)
				for s.Scan() {
					// one ID per line; surrounding whitespace (including the \r of a CRLF line end) is not part of it
					if id := strings.TrimSpace(s.Text()); id != "" {
						ignoreArray = append(ignoreArray, id)
					}
				}
				err = s.Err()
				if err != nil {
					return err
				}
			}

			query, err := openIn(cmd, "query")
			if err != nil {
				return err
			}
			defer query.Close()

			target, err := openIn(cmd, "target")
			if err != nil {
				return err
			}
			defer target.Close()

			var ref io.ReadCloser // only read if a fasta file is given
			if qtype == "fasta" || ttype == "fasta" {
				ref, err = openIn(cmd, "reference")
				if err != nil {
					return err
				}
				defer ref.Close()
			}

			out, err := openOut(cmd, "outfile")
			if err != nil {
				return err
			}
			defer gfio.CloseOut(out, &err)

			err = updown.TopRanking(query, target, ref, out, o.table,
				qtype, ttype, ignoreArray,
				o.sizeTotal, o.sizeUp, o.sizeDown, o.sizeSide, o.sizeSame,
				o.distAll, o.distUp, o.distDown, o.distSide,
				o.thresholdPair, o.thresholdTarget, o.noFill, o.distPush, resolveThreads(o.threads), o.strict)

			return
		},
	}

	c.Flags().StringVarP(&o.query, "query", "", "", "File with sequences to find neighbours for. Either the CSV output of gofasta updown list, or an alignment in fasta format")
	c.Flags().StringVarP(&o.target, "target", "", "", "File of sequences to look for neighbours in. Either the CSV output of gofasta updown list, or an alignment in fasta format")
	c.Flags().StringVarP(&o.outfile, "outfile", "o", "stdout", "CSV-format file of closest neighbours to write")
	c.Flags().BoolVarP(&o.table, "table", "", false, "Write a long-form table of the output")
	c.Flags().StringVarP(&u.reference, "reference", "r", "", "Reference sequence, in fasta format - only required if --query and --target are fasta files")
	c.Flags().StringVarP(&o.ignore, "ignore", "", "", "Optional plain text file of IDs to ignore in the target file when searching for neighbours")

	c.Flags().IntVarP(&o.distAll, "dist-all", "", 0, "Maximum allowed SNP-distance between target and query sequence in any direction. Overrides the settings below")
	c.Flags().IntVarP(&o.distUp, "dist-up", "", 0, "Maximum allowed SNP-distance from query for sequences in the parent bin")
	c.Flags().IntVarP(&o.distDown, "dist-down", "", 0, "Maximum allowed SNP-distance from query for sequences in the child bin")
	c.Flags().IntVarP(&o.distSide, "dist-side", "", 0, "Maximum allowed SNP-distance from query for sequences in the sibling bin")

	c.Flags().IntVarP(&o.sizeTotal, "size-total", "", 0, "Max number of neighbours to find (attempts to split equally between same/up/down/side). A hard limit")
	c.Flags().IntVarP(&o.sizeUp, "size-up", "", 0, "Max number of closest parent sequences to find, if size-total not specified. A soft limit unless --no-fill")
	c.Flags().IntVarP(&o.sizeDown, "size-down", "", 0, "Max number of closest child sequences to find, if size-total not specified. A soft limit unless --no-fill")
	c.Flags().IntVarP(&o.sizeSide, "size-side", "", 0, "Max number of closest sibling sequences to find, if size-total not specified. A soft limit unless --no-fill")
	c.Flags().IntVarP(&o.sizeSame, "size-same", "", 0, "Max number of identical sequences to find, if size-total not specified. A soft limit unless --no-fill")

	c.Flags().Float32VarP(&o.thresholdPair, "threshold-pair", "", 0.1, "Up to this proportion of consequential sites is allowed to be ambiguous in either sequence for each pairwise comparison")
	c.Flags().IntVarP(&o.thresholdTarget, "threshold-target", "", 10000, "Target can have at most this number of ambiguities to be considered")

	c.Flags().IntVarP(&o.distPush, "dist-push", "", 0, "Push the --dist boundaries outwards so that bins have at least these many closest SNP-distances for which there are neighbours, where possible")
	c.Flags().BoolVarP(&o.noFill, "no-fill", "", false, "Don't make up for a shortfall in any of --size-up, -down, -side or -same by increasing the count for other bins")

	c.Flags().Lookup("table").NoOptDefVal = "true"
	c.Flags().Lookup("no-fill").NoOptDefVal = "true"

	addThreadsFlag(c.Flags(), &o.threads)
	addStrictFlag(c, &o.strict)
	c.MarkFlagRequired("query")
	c.MarkFlagRequired("target")
	c.Flags().SortFlags = false
	return c
}

// toprankingInputTypes returns the types ("csv" or "fasta") of --query and --target, from their
// extensions (in any case), and checks that reference (--reference) is given if either is a fasta file
func toprankingInputTypes(o *toprankingOptions, reference string) (string, string, error) {
	inputType := func(path, flag string) (string, error) {
		switch strings.ToLower(filepath.Ext(strings.TrimSuffix(strings.ToLower(path), ".gz"))) {
		case ".csv":
			return "csv", nil
		case ".fasta", ".fa", ".fna", ".fas":
			return "fasta", nil
		}
		return "", errors.New("couldn't tell if --" + flag + " was a .csv or a fasta (.fasta, .fa, .fna, .fas) file (optionally gzipped, .gz)")
	}
	qtype, err := inputType(o.query, "query")
	if err != nil {
		return "", "", err
	}
	ttype, err := inputType(o.target, "target")
	if err != nil {
		return "", "", err
	}
	if (qtype == "fasta" || ttype == "fasta") && len(reference) == 0 {
		return "", "", errors.New("if your either of your input files are fastas, you must provide a --reference")
	}
	return qtype, ttype, nil
}

func toprankingInputTypesErr(o *toprankingOptions, reference string) error {
	_, _, err := toprankingInputTypes(o, reference)
	return err
}
