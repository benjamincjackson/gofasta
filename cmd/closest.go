package cmd

import (
	"errors"
	"runtime"
	"strconv"
	"strings"

	"github.com/spf13/cobra"

	"github.com/virus-evolution/gofasta/pkg/closest"
	"github.com/virus-evolution/gofasta/pkg/gfio"
)

// closestOptions holds the values of the flags of the command newClosestCmd returns
type closestOptions struct {
	threads int
	query   string
	target  string
	outfile string
	n       int
	dist    string
	measure string
	table   bool
	strict  bool
}

func newClosestCmd() *cobra.Command {
	var o closestOptions
	c := &cobra.Command{
		Use:   "closest",
		Short: "Find the closest sequence(s) to a query by genetic distance",
		Long: `Find the closest sequence(s) to a query by genetic distance

The usecase is imagined to be a small number of query sequences, whose nearest neighbours
by genetic distance need to be found amongst a large number of target sequences. The query 
alignment is read into memory and the target alignment is streamed from disk and iterated 
over once, so it can be arbitrarily large.

You can find the single closest neighbour like:

	gofasta closest -t 2 --query query.fasta --target target.fasta -o closest.csv

and the output will be a CSV format file with the headers query, closest, distance, SNPs.

Or you can find the nearest n neighbours:

	gofasta closest -t 2 -n 1000 --query query.fasta --target target.fasta -o closest.n1000.csv

or you can find all the neighbours under a -d distance away:

	gofasta closest -t 2 --measure snp -d 5 --query query.fasta --target target.fasta -o closest.d5snps.csv

and the default output will be a CSV format file with the headers query, closest.  The 'closest' column
is a ";"-delimited list of neighbours, closest first. Ties for distance are broken by genome completeness.

You can combine -d with -n to find the nearest n neighbours less than or equal to a distance, d.

Possible measures of distance are raw number of nucleotide changes per site (the default, raw), raw number
of nucleotide changes in total (snp), or Tamura and Nei's 1993 evolutionary distance (tn93).

Use --table in combination with the -n and/or -d flags to write a long-form output including the distance
between every pair.
`,
		PreRunE: func(cmd *cobra.Command, args []string) error {
			return firstError(
				intAtLeast(cmd, 0, "threads", "number"),
				limitAtLeast(cmd, 0, "max-dist"),
				outputNotInput(cmd, "outfile", "query", "target"),
				measureErr(o.measure),
			)
		},
		RunE: func(cmd *cobra.Command, args []string) (err error) {

			queryIn, err := openIn(cmd, "query")
			if err != nil {
				return err
			}
			defer queryIn.Close()

			targetIn, err := openIn(cmd, "target")
			if err != nil {
				return err
			}
			defer targetIn.Close()

			measure, err := closestMeasureName(o.measure)
			if err != nil {
				return err
			}

			dist := -1.0
			if o.dist != "" {
				dist, err = strconv.ParseFloat(o.dist, 64)
				if err != nil {
					return err
				}
			}

			closestOut, err := openOut(cmd, "outfile")
			if err != nil {
				return err
			}
			defer gfio.CloseOut(closestOut, &err)

			// -t limits the CPUs gofasta uses, not just the number of workers
			if o.threads > 0 && o.threads < runtime.NumCPU() {
				runtime.GOMAXPROCS(o.threads)
			}

			if o.n > 0 || dist != -1.0 {
				err = closest.ClosestN(o.n, dist, queryIn, targetIn, measure, closestOut, o.table, o.threads, o.strict)
			} else {
				err = closest.Closest(queryIn, targetIn, measure, closestOut, o.threads, o.strict)
			}

			return err
		},
	}

	c.Flags().IntVarP(&o.threads, "threads", "t", 0, "Number of CPUs to use (Default: all available CPUs)")
	c.Flags().StringVarP(&o.query, "query", "", "", "Alignment of sequences to find neighbours for, in fasta format")
	c.Flags().StringVarP(&o.target, "target", "", "", "Alignment of sequences to search for neighbours in, in fasta format")
	c.Flags().StringVarP(&o.measure, "measure", "m", "raw", "Which distance measure to use (raw, snp or tn93)")
	c.Flags().IntVarP(&o.n, "number", "n", 0, "(Optional) the closest n sequences to each query will be returned")
	c.Flags().StringVarP(&o.dist, "max-dist", "d", "", "(Optional) return all sequences less than or equal to this distance away (inf for no limit)")
	c.Flags().StringVarP(&o.outfile, "outfile", "o", "stdout", "The output file to write")
	c.Flags().BoolVarP(&o.table, "table", "", false, "Write a long-form table of the output")

	addStrictFlag(c, &o.strict)
	c.MarkFlagRequired("query")
	c.MarkFlagRequired("target")
	c.Flags().SortFlags = false
	return c
}

// closestMeasureName returns the distance measure measure (-m) names: raw, snp or tn93, in any case
func closestMeasureName(measure string) (string, error) {
	switch m := strings.ToLower(measure); m {
	case "raw", "snp", "tn93":
		return m, nil
	}
	return "", errors.New("Couldn't tell which distance --measure / -m to use (choose one of \"raw\", \"snp\" or \"tn93\")")
}

func measureErr(measure string) error {
	_, err := closestMeasureName(measure)
	return err
}
