package sam

import (
	"context"
	"errors"
	"io"
	"os"

	"github.com/virus-evolution/gofasta/pkg/encoding"
	"github.com/virus-evolution/gofasta/pkg/fasta"
	"github.com/virus-evolution/gofasta/pkg/pipeline"
	"github.com/virus-evolution/gofasta/pkg/variants"
)

// Variants annotates amino acid, insertion, deletion, and nucleotide (anything
// outside of codons with an amino acid change) mutations relative to a reference
// sequence from pairwise alignments in sam format. Genome annotations are
// derived from a annotation file in genbank or gff version 3 format
//
// A query whose records can't be used is reported on stderr and skipped, unless strict, when it stops the run.
// threads is the number of goroutines that process records (fewer than 1 means 1).
func Variants(samIn, refIn io.Reader, refFromFile bool, annoIn io.Reader, annoSuffix string, out io.Writer, start, end int, aggregate bool, threshold float64, appendSNP bool, threads int, strict bool) error {

	var ref fasta.EncodedRecord
	if refFromFile {
		refs, err := fasta.LoadEncodeAlignment(refIn, false, false, false)
		if err != nil {
			return err
		}
		if len(refs) > 1 {
			return errors.New("more than one record in --reference")
		}
		ref = refs[0]
	}

	ref, cdsregions, intregions, err := variants.LoadAnnotation(annoIn, annoSuffix, ref)
	if err != nil {
		return err
	}

	// the output header is written before the sam file is read
	var lw *variants.LineWriter
	var agg *variants.Aggregator
	if aggregate {
		agg, err = variants.NewAggregator(out, start, end, appendSNP, ref.ID)
	} else {
		lw, err = variants.NewLineWriter(out, start, end, appendSNP, ref.ID)
	}
	if err != nil {
		return err
	}

	src, err := newSamSource(samIn)
	if err != nil {
		return err
	}
	refSeq := ref.Decode().Seq
	if err := src.checkRefLength(len(refSeq)); err != nil {
		return err
	}

	EA := encoding.MakeEncodingArray()
	refBytes := []byte(refSeq) // shared by the workers, which only read it
	work := func(group samRecords) (variants.AnnoStructs, error) {
		group, err := src.prepareQuery(group, refBytes, len(refBytes))
		if err != nil {
			return variants.AnnoStructs{}, err
		}
		pair, err := blockToPairwiseAlignment(group, refBytes, false)
		if err != nil {
			return variants.AnnoStructs{}, pipeline.Named(group.records[0].Name, err)
		}
		for i, nuc := range pair.query {
			pair.query[i] = EA[nuc]
		}
		for i, nuc := range pair.ref {
			pair.ref[i] = EA[nuc]
		}
		offsetRefCoord := variants.RefToMSAOffsets(pair.ref)
		AS, err := variants.GetVariantsPair(pair.ref, pair.query, pair.refname, pair.queryname, pair.idx, cdsregions, intregions, offsetRefCoord)
		return AS, pipeline.Named(pair.queryname, err)
	}

	write := lw.Write
	if aggregate {
		write = agg.Add
	}
	h := pipeline.NewHandler(strict, os.Stderr, queryNoun, write)
	if err := pipeline.Run(context.Background(), threads, src.next, work, withNotes(src, h.Sink)); err != nil {
		return err
	}
	src.writeTail()
	if err := h.Finish(); err != nil {
		return err
	}
	if aggregate {
		return agg.Write(threshold)
	}
	return nil
}
