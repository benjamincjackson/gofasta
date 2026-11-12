package fasta

import (
	"bytes"
	"errors"
	"fmt"
	"io"
	"reflect"
	"strings"
	"testing"

	"github.com/virus-evolution/gofasta/pkg/encoding"
)

func TestGetAlignmentDims(t *testing.T) {
	alignmentData := []byte(
		`>Target1
ATGATC
>Target2
ATGATG
>Target3
ATTTTC
`)

	alignment := bytes.NewReader(alignmentData)

	n, l, err := getAlignmentDims(alignment)
	if err != nil {
		t.Error(err)
	}

	if n != 3 {
		t.Errorf("wrong number of sequences in TestGetAlignmentDims()")
	}

	if l != 6 {
		t.Errorf("wrong number of sequences in TestGetAlignmentDims()")
	}
}

// readAll returns every record from next, in order, or the first error
func readAll(next func() (Record, error)) ([]Record, error) {
	var records []Record
	for {
		record, err := next()
		if err == io.EOF {
			return records, nil
		}
		if err != nil {
			return records, err
		}
		records = append(records, record)
	}
}

func TestAlignmentSourceRoundTrip(t *testing.T) {
	alignmentData := `>Target1
ATGATC
>Target2
ATGATG
>Target3
ATTTTC
`
	records, err := readAll(NewAlignmentSource(strings.NewReader(alignmentData)))
	if err != nil {
		t.Fatal(err)
	}
	out := new(bytes.Buffer)
	for i, record := range records {
		if record.Idx != i {
			t.Errorf("record %d has Idx %d", i, record.Idx)
		}
		if err := WriteRecord(out, record, 0); err != nil {
			t.Fatal(err)
		}
	}
	if out.String() != alignmentData {
		t.Errorf("got\n%s\nwant\n%s", out.String(), alignmentData)
	}
}

func TestAlignmentSourceShortSeqs(t *testing.T) {
	alignmentData := `>Target1
ATGATC
>TargetShort1
ATTAT
>Target2
ATGATC
`
	_, err := readAll(NewAlignmentSource(strings.NewReader(alignmentData)))
	if !errors.Is(err, ErrDiffLenSeqs) {
		t.Errorf("want ErrDiffLenSeqs, got %v", err)
	}
}

func TestLoadEncodeAlignment(t *testing.T) {
	alignmentData := []byte(
		`>Target1
ATGATC
>Target2
ATGATG
>Target3
ATTTTC
`)

	alignment := bytes.NewReader(alignmentData)
	out := new(bytes.Buffer)

	EFRs, err := LoadEncodeAlignment(alignment, false, false, false)
	if err != nil {
		t.Error(err)
	}

	DA := encoding.MakeDecodingArray()

	for _, FR := range EFRs {
		out.Write([]byte(">" + FR.ID + "\n"))
		for _, nuc := range FR.Seq {
			out.Write([]byte(DA[nuc]))
		}
		out.Write([]byte("\n"))
	}

	if string(out.Bytes()) != `>Target1
ATGATC
>Target2
ATGATG
>Target3
ATTTTC
` {
		t.Errorf("problem in TestLoadEncodeAlignment()")
		fmt.Print(string(out.Bytes()))
	}
}

func TestLoadEncodeAlignmentShortSeqs(t *testing.T) {
	alignmentData := []byte(
		`>Target1
ATGATC
>Target2
ATGAT
>Target3
ATTTTC
`)

	alignment := bytes.NewReader(alignmentData)

	_, err := LoadEncodeAlignment(alignment, false, false, false)
	if !errors.Is(err, ErrDiffLenSeqs) {
		t.Errorf("want ErrDiffLenSeqs, got %v", err)
	}
}

func TestLoadEncodeScoreAlignment(t *testing.T) {
	alignmentData := []byte(
		`>Target1
ATGATN
>Target2
ATGANN
>Target3
ATTR-N
`)

	FRs, err := LoadEncodeAlignment(bytes.NewReader(alignmentData), false, true, true)
	if err != nil {
		t.Fatal(err)
	}

	DA := encoding.MakeDecodingArray()
	alOut := new(bytes.Buffer)
	scoreOut := make([]int64, 0)
	for _, FR := range FRs {
		alOut.Write([]byte(">" + FR.ID + "\n"))
		for _, nuc := range FR.Seq {
			alOut.Write([]byte(DA[nuc]))
		}
		alOut.Write([]byte("\n"))
		scoreOut = append(scoreOut, FR.Score)
	}

	if string(alOut.Bytes()) != `>Target1
ATGATN
>Target2
ATGANN
>Target3
ATTR-N
` {
		t.Errorf("problem in TestLoadEncodeScoreAlignment() (alignment)")
	}

	if !reflect.DeepEqual(scoreOut, []int64{63, 54, 48}) {
		t.Errorf("problem in TestLoadEncodeScoreAlignment() (score): got %v", scoreOut)
	}

	counts := [][]int{{2, 2, 1, 0}, {2, 1, 1, 0}, {1, 2, 0, 0}}
	for i := range FRs {
		got := []int{FRs[i].Count_A, FRs[i].Count_T, FRs[i].Count_G, FRs[i].Count_C}
		if !reflect.DeepEqual(got, counts[i]) {
			t.Errorf("problem in TestLoadEncodeScoreAlignment() (counts): record %d got %v, want %v", i, got, counts[i])
		}
	}
}

func TestWriteRecordWrap(t *testing.T) {
	record := Record{ID: "Seq1", Seq: "ATGATGATG"}
	for _, tc := range []struct {
		wrap int
		want string
	}{
		{0, ">Seq1\nATGATGATG\n"},
		{-1, ">Seq1\nATGATGATG\n"},
		{80, ">Seq1\nATGATGATG\n"},
		{9, ">Seq1\nATGATGATG\n"},
		{5, ">Seq1\nATGAT\nGATG\n"},
		{3, ">Seq1\nATG\nATG\nATG\n"},
	} {
		out := new(bytes.Buffer)
		if err := WriteRecord(out, record, tc.wrap); err != nil {
			t.Fatal(err)
		}
		if out.String() != tc.want {
			t.Errorf("wrap %d: got %q, want %q", tc.wrap, out.String(), tc.want)
		}
	}
}
