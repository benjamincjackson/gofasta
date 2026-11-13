package genbank

import (
	"errors"
	"reflect"
	"strings"
	"testing"
	"testing/iotest"
)

// gbRecord wraps feature table lines in a minimal one-record GenBank file
func gbRecord(featureLines ...string) string {
	return "LOCUS       X 30 bp\nFEATURES             Location/Qualifiers\n" + strings.Join(featureLines, "\n") +
		"\nORIGIN\n        1 atgaaataaa tgaaataaat gaaataaatg\n//\n"
}

func TestReadGenBankFeatureTable(t *testing.T) {
	gb, err := ReadGenBank(strings.NewReader(gbRecord(
		"     CDS             join(1..3,",
		"                     4..9)",
		"                     /gene=\"a\"",
		"                     /note=\"x=y, as in \"\"quoted\"\" text that wraps",
		"                     onto a second line\"",
		"                     /db_xref=\"GeneID:1\"",
		"                     /db_xref=\"GeneID:2\"",
		"                     /transl_except=(pos:complement(4..6),",
		"                     aa:TERM)",
		"                     /pseudo",
		"                     /codon_start=1",
		"                     /translation=\"MK",
		"                     ML\"",
		"     gene            complement(10..18)",
		"                     /gene=\"b\"",
		"                     /ribosomal_slippage",
	)))
	if err != nil {
		t.Fatal(err)
	}
	want := []GenbankFeature{
		{Feature: "CDS", Location: Location{Representation: "join(1..3,4..9)"}, Info: map[string][]string{
			"gene":          {"a"},
			"note":          {`x=y, as in "quoted" text that wraps onto a second line`},
			"db_xref":       {"GeneID:1", "GeneID:2"},
			"transl_except": {"(pos:complement(4..6),aa:TERM)"},
			"pseudo":        {""},
			"codon_start":   {"1"},
			"translation":   {"MKML"},
		}},
		{Feature: "gene", Location: Location{Representation: "complement(10..18)"}, Info: map[string][]string{
			"gene":               {"b"},
			"ribosomal_slippage": {""},
		}},
	}
	if !reflect.DeepEqual(gb.FEATURES, want) {
		t.Errorf("got\n%#v\nwant\n%#v", gb.FEATURES, want)
	}
	f := gb.FEATURES[0]
	if f.Value("db_xref") != "GeneID:1" || f.Value("absent") != "" || !f.HasAttribute("pseudo") {
		t.Errorf("Value/HasAttribute: got %q, %q, %v", f.Value("db_xref"), f.Value("absent"), f.HasAttribute("pseudo"))
	}
	if string(gb.ORIGIN) != "atgaaataaatgaaataaatgaaataaatg" {
		t.Errorf("got ORIGIN %q", gb.ORIGIN)
	}
}

// A file must hold one record
func TestReadGenBankMoreThanOneRecord(t *testing.T) {
	two := gbRecord("     CDS             1..9", "                     /gene=\"a\"")
	if _, err := ReadGenBank(strings.NewReader(two + "\n" + two)); !errors.Is(err, errMultipleRecords) {
		t.Errorf("two records: want errMultipleRecords, got %v", err)
	}
	if _, err := ReadGenBank(strings.NewReader(two + "\n  \n")); err != nil {
		t.Errorf("blank lines after the record: unexpected error %v", err)
	}
}

// Lines longer than bufio.Scanner's default limit are read, and read errors are reported
func TestReadGenBankLongLinesAndReadErrors(t *testing.T) {
	long := strings.Repeat("M", 100_000)
	gb, err := ReadGenBank(strings.NewReader(gbRecord("     CDS             1..9", "                     /translation=\""+long+"\"")))
	if err != nil {
		t.Fatal(err)
	}
	if got := gb.FEATURES[0].Value("translation"); got != long {
		t.Errorf("got a translation of length %d, want %d", len(got), len(long))
	}
	r := iotest.TimeoutReader(strings.NewReader(gbRecord("     CDS             1..9")))
	if _, err := ReadGenBank(iotest.HalfReader(r)); err == nil {
		t.Errorf("read error: want an error, got none")
	}
}
