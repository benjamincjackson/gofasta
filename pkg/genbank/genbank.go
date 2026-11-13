// Package genbank reads annotation files in GenBank format.
package genbank

import (
	"bufio"
	"errors"
	"io"
	"strings"
	"unicode"
	"unicode/utf8"
)

var (
	errFeaturesBeforeFeature = errors.New("Error parsing Genbank FEATURES: qualifier before the first feature")
	errMultipleRecords       = errors.New("Error parsing Genbank file: it has more than one record (only one reference is supported)")
)

// maxLineLength is the longest line ReadGenBank accepts
const maxLineLength = 1 << 30

// qualifierIndent is the column (0-based) where the INSDC feature table puts locations and qualifiers;
// feature keys start at column 5. A line indented less than this that does not start with '/' starts a
// feature.
const qualifierIndent = 21

// Genbank is a master struct containing information from a single genbank record
type Genbank struct {
	LOCUS struct {
		Name     string
		Length   int
		Type     string
		Division string
		Date     string
	} // NOT implemented
	DEFINITION string // NOT implemented
	ACCESSION  string // NOT implemented
	VERSION    string // NOT implemented
	KEYWORDS   string // NOT implemented
	SOURCE     struct {
		Source   string
		Organism string
	} // NOT implemented
	REFERENCE struct {
		Authors string
		Title   string
		Journal string
		Pubmed  string
		Remark  string
	} // NOT implemented
	COMMENT  string           // NOT implemented
	FEATURES []GenbankFeature // implemented
	ORIGIN   []byte           // implemented
}

// GenbankFeature contains information about one feature from a genbank record's
// FEATURES section
type GenbankFeature struct {
	Feature  string
	Location Location
	// Info holds the qualifiers: every value of each, in file order. A qualifier without a value (such
	// as /pseudo) has the value "". Quoted values are unquoted ("" is a quote); wrapped lines of a quoted
	// value are joined with a space, except for /translation, where they are joined directly.
	Info map[string][]string
}

// True/False a feature contains the given attribute
func (F *GenbankFeature) HasAttribute(tag string) bool {
	_, ok := F.Info[tag]
	return ok
}

// Value returns the first value of the given qualifier, or "" if the feature does not have it
func (F *GenbankFeature) Value(tag string) string {
	if v := F.Info[tag]; len(v) > 0 {
		return v[0]
	}
	return ""
}

// parseGenbankFEATURES gets the FEATURES info from a genbank file, following the INSDC feature table
// layout: a feature key, then its location (which may wrap onto further lines), then its qualifiers
// (/key, /key=value or /key="value", where a quoted value may wrap onto further lines)
func parseGenbankFEATURES(lines []string) ([]GenbankFeature, error) {
	features := make([]GenbankFeature, 0)

	var cur *GenbankFeature
	var location strings.Builder
	inLocation := false // the current feature's location may still be continued

	var key string
	var value strings.Builder
	haveQualifier := false
	inQuote := false

	// readQuoted adds the text of a quoted value to value, up to its closing quote if there is one
	readQuoted := func(text string) {
		for i := 0; i < len(text); i++ {
			if text[i] != '"' {
				value.WriteByte(text[i])
				continue
			}
			if i+1 < len(text) && text[i+1] == '"' {
				value.WriteByte('"')
				i++
				continue
			}
			inQuote = false
			return
		}
	}
	endQualifier := func() {
		if haveQualifier {
			cur.Info[key] = append(cur.Info[key], value.String())
		}
		haveQualifier, inQuote = false, false
		key = ""
		value.Reset()
	}
	endFeature := func() {
		endQualifier()
		if cur != nil {
			cur.Location = Location{Representation: location.String()}
			features = append(features, *cur)
		}
		cur = nil
		location.Reset()
	}

	for _, line := range lines {
		text := strings.TrimSpace(line)
		// blank (or whitespace-only) lines carry no information
		if text == "" {
			continue
		}

		if inQuote {
			if key != "translation" && value.Len() > 0 {
				value.WriteByte(' ')
			}
			readQuoted(text)
			continue
		}

		indent := len(line) - len(strings.TrimLeft(line, " \t"))
		switch {
		case text[0] != '/' && indent < qualifierIndent:
			endFeature()
			fields := strings.Fields(text)
			cur = &GenbankFeature{Feature: fields[0], Info: make(map[string][]string)}
			location.WriteString(strings.Join(fields[1:], ""))
			inLocation = true

		case text[0] == '/':
			if cur == nil {
				return nil, errFeaturesBeforeFeature
			}
			endQualifier()
			inLocation = false
			haveQualifier = true
			body := text[1:]
			i := strings.IndexByte(body, '=')
			if i < 0 {
				key = body
				continue
			}
			key = body[:i]
			if v := body[i+1:]; strings.HasPrefix(v, "\"") {
				inQuote = true
				readQuoted(v[1:])
			} else {
				value.WriteString(v)
			}

		default:
			// a continuation line: of the location, or of an unquoted value (e.g. /transl_except)
			if cur == nil {
				return nil, errFeaturesBeforeFeature
			}
			if inLocation {
				location.WriteString(strings.ReplaceAll(text, " ", ""))
			} else if haveQualifier {
				value.WriteString(text)
			}
		}
	}
	endFeature()

	return features, nil
}

// parseGenbankORIGIN gets the sequence from the lines of a genbank file's ORIGIN: its letters
func parseGenbankORIGIN(lines []string) []byte {
	seq := make([]byte, 0)
	for _, line := range lines {
		for _, c := range line {
			if unicode.IsLetter(c) {
				seq = utf8.AppendRune(seq, c)
			}
		}
	}
	return seq
}

// ReadGenBank reads a genbank annotation file and returns a struct that contains
// parsed versions of the fields it contains. Not all fields are currently parsed. The file must hold
// one record: anything after the record's closing "//" line is an error.
func ReadGenBank(r io.Reader) (Genbank, error) {
	gb := Genbank{}

	s := bufio.NewScanner(r)
	s.Buffer(make([]byte, 0, 64*1024), maxLineLength)

	var header string
	var lines []string
	ended := false // the record's closing "//" has been read

	parseField := func() error {
		switch header {
		case "FEATURES":
			var err error
			gb.FEATURES, err = parseGenbankFEATURES(lines)
			return err
		case "ORIGIN":
			gb.ORIGIN = parseGenbankORIGIN(lines)
		}
		return nil
	}

	for s.Scan() {
		line := s.Text()
		trimmed := strings.TrimSpace(line)
		if trimmed == "" {
			continue
		}
		if ended {
			return Genbank{}, errMultipleRecords
		}
		if trimmed == "//" {
			ended = true
			continue
		}
		// a line that starts with a capital letter starts a field
		if 'A' <= line[0] && line[0] <= 'Z' {
			if err := parseField(); err != nil {
				return Genbank{}, err
			}
			header = strings.Fields(line)[0]
			lines = make([]string, 0)
			continue
		}

		lines = append(lines, line)
	}
	if err := s.Err(); err != nil {
		return Genbank{}, errors.New("Error reading Genbank file: " + err.Error())
	}
	if err := parseField(); err != nil {
		return Genbank{}, err
	}

	return gb, nil
}
