package main

import (
	"testing"

	gennotice "github.com/virus-evolution/gofasta/internal/notice"
)

// NOTICE is what the generator writes for the current dependencies: after adding or updating one,
// run go generate. The Go standard library's licence text depends on the Go version that ran the
// generator, so it is only required to be a BSD 3-clause licence.
func TestNoticeIsUpToDate(t *testing.T) {
	want, err := gennotice.Generate(".")
	if err != nil {
		t.Fatal(err)
	}
	got, goLicence := gennotice.WithoutGoLicence(notice)
	want, _ = gennotice.WithoutGoLicence(want)
	if got != want {
		t.Error("NOTICE is out of date: run go generate")
	}
	if l, err := gennotice.Identify(goLicence); err != nil || l != `BSD 3-Clause "New" or "Revised" License` {
		t.Errorf("NOTICE's Go standard library licence: got %q, %v, want the BSD 3-Clause licence", l, err)
	}
}
