package main

import (
	"testing"

	gennotice "github.com/virus-evolution/gofasta/internal/notice"
)

// NOTICE is what the generator writes for the current dependencies: after adding or updating one,
// run go generate
func TestNoticeIsUpToDate(t *testing.T) {
	want, err := gennotice.Generate(".")
	if err != nil {
		t.Fatal(err)
	}
	if notice != want {
		t.Error("NOTICE is out of date: run go generate")
	}
}
