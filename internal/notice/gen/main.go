// Command gen writes gofasta's NOTICE file. Run it with "go generate" from the repository root.
package main

import (
	"fmt"
	"os"

	"github.com/virus-evolution/gofasta/internal/notice"
)

func main() {
	text, err := notice.Generate(".")
	if err != nil {
		fmt.Fprintln(os.Stderr, "generating NOTICE:", err)
		os.Exit(1)
	}
	if err := os.WriteFile("NOTICE", []byte(text), 0o644); err != nil {
		fmt.Fprintln(os.Stderr, "writing NOTICE:", err)
		os.Exit(1)
	}
}
