package main

import (
	_ "embed"

	"github.com/virus-evolution/gofasta/cmd"
)

//go:generate go run ./internal/notice/gen

// notice is the licence notices of the third-party software compiled into gofasta (see NOTICE, which
// go generate writes)
//
//go:embed NOTICE
var notice string

func main() {
	cmd.Notice = notice
	cmd.Execute()
}
