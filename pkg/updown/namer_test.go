package updown

import "github.com/virus-evolution/gofasta/pkg/pipeline"

// records give their name to the warning if processing them panics
var (
	_ pipeline.Namer = csvRow{}
	_ pipeline.Namer = targetIn{}
)
