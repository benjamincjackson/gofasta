package sam

import "github.com/virus-evolution/gofasta/pkg/pipeline"

// queries give their name to the warning if processing them panics
var _ pipeline.Namer = samRecords{}
