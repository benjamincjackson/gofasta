package notice

import (
	"strings"
	"testing"
)

func TestIdentify(t *testing.T) {
	for text, want := range map[string]string{
		"Apache License\n Version 2.0, January 2004\n...": "Apache License 2.0",
		"Redistribution and use in source and binary forms, with or without\nmodification ... Neither the name of X nor the names of its contributors may be used to endorse or promote products derived from this software": `BSD 3-Clause "New" or "Revised" License`,
		"Redistribution and use in source and binary forms, with or without modification, are permitted":                                                                                                                     `BSD 2-Clause "Simplified" License`,
		"Permission is hereby granted, free of charge, to any person obtaining a copy of this software":                                                                                                                      "MIT License",
	} {
		if got, err := identify(text); err != nil || got != want {
			t.Errorf("%.30q: got %q, %v, want %q", text, got, err, want)
		}
	}
	// anything else needs a person to look at it
	for _, text := range []string{"GNU GENERAL PUBLIC LICENSE Version 3", "All rights reserved.", ""} {
		if got, err := identify(text); err == nil {
			t.Errorf("%q: got %q, want an error", text, got)
		}
	}
}

// Components with the same licence text share one block; the summary lists each licence once
func TestFormatGroupsIdenticalTexts(t *testing.T) {
	out := format([]component{
		{name: "a v1", licence: "MIT License", text: "text A"},
		{name: "b v1", licence: "Apache License 2.0", text: "text B"},
		{name: "c v1", licence: "MIT License", text: "text A"},
	})
	if strings.Count(out, "text A") != 1 || !strings.Contains(out, "MIT License, used by:\n  a v1\n  c v1\n") {
		t.Errorf("identical texts not grouped:\n%s", out)
	}
	if strings.Count(out, "  - MIT License\n") != 1 || strings.Index(out, "  - Apache") > strings.Index(out, "  - MIT") {
		t.Errorf("summary should list each licence once, sorted:\n%s", out)
	}
	if i, j := strings.Index(out, "Apache License 2.0, used by"), strings.Index(out, "MIT License, used by"); i > j {
		t.Errorf("blocks should be sorted by licence:\n%s", out)
	}
}
