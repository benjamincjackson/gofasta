// Package notice generates gofasta's NOTICE file: the licence and notice texts of the third-party Go
// modules compiled into gofasta, and of the Go standard library, copied verbatim from their sources.
package notice

import (
	"bytes"
	"encoding/json"
	"errors"
	"fmt"
	"io"
	"os"
	"os/exec"
	"path/filepath"
	"regexp"
	"sort"
	"strings"
)

// Platforms are the GOOS/GOARCH pairs gofasta is built for: modules compiled in on any of them are
// included
var Platforms = [][2]string{{"linux", "amd64"}, {"linux", "arm64"}, {"darwin", "amd64"}, {"darwin", "arm64"}}

// licenceFiles and noticeFiles are the names (case-insensitive) a module's licence and notice files
// may have, in order of preference
var (
	licenceFiles = []string{"LICENSE", "LICENSE.txt", "LICENSE.md", "LICENCE", "LICENCE.txt", "LICENCE.md", "COPYING"}
	noticeFiles  = []string{"NOTICE", "NOTICE.txt", "NOTICE.md"}
)

// component is one thing compiled into gofasta: a module, or the Go standard library
type component struct {
	name    string // e.g. "github.com/biogo/hts v1.4.5 (https://github.com/biogo/hts)"
	licence string // SPDX full name
	text    string // the licence text, followed by the upstream notice, if any
}

type module struct {
	Path, Version, Dir string
	Main               bool
	Replace            *module
}

type pkg struct {
	Module *module
}

// Generate returns the NOTICE text for the main module in dir
func Generate(dir string) (string, error) {
	modules := make(map[string]module)
	for _, p := range Platforms {
		if err := listModules(dir, p[0], p[1], modules); err != nil {
			return "", err
		}
	}

	components := make([]component, 0, len(modules)+1)
	paths := make([]string, 0, len(modules))
	for path := range modules {
		paths = append(paths, path)
	}
	sort.Strings(paths)
	for _, path := range paths {
		m := modules[path]
		src := m
		if m.Replace != nil {
			src = *m.Replace
		}
		text, licence, err := licenceOf(src.Dir)
		if err != nil {
			return "", fmt.Errorf("%s %s: %w", m.Path, m.Version, err)
		}
		components = append(components, component{name: m.Path + " " + m.Version + " (https://" + m.Path + ")", licence: licence, text: text})
	}

	goroot, err := goEnv(dir, "GOROOT")
	if err != nil {
		return "", err
	}
	text, licence, err := licenceOf(goroot)
	if err != nil {
		return "", fmt.Errorf("the Go standard library: %w", err)
	}
	components = append(components, component{name: "the Go standard library and runtime (https://go.dev)", licence: licence, text: text})

	return format(components), nil
}

// listModules adds the modules other than the main one that are compiled in for goos/goarch
func listModules(dir, goos, goarch string, modules map[string]module) error {
	c := exec.Command("go", "list", "-deps", "-json", ".")
	c.Dir = dir
	c.Env = append(os.Environ(), "GOOS="+goos, "GOARCH="+goarch, "CGO_ENABLED=0")
	out, err := c.Output()
	if err != nil {
		return fmt.Errorf("go list for %s/%s: %w", goos, goarch, err)
	}
	dec := json.NewDecoder(bytes.NewReader(out))
	for {
		var p pkg
		if err := dec.Decode(&p); err == io.EOF {
			return nil
		} else if err != nil {
			return err
		}
		if p.Module != nil && !p.Module.Main {
			modules[p.Module.Path] = *p.Module
		}
	}
}

func goEnv(dir, name string) (string, error) {
	c := exec.Command("go", "env", name)
	c.Dir = dir
	out, err := c.Output()
	return strings.TrimSpace(string(out)), err
}

// licenceOf returns the licence text in dir (followed by its notice file, if any) and the licence's
// SPDX full name
func licenceOf(dir string) (string, string, error) {
	text, err := readFirst(dir, licenceFiles)
	if err != nil {
		return "", "", err
	}
	if text == "" {
		return "", "", errors.New("no licence file found in " + dir)
	}
	licence, err := identify(text)
	if err != nil {
		return "", "", err
	}
	notice, err := readFirst(dir, noticeFiles)
	if err != nil {
		return "", "", err
	}
	if notice != "" {
		text += "\n\n" + notice
	}
	return text, licence, nil
}

// readFirst returns the trimmed content of the first file in dir whose name matches one of names
// (case-insensitive), or "" if there is none
func readFirst(dir string, names []string) (string, error) {
	entries, err := os.ReadDir(dir)
	if err != nil {
		return "", err
	}
	for _, name := range names {
		for _, e := range entries {
			if !e.IsDir() && strings.EqualFold(e.Name(), name) {
				b, err := os.ReadFile(filepath.Join(dir, e.Name()))
				if err != nil {
					return "", err
				}
				return strings.Trim(string(b), "\n"), nil
			}
		}
	}
	return "", nil
}

var space = regexp.MustCompile(`\s+`)

// identify names the licence of text. Only licences known to be acceptable are recognised: anything
// else is an error, to be looked at by a person.
func identify(text string) (string, error) {
	t := space.ReplaceAllString(text, " ")
	has := func(s string) bool { return strings.Contains(t, s) }
	switch {
	case has("Apache License") && has("Version 2.0, January 2004"):
		return "Apache License 2.0", nil
	case has("Redistribution and use in source and binary forms") && has("may be used to endorse or promote products derived"):
		return `BSD 3-Clause "New" or "Revised" License`, nil
	case has("Redistribution and use in source and binary forms"):
		return `BSD 2-Clause "Simplified" License`, nil
	case has("Permission is hereby granted, free of charge, to any person obtaining a copy"):
		return "MIT License", nil
	case has("Permission to use, copy, modify, and/or distribute this software for any purpose"):
		return "ISC License", nil
	}
	return "", errors.New("unrecognised licence")
}

// format writes the NOTICE: a summary, then each distinct licence text with the components that use it
func format(components []component) string {
	type group struct {
		licence, text string
		users         []string
	}
	groups := make([]*group, 0)
	byText := make(map[string]*group)
	licences := make(map[string]bool)
	for _, c := range components {
		licences[c.licence] = true
		g, ok := byText[c.text]
		if !ok {
			g = &group{licence: c.licence, text: c.text}
			byText[c.text] = g
			groups = append(groups, g)
		}
		g.users = append(g.users, c.name)
	}
	sort.SliceStable(groups, func(i, j int) bool { return groups[i].licence < groups[j].licence })
	names := make([]string, 0, len(licences))
	for l := range licences {
		names = append(names, l)
	}
	sort.Strings(names)

	var b strings.Builder
	b.WriteString("gofasta is licensed under the MIT License; see LICENSE.\n\n")
	b.WriteString("gofasta incorporates third-party Go modules, and is built with the Go standard\n")
	b.WriteString("library, which are distributed under the following licences:\n\n")
	for _, n := range names {
		b.WriteString("  - " + n + "\n")
	}
	b.WriteString("\nTheir copyright and licence notices are reproduced in full in the gofasta\n")
	b.WriteString("NOTICE file, as those licences require.\n")
	rule := strings.Repeat("-", 80)
	for _, g := range groups {
		b.WriteString("\n" + rule + "\n" + g.licence + ", used by:\n")
		for _, u := range g.users {
			b.WriteString("  " + u + "\n")
		}
		b.WriteString(rule + "\n\n" + g.text + "\n")
	}
	return b.String()
}
