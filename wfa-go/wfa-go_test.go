package main

import (
	"context"
	"flag"
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"testing"
	"time"
)

// TestCLIHelper runs the real entry point in an isolated process.
func TestCLIHelper(t *testing.T) {
	if os.Getenv("WFA_CLI_TEST_HELPER") != "1" {
		return
	}
	for i, arg := range os.Args {
		if arg == "--" {
			os.Args = append([]string{"wfa-go"}, os.Args[i+1:]...)
			break
		}
	}
	flag.CommandLine = flag.NewFlagSet("wfa-go", flag.ExitOnError)
	main()
	os.Exit(0)
}

// cliCommand keeps malformed input and output failures from stopping the suite.
func cliCommand(t *testing.T, args ...string) *exec.Cmd {
	t.Helper()
	ctx, cancel := context.WithTimeout(context.Background(), 10*time.Second)
	t.Cleanup(cancel)
	cmd := exec.CommandContext(ctx, os.Args[0], append([]string{"-test.run=^TestCLIHelper$", "--"}, args...)...)
	cmd.Env = append(os.Environ(), "WFA_CLI_TEST_HELPER=1")
	return cmd
}

// TestCLIInputPairs requires complete, correctly prefixed sequence pairs.
func TestCLIInputPairs(t *testing.T) {
	for _, tc := range []struct {
		name, input string
		valid       bool
	}{
		{"valid", ">AC\n<AC\n", true},
		{"crlf", ">AC\r\n<AC\r\n", true},
		{"no-final-newline", ">AC\n<AC", true},
		{"missing-target", ">AC\n<AC\n>GT\n", false},
		{"blank-query", "\n<AC\n", false},
		{"empty-query", ">\n<AC\n", false},
		{"empty-target", ">AC\n<\n", false},
		{"query-prefix", "AC\n<AC\n", false},
		{"target-prefix", ">AC\n>AC\n", false},
	} {
		t.Run(tc.name, func(t *testing.T) {
			p := filepath.Join(t.TempDir(), "pairs.txt")
			if err := os.WriteFile(p, []byte(tc.input), 0600); err != nil {
				t.Fatal(err)
			}
			output, err := cliCommand(t, "-i", p).CombinedOutput()
			if tc.valid {
				if err != nil || !strings.Contains(string(output), "cigar   2M") {
					t.Fatalf("%v: %s", err, output)
				}
			} else {
				if err == nil || strings.Contains(string(output), "panic:") {
					t.Fatalf("expected input error, got %v: %s", err, output)
				}
			}
		})
	}
}

// TestCLITrimWithoutMatches prints a finite identity and an empty region.
func TestCLITrimWithoutMatches(t *testing.T) {
	output, err := cliCommand(t, "-t", "A", "C").CombinedOutput()
	if err != nil {
		t.Fatalf("%v: %s", err, output)
	}
	if strings.Contains(string(output), "NaN") {
		t.Fatalf("invalid identity: %s", output)
	}
}

// TestCLILongInput accepts valid sequence lines beyond Scanner's default limit.
func TestCLILongInput(t *testing.T) {
	seq := strings.Repeat("A", 70000)
	p := filepath.Join(t.TempDir(), "long.txt")
	if err := os.WriteFile(p, []byte(">"+seq+"\n<"+seq+"\n"), 0600); err != nil {
		t.Fatal(err)
	}
	output, err := cliCommand(t, "-N", "-i", p).CombinedOutput()
	if err != nil {
		t.Fatalf("%v: %s", err, output)
	}
}

// TestCLISemiGlobalReadme preserves the original -g example and its long match.
func TestCLISemiGlobalReadme(t *testing.T) {
	output, err := cliCommand(t, "-g", "Bioinformatics helps Biology", "We learn bioinformatics to help biologists").CombinedOutput()
	if err != nil {
		t.Fatalf("%v: %s", err, output)
	}
	for _, want := range []string{"cigar   9I1X14M3I4M1D1M1X5M1X3I", "align-score : 32", "match-region: q[2, 27]/28 vs t[11, 38]/42", "align-length: 29, matches: 24 (82.76%)"} {
		if !strings.Contains(string(output), want) {
			t.Errorf("missing %q in output:\n%s", want, output)
		}
	}
}
