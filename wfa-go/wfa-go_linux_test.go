package main

import (
	"bytes"
	"os"
	"strings"
	"testing"
)

// TestCLIOutputFailure checks both buffered flushes and writes during formatting.
func TestCLIOutputFailure(t *testing.T) {
	for _, n := range []int{2, 8192} {
		f, err := os.OpenFile("/dev/full", os.O_WRONLY, 0)
		if err != nil {
			t.Fatal(err)
		}
		cmd := cliCommand(t, strings.Repeat("A", n), strings.Repeat("A", n))
		cmd.Stdout = f
		var stderr bytes.Buffer
		cmd.Stderr = &stderr
		err = cmd.Run()
		f.Close()
		if err == nil || stderr.Len() == 0 {
			t.Fatalf("expected output error, got %v: %s", err, stderr.String())
		}
	}
}
