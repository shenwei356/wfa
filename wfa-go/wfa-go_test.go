package main

import (
	"bytes"
	"flag"
	"fmt"
	"os"
	"os/exec"
	"path/filepath"
	"runtime"
	"testing"
)

// TestCLIHelper runs the CLI in a subprocess to isolate flags and profiling.
func TestCLIHelper(t *testing.T) {
	if os.Getenv("WFA_TEST_CLI_HELPER") != "1" {
		return
	}
	for i, arg := range os.Args {
		if arg == "--" {
			os.Args = append([]string{"wfa-go"}, os.Args[i+1:]...)
			break
		}
	}
	flag.CommandLine = flag.NewFlagSet("wfa-go", flag.ContinueOnError)
	if err := run(); err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}
	os.Exit(0)
}

// TestProfileInputAliases ensures profiling cannot truncate its sequence input,
// including file aliases, and still works with a separate input file.
func TestProfileInputAliases(t *testing.T) {
	aliases := []string{"same", "absolute", "hardlink", "separate"}
	if runtime.GOOS != "windows" { // Windows symlinks require extra privileges.
		aliases = append(aliases, "symlink")
	}
	for _, mode := range []struct{ flag, file string }{{"-p", "cpu.pprof"}, {"-m", "mem.pprof"}} {
		for _, alias := range aliases {
			t.Run(mode.file+"/"+alias, func(t *testing.T) {
				dir := t.TempDir()
				profile := filepath.Join(dir, mode.file)
				data := []byte(">ACGT\n<ACGT\n")
				if err := os.WriteFile(profile, data, 0600); err != nil {
					t.Fatal(err)
				}
				input := mode.file
				if alias == "absolute" {
					input = profile
				} else if alias != "same" {
					input = "query.txt"
					path := filepath.Join(dir, input)
					var err error
					switch alias {
					case "hardlink":
						err = os.Link(profile, path)
					case "symlink":
						err = os.Symlink(profile, path)
					case "separate":
						err = os.WriteFile(path, data, 0600)
					}
					if err != nil {
						t.Fatal(err)
					}
				}
				cmd := exec.Command(os.Args[0], "-test.run=^TestCLIHelper$", "--", mode.flag, "-i", input)
				cmd.Dir = dir
				cmd.Env = append(os.Environ(), "WFA_TEST_CLI_HELPER=1")
				out, err := cmd.CombinedOutput()
				if alias == "separate" {
					if err != nil || !bytes.Contains(out, []byte("4M")) {
						t.Fatalf("profiling a separate input failed: %v\n%s", err, out)
					}
				} else if err == nil || !bytes.Contains(out, []byte("aliases profile output")) {
					t.Errorf("expected alias error, got %v:\n%s", err, out)
				}
				inputPath := input
				if !filepath.IsAbs(inputPath) {
					inputPath = filepath.Join(dir, input)
				}
				got, err := os.ReadFile(inputPath)
				if err != nil || !bytes.Equal(got, data) {
					t.Fatalf("input changed after profiling: %v, got %d bytes, want %d", err, len(got), len(data))
				}
			})
		}
	}
}
