#!/usr/bin/env python3
"""Guard the API surface: the shipped iteration loops must call nothing but public methods.

The shared-library API is built on the idea that the tools' own loops are the first client
of the same public building blocks an external caller gets. If the built-in loop can be
written using only public methods, then any caller can write any loop. The moment a loop
reaches for a protected or private member, that guarantee is gone - silently, because it
still compiles.

This checks the property directly: parse every member call made inside the loop bodies, and
fail if any of them resolves to a non-public member of the owning class.

Deliberately conservative: it only complains about names it can positively place in the
class's protected/private section, so unknown or free functions are ignored rather than
guessed at. False negatives are acceptable here; false positives would make it noise.
"""
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "src" / "libs" / "pestpp_common"

# (file, function signature fragment, header, class) for each shipped loop
LOOPS = [
    ("EnsembleSmoother.cpp", "IterEnsembleSmoother::iterate_2_solution",
     "EnsembleMethodUtils.h", "EnsembleMethod"),
    ("DataAssimilator.cpp", "DataAssimilator::da_update",
     "EnsembleMethodUtils.h", "EnsembleMethod"),
    ("MOEA.cpp", "MOEA::iterate_to_solution", "MOEA.h", "MOEA"),
    ("SQP.cpp", "SeqQuadProgram::iterate_2_solution", "SQP.h", "SeqQuadProgram"),
]

# language keywords and common locals that look like calls but are not member calls
IGNORE = {
    "if", "for", "while", "switch", "return", "catch", "sizeof", "throw",
    "abs", "max", "min", "break", "else", "do", "new", "delete", "static_cast",
    "dynamic_cast", "const_cast", "reinterpret_cast", "string", "int", "double",
    "bool", "vector", "map", "set", "pair", "stringstream", "ofstream", "ifstream",
    "size_t", "find", "make_pair", "to_string", "move", "sort", "count",
}


def function_body(path, signature):
    """Return the text of the function whose definition line contains `signature`."""
    lines = path.read_text(errors="ignore").split("\n")
    start = next((i for i, l in enumerate(lines) if signature in l and "::" in l), None)
    if start is None:
        raise SystemExit(f"FAIL: could not find {signature} in {path.name}")
    depth, started, out = 0, False, []
    for line in lines[start:]:
        out.append(line)
        depth += line.count("{") - line.count("}")
        if "{" in line:
            started = True
        if started and depth <= 0:
            break
    return "\n".join(out)


def sections(path, cls):
    """Return (public_names, nonpublic_names) for `cls` in header `path`."""
    text = path.read_text(errors="ignore")
    i = text.index(f"class {cls}")
    body = text[i:]
    # stop at the end of this class: the next top-level 'class X' declaration
    nxt = re.search(r"\nclass [A-Za-z_]", body[10:])
    if nxt:
        body = body[: nxt.start() + 10]

    public, nonpublic, current = set(), set(), None
    for line in body.split("\n"):
        stripped = line.strip()
        if stripped.startswith("public:"):
            current = public
        elif stripped.startswith("protected:") or stripped.startswith("private:"):
            current = nonpublic
        elif current is not None:
            for m in re.finditer(r"(?<![\w:.>])([a-z_][\w]*)\s*\(", line):
                current.add(m.group(1))
    return public, nonpublic


def main():
    violations = []
    for cpp, signature, header, cls in LOOPS:
        body = function_body(SRC / cpp, signature)
        public, nonpublic = sections(SRC / header, cls)
        called = {
            m.group(1)
            for m in re.finditer(r"(?<![\w.>:])([a-z_][\w]*)\s*\(", body)
        }
        for name in sorted(called - IGNORE):
            # only flag names we can positively place as non-public
            if name in nonpublic and name not in public:
                violations.append(f"  {cpp}:{signature} calls non-public '{cls}::{name}()'")

    if violations:
        print("FAIL: shipped loops reach for non-public members.")
        print("The API guarantee is that a caller can write any loop the tools can write,")
        print("so these must either become public or be replaced by a public equivalent:")
        print("\n".join(violations))
        return 1

    print(f"OK: all {len(LOOPS)} shipped loops use public methods only")
    return 0


if __name__ == "__main__":
    sys.exit(main())
