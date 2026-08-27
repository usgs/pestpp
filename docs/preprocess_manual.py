#!/usr/bin/env python3
"""
Preprocess pestpp_users_manual.md for Sphinx/MyST compatibility.

Converts HTML anchor tags in headings to MyST named targets so that
internal links like [Section](#s5-1) resolve correctly in Sphinx.

  Before: ## <a id='s5-1' />1.1 PEST++ and PEST
  After:  (s5-1)=
          ## 1.1 PEST++ and PEST

Also rewrites the version on the title page from config_os.h, so the published manual cannot
disagree with the code it documents. The manual had been sitting at 5.2.28 while the code moved
on, and nothing pointed that out.
"""
import re
import sys
from pathlib import Path

# the one place the version is defined - the same line cmake, cpack and release.yml parse
CONFIG_OS_H = Path(__file__).resolve().parent.parent / "src" / "libs" / "common" / "config_os.h"


def read_version():
    """The version from config_os.h, or None if it cannot be read.

    None rather than an exception on purpose: a docs build that cannot find the header should
    still produce a manual, just with whatever version the file already carries.
    """
    try:
        text = CONFIG_OS_H.read_text(encoding="utf-8")
    except OSError:
        return None
    m = re.search(r'#define\s+PESTPP_VERSION\s+"([^"]+)"', text)
    return m.group(1) if m else None


def preprocess(input_path, output_path):
    content = Path(input_path).read_text(encoding="utf-8")

    # Title page version. It sits on a line of its own, deliberately - it used to BE the first
    # heading, which meant the table of contents mirrored it and the version had to be patched
    # in two places and two syntaxes. Now the heading is stable title text and this is the only
    # line carrying a version, so only the version characters ever change.
    version = read_version()
    if version:
        content, n = re.subn(
            r'^Version[ \t]+\S+[ \t]*$',
            f"Version {version}",
            content,
            count=1,
            flags=re.MULTILINE,
        )
        if n == 0:
            print(f"warning: no 'Version ...' line to update in {input_path}",
                  file=sys.stderr)

    # Convert:  # <a id='s5' />1. Introduction
    # To:       (s5)=
    #           # 1. Introduction
    def replace_heading_anchor(m):
        hashes = m.group(1)
        anchor_id = m.group(2)
        heading_text = m.group(3).strip()
        return f"({anchor_id})=\n{hashes} {heading_text}"

    content = re.sub(
        r'^(#{1,6})[ \t]*<a\s+id=[\'"]([^\'"]+)[\'"][ \t]*/?>[ \t]*(.*)$',
        replace_heading_anchor,
        content,
        flags=re.MULTILINE,
    )

    Path(output_path).write_text(content, encoding="utf-8")


if __name__ == "__main__":
    preprocess(sys.argv[1], sys.argv[2])
