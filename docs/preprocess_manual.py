#!/usr/bin/env python3
"""
Preprocess pestpp_users_manual.md for Sphinx/MyST compatibility.

Converts HTML anchor tags in headings to MyST named targets so that
internal links like [Section](#s5-1) resolve correctly in Sphinx.

  Before: ## <a id='s5-1' />1.1 PEST++ and PEST
  After:  (s5-1)=
          ## 1.1 PEST++ and PEST
"""
import re
import sys
from pathlib import Path


def preprocess(input_path, output_path):
    content = Path(input_path).read_text(encoding="utf-8")

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
