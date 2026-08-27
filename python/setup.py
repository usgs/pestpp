"""Version lookup for the pestpp distribution.

Everything else about the package is declared in pyproject.toml. This file exists only to read
the version from src/libs/common/config_os.h - the one place it is defined, and the same line
cmake, cpack, release.yml and the docs build all parse. Restating it here would make a second
source of truth, and the way that fails is publishing a wheel whose version disagrees with the
library inside it.
"""
import re
from pathlib import Path

from setuptools import setup

CONFIG_OS_H = Path(__file__).resolve().parent.parent / "src" / "libs" / "common" / "config_os.h"


def read_version():
    text = CONFIG_OS_H.read_text(encoding="utf-8")
    m = re.search(r'#define\s+PESTPP_VERSION\s+"([^"]+)"', text)
    if not m:
        raise RuntimeError("no PESTPP_VERSION in {0}".format(CONFIG_OS_H))
    # config_os.h carries release candidates as "6.0.0rc1"; PEP 440 wants "6.0.0rc1" too, so
    # the common cases need no translation. A bare trailing letter ("6.0.0a") is NOT valid
    # PEP 440 on its own - it means alpha, and needs a number after it.
    version = m.group(1)
    version = re.sub(r"(\d)([ab])$", r"\1\g<2>0", version)
    return version


setup(version=read_version())
