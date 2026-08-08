import argparse
import pathlib as pl
import platform
import shlex
import subprocess


parser = argparse.ArgumentParser()
# parser.add_argument("--compiler", type=str, default="g++")
parser.add_argument("--buildtype", type=str, default="release")
parser.add_argument("action")
args = parser.parse_args()

print(f"os platform: {platform.system().lower()}")

# Remove all files from bin folder. This used to clear bin/<ostag> and its parent; the
# per-platform sub-directory no longer exists (see the note on local_install in CMakeLists.txt),
# and looking for it here would have made this quietly stop cleaning anything at all.
bin_dir = pl.Path.cwd() / "bin"
print(f"bin_dir = {bin_dir}")
if bin_dir.is_dir():
    for path in bin_dir.iterdir():
        if path.is_file():
            print(f"removing...'{path}'")
            path.unlink()

command = [
    "pixi",
    "run",
]
if args.action == "rebuild":
    if args.buildtype == "release":
        command += ["build-release"]
    elif args.buildtype == "debug":
        command += ["build-debug"]
    else:
        raise ValueError(f"action '{args.action}' not supported")
else:
    command += ["build"]

print("Run:", shlex.join(command))
subprocess.run(command, check=True)

