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

os_platform = platform.system().lower()
print(f"os platform: {os_platform}")
if "linux" in os_platform:
    _ostag = "linux"
elif "darwin" in os_platform:
    _ostag = "mac"
elif "windows" in os_platform:
    _ostag = "win"
else:
    _ostag = "unknown"

# Remove all files from bin folder
bin_dir = pl.Path.cwd() / f"bin/{_ostag}"
print(f"bin_dir = {bin_dir}")
if bin_dir.is_dir():
    for on_dir in [bin_dir, bin_dir.parent]:
        for path in on_dir.iterdir():
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

