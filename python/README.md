# the pest++ python api

run pestpp-ies, -da, -glm, -mou, -sqp and -opt from python, in process, instead of writing a
control file, launching an exe and reading what it left behind.

this is a 0.x preview - see [stability](#stability) before you build anything on it.

## two layers

| module | what it is | needs |
|---|---|---|
| `pestpp_lib.py` | thin binding - one method per c function, same names | numpy |
| `pestpp.py` | the friendly layer - pandas/pyemu objects, tool classes | numpy, pandas, pyemu |

use `pestpp.py` unless you have a reason not to. `pestpp_lib.py` doesnt import pyemu, so if you
are building something directly on the c api you dont have to drag pyemu along with it.

## install

there is no pip package yet - see [stability](#stability). two ways to set it up:

### from a release bundle (no build)

the release tarball has everything - the library, the header and these python files:

```
bin/pestpp-api.dylib   (or .so / .dll)
bin/pestpp-ies, -da, -glm, ...
include/pestpp-api.h
python/pestpp.py, pestpp_lib.py, README.md
```

so just extract it and put the python directory on your path:

```bash
tar xzf pestpp-5.2.x-mac.tar.gz          # or unzip the windows bundle
```

```python
import sys
sys.path.insert(0, "/path/to/pestpp/python")
from pestpp import Ies
```

nothing else to set. `find_library()` looks in `../bin` relative to itself, which is where the
library is in this layout.

dont mix and match a python directory from one version with a library from another - they move
together right now, and if they dont match you get a ctypes error about a missing symbol, which
sends you looking in the wrong place.

### from a build tree

if you build pest++ yourself there is nothing else to do - `find_library()` looks under
`<repo>/build` and takes the newest one it finds.

### pointing it at a specific library

```bash
export PESTPP_LIB=/somewhere/pestpp-api.dylib
```

worth setting any time more than one library could match. the search picks the newest by
modification time, so an installed copy next to the exes and a fresh build tree will compete.

## checking what you have

```bash
python python/pestpp_lib.py
```

```
pestpp API environment:
  library_path     /Users/you/pestpp/build/src/libs/pestpp_capi/pestpp-api.dylib
  library_mtime    2026-08-12 09:40:09
  library_size     9641152
  api_version      0.4.0
  pestpp_version   5.2.30
  python           3.13.12
  platform         macOS-26.4.1-arm64-arm-64bit-Mach-O
```

please paste this into any bug report. the `library_path` line is the reason it exists - since
the search picks the newest match, two people with the same checkout can be running different
libraries and nothing in the behavior tells you which. it also works when the library wont load
at all, which is the case you most need it for - it puts the problem in an `error` field
instead of raising.

also available as `pestpp_lib.api_info()` (a dict) and `pestpp_lib.format_api_info()` (the block
above). both can be imported from `pestpp` too.

## hello world

```python
from pestpp import Ies

with Ies.from_pst("pest.pst", workdir="my_case") as ies:
    ies.initialize()
    while not ies.should_terminate:
        ies.solve()
        print(ies.iteration, ies.phi.mean())
    ies.finalize()
```

`examples/pestpp_api_demo.ipynb` is the longer version. `examples/pestpp_progress_demo.ipynb`
shows how to watch a run, and `examples/pestpp_staged_ies_demo.ipynb` shows driving one in
stages.

## reading a finished run

dont use this api for that - use `pyemu.Results(master_dir, case)`. helpers that read the output
files were written here and then taken back out, because there is no reason to have two ways to
read the same csv. this api is for the stuff that only exists while a tool is running, or that
changes what it does next.

## what each tool can do

the tools arent the same shape and the api says so instead of pretending. calls that dont apply
give you a `PestppError` with a reason rather than something empty - an empty dataframe looks
like "nothing yet" when what you really did was call the wrong tool.

| | ies | da | glm | mou | sqp | opt |
|---|:-:|:-:|:-:|:-:|:-:|:-:|
| `par_df()` / `obs_df()` | y | y | | y | y | |
| `phi` | y | y | | | | |
| `par_vector()` | | | y | | | y |
| `jco()` / jacobian steps | | | y | | | |
| `regularization()` | | | y | | | |
| `obs_vector()` | | | y | | | y |
| stacks / risk / chance | | | | y | y | y |
| constraints / feasibility | | | | y | y | y |
| objective values | | | | | | y |
| decision variable names | | | | | | y |
| assimilation cycles | | y | | | | |

the pattern behind it:

- **glm and opt carry one parameter vector, not an ensemble**, so the ensemble calls and the
  phi-over-realizations calls error out. use `par_vector()` for glm and `dec_var_vector()` for
  opt.
- **mou and sqp have populations but no phi** - mou has objectives and a pareto front, sqp has
  an objective function. phi over realizations is an ies/da thing.
- **mou, sqp and opt have the chance/stack machinery**; ies, da and glm dont even have the
  methods on the class.
- **mou, sqp and opt have constraints**, so the constraint names, senses and violation total
  are shared across those three rather than being opt-only.
- **opt is the only one with a single objective value that moves per iteration**, which is why
  it gets its own calls instead of going through phi.

this table came from actually calling everything on every tool, not from reading the source.

## stability

**version 0.x. expect breaking changes. pin to a commit.**

this is out for feedback from people who know what they are doing, and that feedback is going to
change its shape. until 1.0:

- the abi minor number goes up for additions, and might also change for breaking changes
- `pestpp_get_api_version()` tells you what the library you loaded provides - compare it against
  the `PESTPP_API_VERSION_*` values in `include/pestpp-api.h` that you built against
- the python files move with the library, so use the ones out of the same bundle
  (or the same commit, if you are building it yourself)

if you need something that wont move on you, run the exes.

if something surprises you, or you had to write a workaround, please say so - and include the
`python pestpp_lib.py` block.
