"""Progress rendering for long-running pest++ work.

`solve()` is silent for as long as the model runs take, which on a real problem is all of the
time. This turns the run manager's counters and the per-iteration phi into something you can
watch, in a terminal or a notebook, without either one needing to know which it is.

Two things shape the implementation, and both are worth knowing before changing it:

**Nothing is imported that pest++ does not already require.** The rest of the python layer
depends on numpy and pandas and nothing else, and a progress bar is not a good reason to add
tqdm or ipywidgets to a modelling stack. A terminal bar is a carriage return, and a notebook
bar is IPython's display protocol, which is present by definition when you are in a notebook.

**Progress goes to STDERR, never stdout.** A session opened with ``quiet=True`` has the
library dup2 its own output over file descriptor 1, and python's stdout is that same
descriptor - so anything printed while a run is in flight ends up in ``pestpp.stdout.log``
rather than on screen. stderr is left alone. In a notebook the question does not arise: the
display protocol is a kernel message, not a write to a descriptor.
"""
from __future__ import annotations

import os
import sys
import time

__all__ = ["Progress", "TextProgress", "NotebookProgress", "auto", "in_notebook"]


def in_notebook() -> bool:
    """Whether we are running under a notebook frontend that can render a display update.

    A terminal IPython session is deliberately NOT a notebook here: it has ``get_ipython()``
    but renders HTML as a repr string, which is worse than the text bar it would otherwise get.
    """
    try:
        from IPython import get_ipython
    except Exception:
        return False
    try:
        shell = get_ipython()
    except Exception:
        return False
    return shell is not None and shell.__class__.__name__ == "ZMQInteractiveShell"


def _fmt(value) -> str:
    """Numbers short enough to sit in a status line, everything else as-is."""
    if isinstance(value, float):
        if value != value:                 # NaN
            return "-"
        if value == 0 or 1e-3 <= abs(value) < 1e5:
            return "{0:.4g}".format(value)
        return "{0:.3e}".format(value)
    return str(value)


class Progress:
    """The no-op renderer, and the interface the others implement.

    Every hook is safe to call at any time and in any order, so calling code never has to
    branch on whether progress is switched on.
    """

    def start(self, label: str = "", total: int | None = None) -> None:
        """Begin a phase. Calling it again starts a new one."""

    def update(self, done: int | None = None, total: int | None = None, **fields) -> None:
        """Advance the current phase. ``fields`` are shown as ``name=value``."""

    def note(self, text: str) -> None:
        """A line that stays, rather than being overwritten by the next update."""

    def close(self) -> None:
        """Finish: leave the final state visible and stop drawing."""

    # a Progress is usable as a context manager so close() cannot be forgotten
    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()
        return False


class _Rendered(Progress):
    """Shared state and formatting for the renderers that actually draw something."""

    #: seconds between redraws. Fast enough to look live, slow enough that a serial run
    #: manager finishing thousands of quick runs does not spend its time formatting strings.
    min_interval = 0.1

    def __init__(self):
        self._label = ""
        self._total = None
        self._done = 0
        self._fields = {}
        self._started = time.time()
        self._last_draw = 0.0
        self._closed = False

    def start(self, label="", total=None):
        self._label = label
        self._total = total
        self._done = 0
        self._fields = {}
        self._started = time.time()
        self._last_draw = 0.0
        self._closed = False
        self._draw(force=True)

    def update(self, done=None, total=None, **fields):
        if done is not None:
            self._done = done
        if total is not None:
            self._total = total
        self._fields.update(fields)
        self._draw()

    def close(self):
        if self._closed:
            return
        self._draw(force=True, final=True)
        self._closed = True

    # -- formatting ---------------------------------------------------------------------

    def _fraction(self):
        if not self._total:
            return None
        return max(0.0, min(1.0, float(self._done) / float(self._total)))

    def _status(self) -> str:
        bits = []
        if self._label:
            bits.append(self._label)
        if self._total:
            bits.append("{0}/{1}".format(self._done, self._total))
        elif self._done:
            bits.append(str(self._done))
        for k, v in self._fields.items():
            bits.append("{0}={1}".format(k, _fmt(v)))
        bits.append("{0:.1f}s".format(time.time() - self._started))
        return "  ".join(bits)

    def _draw(self, force=False, final=False):
        now = time.time()
        if (not force) and (now - self._last_draw < self.min_interval):
            return
        self._last_draw = now
        self._paint(final)

    def _paint(self, final):
        raise NotImplementedError


class TextProgress(_Rendered):
    """A single rewriting line on a terminal; plain lines when the output is not a terminal.

    The distinction matters more than it looks: carriage returns in a redirected log turn it
    into one unreadable line, so a non-tty gets whole lines, and only when something changed
    enough to be worth a line.
    """

    bar_width = 24

    def __init__(self, stream=None):
        super().__init__()
        # stderr, not stdout - see the module docstring
        self._stream = stream if stream is not None else sys.stderr
        try:
            self._tty = self._stream.isatty()
        except Exception:
            self._tty = False
        self._last_logged = -1.0
        self._width = 0

    def note(self, text):
        if self._tty and self._width:
            self._stream.write("\r" + " " * self._width + "\r")
            self._width = 0
        self._stream.write(text.rstrip() + "\n")
        self._stream.flush()

    def _bar(self) -> str:
        frac = self._fraction()
        if frac is None:
            return ""
        filled = int(round(frac * self.bar_width))
        return "[{0}{1}] {2:3.0f}%  ".format(
            "#" * filled, "." * (self.bar_width - filled), frac * 100.0)

    def _paint(self, final):
        line = self._bar() + self._status()
        if self._tty:
            pad = max(0, self._width - len(line))
            self._stream.write("\r" + line + " " * pad)
            self._width = len(line)
            if final:
                self._stream.write("\n")
                self._width = 0
            self._stream.flush()
            return
        # not a terminal: whole lines, and only on real movement, so a log stays readable
        frac = self._fraction()
        step = 0.0 if frac is None else frac
        if final or (self._last_logged < 0) or (step - self._last_logged >= 0.25):
            self._last_logged = step
            self._stream.write(line + "\n")
            self._stream.flush()


class NotebookProgress(_Rendered):
    """One output area, updated in place, using IPython's display protocol.

    HTML rather than ipywidgets: it needs no extra package, survives notebook reload as static
    output, and renders in nbconvert. The cost is that it cannot be interacted with, which a
    progress bar does not need.
    """

    def __init__(self):
        super().__init__()
        from IPython.display import display, HTML
        self._display, self._HTML = display, HTML
        self._handle_id = "pestpp-progress-{0}".format(id(self))
        self._shown = False
        self._notes = []

    def note(self, text):
        self._notes.append(text)
        self._draw(force=True)

    def _paint(self, final):
        frac = self._fraction()
        pct = 0.0 if frac is None else frac * 100.0
        # an indeterminate phase gets a full-width neutral bar rather than an empty one
        width = 100.0 if frac is None else pct
        colour = "#3d7ea6" if not final else "#2e7d4f"
        bar = (
            '<div style="background:#e9e9e9;border-radius:3px;height:10px;width:100%;'
            'overflow:hidden;margin:2px 0">'
            '<div style="background:{0};height:10px;width:{1:.1f}%"></div></div>'
        ).format(colour, width)
        notes = "".join(
            '<div style="color:#666;font-size:90%">{0}</div>'.format(_escape(n))
            for n in self._notes[-4:])
        html = (
            '<div style="font-family:monospace;font-size:90%;max-width:46em">'
            '{0}<div>{1}{2}</div>{3}</div>'
        ).format(notes, "" if frac is None else "{0:3.0f}%  ".format(pct),
                 _escape(self._status()), bar)
        payload = self._HTML(html)
        if not self._shown:
            self._display(payload, display_id=self._handle_id)
            self._shown = True
            return
        from IPython.display import update_display
        update_display(payload, display_id=self._handle_id)


def _escape(text: str) -> str:
    return (str(text).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"))


def auto(enabled=True, stream=None) -> Progress:
    """The right renderer for wherever this is running.

    ``enabled`` may be a bool, or a :class:`Progress` to use as-is, which is what lets an
    argument spelled ``progress=True`` and one spelled ``progress=my_renderer`` share a path.
    Disabled - and any environment that cannot render - gives the no-op, so calling code never
    branches.
    """
    if isinstance(enabled, Progress):
        return enabled
    if not enabled:
        return Progress()
    if os.environ.get("PESTPP_NO_PROGRESS"):
        return Progress()
    if in_notebook():
        try:
            return NotebookProgress()
        except Exception:
            pass          # IPython present but display unavailable; fall back to text
    return TextProgress(stream=stream)
