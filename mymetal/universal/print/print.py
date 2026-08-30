import shutil
import sys
from pathlib import Path


def pr(text=""):
    """Print a plain message (alias of ``print``)."""
    print(text)

def er(text=""):
    """Print an error-line message without exiting (see :func:`fail`)."""
    print(text)

def warn(msg):
    """Print a warning message prefixed with a warning emoji."""
    print(f"⚠️  {msg}")

def fail(msg):
    """Print an error message prefixed with a cross emoji, then raise ``SystemExit(1)``.

    Use this for unrecoverable errors where continuing would produce
    misleading output. The message is printed to stdout before the exit.
    """
    print(f"❌ ERROR: {msg}")
    raise SystemExit(1)


def confirm_prepare_outdir(path_out, force: bool = False) -> None:
    """Ask before deleting an existing output directory, then leave it absent.

    The python and bash workflows share one policy for an existing ``y_*``
    output, mirroring the bash ``rm -I`` prompt:

    - missing         -> return; the caller creates it.
    - exists + force   -> delete without asking (automation / ``-force``).
    - exists + a tty   -> ask ``[y/N]`` (blank = No); delete on yes, abort on no.
    - exists + no tty  -> abort; there is no terminal to confirm (piped/sbatch).

    Args:
        path_out (str | Path): Directory the workflow is about to (re)create.
        force (bool): Skip the prompt and delete unconditionally.
    """
    path_out = Path(path_out)
    if not path_out.exists():
        return
    # forced: matches the old python default of deleting silently, but opt-in
    if force:
        warn("deleting existing output (forced): %s" % path_out)
        shutil.rmtree(path_out)
        return
    # no terminal to ask on -> never delete data behind the user's back
    if not sys.stdin.isatty():
        fail("output already exists and there is no terminal to confirm "
             "deletion: %s\n   re-run in a terminal, or pass -force to "
             "overwrite" % path_out)
    answer = input("⚠️  %s already exists. Delete and recreate? [y/N] "
                   % path_out)
    if answer.strip().lower() not in ("y", "yes"):
        fail("aborted by user; existing output kept: %s" % path_out)
    shutil.rmtree(path_out)