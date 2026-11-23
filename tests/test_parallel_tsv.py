import os
import subprocess
import sys
import threading
from pathlib import Path

from app.services import parallel_tsv

REPO = Path(__file__).resolve().parents[1]


def test_fork_map_keeps_order_and_runs_in_workers_when_single_threaded():
    code = (
        "import os, sys; sys.path.insert(0, 'backend');"
        "from app.services.parallel_tsv import fork_map;"
        "parent = os.getpid();"
        "out = fork_map(lambda x: (x * x, os.getpid() != parent), range(6), 3);"
        "print([v for v, _ in out], any(child for _, child in out))"
    )
    proc = subprocess.run([sys.executable, "-c", code], cwd=REPO, capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr
    assert proc.stdout.strip() == "[0, 1, 4, 9, 16, 25] True"


def test_fork_map_falls_back_to_serial_with_other_threads():
    stop = threading.Event()
    thread = threading.Thread(target=stop.wait)
    thread.start()
    try:
        assert parallel_tsv.fork_map(lambda _x: os.getpid(), range(3), 3) == [os.getpid()] * 3
    finally:
        stop.set()
        thread.join()


def test_fork_map_raises_instead_of_hanging_when_a_worker_is_killed():
    # multiprocessing.Pool waits forever when a worker is SIGKILLed (e.g. by
    # the OOM killer on a WGS batch); fork_map must fail fast instead.
    code = (
        "import os, signal, sys, time; sys.path.insert(0, 'backend');"
        "from app.services.parallel_tsv import fork_map, WorkerDied\n"
        "def work(x):\n"
        "    if x == 2: os.kill(os.getpid(), signal.SIGKILL)\n"
        "    time.sleep(0.1); return x\n"
        "try:\n"
        "    fork_map(work, range(6), 3)\n"
        "except WorkerDied as exc:\n"
        "    print('died', 'out of memory' in str(exc))\n"
    )
    proc = subprocess.run([sys.executable, "-c", code], cwd=REPO, capture_output=True, text=True,
                          timeout=60)
    assert proc.returncode == 0, proc.stderr
    assert proc.stdout.strip() == "died True"
