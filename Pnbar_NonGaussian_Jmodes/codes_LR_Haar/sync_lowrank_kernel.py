#!/usr/bin/env python3
"""Copy fast_loop_hafnian.py from probability_nbar_J2 into codes_LR_Haar/lowrank/."""

from __future__ import annotations

import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
HAAR_ROOT = Path(__file__).resolve().parent
src = ROOT / "probability_nbar_J2" / "lowrank" / "fast_loop_hafnian.py"
dst = HAAR_ROOT / "lowrank" / "fast_loop_hafnian.py"
shutil.copy2(src, dst)
print(f"Copied {src.name} -> {dst}")
