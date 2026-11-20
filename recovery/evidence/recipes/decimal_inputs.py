"""Extract the nine versioned literal inputs; do not embed their coordinates here."""
import sys
from pathlib import Path

source, destination = map(Path, sys.argv[1:])
driver = source.read_text()
block = driver[driver.index("struct wkt_input"):driver.index("std::vector<input> make_inputs")]
if block.count('{"') != 9:
    raise ValueError("unexpected benchmark revision")
with destination.open("x") as f:
    f.write(block)
