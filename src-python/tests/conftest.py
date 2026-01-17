"""Configuration file for pytest."""
import sys
from pathlib import Path

# Add src-python to the path
src_path = Path(__file__).parent.parent
if str(src_path) not in sys.path:
    sys.path.insert(0, str(src_path))

# Debug: print sys.path to see what's happening
print(f"DEBUG conftest.py: Added {src_path} to sys.path")
print(f"DEBUG conftest.py: sys.path = {sys.path}")
