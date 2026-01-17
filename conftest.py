"""Configuration file for pytest at project root."""
import sys
from pathlib import Path

# Add src-python to the path  
src_path = Path(__file__).parent / "src-python"
if str(src_path) not in sys.path:
    sys.path.insert(0, str(src_path))
