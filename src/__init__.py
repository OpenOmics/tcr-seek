from pathlib import Path

version = (Path(__file__).resolve().parents[1] / "VERSION").read_text().strip()
