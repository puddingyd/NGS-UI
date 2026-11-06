"""Legacy entry point for the single shared daily file index."""


def main():
    from .file_index_update import main as update
    return update()


if __name__ == "__main__":
    # Old installed systemd units still invoke this module. Read the running
    # web service's path environment before importing application config.
    import runpy
    from pathlib import Path
    runpy.run_path(str(Path(__file__).resolve().parents[3] / "scripts/update_file_index.py"), run_name="__main__")
