"""TrainCraft command-line interface (thin shell over the public API)."""

from __future__ import annotations

import logging
from pathlib import Path

import typer

from .config import dump_starter_config, load_config
from .core import available
from .orchestration import STAGE_ORDER, run_pipeline, run_stage

app = typer.Typer(
    add_completion=False,
    help="TrainCraft: modular MLIP dataset generation & active learning.",
)


def _setup_logging() -> None:
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s"
    )


@app.command()
def run(config: Path = typer.Argument(..., exists=True, readable=True)) -> None:
    """Run the entire workflow declared in CONFIG (a single TOML)."""
    _setup_logging()
    cfg = load_config(config)
    summary = run_pipeline(cfg)
    typer.echo("Done:")
    for key, value in summary.items():
        typer.echo(f"  {key}: {value}")


@app.command()
def stage(
    name: str = typer.Argument(..., help=f"one of: {', '.join(STAGE_ORDER)}"),
    config: Path = typer.Argument(..., exists=True, readable=True),
    force: bool = typer.Option(False, "--force", help="recompute even if the artifact exists"),
) -> None:
    """Run a single pipeline stage standalone (reads/writes workspace artifacts).

    This is what the Slurm executor dispatches as separate jobs; each stage picks
    up the previous stage's artifact from the run workspace.
    """
    _setup_logging()
    cfg = load_config(config)
    frames = run_stage(name, cfg, force=force)
    typer.echo(f"stage '{name}': {len(frames)} frames")


@app.command()
def validate(config: Path = typer.Argument(..., exists=True, readable=True)) -> None:
    """Validate CONFIG and print the resolved stages."""
    cfg = load_config(config)
    stages = [
        name
        for name in ("geometry", "calculator", "sampling", "selection", "labeling", "dataset")
        if getattr(cfg, name) is not None
    ]
    typer.echo("OK: config is valid")
    typer.echo(f"  run: {cfg.run.name}")
    typer.echo(f"  stages: {', '.join(stages) or '(none)'}")


@app.command()
def new(path: Path = typer.Argument(...)) -> None:
    """Write a starter config to PATH."""
    if path.exists():
        typer.echo(f"refusing to overwrite existing file: {path}")
        raise typer.Exit(code=1)
    path.write_text(dump_starter_config())
    typer.echo(f"wrote starter config: {path}")


@app.command()
def plugins() -> None:
    """List registered plugins by kind."""
    for kind in ("source", "builder", "transform", "calculator", "sampler", "selector"):
        typer.echo(f"{kind}: {', '.join(available(kind)) or '(none)'}")


if __name__ == "__main__":  # pragma: no cover
    app()
