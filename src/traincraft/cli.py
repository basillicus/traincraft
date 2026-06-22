"""TrainCraft command-line interface (thin shell over the public API)."""

from __future__ import annotations

import logging
from pathlib import Path

import typer

from .config import dump_starter_config, load_config
from .core import available
from .orchestration import STAGE_ORDER, run_pipeline, run_stage, submit_slurm

app = typer.Typer(
    add_completion=False,
    help="TrainCraft: modular MLIP dataset generation & active learning.",
)


def _setup_logging() -> None:
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s"
    )


@app.command()
def run(
    config: Path = typer.Argument(..., exists=True, readable=True),
    force: bool = typer.Option(
        False, "--force", "-f",
        help="recompute every stage, ignoring cached artifacts (a clean slate)",
    ),
) -> None:
    """Run the workflow in CONFIG. Routes to Slurm if [orchestration].engine="slurm".

    Stages are cached: a stage reruns only if its config (or an upstream stage)
    changed since the last run. Pass --force to recompute everything regardless.
    """
    _setup_logging()
    cfg = load_config(config)
    if cfg.orchestration is not None and cfg.orchestration.engine == "slurm":
        jobs = submit_slurm(cfg, str(config))
        typer.echo("Submitted Slurm pipeline:")
        for job in jobs:
            jid = f" (job {job.job_id})" if job.job_id else ""
            typer.echo(f"  {job.stage}: {job.script}{jid}")
        return
    summary = run_pipeline(cfg, force=force)
    typer.echo("Done:")
    for key, value in summary.items():
        typer.echo(f"  {key}: {value}")


@app.command()
def submit(
    config: Path = typer.Argument(..., exists=True, readable=True),
    dry_run: bool = typer.Option(False, "--dry-run", help="write sbatch scripts but don't submit"),
) -> None:
    """Render + submit the workflow as dependency-chained Slurm jobs (Apptainer)."""
    _setup_logging()
    cfg = load_config(config)
    if cfg.orchestration is None or cfg.orchestration.slurm is None:
        typer.echo("CONFIG has no [orchestration.slurm] section")
        raise typer.Exit(code=1)
    jobs = submit_slurm(cfg, str(config), dry_run=dry_run)
    verb = "Rendered" if dry_run else "Submitted"
    typer.echo(f"{verb} {len(jobs)} stage(s):")
    for job in jobs:
        dep = f" after {job.depends_on}" if job.depends_on else ""
        jid = f" (job {job.job_id})" if job.job_id else ""
        typer.echo(f"  {job.stage}{dep}: {job.script}{jid}")


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
        for name in (
            "geometry", "calculator", "sampling", "selection",
            "labeling", "dataset", "training",
        )
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
    for kind in (
        "source", "builder", "transform", "calculator", "sampler", "selector", "trainer",
    ):
        typer.echo(f"{kind}: {', '.join(available(kind)) or '(none)'}")


if __name__ == "__main__":  # pragma: no cover
    app()
