# Concepts

These pages explain the *why* behind TrainCraft's design — the abstractions,
conventions, and invariants that govern every part of the system.

You don't need to read these to use TrainCraft, but they'll help you:

- Write custom plugins (builders, calculators, samplers)
- Understand why your configs behave the way they do
- Debug unexpected results
- Design extensions that fit naturally into the architecture

| Concept | What it covers |
|---|---|
| [Architecture](architecture.md) | The full pipeline and plugin system |
| [Geometry System](geometry.md) | Sources, builders, and transform composition |
| [Fragment Identity](fragments.md) | How atoms are grouped for MC sampling |
| [Provenance](provenance.md) | How every frame records its history |
| [Selection Funnel](selection.md) | The stages and their ordering |
