# Offline laptop app

The same container that runs on a cluster is also a **zero-code app you can use with
no internet** -- on a plane, in a seminar room with flaky wifi, anywhere. Pull the
image once while you have a connection; after that everything runs locally in your
browser, computing exactly, with **your machine's cores/RAM and the compute limits
shown right in the GUI**.

---

## Run it

Pull once (with internet):

```bash
apptainer pull quiverlab.sif docker://ghcr.io/MarcoArmenta/quiverlab:latest
# or, with Docker:
docker pull ghcr.io/MarcoArmenta/quiverlab:latest
```

Then, offline, start the GUI and open the printed URL:

```bash
# Apptainer:
apptainer run quiverlab.sif gui
# Docker (publish the port):
docker run --rm -p 8000:8000 ghcr.io/MarcoArmenta/quiverlab:latest gui
```

Open **<http://localhost:8000>**. Draw or pick an algebra, choose a field and what
to compute, and read exact results with rendered mathematics -- no code, no account,
no network. Math renders from vendored KaTeX, and the print-ready worked-steps HTML
report (export to PDF from your browser) is self-contained, so nothing is fetched at
runtime.

---

## Memory and time, always visible

Because you are computing on your own hardware, the app is explicit about what that
hardware can do. The GUI detects and shows:

- your machine's **cores and RAM** (from the container's cgroup-/host-aware probe),
- an **estimate** of the size (cells) and **time** for the example you set up, and
- an **estimate of the memory** it will need, next to the worker's memory cap
  (`RLIMIT_AS`) -- so you can see the limits you are computing under before you hit
  them, rather than after.

An example that would blow past your laptop's memory is flagged up front (the same
sizing that routes oversized jobs off the instant tier on the web service), so you
can shrink the degree range instead of waiting for an out-of-memory stop.

GPUs are irrelevant here: quiverlab is exact CPU arithmetic. If your machine has a
GPU it is simply ignored.

---

## Precomputed examples

The image ships a **seeded example cache**: a set of examples computed at build time,
so opening them is instant and needs no compute at all. (The curated list is an open
decision; the mechanism and a placeholder manifest ship now.) When you open one, the
GUI notes it was served from the cache.

---

## Sending an example to a cluster

Anything you build in the offline GUI you can also run big on a cluster: use the
**"Config (YAML)"** export button to download the exact config, then follow
[Run on your HPC cluster](hpc.md). The batch container and the offline app are the
same image, so a config that renders here runs there unchanged.

---

## Platform note (honest)

- **Linux:** Apptainer or Docker, both work.
- **Windows:** Docker Desktop works (`docker run -p 8000:8000 ... gui`).
- **macOS:** use **Docker** (`docker run -p 8000:8000 ... gui`) or the plain pip
  CLI (`pip install "quiverlab[hpc,web]"; quiverlab-hpc gui`). Apptainer has no
  native macOS build -- it needs a Linux VM -- so on a Mac the Docker or pip paths
  are the offline app.

The pip path is fully offline too once installed: `quiverlab-hpc gui` serves the same
local app with no network calls.
