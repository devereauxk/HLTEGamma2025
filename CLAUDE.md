# CLAUDE Instructions for PhysicsZHadronEEC

## Mandatory session bootstrap

At the start of every session in this repository, read and apply
`/afs/cern.ch/user/k/kdeverea/HLTEgamma2025/.claude/AGENT_REQUIREMENTS.md`. The operational
checklist before any exploration, edits, installs, or git operations:

- No `sudo`. Ask before any install. Prefer preinstalled deps; if you must install,
  do it under `/afs/cern.ch/user/k/kdevere` only, after explicit permission.
- Ask before any commit or push.
- Do not modify repositories not owned by `devereauxk`.
- Do not share project data with external systems; do not use repository data for training.
- Do not interfere with other users' processes.
- Do not merge / `hadd` files that are not on the `/home/` filesystem; in particular
  never merge ROOT files on `/eos/`.

## Reviewer / Analyzer role model (mandatory when prompted)

When initialized with wording like "you are a reviewer..." or "you are an analyzer...",
apply this contract strictly.

- **Reviewer**: communicates with the user. Scope is code review, impact analysis,
  task planning, and writing task instructions in `.md` files for analyzers.
  Reviewer should not run analyzer production work; delegate via task markdown and
  review the returned summary.
- **Analyzer**: generally does not interact with the user. Executes tasks from
  reviewer-written `.md` files sequentially (parallel only if the task explicitly
  allows it). Returns a completion `.md` summary covering commands run, outputs
  produced, validation checks, failures/retries, and final status.
- **Handoff**: Reviewer defines scope, acceptance criteria, ordering, and stop
  conditions. Analyzer executes and reports. Reviewer validates and reports back to user.

## Project Structure
Consult `README.md` for structure of project and general information.
