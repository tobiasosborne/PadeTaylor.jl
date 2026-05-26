# Memory index

- [Scale covariance is a core principle](scale-covariance-core-principle.md) — never hardwire an absolute length scale; h derives from local pole spacing.
- [Schwarz reflection out of scope](schwarz-reflection-out-of-scope.md) — user uncomfortable with it; don't propose it.
- [Long Julia agents: verbose + eager flush](long-julia-agents-verbose-flush.md) — long renders / Pkg.test runs need `flush(stdout)` phase markers + `stdbuf -oL` + log tail; never poll exit status at short intervals.
