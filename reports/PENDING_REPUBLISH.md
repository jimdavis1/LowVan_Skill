# Pending republication

The Crinivirus module was corrected twice on 29 September 2026 — 23 features
to 15, then back to 17 when the duplicate-call check was found to be missing a
containment floor. **The committed HTML in `reports/` is current; the live
Artifact pages are not.**

The artifact service began refusing content fetches at about 14:38 (listing
works, `action: "read"` and in-place publish both fail with
`artifact content fetch failed (network error)`). These four pages need
republishing to their existing URLs when it recovers:

| page | URL | live shows | should show |
|---|---|---|---|
| Crinivirus String Collapse | https://claude.ai/artifact/9vSCiWfRo6gZzrfYcvWzYQ | 15 targets | 17 targets, 204 → 17 (12.0×) |
| Crinivirus PSSM Registry | https://claude.ai/artifact/NvZDJv9BTWvY2B4dz1F1dQ | 15 features, 95 PSSMs | 17 features, 107 PSSMs |
| Crinivirus Vocabulary Saturation | https://claude.ai/artifact/FNs7pvEmyyGY3KjHdooNV1 | 14-string module curve | 16-string module curve |
| Crinivirus Coverage Audit | https://claude.ai/artifact/8sUQEueUJmA6CkFU4nziGF | 15 features, 7 essential | 17 features, 8 essential |

Republish with:

    Artifact(file_path="reports/crinivirus-<page>.html", url="<the URL above>")

Nothing else is affected — the other modules' pages match their committed
sources.
