# 28 pages regenerated but NOT yet live (7 October 2026)

The 3D annotation fix (`RNA-directed RNA polymerase 3D` ->
`RNA-dependent RNA polymerase`) is complete in the repo, the runtime and S1
Table, and all 42 affected pages are regenerated on disk. 14 reached
claude.ai before publishing started failing.

## Why it stopped

`*.frame.claudeusercontent.com` does not resolve to Anthropic from this
machine:

```
frame.claudeusercontent.com is an alias for honeypot-malware.anl.gov
honeypot-malware.anl.gov has address 146.139.125.24
```

Argonne is DNS-sinkholing that domain to a malware honeypot. The Artifact
tool reads a page over that host to confirm the target is not a review page
before overwriting it, so the publish is refused. Each attempt also lands in
ANL's security logs as a suspected malware callback, which is why this was
stopped rather than retried.

**To finish:** ask ANL networking to allow `*.frame.claudeusercontent.com`,
or republish these from a network that resolves it normally. Nothing needs
rebuilding -- the files on disk are correct and committed.

## Already live (14)

enterovirus x3, aphthovirus x3, teschovirus x3, parechovirus x3,
cardiovirus-pssm-registry, cardiovirus-string-collapse

## Still stale on claude.ai (28)

| page | artifact id |
|---|---|
| cardiovirus-coverage-audit | CWpYjvdg3sbjLaRp68LFti |
| senecavirus-coverage-audit | Xr21mDLyZzcSebbLyeU1T7 |
| senecavirus-pssm-registry | XJMxY4EvPdae8nU5LNiJfq |
| senecavirus-string-collapse | Tbtxi52p4gQ7wdFwbsTPLq |
| avihepatovirus-coverage-audit | AUmXPWzEr4vCreQShkZsrz |
| avihepatovirus-pssm-registry | JFVvSbRmU1Fv6nq7Wqi7gC |
| avihepatovirus-string-collapse | 6RHPEZykKvhFBP9wrrH1qU |
| hepatovirus-coverage-audit | 4EgbaxJcXnos8VjP5o5qwq |
| hepatovirus-pssm-registry | HPeVWDKLb1684L7Lksqrug |
| hepatovirus-string-collapse | 93JNSJkTL6tnxWAqqh1Ehj |
| kobuvirus-coverage-audit | 644ijYfi9wT8E9YRUkxu31 |
| kobuvirus-pssm-registry | HXpLhyF1AdFbSbKPFcLgEp |
| kobuvirus-string-collapse | G5BNV47HTuCKGDHstSCYkD |
| sapelovirus-coverage-audit | 88yMxbSKFycQvTt6m6NvZ2 |
| sapelovirus-pssm-registry | R3NPTkzAW38A1af9UtGrT1 |
| sapelovirus-string-collapse | 57VNf18a2R6sdpqC3k3JAR |
| pico_ldr_vp4-coverage-audit | B2Hy5H5UF7VbawdnjMMQCU |
| pico_ldr_vp4-pssm-registry | Pfmpmk7RFRSj9pYXsL4qPZ |
| pico_ldr_vp4-string-collapse | LVc4VEW9iCpooM6ou4iyu7 |
| pico_ldr_vp0-coverage-audit | RYMzkjW5JdUu8pHae7eBUz |
| pico_ldr_vp0-pssm-registry | D3DZs4rxzkWXYR8PgSi7Lf |
| pico_ldr_vp0-string-collapse | DsaQvYVYzTqnmZvFZ93BUc |
| pico_noldr_vp4-coverage-audit | 3Ky6NnngVVgqbc2vTJxL1g |
| pico_noldr_vp4-pssm-registry | CgHakkLnCawRjWZdawJwWB |
| pico_noldr_vp4-string-collapse | D4ncL6xGyYttgDBcpqDJzF |
| pico_noldr_vp0-coverage-audit | 3sFkz4f99xuw9FGSzgxE8Z |
| pico_noldr_vp0-pssm-registry | Lr9uABqCQZxhBqiw83BJSS |
| pico_noldr_vp0-string-collapse | FPxCxJr26eRtF8PJY8wH4p |

Only the label differs; no measurement changed. The pages are correct in the
repo, so this is a publishing backlog, not a data problem.
