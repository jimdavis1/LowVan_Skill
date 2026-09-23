# Reports

The source HTML for every published page. **Artifacts are published to
claude.ai, not to GitHub** — opening a file here in a browser renders it
locally, and the URL column is where the shared copy lives. Keeping the source
here means any page can be republished without rebuilding its measurements.

Each page is generated from a facts file by a script in `generators/` or
`../skill/scripts/gen_*.py`. Edit the facts and regenerate; do not hand-edit
the HTML.

## Published pages

| taxon | coverage audit | PSSM registry | string collapse | vocabulary saturation |
|---|---|---|---|---|
| Orthoflavivirus | [audit](orthoflavivirus-coverage-audit.html) · [live](https://claude.ai/artifact/RokLyf9pujJnw37CAri2DR) | [registry](orthoflavivirus-pssm-registry.html) · [live](https://claude.ai/artifact/AkjDvmo28d4bVZCndjymXP) | [strings](orthoflavivirus-string-collapse.html) · [live](https://claude.ai/artifact/MFXjT53roYWz1ymQkrh1ur) | [saturation](orthoflavivirus-vocabulary-saturation.html) · [live](https://claude.ai/artifact/9c5gUBxFf4jHNfYizWHw7m) |
| Hepaciviridae | — | — | [strings](hepaciviridae-string-collapse.html) · [live](https://claude.ai/artifact/SuSzzD4Lgid5Nqh4wsAUsw) | — |
| Alphaflexiviridae | [audit](alphaflexiviridae-coverage-audit.html) · [live](https://claude.ai/artifact/MzC77du44ZP2iBnPogQPx1) | [registry](alphaflexiviridae-pssm-registry.html) · [live](https://claude.ai/artifact/TseN29Pu2NS87fMCJRTP1a) | [strings](alphaflexiviridae-string-collapse.html) · [live](https://claude.ai/artifact/1opeHxTHrAVEx7nsy8DMc1) | [saturation](alphaflexiviridae-vocabulary-saturation.html) · [live](https://claude.ai/artifact/DrNtpzKawMfnQvEFgdvJCP) |
| Allexivirus | [audit](allexivirus-coverage-audit.html) · [live](https://claude.ai/artifact/KAyMy13QTKUsFASLRsoTwo) | [registry](allexivirus-pssm-registry.html) · [live](https://claude.ai/artifact/PR148JCFeHCx5rsiqSPe8d) | [strings](allexivirus-string-collapse.html) · [live](https://claude.ai/artifact/Nod4NnVCU3arrmBqiTPh1J) | [saturation](allexivirus-vocabulary-saturation.html) · [live](https://claude.ai/artifact/S5H57Gg1N2DECvutD8JBTS) |
| Hepeviridae | [audit](hepeviridae-coverage-audit.html) | [registry](hepeviridae-pssm-registry.html) | [strings](hepeviridae-string-collapse.html) | [saturation](hepeviridae-vocabulary-saturation.html) |
| Tobamovirus | [audit](tobamovirus-coverage-audit.html) | [registry](tobamovirus-pssm-registry.html) | [strings](tobamovirus-string-collapse.html) | [saturation](tobamovirus-vocabulary-saturation.html) |
| Trivirinae | [audit](trivirinae-coverage-audit.html) | [registry](trivirinae-pssm-registry.html) | [strings](trivirinae-string-collapse.html) | [saturation](trivirinae-vocabulary-saturation.html) |
| Matonaviridae | [audit](matonaviridae-coverage-audit.html) · [live](https://claude.ai/artifact/AxuNuwq21AAiGws5JvQoWo) | [registry](matonaviridae-pssm-registry.html) | [strings](matonaviridae-string-collapse.html) | [saturation](matonaviridae-vocabulary-saturation.html) |
| Togaviridae | [audit](togaviridae-coverage-audit.html) | [registry](togaviridae-pssm-registry.html) | [strings](togaviridae-string-collapse.html) | — |
| Rhabdoviridae | [audit](rhabdoviridae-coverage-audit.html) | [registry](rhabdoviridae-pssm-registry.html) | — | — |

Not per-taxon: [the module handbook](lowvan-module-handbook.html),
[the BV-BRC string collapse](bvbrc-string-collapse.html) and
[the Togaviridae build anatomy](togaviridae-build-anatomy.html).

Hepaciviridae's other three pages are in progress. A taxon with four pages has
had the whole series done; a dash is a page that was never made, not one that
was made and lost.

## Where the live copies are, if a URL here is stale

`/artifacts` in the Claude Code terminal, or the gallery at
`claude.ai/code/artifacts`. Artifacts are private until shared from the page's
own share menu, so a URL here will not open for someone who has not been
given access.
