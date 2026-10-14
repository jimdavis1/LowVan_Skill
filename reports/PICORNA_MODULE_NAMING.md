# Why the four grouped modules are not named after a taxon

Asked and checked 14 October 2026. **Decision: keep `PICO_LDR_VP0`,
`PICO_NOLDR_VP0`, `PICO_LDR_VP4`, `PICO_NOLDR_VP4` as they are.**

Ten of the fourteen picornavirus modules already follow the convention and are
named for their genus -- Enterovirus, Kobuvirus, Cardiovirus, Parechovirus,
Teschovirus, Sapelovirus, Hepatovirus, Senecavirus, Avihepatovirus. The four
above are the grouped ones: 3-6 genera each, pooled because no member holds
enough sequence to model alone (Tottorivirus has 10 genomes, Anativirus 38).
The names encode what separates them -- leader protein present or absent, and
whether VP0 is cleaved to VP4+VP2.

## The check

ICTV subfamily for every genus in the four modules, from the BV-BRC taxonomy
lineage:

| module | genera | subfamilies |
|---|---|---|
| `PICO_LDR_VP0` | Sicinivirus, Salivirus, Oscivirus, Passerivirus, **Shanbavirus** | Kodimesavirinae x4, **Paavivirinae x1** |
| `PICO_NOLDR_VP0` | Pasivirus, Avisivirus, Limnipivirus, Potamipivirus, Aquamavirus | Paavivirinae x5 -- **pure** |
| `PICO_LDR_VP4` | Bopivirus, Erbovirus, Mischivirus, Hunnivirus, Tottorivirus, **Anativirus** | Caphthovirinae x5, **Ensavirinae x1** |
| `PICO_NOLDR_VP4` | Cosavirus, Dicipivirus, Tremovirus | Caphthovirinae, Kodimesavirinae, Heptrevirinae -- **three for three genera** |

One of four is subfamily-pure. That alone would be reason enough, but the
converse is the decisive part: **every subfamily is split across modules**, and
the splits include the single-genus modules that are already correctly named.

    Caphthovirinae   -> PICO_LDR_VP4 (5 genera) + Cosavirus + Cardiovirus
    Kodimesavirinae  -> PICO_LDR_VP0 (4 genera) + Dicipivirus + Kobuvirus
    Paavivirinae     -> PICO_NOLDR_VP0 (5 genera) + Shanbavirus + Parechovirus
    Heptrevirinae    -> Tremovirus + Hepatovirus
    Ensavirinae      -> Anativirus + Enterovirus

So even `PICO_NOLDR_VP0`, the pure one, cannot take its subfamily name:
"Paavivirinae" would claim Parechovirus and Shanbavirus, which this module does
not annotate. **A taxon name here would assert coverage the module does not
have**, which is worse than a name that needs explaining -- an operator who
sees `Paavivirinae` on a module and feeds it a parechovirus gets nothing, with
no hint why.

## Why this was always going to happen

The project partitions on **genome organisation, not phylogeny**. Leader
presence and VP0 cleavage are exactly the characters that cut across the
picornavirus tree, which is why these four groups exist at all. A grouping
formed on one axis will not generally align with a classification built on
another, and when it does not, naming it as though it did is a false claim
rather than a convenience.

## If this is revisited

The useful change would not be renaming. It would be **splitting a grouped
module once a genus has enough sequence to stand alone** -- Salivirus is the
nearest candidate, being the only human pathogen in `PICO_LDR_VP0` -- at which
point it becomes a genus-named module like the other ten and the residue keeps
the descriptive name. That is the same lever as the Betaflexiviridae split,
and it improves routing as well as naming.
