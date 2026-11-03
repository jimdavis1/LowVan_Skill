# Adenoviridae — progress log

Started 3 Nov 2026 01:46 UTC. Append-only; newest at the bottom. Every entry
is timestamped so a crash can be resumed from the last line.

Scope settled by the curator 3 Nov: build **Mastadenovirus** and
**Aviadenovirus**; measure the next-largest genus and decide. Feature depth:
**every known protein and mature peptide, including signal peptides.**

| | |
|---|---|
| cost meter baseline | transcript record 8607, 2026-11-03T01:46Z |
| tally | `scripts/tally_cost.sh` |

---
- 01:48Z Measured the family holding. 43,636 records; 5,351 >=5 kb; 5,015 >=20 kb.
- 01:48Z Genus holdings >=20 kb: Mastadenovirus 4,362 | Aviadenovirus 450 | blank 80 | Siadenovirus 66 | Barthadenovirus 55 | Ichtadenovirus 2 | Testadenovirus 0.
- 01:48Z CORRECTION to the scoping report: Atadenovirus is NOT absent. BV-BRC files it as Barthadenovirus (duck/ovine/bovine atadenovirus A-D, snake, lizard, bearded dragon), 55 near-complete. Four candidate genera, not three.
- 01:48Z Splicing feasibility settled on NC_001405 (human adenovirus 2): 30 single-exon CDS, 8 two-exon, ZERO with 3+ exons. The tripartite leader is 5' UTR, not CDS. One intron is what the splice caller handles.
- 01:48Z min-len set to 20000, measured: 25000 would cut Siadenovirus 66 -> 54 and that band (23,811-27,209 nt) is the genus's own size range. 20000 costs 40 partial records family-wide.
- 01:48Z DUMP LAUNCHED: dump.sh -> dump.log, prefix A, ~5,015 genomes.
- 01:52Z DUMP COMPLETE in 4 min and check_dump.py PASS. 5,015 genomes, 155,187 feature rows (139,560 protein-coding, 15,627 blank-md5 non-coding), 21,947 distinct sequences. Index alignment clean.
- 01:52Z INTAKE RULE (SKILL §1) applied: dropped 139 'sp.' genomes -> intake/DROPPED_sp.tsv. NOTE the character differs from the parvo case: 118 are 'Human adenovirus sp.' HAdV-B3 clinical isolates, i.e. a submitter naming artefact on an abundantly-represented type, not an unclassified bin. 5 'Siadenovirus sp.' and 3 'Aviadenovirus sp.' ARE the DensovirinaeAmbi case. Dropping all of them costs 3% of Mastadenovirus and nothing in diversity.
- 01:52Z 67 genomes carry NO genus (bat/raccoon/porcupine/Adenovirus capreolus etc.) -> intake/NO_GENUS.tsv. Named viruses but unassignable to a module; held out of training data, kept as a routing-evaluation set.
- 01:52Z IN SCOPE: Mastadenovirus 4,244 | Aviadenovirus 447 | Siadenovirus 61 | Barthadenovirus 55 | Ichtadenovirus 2.
- 02:01Z TRIAGE PASS 1 (build_cols.py, cols1.log). Design: masses and ORF numbers do NOT bind by name -- both go to a pool resolved by homology (blastp >=80% id / >=85% qcov) then mmseqs clustering. Result: Mastadeno 3,356 genomes / 42 named features; Avia 430 / 24; Sia 59 / 22; Bartha 54 / 29. Homology rescued 27,237 Mastadenovirus features out of mass-named strings.
- 02:01Z PASS 1 DEFECTS to fix: (a) Aviadenovirus leaves 7,163 features / 1,810 seqs unplaced = 46% of the genus -- the ORF pool is not being resolved. (b) Mastadenovirus bound only 3,356 of 4,244 in-scope genomes; 888 bound nothing. (c) E1A_243R emptied by consolidation -- 289R/243R are nested splice forms and consolidation must not arbitrate between them, same exemption parvo needed for VP1/VP2. (d) UNC floor too loose: 152 UNC clusters in Mastadeno, most are 1-2 sequences.
- 02:10Z Wrote RESUME_HERE.md (checklist + the 5 fixed defects) and MATURE_PEPTIDES.md.
- 02:10Z MATURE PEPTIDES: measured, not assumed -- ZERO mat_peptide features in all four reference genomes and ZERO mature-peptide deposits among the 21,947 distinct sequences (every 'protein VIII' deposit is 227 aa, the precursor length). They must be built from published adenain cleavage sites, not binned. The one signal peptide in the family is E3 gp19K (Mastadenovirus); SignalP 6.0 is not installed and the published site is the authority anyway.
- 02:10Z Wrote probes.json: 4 modules, every PMID resolved against PubMed first. Probe pass already paid for itself -- it found E3 6.7K (PMID 15507617), which no annotation string in the dump would have surfaced.
- 02:20Z PUSHED 345d0dfc to jimdavis1/LowVan_Skill: reports/adenoviridae/ with the resume doc, MATURE_PEPTIDES.md, probes.json, build_cols.py, dump.sh and the four vocabulary tables.
- 02:20Z PASS 3 -> PASS 4 fixes: (1) normalise() orphaned 40 rule alternatives written in deposit form ('maturation protein' -> 'maturation'); I_LEADER and SHUTOFF both shipped EMPTY because of it. Fixed structurally with expand(), which adds each literal alternative's normalised form at compile time, plus check_rules() which now ABORTS the run if any alternative can never match. (2) 909 Mastadenovirus features were going to NOT_MODELLED on bare region labels ('E3','E4','L2') and 302 on 'splicing factor' (which is 33K's activity, not a locus) -- all now go to the homology pool instead. (3) FIBER routed by what the GENOME carries: one fiber-family deposit -> FIBER, two or more -> longest FIBER1, rest FIBER2.
- 02:25Z PASS 4 EXPOSED AN ERROR I INTRODUCED. The fiber length rule ('longest = fiber-1') is wrong. MEASURED: of 215 Aviadenovirus genomes depositing both, fiber-1 is the longer in only 39 (18%) -- Fowl aviadenovirus 4 deposits f1=431/f2=479 across 99 genomes and only CELO (793/410) matches the intuition. 'fiber-1'/'fiber-2' are POSITIONAL (gene order), not length. The rule swapped 181 deposits each way. 'long'/'short' fibre DO mean length (247 of 252 Mastadenovirus genomes). Pass 5: explicit labels are trusted as deposited; only GENERIC 'fiber' is routed, by homology to the two labelled collections.

## 2026-11-03 13:2x UTC -- Siadenovirus and Barthadenovirus taken up; literature check

CURATOR DECISION: build Siadenovirus and Atadenovirus too. All four modules
are in scope. Curator also asked for the groups to be literature-checked for
(1) completeness, (2) coordinates, (3) FUNCTIONS -- "there are a lot without
functions here", meaning the 224 UNC clusters.

GENUS NAME: ICTV renamed Atadenovirus -> Barthadenovirus in the 2024
ratification (proposal 2023.008D). BV-BRC uses the new name; all the cited
literature predates it and says Atadenovirus. Module keeps Barthadenovirus.

Coordinates taken from annotated reference genomes, not assumed:
  TAdV-3  NC_001958  26,263 bp   Siadenovirus
  OAdV-7  NC_004037  29,576 bp   Barthadenovirus prototype
  SnAdV-1 NC_009989  27,751 bp   Barthadenovirus, reptilian branch
  CELO    reannotation PMID_claude_generated 14604445, Table 1  Aviadenovirus

DEFECTS FOUND AND FIXED THIS PASS

1. LH3 SHIPPED WITH 4 SEQUENCES WHILE 13 MORE SAT UNNAMED. Both atadenovirus
   reference genomes annotate the LH3 locus as "E1B 55K" -- NC_004037
   1916..3064 382 aa, NC_009989 2169..3290 373 aa -- and atadenoviruses have
   NO E1B region. All 13 deposits named E1B 55K / E1B-2 / "E1B protein, large
   T-antigen" hit LH3 at 42-77% identity over full length at 370-393 aa. The
   structure paper predicts exactly this: mastadenovirus E1B-55K was EXAPTED
   from atadenovirus LH3 (PMID_claude_generated 33789897, 28943338).

2. RH3 SHIPPED EMPTY AND IT IS NOT AN ABSENT GENE. Consolidation moved all 4
   RH3 deposits into RH4, plus 7 of RH2's 22. RH1/RH2/RH4/RH6 are a related
   PARALOGUE family -- they are F-box proteins and the paper says they are
   related (PMID_claude_generated 12504581). Added as a PARALOG set, with
   LH1/LH2/LH3 and E4.1/E4.2/E4.3.

3. RH0, RH1 AND RH6 WERE MISSING FROM THE FEATURE SET. The right-hand region
   is numbered RH0..RH6; the module declared RH2..RH5. 7 RH1 deposits went to
   UNASSIGNED and the rest were swept into ONE accessory cluster that mixed
   RH6 + RH1 + RH0 -- three genes, one would-be profile.

4. SHUTOFF SHOULD NEVER HAVE BEEN A FEATURE. Its 8 deposits are 703 aa and
   consolidation put every one into K100, correctly: the late 100K protein IS
   the host shut-off protein. A synonym, not a gene. Two profiles would have
   fired on one locus. Removed; strings moved to K100.

   3 and 4 TOGETHER ARE THE LESSON: an empty collection can be a real protein
   hiding under a wrong name OR a duplicate of a feature already declared.
   Both occurred here, in the same family, in the same pass. Neither can be
   assumed; each has to be chased.

5. I CLAIMED BARTHADENOVIRUS HAS pIX. IT DOES NOT. The old probe note said
   "the atadenovirus capsid structure resolved it together with p32K". The
   structure resolves the GENUS-SPECIFIC proteins LH2, p32K and LH3, and its
   point is that LH3 occupies pIX's STRUCTURAL ROLE. Neither reference genome
   annotates a pIX. Corrected to ABSENT.

6. PERFORMANCE BUG I INTRODUCED: the cluster-adoption pass called rescue()
   once per cluster, rebuilding the named-sequence BLAST database each time --
   365 rebuilds against 81,000 Mastadenovirus sequences. Killed at 25 min.
   Now ONE blastp over the whole pool with votes tallied per cluster.

FUNCTIONS ESTABLISHED (the curator's point 3)
  RH1, RH2, RH4, RH6  F-box proteins; drive degradation of cellular
                      regulatory proteins via SCF ubiquitin-ligase complexes.
                      Rare among viruses. PMID_claude_generated 12504581.
  LH3                 capsid cementing protein, trimeric beta-helix phage
                      tailspike fold, pIX-like triskelion. 33789897/28943338.
  p32K, LH2           genus-specific minor coat / structural proteins. same.
  SIALIDASE (ORF1)    genus-defining; bacterial-sialidase-like domain,
                      virion component; catalytic activity OPEN. 36016389.
  FIBER (Sia)         binds sialic acid on N-linked glycoproteins of B cells,
                      required for infection. 32459612.
  Avia ORF2/12/13     parvovirus NS1/Rep homologues with an INACTIVE
                      ATPase/helicase domain. 14604445.
  Avia ORF16          mono-ADP-ribosyltransferase (Pfam PF01129). 14604445.
  Avia ORF18/19       lipase domain (PF00151) + O-glycosylated S/T region.
  Avia ORF22          binds pRb, activates E2F. 14604445 / 10400747.
  Avia ORF1           functional dUTPase. 14604445.

STILL WITHOUT A FUNCTION, AND THE NUMBER IS KEPT ON PURPOSE: Sia ORF7, ORF8,
hyd, E3; Bartha RH0, RH3, RH5, LH1, E4.1; Avia ORF9/10/11/14/17/20. A
functional string the literature does not support is worse than a number.

PARTITIONING QUESTION OPEN FOR BARTHADENOVIRUS: 55 genomes split 24 ruminant
/ 23 bird / 7 reptile / 1 marsupial. SnAdV-1 carries a DIFFERENT right-hand
set -- ORF1 (112 aa) and 105R (178 aa) in place of most RH genes, only RH1
shared. Building as one module and measuring routing before deciding; a
3-way split would leave the reptilian module at 7 genomes.

## 2026-11-03 -- CURATOR CORRECTION: the database is not the boundary authority

Curator: "The genbank and bvbrc genomes are the same in that they could have
errors in their gene names and boundaries. Making sure our calls are
consistent is okay, but the literature is authoritative if you can find the
actual gene boundaries."

I had used the RefSeq feature tables of NC_001958, NC_004037 and NC_009989 as
the coordinate authority and said so in both inventories. That was wrong, and
I had the counter-example in hand in the same session: NC_004037 annotates the
LH3 locus as "E1B 55K" and atadenoviruses have no E1B region at all. I caught
that error by homology and literature and then went on treating the same
record as ground truth for every other boundary.

WHAT CHANGED
  refcoords.json -> dbcoords.json, with a _WARNING and a _provenance note per
  module. build_cols.py and check_coords.py both restate that these are
  CONSISTENCY checks and never a boundary authority. Both inventories gained a
  section 0 grading every source as Experimental / Computational / Database.

WHAT THE PROVENANCE TURNED OUT TO BE
  TAdV-3 NC_001958 is a PROVISIONAL RefSeq -- "has not yet been subject to
  final NCBI review" -- identical to AF074946. Its genes were called by
  "comparison of the DNA and predicted amino acid sequences with published
  sequences of other Ads" (PMID_claude_generated 9791022). Every Siadenovirus
  boundary in this module is a 1998 homology PREDICTION in an unreviewed
  record, except the 14 transcripts in the 2026 RACE map (41960683).
  The weak E3 call traces to that paper's "very low similarity to ovine E3".

  OAdV-7 has real experimental boundary work: transcription start sites and
  SPLICE JUNCTIONS for the lefthand end, IVa2, E2, P32K and E4 regions
  (PMID_claude_generated 9614874), the right-end RH transcripts (12504581),
  and proteinase cleavage sites measured on purified virions (8659111).

  PROVENANCE OF THE "E1B 55K" ERROR, TRACED: the 1996 genome paper proposed
  only TENTATIVELY that the lefthand group "probably represents the E1A/E1B
  regions". That hardened into /product="E1B 55K" in the database. The 2017
  and 2021 structure papers identified the same protein as LH3, a capsid
  cementing protein. The literature corrected itself; the database did not.

MATURE PEPTIDES -- AN ASSUMPTION I HAD TO WITHDRAW
  MATURE_PEPTIDES.md derived cleavage from the HUMAN adenain consensus
  (M/I/L)XGX/G for the whole family. 8659111 measured the sites on purified
  OAV287 virions and reports that some "did not conform to the known
  specificity of the human Ad proteinase". So that consensus cannot be used
  for Barthadenovirus. Also: p32K is itself a precursor -- the paper reports a
  processed 28 kDa virion polypeptide -- which this module does not model.
  Barthadenovirus mature peptides are BLOCKED on 8659111, not merely unbuilt.

BLOCKED ON PAYWALL -- NOT OPEN ACCESS, ABSTRACTS ONLY
  8659111  OAdV-7 genome + measured proteinase cleavage sites
  9614874  OAdV-7 transcription start sites and splice junctions
  12504581 OAdV-7 right-end transcript map + F-box identification
  18166240 SnAdV-1 genome completion
  These are the authoritative boundaries for Barthadenovirus and they are NOT
  in the module. Named gap, not a silent one. Ask the curator for PDFs.

## 2026-11-03 -- THE FOUR PAPERS ARE IN. Boundaries now from the literature.

Curator supplied all four PDFs. Read in full (two-column extraction by line
x-position; block-level extraction interleaved the columns and WOULD have
produced misread cleavage sites).

PMID_claude_generated 8659111 Vrati 1996 -- OAdV-7 genome + MEASURED cleavage
PMID_claude_generated 9614874 Khatri & Both 1998 -- transcripts and promoters
PMID_claude_generated 12504581 Both 2002 -- right end transcripts, F-box
PMID_claude_generated 18166240 Farkas 2008 -- SnAdV-1

### THE "E1B 55K" NAME IS SETTLED, FROM THE SOURCE
Khatri & Both 1998 Table 1 lists the lefthand ORFs as LH1 9.7 kDa, LH2 14.7
kDa, LH3 42.8 kDa, p32 32.1 kDa, IVa2 37.5 kDa. "(E1B 19 kDa)" and "(E1B 55
kDa)" appear in that table in PARENTHESES with footnote e: "Shows homology
with the Ad5 protein". They are homology NOTES, not identities. The database
turned a parenthetical note into a /product=. LH3 42.8 kDa ~ 382 aa and the
module's LH3 collection has median 382. Settled.

### MATURE PEPTIDES -- NOW MEASURED, NOT ASSUMED
My earlier plan applied the HUMAN adenain consensus family-wide. The paper
shows that would have been wrong FOUR separate ways, each measured:
  - pX is NOT cleaved in OAV. Ad2 res32 / Ad12 res27 sites absent; Ad12 res42
    site not used (P1' is Phe). Table 2 lists no mature size for X.
  - pVI N-terminal type I site is MASKED -- Phe not Ala at P1', residue 29.
  - pVII uses TWO NON-CONFORMING sites, NTGW|G and MRAT|G. Its ORF does not
    even predict a precursor, and there is no site at residue 25 where Ad12
    uses MYGG|A.
  - p32K has no human counterpart at all.
The paper gives the extension: OAV also recognises (N/M/L/I).X.G.X|G and
(M/L/I).X.A.X|G. N at P4 and A at P2 are outside human specificity.

CONFIRMED AGAINST OUR OWN COLLECTIONS, EXACTLY:
  p32K  LAGG|AYR at res 12; mature N-term AYRKRKKKFQRPKPRK == Table 1 peak 7.
        286 -> 274 aa. Genome 114430 deposits BOTH forms, the shorter an exact
        C-terminal substring of the longer -- one collection, two forms, the
        Parvoviridae VP1/VP2 defect again.
  pVII  NTGW at position 9 in 16 of 16 deposits. OAdV-7 111 aa splits
        1-13 / 14-23 / 24-111. Fragment 14-23 is GIGTSSMRAT -- character for
        character the paper's measured peptide II. Mature VII starts
        GLKFSKKQPVRVRPYYR -- character for character Table 1 peak 9. 88 aa.
  pVIII measured peptide GSNRLTSYPDIPPILK starts at residue 118, preceded by
        SELNGW, so the site is LNGW|G. Fragments 117 aa and 100 aa against the
        paper's predicted 13.2 and 11.4 kDa.
  pVI   C-terminal type II site gives the 11-residue proteinase cofactor
        GVNFNTRRYCY; active consensus G.XXXXX.R.R.[yr].C.[fy].

DO NOT MODEL, on the paper's own words: pX (not cleaved), pVI N-terminus
(masked), pIIIa (a type I site LIGH|G sits 20 res from the C-terminus but no
peptide was isolated -- "unclear whether or not the site is used").

### A STRUCTURAL SCAN, REPORTED BUT NOT MODELLED
Searching every Barthadenovirus collection for a shorter deposit that is an
exact C-terminal substring of a longer one FROM THE SAME GENOME found six
pairs: P32K 286->274 (confirmed), PIIIA 573->516, DBP 381->363, RH2 186->103,
RH4 157->144, E4_3 225->217. Only P32K has literature behind it. The other
five are recorded as OBSERVATIONS. A nested shorter deposit is equally
consistent with a submitter annotating a downstream start codon, and nothing
in the data separates the two. None is being built.

### ALSO FROM THE 1996 PAPER
  Penton lacks the RGD motif human adenoviruses use to engage integrins.
  Proteins V and IX are absent from the genome -- confirms the pIX correction.
  Authors' own section accessions: U40839, U31557, U40837, U18755, U40838.
