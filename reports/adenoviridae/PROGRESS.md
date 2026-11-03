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

## 2026-11-03 -- curator corrections 1-3 applied; two silent defects found

CURATOR: "1. split them and be clear about what is the precursor and the mat
peptide. 2. Don't build them unless you can find literature, or base them on
reasonable homology. 3. go ahead with barthadenovirus."

### 1. PRECURSOR / MATURE SPLIT -- DONE
P32K -> P32K_MAT (274 of 286 aa, 12 removed). The test needs TWO signals and
both must hold: the shorter sequence is an EXACT C-terminal substring of a
longer deposit FROM THE SAME GENOME, and it does NOT start with Met, because a
cleaved N-terminus cannot. Nesting alone also matches a submitter annotating a
downstream start codon, which still begins with Met; the Met test alone also
matches a partial deposit. Naming follows the kit: P32K is the CDS precursor,
P32K_MAT the mat_peptide with upstream_ext 0.

### 2. MATURE PEPTIDES -- ONE TRANSFERS, THE REST DO NOT
pVII's first cleavage is MEASURED IN TWO GENERA INDEPENDENTLY. Table 1 of
PMID_claude_generated 8659111 sequenced peptide II from purified virions of
BOTH OAV287 (GIGTSSMRA) and human Ad5 (GLRFPSKMF), and both sit immediately
after .GW|G. The paper also states pVII's N-terminal 16 residues are "highly
conserved between viruses of human, animal, and avian origin".
MEASURED IN OUR OWN COLLECTIONS, all four modules:
  Mastadenovirus  MSILISPSNNTGWG  cut after 13   281/302
  Aviadenovirus   MSILISPNDNRGWG  cut after 13    38/38
  Barthadenovirus MSILVSPSDNTGWG  cut after 13    15/16
  Siadenovirus    MHSVVYSPGDSRGWG cut after 14    11/11
So pVII is built everywhere. NOT transferred: the MRAT|G second cut (the paper
says it does not conform to human specificity, and it is essentially absent
outside Barthadenovirus), p32K (genus-specific), pVIII and pVI (no TAdV-3
evidence). Decisions per candidate in MATURE_PEPTIDE_DECISIONS.md.

### 3. BARTHADENOVIRUS -- building as one module, routing measured later.

### DEFECT: 143 STALE COLLECTIONS FROM EARLIER PASSES
build_cols.py wrote collections but never REMOVED ones it no longer produces.
Pass 1 declared 152 Mastadenovirus UNC clusters, pass 11 declares 55, and
UNC56..UNC152 were still on disk holding 1-2 sequences each -- below the
floor, belonging to a clustering that no longer exists, and absent from every
summary the script prints. SHUTOFF was still there too, two passes after it
was deleted as a duplicate of K100. The next workflow step globs *.fasta and
would have built PSSMs from all of them. Nothing would have errored.
Removed 97 + 29 + 8 + 9 = 143. The removal is scoped: .fasta only, inside the
module's own collections directory, nothing recursive, synonyms.tsv untouched.

### DEFECT: find_misfiled ACCUSED THE WRONG GENOME, AND I ALMOST BELIEVED IT
First run reported `Penguin siadenovirus A` (2005055.4) belonged in
Mastadenovirus on the strength of three >=98% HEXON partners. A >=90%
cross-genus hexon identity is extraordinary, so it was checked instead of
applied. The three partners are 10535.510/511/512, all "unidentified
adenovirus" -- BV-BRC taxon 10535, a single bucket of 111 genomes of every
genus, filed under Mastadenovirus by default. THE DIRECTION WAS BACKWARDS:
those three are siadenoviruses in the wrong module, and the one NAMED genome
in the clade was accused by the three unnamed ones.

Two fixes to the tool, not the output:
  - an unnamed genome ("unidentified", "unclassified", "untyped", " sp.") is
    never a VOUCHER, only ever a query. It carries no taxonomic claim.
  - a single >=95% full-length HEXON match is conclusive. The two-voucher rule
    assumes the right clade has more than one named member in the dump;
    here it has exactly one, so own=0 foreign=1 and all three were missed.
The bad entry was deleted from misfiled.json by hand -- the file merges rather
than replaces, which is right for preserving applied fixes but wrong for an
artefact.
RESULT: 10535.510, 10535.511, 10535.512  Mastadenovirus -> Siadenovirus.

### ALL FOUR MODULES NOW PASS THE LITERATURE CHECK
Mastadenovirus and Aviadenovirus had NO probe patterns, so their check could
not run at all. Both lists completed; 13 new PMIDs resolved against PubMed and
each title read. Added probes for the ten CELO ORFs now named as features.
  Mastadenovirus  29/30 PRESENT. ABSENT: VA RNA I/II -- correct, they are
                  non-coding RNAs and cannot appear in a protein dump.
  Aviadenovirus   17/19 PRESENT. ABSENT: E1 region and protein IX -- both
                  real biology, aviadenoviruses have neither.
  Siadenovirus     9/9,  Barthadenovirus 16/16, except the decided pIX ABSENT.

## 2026-11-03 -- PSSMs built for all four modules, in parallel

PARALLELISM USED: 4 concurrent PSSM builds, then 4 concurrent retries, then 4
concurrent leftover passes, plus rep contigs alongside. NOT parallelised:
build_cols.py (one process writes all four modules, two copies would race) and
coverage runs (three kernel panics on Parvoviridae at ~515 vouchers/genome
against a ~245,000 ceiling; Mastadenovirus is 3,353 genomes).

MODULE JSONs: 292 features, every one with an annotation string. The
manuscript rule is that the annotation says what the protein DOES and the
symbolic name goes in gene_symbol -- "Major capsid protein hexon" + symbol
"hexon", never "hexon protein". Adenoviridae contributes nothing to the
controlled vocabulary yet (all 35 taxa in it are RNA viruses), so these are
coined in house style. Discovered clusters get "Uncharacterized
lineage-specific protein N", the canonical string.
NUMBERS KEPT ON PURPOSE, because no paper supports a function: Avia
ORF9/10/11/14/17/20, Sia ORF7/ORF8, Bartha RH0/RH3/RH5/LH1/E4.1.
P32K_MAT carries upstream_ext 0 -- with the default of 1 the annotator scans
upstream from the cleaved N-terminus, finds the PRECURSOR's start codon and
emits the mature peptide on the precursor's own coordinates. That shipped on
rabies G and on six norovirus products at once.

### BARTHADENOVIRUS: 17 of 43 FEATURES BUILT NOTHING, INCLUDING CORE GENES
DBP, IVa2, 52K, IIIa, protease, pVI, pX and UXP all produced no profile. Not a
count problem -- they have 12-22 sequences each. TRACED ONE INSTEAD OF
GUESSING: E4_1 has 13 input sequences and ALL 13 went to Leftover_Seqs.aa.
Zero clusters formed, because at -mi 0.8 -mc 0.8 no two of the thirteen are
80% identical over 80% coverage. The genus spans ruminant, avian and
reptilian hosts and 0.8 is the wrong identity floor for it. This is the
documented Betarhabdovirinae case in references/pipeline.md.

-mi is NOT a coverage dial -- it is the identity at which two sequences are
declared the same protein. So the retry climbs DOWN one step at a time
(0.6, 0.5, 0.4, 0.3) and stops at the FIRST setting that yields a profile, so
nothing is clustered looser than it needs. Every rescue recorded in
BUILD_PARAMS with the setting that did it.
RESULT: all 10 core-gene failures rescued at mi=0.6 mc=0.7. Siadenovirus UXP
rescued at 0.6. Mastadenovirus UNC14 rescued.

THIS IS ALSO EVIDENCE ON THE PARTITIONING QUESTION: a genus that needs
mi ~0.6 on its CORE genes is one whose members barely cluster together.
Routing numbers to follow before any split is proposed.

### STILL WITHOUT A PROFILE AFTER THE LADDER
Bartha  LH1 (n=11), RH4 (n=8), UNC3 (n=7), RH6 (n=3), RH0 (n=2), RH3 (n=1),
        P32K_MAT (n=1)
Mast    UNC45 (n=4), UNC50 (n=3), UNC51 (n=3)
Avia    UNC81 (n=3), UNC85 (n=3)
Three or four sequences too divergent to cluster even at 30% identity probably
are not one homology group. The leftover pass (build_leftover_pssms.py) is the
designed route for sub-threshold sequences and is running now. Anything still
unprofiled will be REPORTED as declared-but-unprofiled, not quietly dropped: a
feature declared with nothing behind it makes the annotator report it missing
on every genome.

### REP CONTIGS
Siadenovirus      8 references cover 64/64 genomes = 100.0%
Barthadenovirus  14 references cover 55/55 genomes = 100.0%
Mastadenovirus and Aviadenovirus running.
