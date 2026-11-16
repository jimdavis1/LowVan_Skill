#!/usr/bin/env python3
"""Coverage-audit facts + registry notes for the five Adenoviridae modules.
Numbers: pages/<M>.measure.json, <M>/ASSIGN.tsv. Prose: below, written by hand."""
import json,os,csv,collections
W=os.path.expanduser('~/lowvan_scratch/adeno_v3')
DATE='16 November 2026'
NAMES={'Siadenovirus':'Siadenovirus','Barthadenovirus':'Barthadenovirus (formerly Atadenovirus)',
 'Aviadenovirus':'Aviadenovirus','Mastadenovirus_MultiFiber':'Mastadenovirus, multi-fiber genomes',
 'Mastadenovirus':'Mastadenovirus, one-fiber genomes'}
SCOPE={
 'Siadenovirus':'Siadenoviruses of birds, reptiles and amphibians (turkey adenovirus 3, frog adenovirus 1, penguin and raptor siadenoviruses).',
 'Barthadenovirus':'Barthadenoviruses of ruminants, reptiles, birds and marsupials, with their lineage-specific RH and LH gene families.',
 'Aviadenovirus':'Aviadenoviruses of poultry and wild birds, including fowl adenoviruses A&ndash;E with one or two fiber genes.',
 'Mastadenovirus_MultiFiber':'Mastadenoviruses that carry more than one fiber gene: human adenovirus F (types 40 and 41) and the simian and rhesus adenoviruses whose fiber census places them here.',
 'Mastadenovirus':'Mastadenoviruses with one fiber gene: human adenovirus species A&ndash;E, the great-ape and Old World monkey adenoviruses, and the canine, bovine, porcine, bat, rodent and other mammalian lineages.'}
SHARED_GAPS=[
 ['Cleaved mature-peptide ends','A cleaved end lands where the tBLASTn HSP lands; nothing moves it onto the AVP motif','No setting does this: it needs an annotator step that snaps the end to the motif'],
 ['UXP three-exon product','get_splice_variant_features.pl supports one junction; UXP is called as the U-exon ORF only','Needs multi-junction support in the splice caller'],
 ['Unspliced calls of spliced genes start off-Met','pol, pTP and IVa2 are called on the body exon','Expected for body-exon calls; the leader exons are not modelled'],
 ['Products not modelled','pVIc, mu, pVIIn, the propeptides, iTP and the 52K products','11&ndash;19 aa peptides are too short for a profile']]
GAPS={
 'Siadenovirus':[['10 genomes have no cytosine','Deposits with 0.0% C (C&rarr;T converted); every essential gene is &ldquo;missing&rdquo;','Not a module problem; quality scoring flags them correctly'],
   ['VIII-N, VIII-C, TP and IIIa mature peptides','The AVP site is in fewer than 60% of siadenovirus precursors','Not installed; recorded as gaps']],
 'Barthadenovirus':[['No 22K feature','21 deposits named &ldquo;22K&rdquo; and 12 &ldquo;Pac22k&rdquo; are unassigned; 22K is a conserved central gene','A missing declaration, found by this audit; a build fix, not a cutoff'],
   ['VIII-C ends','7 of 73 calls start on the AVP site; most start 5&ndash;8 aa past a site that is present','HSP-bound cleaved end (see above)'],
   ['IIIa ends','59% of calls end on the site (51% in the subset)','HSP-bound cleaved end'],
   ['RH family','RH_A&ndash;RH_R and UNC_RHP_1&ndash;3 are positional orthogroups present in 0&ndash;36% of genomes','Lineage-specific by nature; never essential']],
 'Aviadenovirus':[['Unmodelled avian ORFs','ORF0 (280), ORF1B (261), ORF20A (257), ORF19 (157), ORF1C (142), ORF29 (127) deposits are unassigned','Features not yet declared; needs orthogroup curation, not tuning'],
   ['IIIa ends','72% of calls end on the AVP site (57% in the subset)','HSP-bound cleaved end'],
   ['Rare fiber classes','FIBER_A17&ndash;A33 are each called in 0&ndash;3% of genomes','Real but rare classes; never essential']],
 'Mastadenovirus_MultiFiber':[['Partial and interrupted genomes','Hexon is called in 85% and pol in 84% of near-complete genomes; checked genomes carry internal stops (e.g. 130309.457: DBP 775 bits and E1B 55K 954 bits, both with stops)','The sequence, not the profile; such genomes are scored poor'],
   ['FIBER_CL2, CL34, CL1B','Called in 5&ndash;13% of genomes','Real but lineage-restricted classes; never essential']],
 'Mastadenovirus':[['HAdV-D IIIa start','In 316 genomes the pIIIa precursor starts ~29 aa upstream of the IIIa product, which starts at the Met homologous to HAdV-5 MQR','Precursor and product disagree; unresolved'],
   ['VII in divergent one-fiber lineages','8 lemur/bat/tree-shrew genomes: the VII profile scores 73.6 bits but no stop-free HSP','Lineages not represented in the VII profiles'],
   ['3 genomes still route to MultiFiber','47680.6, 3075435.12, 2050578.14 after the rep-contig swap','25-reference budget is used'],
   ['14 genomes have no cytosine','e.g. HAdV-E 130308.18&ndash;21 (0.0% C)','Not a module problem; flagged poor correctly'],
   ['Lineage-specific E3/E4 genes','gp19K 87%, RID-&alpha;/&beta; 94&ndash;95%, E4orf2&ndash;4 94&ndash;96% of genomes but 59&ndash;72% of lineages','Primate-lineage genes; deliberately not essential']]}
MAPNOTE={
 'Siadenovirus':'Turkey adenovirus 3 carries the siadenovirus core and its few genus genes (sialidase, E3, hydrophobic protein, ORF7, ORF8).',
 'Barthadenovirus':'The map shows the 22K gap: 33K is called after 100K, and no 22K is called, because the module never declared it (section 07). The RH genes at the right end are the positional orthogroups.',
 'Aviadenovirus':'Fowl adenovirus 1 (Phelps) carries two fiber genes and lineage-specific ORFs at both ends; several deposited ORFs at the ends are not yet declared (section 07).',
 'Mastadenovirus_MultiFiber':'Human adenovirus 41 with its two fiber genes, the long and short fibers that put it in this module.',
 'Mastadenovirus':'A species C human adenovirus: E1 at the left end, the conserved core in the middle, E3 and E4 at the right end.'}
REF={'Siadenovirus':'41678.98','Barthadenovirus':'70333.3','Aviadenovirus':'190061.121','Mastadenovirus_MultiFiber':'10524.1021','Mastadenovirus':'129951.1438'}
STAND={
 'Siadenovirus':'Siadenoviruses carry the smallest adenovirus genomes, about 26 kb, and the core gene set with very few extras. The module calls every core gene in every genome that has sequence to call; the 10 poor genomes are deposits with no cytosine at all.',
 'Barthadenovirus':'Barthadenoviruses add two lineage-specific gene families, RH and LH, to the conserved core. The module calls the core in every near-complete genome and splits the RH family into positional orthogroups; the audit also found that 22K was never declared.',
 'Aviadenovirus':'Aviadenoviruses have the largest adenovirus genomes, about 44 kb, with dozens of lineage-specific ORFs at both ends. The core is called in essentially every near-complete genome; the open work is the avian ORFs the module does not yet declare.',
 'Mastadenovirus_MultiFiber':'Human adenovirus F and its simian relatives carry two or more fiber genes and are kept apart from the one-fiber mastadenoviruses. Core calling is high; most poor genomes here are partial or carry interrupted genes, which the checked examples confirm is the sequence and not the profile.',
 'Mastadenovirus':'The largest module: 3,821 genomes, about 85% of them human adenoviruses by deposited name, plus every one-fiber mammalian lineage. Core genes are called in 98&ndash;99% of near-complete genomes; the routing fix moved 36 divergent genomes back from the multi-fiber module.'}
def assign_counts(m):
    c=collections.Counter(r['verdict'] for r in csv.DictReader(open(f'{W}/{m}/ASSIGN.tsv'),delimiter='\t'))
    return c
for m in NAMES:
    M=json.load(open(f'{W}/pages/{m}.measure.json')); s=m.lower().replace('_','-')
    A=assign_counts(m); col=M['collapse']; q=M['quality']; tot=q.get('good',0)+q.get('poor',0)
    fold=[l for l in M['rarefaction'] if 'fold' in l][0].split('->')[-1].strip()
    cr=M['call_rate']; ess=M['essential']
    ro=M['routing']
    tiles=[[str(M['features']),'features declared'],[str(M['pssms']),'PSSMs'],
           [f"{ro['own_after']}/{ro['panel']}",'genomes routed to this module'],
           [f"{100*q.get('good',0)/tot:.1f}%",'subset genomes good quality'],
           [str(len(ess)),'essential features'],[fold.replace('-fold','&times;'),'fewer distinct strings at equal depth']]
    S=[]
    S.append(dict(num='01',title='What the module covers',lede=SCOPE[m],blocks=[
      dict(type='table',headers=['measure','value'],rows=[['genomes assigned (by sequence)',f"{M['genomes']:,}"],['70% nucleotide clusters',str(M['clusters70'])],
        ['near-complete genomes (&ge;90% of median length)',f"{M['nearcomplete']:,}"],['median genome length',f"{M['median_len']:,} nt"]]),
      dict(type='figure',**{'from':f'{s}-genome-map.html'}),
      dict(type='prose',html=MAPNOTE[m]),
      dict(type='prose',html='Genomes were assigned to modules from sequence alone: the genus by best hexon match, and one-fiber versus multi-fiber mastadenoviruses by a census of intact fiber genes classed by literature orthology (Podgorski 2023, PMID 36610683). Deposited names were never used to assign a genome.')]))
    S.append(dict(num='02',title='How it was built',lede='Membership by homology to curated seed collections; the deposited name is recorded but never decides.',blocks=[
      dict(type='table',headers=['step','effect'],rows=[
        ['homology assignment',f"{A['ASSIGNED']:,} unique proteins assigned (&ge;80% coverage both ways, best feature &ge;1.2&times; the runner-up)"],
        ['left out',f"{A['NO_HIT']:,} no hit, {A['AMBIGUOUS']:,} ambiguous, {A['LOW_COVERAGE']:,} low coverage, {A['PARTIAL_BY_NAME']:,} submitter-marked partial"],
        ['profiles',f"{M['pssms']} PSSMs over {M['features']-len(M['special'])} profile features"],
        ['special features',', '.join(M['special']) or 'none'],['mature peptides',', '.join(M['mat']) or 'none']]),
      dict(type='prose',html='Paralogue families were split into orthogroup features by genome position, and each is called by sequence, never by location. Features with no literature function are named &ldquo;Uncharacterized lineage-specific protein&rdquo;. Splice junctions come from the literature and were extended to other lineages by aligning donor and acceptor sites; every spliced product was translated and checked in frame.'),
      dict(type='prose',html='Self-recall (each feature&rsquo;s profiles searched against its own collection) is 100% for nearly every feature; the registry page lists every feature.')]))
    cut=[[k,str(b),src] for k,b,src in M['cutoffs']][:14]
    S.append(dict(num='03',title='Where the cutoffs sit',lede='Cutoffs left at the default are not listed; every one below was set from a measurement.',blocks=[
      dict(type='table',headers=['feature','bit_cutoff','how it was set'],rows=cut) if cut else dict(type='prose',html='No feature needed a cutoff away from the default.'),
      dict(type='prose',html='Every CDS scans to its stop codon (<code>downstream_ext 1</code>); mature peptides set both extension flags from their own termini. In the 916-genome check, every unspliced CDS call ends on a stop codon except calls that stop at an N gap or an IUPAC base in the assembly.')]))
    rb=[['25' if M['reps']==25 else str(M['reps']),'rep contigs'],[f"{ro['own_after']}/{ro['panel']}",'routed to this module']]
    rp=('All genomes assigned to this module route to it.' if ro['own_after']==ro['panel'] else
        f"Before the fix {ro['own_before']:,} of {ro['panel']:,} routed here. Rep contigs 12&ndash;14, which were the best match for 0, 2 and 0 genomes, were replaced with Lemur mastadenovirus 204, SAdV GZ2-7 and Murine adenovirus 3; a simulation with the annotator&rsquo;s BLASTn settings moved none of 1,188 correctly routed genomes (smallest margin 837 bits). {ro['own_after']:,} now route here; the 11 that do not are 7 HAdV-F genomes with broken fibers (correctly multi-fiber), one HAdV-F-like genome and 3 divergent genomes the 25-reference budget does not reach.")
    S.append(dict(num='04',title='Does it route?',lede='A genome no rep contig claims is never annotated at all.',blocks=[dict(type='tiles',items=rb),dict(type='prose',html=rp)]))
    show=ess+[k for k in ('HEXON','PVII','FIBER','DBP','PX','PVI','E1B19K','E4ORF6') if k in cr and k not in ess]
    bars=[[k,str(cr[k]),cr[k],'essential' if k in ess else ''] for k in show]
    low=min(cr[k] for k in ess)
    S.append(dict(num='05',title='Does it call the proteins?',lede=f"Calls per near-complete routed genome, full panel ({M['nearcomplete']:,} genomes).",blocks=[
      dict(type='bars',unit='%',rows=bars),
      dict(type='verdict',kind='ok' if low>=85 else 'warn',html=f"<b>Lowest essential feature: {low}%.</b> " + ('Every essential feature is called in at least 85% of near-complete genomes.' if low>=85 else 'Below 85% in this module; section 08 shows the genomes are partial or interrupted rather than missed by the profile.'))]))
    ce=[l for l in M['calls_audit']]
    S.append(dict(num='06',title='Are the calls right?',lede='Duplicates, stop codons and cleavage sites, measured on every call.',blocks=[
      dict(type='pre',text='Features with extra calls in one genome (full panel):\n'+('\n'.join('  '+l for l in ce) if ce else '  none')),
      dict(type='pre',text='Mature-peptide ends (916-genome subset): product inside its own precursor, cut ends on the AVP motif\n'+'\n'.join('  '+l for l in M['matpep_ends_dext'])),
      dict(type='prose',html='Extra calls listed as INSEPARABLE are documented multi-copy cases or single broken genomes; none is a second feature stacked on another feature&rsquo;s locus. Same-frame overlaps between different features are limited to single genomes.')]))
    S.append(dict(num='07',title='What it cannot do',lede='Named as prominently as what works.',blocks=[
      dict(type='table',headers=['gap','evidence','why tuning will not close it'],rows=GAPS[m]+SHARED_GAPS)]))
    fl=[[f.replace('Genome is','genome is').replace('Contig has','contig has'),str(n)] for f,n in M['quality_flags'][:6]]
    S.append(dict(num='08',title='What the run found',lede='Genome quality from the GTO pipeline on the 916-genome subset, with the essential genes set by the curator&rsquo;s rules.',blocks=[
      dict(type='tiles',items=[[str(q.get('good',0)),f"good ({100*q.get('good',0)/tot:.1f}%)"],[str(q.get('poor',0)),f"poor ({100*q.get('poor',0)/tot:.1f}%)"],[str(len(ess)),'essential features']]),
      dict(type='table',headers=['flag','n'],rows=fl),
      dict(type='prose',html='Essential genes (copy_num) follow the curator&rsquo;s rules: never on uncharacterized proteins, special features or lineage-specific genes, and excluded wherever a cutoff, copy-number or coverage problem is on record. Essential here: '+', '.join(ess)+'.'),
      dict(type='prose',html='The subset deliberately over-samples unusual genomes, so its quality rate is a floor for the module, not an estimate of it. Poor genomes are dominated by sequence problems: deposits with no cytosine, N gaps, partial genomes and genes interrupted by stop codons.')]))
    S.append(dict(num='09',title='What the vocabulary bought',lede='The point of LowVan is a controlled vocabulary. This is what it cost the free text.',blocks=[
      dict(type='tiles',items=[[f"{col['strings']:,}",'BV-BRC strings bound'],[str(col['targets']),'annotation targets'],[fold.replace('-fold','&times;'),'fewer distinct strings at equal depth'],[f"{M['synonyms'].get('UNASSIGNED',0):,}",'deposited features unbound']]),
      dict(type='figure',**{'from':f'{s}-vocabulary-saturation.html'}),
      dict(type='prose',html='Unbound deposits are mostly &ldquo;hypothetical protein&rdquo; and deposits whose sequence was ambiguous between two features; several name proteins the module does call (for example 22K and E4 ORF1 in the one-fiber mastadenoviruses), so the unbound count is not a count of missed proteins.')]))
    F=dict(title=f"{m.replace('_',' ')} Coverage Audit",eyebrow=f"LowVan &middot; Adenoviridae &middot; {m} module &middot; {DATE}",
           h1=f"Does the {NAMES[m]} module work?",standfirst=STAND[m],headline_tiles=tiles,sections=S,
           source_note='Built with the LowVan module kit at <code>github.com/jimdavis1/LowVan_Skill</code> (modules/'+m+'). Every figure comes from a measurement file in the build; the known-issues log is in BUILD_NOTES.md.')
    json.dump(F,open(f'{W}/pages/facts/{s}_audit.json','w'),indent=1)
    notes=dict(standfirst=STAND[m],modules={m:SCOPE[m]},gaps='Features not modelled and known limits are listed in section 07 of the coverage audit and in BUILD_NOTES.md.')
    json.dump(notes,open(f'{W}/pages/facts/{s}_notes.json','w'),indent=1)
    print(m,'facts written')
