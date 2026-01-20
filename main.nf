#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-

==============================
IPAW: HiRIEF II varDB pipeline
==============================
@Authors
Jorrit Boekel @glormph
Yafeng Zhu @yafeng

https://github.com/lehtiolab/proteogenomics-analysis-workflow
*/

nf_required_version = '19.04.0'
if( ! nextflow.version.matches(">= ${nf_required_version}") ){
  println("Nextflow version too old, ${nf_required_version} required")
  exit(1)
}


mods = file(params.mods)
knownproteins = file(params.knownproteins)
blastdb = file(params.blastdb)
gtffile = file(params.gtf)
snpfa = file(params.snpfa)
dbsnp = params.dbsnp ? file(params.dbsnp) : false
cosmic = params.cosmic ? file(params.cosmic) : false
genomefa = file(params.genome)
tdb = file(params.tdb)
normalpsms = params.normalpsms ? file(params.normalpsms) : false
timstof = params.timstof ?: false
file_extension = params.input_format == 'timstof_d' ? '*.d' : '*.mzML'

// Search engine selection (default: msgf)
search_engine = params.search_engine ?: 'msgf'
if (!(search_engine in ['msgf', 'fragpipe'])) {
  println("Invalid search_engine: ${search_engine}. Must be 'msgf' or 'fragpipe'")
  exit(1)
}

// Database splits for FragPipe nonspecific search (default: 8)
db_splits = params.db_splits ?: 8

// Number of top PSMs to keep per spectrum (default: 1)
output_report_topN = params.output_report_topN ?: 1

/* PIPELINE START */

// Either feed an mzmldef file (tab separated lines with filepath\tsetname), or /path/to/files
if (!params.mzmldef && !params.input) {
  Channel
    .fromPath(params.mzmls)
    .map { it -> [it, 'NA', 'NA', 'NA'] }
    .set { mzml_raw }
} else {
  header = ['mzmlfile', 'setname', 'plate', 'fraction']
  mzmldef = params.mzmldef ?: params.input
  mzmllines = file(mzmldef).readLines().collect { it.tokenize('\t') }
  if (mzmllines[0] == header) {
    mzmllines.remove(0)
  }
  Channel
    .from(mzmllines)
    .map { it -> [file(it[0]), it[1], it[2] ? it[2] : 'NA', it[3] ? it[3].toInteger() : 'NA' ]}
    .set { mzml_raw }
}


// Route files based on search engine
// For FragPipe: ALL files bypass filterMzML (go to passthrough)
// For MSGF+: Only .mzML files go to filterMzML, .d files bypass
if (search_engine == 'fragpipe') {
  // FragPipe mode: send everything to passthrough, nothing to filter
  Channel.empty().set { mzml_to_filter }
  mzml_raw.set { mzml_to_passthrough }
} else {
  // MSGF+ mode: use choice to split by file extension
  mzml_to_filter = Channel.create()
  mzml_to_passthrough = Channel.create()
  mzml_raw.choice(mzml_to_filter, mzml_to_passthrough) { item ->
    item[0].name.endsWith('.d') ? 1 : 0
  }
}

process filterMzML {
  
  tag "${mzmlfile.baseName}"
  
  input:
  set file(mzmlfile), val(setname), val(plate), val(fraction) from mzml_to_filter
  
  output:
  set file("${mzmlfile.baseName}_filtered.mzML"), val(setname), val(plate), val(fraction) into mzml_filtered
  
  script:
  """
  ${workflow.projectDir}/bin/extract_spectra.py $mzmlfile ${mzmlfile.baseName}_filtered.mzML 10000000
  """
}

process passThroughRaw {
  
  tag "${infile.baseName}"
  
  input:
  set file(infile), val(setname), val(plate), val(fraction) from mzml_to_passthrough
  
  output:
  set file(infile), val(setname), val(plate), val(fraction) into mzml_passthrough
  
  script:
  """
  echo "Passing through ${infile.baseName} for FragPipe"
  """
}

// Merge filtered and passthrough files
// When filterMzML doesn't run (FragPipe mode), mzml_filtered is empty but still exists
// because Channel.empty() goes to mzml_to_filter, producing an empty mzml_filtered
mzml_filtered
  .mix(mzml_passthrough)
  .tap { sets; mzmlcounter }
  .tap { strips }
  .map { it -> [it[1], it[0].baseName.replaceFirst(/_filtered$/, ""), it[0], it[2], it[3]] }
  .into { mzmlfiles; groupset_mzmls; mzml_isobaric; mzml_presearch }

mzmlcounter
  .count()
  .subscribe { println "$it files in analysis" }
  .set { mzmlcount_psm }
sets
  .map{ it -> it[1] }
  .unique()
  .tap { sets_for_emtpybam; sets_for_denoms; sets_for_six }
  .collect()
  .subscribe { println "Detected setnames: ${it.join(', ')}" }

strips
  .map { it -> [it[1], it[2]] }
  .unique()
  .groupTuple()
  .set { strips_for_six }

if (params.pisepdb) {
  sets_for_six
    .toList()
    .map { it -> [it, normalpsms]}
    .set { normpsms }
} else {
  sets_for_six.set{ normpsms }
}

process splitSetNormalSearchPsms {
  //  normal search psm table, split on set col, collect files

  when: params.pisepdb
  input:
  set val(setnames), file('normalpsms') from normpsms
  output:
  set val(setnames), file({setnames.collect() { it + '.tsv' }}) into setnormpsms
  """
  msstitch split -i normalpsms --splitcol bioset
  """
}

setnormpsms
  .transpose()
  .join(strips_for_six)
  .set { setnormpsmtable } 

process splitPlateNormalSearchPsms {
  // create pep tables, split on plate, collect

  when: params.pisepdb
  input:
  set val(setname), file(normpsm), val(stripnames) from setnormpsmtable
  output:
  set val(setname), val(stripnames), file({stripnames.collect() { it + '.tsv' }}) into setplatepsms
  """
  msstitch split -i $normpsm --splitcol `python -c 'with open("$normpsm") as fp: h=next(fp).strip().split("\\t");print(h.index("Strip")+1)'`
  """
}

setplatepsms
  .transpose()
  .set { setplatepsmtable }

process normalSearchPsmsToPeptides {
  // create pep tables, split on plate, collect

  when: params.pisepdb
  input:
  set val(setname), val(strip), file(normpsm) from setplatepsmtable
  output:
  set val(setname), val(strip), file('peptides') into setplatepeptides
  """
  msstitch peptides -i $normpsm -o peptides --scorecolpattern area --spectracol 1 
  """
}

pipep = params.pisepdb ? Channel.fromPath(params.pisepdb) : Channel.empty()
varnov_peptides = params.pisepdb ? Channel.fromPath(params.pisepdb) : Channel.from(tdb)
setplatepeptides
  .combine(pipep)
  .set { sixftcreation_in }

process create6FTDB {
  // create 6FT DB per peptable-plate, collect fr 

  when: params.pisepdb

  input:
  set val(setname), val(stripname), file(peptides), file(pipeptides) from sixftcreation_in

  output:
  set val(setname), val(stripname), file('target_fr*.fasta') into t_splitdb

  script:
  strip = params.strips[stripname]
  """
  echo \'${groovy.json.JsonOutput.toJson(strip)}\' >> strip.json
  pi_database_splitter.py -i $pipeptides -p $peptides --stripdef strip.json --deltacolpattern Delta --fraccolpattern Fraction --fdrcolpattern '^q-value' --picutoff 0.2 --fdrcutoff 0.0 --maxlen $params.maxlen --minlen $params.minlen
  """
}


// channel match plate/fr/mzML
if (params.pisepdb) {
  t_splitdb
    .transpose()
    .map { it -> ["${it[0]}_${it[1]}_${it[2].baseName.replaceFirst(/.*_fr[0]*/, "")}", it[2]]}
    .set { db_w_id }
  mzml_presearch
    .map { it -> ["${it[0]}_${it[3]}_${it[4]}", it[0], it[1], it[2]] }  // add set_strip_fr identifier
    .into { mzml_dbid; mzml_dbfilter }
} else {
  Channel.from([['NA', tdb]]).set { db_w_id }
  mzml_presearch
    .map { it -> ["NA", it[0], it[1], it[2]] }
    .into { mzml_dbid; mzml_dbfilter }
}

process makeTargetSeqLookup {

  input:
  file(tdb) from varnov_peptides
  file(knownproteins)

  output:
  set file('mslookup_db.sqlite'), file('decoy_known.fa') into target_seq_lookup

  script:
  """
  # create text file of pi sep (much faster to import to SQLite than fasta)
  ${params.pisepdb ? "cut -f2 $tdb | sort -u > targetseq.txt" : "grep -v '^>' $tdb > targetseq.txt"}
  # Add trypsinized known proteins to txt file
  msstitch trypsinize -i $knownproteins -o knowntryp
  grep -v '^>' knowntryp >> targetseq.txt

  # TODO parametrize notrypsin?
  msstitch storeseq -i targetseq.txt --minlen $params.minlen ${params.pisepdb ? '--notrypsin': ''}
  msstitch makedecoy -i $knownproteins --dbfile mslookup_db.sqlite -o decoy_known.fa --scramble tryp_rev --minlen $params.minlen
  """
}

// Join DB with mzmls, use join filters out DBs without a match (in case of missing mzML fractions) 
// or duplicates (when having reruns, or when not running pI separared DBs in which case you only need 
// this process once). So we dont generate more decoys than necessary
db_w_id
  .join(mzml_dbfilter)
  .map { it -> [it[0], it[1]] }
  .combine(target_seq_lookup)
  .set { db_filtered }

process concatFasta {
  input:
  set val(dbid), file(db), file(targetlookup), file('decoy_known.fa') from db_filtered
  file knownproteins

  output:
  set val(dbid), file('td_concat.fa') into db_concatdecoy

  script:
  """
  # copy DB for faster access on network FS
  cp ${targetlookup} localdb.sql
  
  # Generate decoys from VarDB
  msstitch makedecoy -i $db --dbfile localdb.sql -o decoy_db.fa --scramble tryp_rev --minlen $params.minlen ${params.pisepdb ? '--notrypsin': ''}
  
  # Combine everything into one file
  cat $db $knownproteins decoy_db.fa decoy_known.fa > td_concat_all.fa
  
  # Deduplicate by ID (keep first occurrence) - applied to EVERYTHING
  awk '/^>/ {id=\$1; if (!seen[id]++) {keep=1; print} else {keep=0}} !/^>/ {if (keep) print}' \
    td_concat_all.fa > td_concat.fa
  
  # Report deduplication stats
  BEFORE=\$(grep -c "^>" td_concat_all.fa)
  AFTER=\$(grep -c "^>" td_concat.fa)
  echo "Deduplication: \$BEFORE -> \$AFTER sequences (removed \$((BEFORE - AFTER)) duplicates)"
  
  # Cleanup
  rm decoy_db.fa localdb.sql td_concat_all.fa
  """
}

// Now re-match the DB (now with decoy) with mzML, this is needed to fan out if more mzMLs than
// DBs have been used so we use the cross operator
db_concatdecoy
  .cross(mzml_dbid) // gives two Arrays so unfold them in next map step
  .map { it -> [it[0][0], it[0][1], it[1][1], it[1][2], it[1][3]] } // dbid, db, set, sample, file
  .set { mzml_search }

process prepareFilterDB {

  input:
  path knownproteins
  path snpfa

  output:
  file('knownprot.sqlite') into protseqdb
  file('snpprot.sqlite') into snpseqdb
  file('mslookup_db.sqlite') into trypseqdb

  """
  msstitch storeseq -i $knownproteins --minlen ${params.minlen} --fullprotein --minlen 7 -o knownprot.sqlite
  msstitch storeseq -i $snpfa --minlen ${params.minlen} -o snpprot.sqlite --fullprotein --minlen 7
  msstitch storeseq -i $knownproteins --insourcefrag --minlen ${params.minlen}
  """
}


process IsobaricQuant {

  when: !params.quantlookup && params.isobaric && search_engine != 'fragpipe'

  input:
  set val(setname), val(sample), file(infile), val(strip), val(fraction) from mzml_isobaric

  output:
  set val(sample), file("${infile}.consensusXML") into isobaricxml

  script:
  activationtype = [hcd:'High-energy collision-induced dissociation', cid:'Collision-induced dissociation', etd:'Electron transfer dissociation'][params.activation]
  isobtype = setisobaric && setisobaric[setname] ? setisobaric[setname] : false
  isobtype = isobtype == 'tmtpro' ? 'tmt16plex' : isobtype
  plextype = isobtype ? isobtype.replaceFirst(/[0-9]+plex/, "") : 'false'
  massshift = [tmt:0.0013, itraq:0.00125, false:0][plextype]
  """
  IsobaricAnalyzer  -type $isobtype -in $infile -out "${infile}.consensusXML" -extraction:select_activation "$activationtype" -extraction:reporter_mass_shift $massshift -extraction:min_precursor_intensity 1.0 -extraction:keep_unannotated_precursor true -quantification:isotope_correction true -t 10ppm -ti 0,1 
  """
}

isobaricxml
  .ifEmpty(['NA', 'NA'])
  .toList()
  .flatMap { it.sort({a, b -> a[0] <=> b[0]}) }
  .map { it -> it[1] }
  .collect()
  .set { sorted_isoxml }


mzmlfiles
  .toList()
  .map { it.sort( {a, b -> a[1] <=> b[1]}) }
  .map { it -> [it.collect() { it[0] }, it.collect() { it[2] }] }
  .set{ mzmlfiles_all }


process createNewSpectraLookup {

  publishDir "${params.outdir}", mode: 'copy', overwrite: true, saveAs: {it == 'mslookup_db.sqlite' ? 'quant_lookup.sql' : null }

  when: !params.quantlookup && search_engine != 'fragpipe'

  input:
  file(isobxmls) from sorted_isoxml 
  set val(setnames), file(mzmlfiles) from mzmlfiles_all
  
  output:
  file('mslookup_db.sqlite') into newspeclookup

  script:
  if(params.isobaric)
  """
  msstitch storespectra --spectra ${mzmlfiles.join(' ')} --setnames ${setnames.join(' ')}
  msstitch storequant --dbfile mslookup_db.sqlite --isobaric ${isobxmls.join(' ')} --spectra ${mzmlfiles.join(' ')}
  """
  else
  """
  msstitch storespectra --spectra ${mzmlfiles.join(' ')} --setnames ${setnames.join(' ')}
  """
}


// Set up spec_lookup channel based on mode
if (params.quantlookup) {
  Channel.fromPath(params.quantlookup).set { spec_lookup }
} else if (search_engine == 'fragpipe') {
  // FragPipe/timsTOF mode: create dummy file channel (createPSMTable handles this case)
  Channel.from(file('NO_LOOKUP')).set { spec_lookup }
} else {
  newspeclookup.set { spec_lookup }
} 


// ============================================================================
// SEARCH ENGINE PROCESSES
// ============================================================================

if (search_engine == 'msgf') {
  mzml_search.set { mzml_msgf }
  Channel.empty().set { mzml_fragpipe }
} else {
  mzml_search.set { mzml_fragpipe }
  Channel.empty().set { mzml_msgf }
}

process msgfPlus {

  when: search_engine == 'msgf'

  input:
  set val(setfr_id), file(db), val(setname), val(sample), file(mzml) from mzml_msgf
  file mods
  
  output:
  set val(setname), val(sample), file("${sample}.mzid") into mzids_msgf
  set val(setname), file("${sample}.mzid"), file('out.mzid.tsv') into mzidtsvs_msgf
  
  script:
  msgfprotocol = 0
  """
  msgf_plus -Xmx${task.memory.toMega()}M -d $db -s $mzml -o "${sample}.mzid" -thread ${task.cpus * params.threadspercore} -mod $mods -tda 0 -t 15ppm -ti -1,2 -m 0 -inst 3 -e 0 -protocol ${msgfprotocol} -ntt 2 -minLength 8 -maxLength 15 -minCharge 2 -maxCharge 4 -n 1 -addFeatures 1
  msgf_plus -Xmx${task.memory.toMega()}M edu.ucsd.msjava.ui.MzIDToTsv -i "${sample}.mzid" -o out.mzid.tsv
  rm -f td_concat.c*
  """
}

process fragPipeSearch {

  when: search_engine == 'fragpipe'
  
  conda "/home/il364/.conda/envs/fragpipe"

  input:
  set val(setfr_id), path(db), val(setname), val(sample), path(mzml) from mzml_fragpipe
  path mods
  
  output:
  set val(setname), val(sample), file("${sample}.pepXML") into pepxml_files
  set val(setname), val(sample), file("${sample}.pin") into pin_files
  
  script:
  """
  # Parse MSGF+ mods file to MSFragger format
  parse_msgf_mods_to_fragger.py $mods fragger_mods.txt
  
  # Split database into ${db_splits} parts for nonspecific search
  mkdir -p split_dbs
  
  # Count proteins and calculate split size
  TOTAL_PROTS=\$(grep -c "^>" $db)
  SPLIT_SIZE=\$(( (TOTAL_PROTS + ${db_splits} - 1) / ${db_splits} ))
  
  # Split FASTA file
  awk -v size="\$SPLIT_SIZE" -v prefix="split_dbs/part_" '
    BEGIN { filenum = 1; count = 0 }
    /^>/ { 
      if (count >= size) { filenum++; count = 0 }
      count++
    }
    { print > (prefix filenum ".fa") }
  ' $db
  
  echo "Split database into \$(ls split_dbs/*.fa | wc -l) parts"
  
  # Create base MSFragger parameter file
  cat > fragger_base.params <<EOF
num_threads = ${task.cpus * params.threadspercore}
decoy_prefix = decoy_

precursor_mass_lower = -20.0
precursor_mass_upper = 20.0
precursor_mass_units = 1
precursor_true_tolerance = 20.0
precursor_true_units = 1
fragment_mass_tolerance = 20.0
fragment_mass_units = 1
calibrate_mass = 2

isotope_error = 0
delta_mass_exclude_ranges = (-1.5,3.5)
fragment_ion_series = b,y

search_enzyme_name = nonspecific
search_enzyme_sense_1 = C
search_enzyme_cut_1 = @
search_enzyme_nocut_1 = 
allowed_missed_cleavage_1 = ${params.maxmiscleav}
search_enzyme_sense_2 = C
search_enzyme_cut_2 = @
search_enzyme_nocut_2 = 
allowed_missed_cleavage_2 = 2
num_enzyme_termini = 0

clip_nTerm_M = 1

output_format = pepxml_pin
output_report_topN = 1
output_max_expect = 50.0

precursor_charge = 2 3
digest_min_length = ${params.minlen}
digest_max_length = ${params.maxlen}
digest_mass_range = 300.0 5000.0

deisotope = 1
deneutralloss = 1
minimum_peaks = 15
use_topN_peaks = 300
min_matched_fragments = 4
intensity_transform = 1
activation_types = all
analyzer_types = all
remove_precursor_peak = 1
remove_precursor_range = -1.5,1.5

\$(cat fragger_mods.txt)
EOF

  # Run MSFragger on each database split
  mkdir -p results
  for split_db in split_dbs/*.fa; do
    split_name=\$(basename \$split_db .fa)
    echo "Searching against \$split_name..."
    
    # Create params file with this database
    cp fragger_base.params fragger_\${split_name}.params
    echo "database_name = \$split_db" >> fragger_\${split_name}.params
    
    # Run MSFragger
    java -Xmx${task.memory.toMega()}M -jar ${params.msfragger_jar} fragger_\${split_name}.params $mzml
    
    # Get output basename
    BASENAME=\$(basename $mzml .mzML)
    BASENAME=\$(basename \$BASENAME .d)
    
    # Move outputs to results folder with split name
    mv \${BASENAME}.pepXML results/\${split_name}.pepXML
    mv \${BASENAME}.pin results/\${split_name}.pin
  done
  
import glob
import csv
import collections
import os

pin_files = sorted(glob.glob('results/*.pin'))
print(f"Merging {len(pin_files)} PIN files from database splits")

# Read header
with open(pin_files[0]) as f:
    header = f.readline().strip().split('\t')

specId_idx = header.index('SpecId')
hyperscore_idx = header.index('hyperscore') if 'hyperscore' in header else None
rank_idx = header.index('rank') if 'rank' in header else None

# Group PSMs by spectrum
spectrum_psms = collections.defaultdict(list)
csv.field_size_limit(2**31 - 1)

for pin_file in pin_files:
    with open(pin_file) as f:
        reader = csv.reader(f, delimiter='\t')
        next(reader)  # skip header
        for row in reader:
            if len(row) < len(header):
                continue
            spec_id = row[specId_idx]
            # Remove rank suffix to get spectrum key
            spec_key = spec_id.rsplit('_', 1)[0] if '_' in spec_id else spec_id.rsplit('.', 1)[0]
            
            # Get hyperscore for sorting
            score = float(row[hyperscore_idx]) if hyperscore_idx else 0
            spectrum_psms[spec_key].append((score, row))

print(f"Found {len(spectrum_psms)} unique spectra")

# Write merged PIN with only top hit per spectrum
total_written = 0
with open('${sample}.pin', 'w') as out:
    out.write('\t'.join(header) + '\n')
    for spec_key in sorted(spectrum_psms.keys()):
        psms = spectrum_psms[spec_key]
        psms.sort(key=lambda x: x[0], reverse=True)  # Sort by hyperscore desc
        
        # Keep only top hit
        score, row = psms[0]
        if rank_idx is not None:
            row[rank_idx] = '1'
        row[specId_idx] = f"{spec_key}_1"
        out.write('\t'.join(row) + '\n')
        total_written += 1

print(f"Wrote {total_written} PSMs (top hit per spectrum)")
print(f"PIN file size: {os.path.getsize('${sample}.pin') / 1e6:.1f} MB")
MERGE_PIN_EOF
  
  # Merge pepXML files using shell
  # Get header from first file (everything up to first spectrum_query)
  FIRST_PEPXML=\$(ls results/*.pepXML | head -1)
  
  # Extract header (everything before <spectrum_query)
  awk '/<spectrum_query/{exit} {print}' \$FIRST_PEPXML > ${sample}.pepXML
  
  # Extract and append all spectrum_query elements from all files
  for f in results/*.pepXML; do
    # Extract spectrum_query blocks
    awk '/<spectrum_query/,/<\\/spectrum_query>/' "\$f" >> ${sample}.pepXML
  done
  
  # Add closing tags
  echo '</msms_run_summary>' >> ${sample}.pepXML
  echo '</msms_pipeline_analysis>' >> ${sample}.pepXML
  
  echo "Merged pepXML complete"
  """
}

// For FragPipe: create TSV directly from pepXML (skip mzid conversion)
// PIN files go directly to percolator
pepxml_files.set { pepxml_for_tsv }

// Route PIN files to percolator
pin_files
  .groupTuple()
  .set { pin_for_percolator }

process pepXMLToTSV {

  when: search_engine == 'fragpipe'
  
  input:
  set val(setname), val(sample), path(pepxml) from pepxml_for_tsv
  
  output:
  set val(setname), val(sample), file("${sample}.tsv") into tsvs_fragpipe
  
  script:
  """
  #!/usr/bin/env python3
  import xml.etree.ElementTree as ET
  import csv
  import re
  
  # Parse pepXML and create TSV
  tree = ET.parse('${pepxml}')
  root = tree.getroot()
  
  # Handle namespace
  ns_match = re.match(r'\\{(.*)\\}', root.tag)
  ns = {'pep': ns_match.group(1)} if ns_match else {}
  
  rows = []
  
  # Find spectrum queries
  if ns:
      queries = root.findall('.//pep:spectrum_query', ns)
  else:
      queries = root.findall('.//spectrum_query')
  
  for sq in queries:
      spectrum = sq.get('spectrum', '')
      start_scan = sq.get('start_scan', '')
      precursor_mass = sq.get('precursor_neutral_mass', '')
      charge = sq.get('assumed_charge', '')
      
      # Get search hits
      if ns:
          hits = sq.findall('.//pep:search_hit', ns)
      else:
          hits = sq.findall('.//search_hit')
      
      for hit in hits:
          peptide = hit.get('peptide', '')
          protein = hit.get('protein', '')
          calc_mass = hit.get('calc_neutral_pep_mass', '')
          num_matched = hit.get('num_matched_ions', '')
          tot_ions = hit.get('tot_num_ions', '')
          
          # Get scores
          scores = {}
          if ns:
              score_elems = hit.findall('.//pep:search_score', ns)
          else:
              score_elems = hit.findall('.//search_score')
          for s in score_elems:
              scores[s.get('name', '')] = s.get('value', '')
          
          rows.append({
              '#SpecFile': '${sample}',
              'SpecID': spectrum,
              'ScanNum': start_scan,
              'Peptide': peptide,
              'Protein': protein,
              'Charge': charge,
              'PrecursorMass': precursor_mass,
              'CalcMass': calc_mass,
              'MatchedIons': num_matched,
              'TotalIons': tot_ions,
              'Hyperscore': scores.get('hyperscore', ''),
              'Expect': scores.get('expect', ''),
              'Nextscore': scores.get('nextscore', ''),
          })
  
  # Write TSV
  if rows:
      with open('${sample}.tsv', 'w', newline='') as f:
          writer = csv.DictWriter(f, fieldnames=rows[0].keys(), delimiter='\\t')
          writer.writeheader()
          writer.writerows(rows)
      print(f"Created TSV with {len(rows)} PSMs")
  else:
      # Create empty file with header
      with open('${sample}.tsv', 'w') as f:
          f.write('#SpecFile\\tSpecID\\tScanNum\\tPeptide\\tProtein\\tCharge\\n')
      print("Warning: No PSMs found in pepXML")
  """
}

// Group TSVs by setname for downstream
tsvs_fragpipe
  .map { setname, sample, tsv -> [setname, tsv] }
  .groupTuple()
  .set { grouped_tsvs_fragpipe }

// For FragPipe: create mzidtsvs-compatible channel (TSV only, no mzid needed for timsTOF path)
// The svmToTSV_timsTOF process will be replaced with a FragPipe-specific version
if (search_engine == 'fragpipe') {
  // Create dummy structure to match expected format [setname, mzid_files, tsv_files]
  grouped_tsvs_fragpipe
    .map { setname, tsvs -> [setname, tsvs, tsvs] }  // Use TSVs as placeholder for mzids
    .set { mzidtsvs }
} else {
  mzidtsvs_msgf.set { mzidtsvs }
}

// ============================================================================
// PERCOLATOR SECTION - Separate paths for MSGF+ and FragPipe
// ============================================================================

// MSGF+ path: mzid files work with msgf2pin
mzids_msgf
  .groupTuple()
  .set { mzids_2pin }

process percolator_msgf {
  when: search_engine == 'msgf'
  
  input:
  set val(setname), val(samples), file('mzid?') from mzids_2pin
  
  output:
  set val(setname), file('perco.xml') into percolated_msgf
  
  script:
  """
  mkdir mzids
  count=1;for sam in ${samples.join(' ')}; do ln -s `pwd`/mzid\$count mzids/\${sam}.mzid; echo mzids/\${sam}.mzid >> metafile; ((count++));done
  msgf2pin -o percoin.xml -P "decoy_" metafile
  percolator -j percoin.xml -X perco.xml -N 500000 --decoy-xml-output -y --subset-max-train 300000
  """
}

// FragPipe path: use PIN files directly from MSFragger (output_format = pepxml_pin)
process percolator_fragpipe {
  when: search_engine == 'fragpipe'
  
  input:
  set val(setname), val(samples), file('*.pin') from pin_for_percolator
  
  output:
  set val(setname), file('perco.xml') into percolated_fragpipe
  
  script:
  def topN = params.output_report_topN ?: 1
  """
  #!/usr/bin/env python3
  import glob
  import os
  import csv
  import collections
  
  pin_files = sorted(glob.glob('*.pin'))
  print(f"Processing {len(pin_files)} PIN files", flush=True)
  
  # Read header from first file
  with open(pin_files[0]) as f:
      header = f.readline().strip().split('\\t')
  
  # Find important column indices
  specId_idx = header.index('SpecId')
  hyperscore_idx = header.index('hyperscore') if 'hyperscore' in header else None
  log10_evalue_idx = header.index('log10_evalue') if 'log10_evalue' in header else None
  rank_idx = header.index('rank') if 'rank' in header else None
  
  # Determine score column for sorting (prefer hyperscore, fallback to log10_evalue)
  if hyperscore_idx is not None:
      score_idx = hyperscore_idx
      higher_is_better = True
  elif log10_evalue_idx is not None:
      score_idx = log10_evalue_idx
      higher_is_better = False  # lower log10_evalue = better
  else:
      # Fallback: find any score-like column
      score_idx = None
      for i, col in enumerate(header):
          if 'score' in col.lower() and i > 2:
              score_idx = i
              higher_is_better = True
              break
  
  print(f"Using column {header[score_idx] if score_idx else 'none'} for scoring", flush=True)
  
  # Group all PSMs by spectrum (without rank suffix)
  # SpecId format: filename.scan.scan.charge_rank -> group by filename.scan.scan.charge
  spectrum_psms = collections.defaultdict(list)
  total_psms = 0
  
  csv.field_size_limit(2**31 - 1)
  
  for pin_file in pin_files:
      with open(pin_file) as f:
          reader = csv.reader(f, delimiter='\\t')
          next(reader)  # skip header
          for row in reader:
              if len(row) < len(header):
                  continue
              # Extract spectrum key (remove rank suffix if present)
              spec_id = row[specId_idx]
              # Handle format: name.scan.scan.charge_rank or name_rank
              if '_' in spec_id:
                  spec_key = spec_id.rsplit('_', 1)[0]
              else:
                  spec_key = spec_id.rsplit('.', 1)[0]
              
              # Get score for sorting
              if score_idx is not None:
                  try:
                      score = float(row[score_idx])
                  except (ValueError, IndexError):
                      score = 0
              else:
                  score = 0
              
              spectrum_psms[spec_key].append((score, row))
              total_psms += 1
      
      print(f"  Read {pin_file}: {total_psms} total PSMs so far", flush=True)
  
  print(f"Total PSMs: {total_psms}", flush=True)
  print(f"Unique spectra: {len(spectrum_psms)}", flush=True)
  
  # Sort each spectrum's PSMs and keep top N
  topN = ${topN}
  kept_psms = 0
  
  with open('merged.pin', 'w') as out:
      out.write('\\t'.join(header) + '\\n')
      
      for spec_key in sorted(spectrum_psms.keys()):
          psms = spectrum_psms[spec_key]
          # Sort by score (higher hyperscore = better, lower log10_evalue = better)
          if hyperscore_idx is not None:
              psms.sort(key=lambda x: x[0], reverse=True)  # higher is better
          else:
              psms.sort(key=lambda x: x[0], reverse=False)  # lower is better
          
          # Keep top N and update ranks
          for new_rank, (score, row) in enumerate(psms[:topN], 1):
              # Update rank column if it exists
              if rank_idx is not None:
                  row[rank_idx] = str(new_rank)
              # Update SpecId to include new rank
              row[specId_idx] = f"{spec_key}_{new_rank}"
              out.write('\\t'.join(row) + '\\n')
              kept_psms += 1
  
  print(f"Kept {kept_psms} PSMs ({100*kept_psms/max(1,total_psms):.1f}%)", flush=True)
  print(f"Merged PIN file size: {os.path.getsize('merged.pin') / 1e9:.2f} GB", flush=True)
  
  # Run Percolator
  import subprocess
  print("Running Percolator...", flush=True)
  subprocess.run([
      'percolator', '-j', 'merged.pin', '-X', 'perco.xml',
      '-N', '500000', '--decoy-xml-output', '-y', '--subset-max-train', '300000'
  ], check=True)
  print("Percolator complete", flush=True)
  """
}

// Merge percolator outputs from both search engines
percolated_msgf
  .mix(percolated_fragpipe)
  .set { percolated }


percolated
  .tap { var_percolated }
  .set { nov_percolated }


process getVariantPercolator {

  when: params.varheaders

  input:
  set val(setname), file(x) from var_percolated

  output:
  set val(setname), val('variant'), file("${x}_h0.xml") into var_perco

  script:
  """
  msstitch splitperco -i $x --protheaders "known:${params.knownheaders}|novel:${params.varheaders}"
  """
}


process getNovelPercolator {

  when: params.novheaders

  input:
  set val(setname), file(x) from nov_percolated

  output:
  set val(setname), val('novel'), file("${x}_h0.xml") into nov_perco

  script:
  """
  msstitch splitperco -i $x --protheaders "known:${params.knownheaders}|novel:${params.novheaders}"
  """
}


nov_perco
  .concat(var_perco)
  .set { splitperco }


process filterPercolator {

  input:
  set val(setname), val(peptype), file(perco) from splitperco
  file 'trypseqdb' from trypseqdb
  file 'protseqdb' from protseqdb
  file knownproteins

  output:
  set val(setname), val(peptype), file('filtprot') into filtered_perco
  
  script:
  if (params.noclassfdr)
  """
  mv $perco filtprot
  """
  else
  """
  msstitch filterperco -i $perco -o filtprot --fullprotein --dbfile protseqdb --minlen $params.minlen --deamidate
  """
}

nov_filtered_perco = Channel.create()
var_filtered_perco = Channel.create()
filtered_perco
  .choice( var_filtered_perco, nov_filtered_perco) { it -> it[1] == 'variant' ? 0 : 1 }

mzidtsvs
  .groupTuple()
  .tap { variantmzidtsv }
  .join(nov_filtered_perco)
  .set { nov_mzperco }

variantmzidtsv
  .join(var_filtered_perco)
  .concat(nov_mzperco)
  .set { allmzperco }
  
if (params.timstof) {
  allmzperco.set { allmzperco_timstof }
  Channel.empty().set { allmzperco_standard }
} else {
  allmzperco.set { allmzperco_standard }
  Channel.empty().set { allmzperco_timstof }
}

process svmToTSV {

  when: !params.timstof

  input:
  set val(setname), file('mzid?'), file('tsv?'), val(peptype), file(perco) from allmzperco_standard

  output:
  set val(setname), val(peptype), file('target.tsv') into mzidtsv_perco_standard

  script:
  """
  tsvs=""
  mzids=""
  count=1; for tsvfn in \$(ls tsv*)
    do 
    tsvs="\${tsvs} tsv\${count}"
    mzids="\${mzids} mzid\${count}"
    ((count++))
    done
  mkdir outtables
  msstitch perco2psm --perco $perco -d outtables -i \$tsvs --mzids \$mzids
  msstitch concat -i outtables/* -o psms
  msstitch split -i psms --splitcol TD
  """
}

process svmToTSV_timsTOF {

  when: params.timstof

  input:
  set val(setname), file('mzid?'), file('tsv?'), val(peptype), file(perco) from allmzperco_timstof

  output:
  set val(setname), val(peptype), file('target.tsv') into mzidtsv_perco_timstof

  script:
  """
  #!/usr/bin/env python3
  import xml.etree.ElementTree as ET
  import glob
  import csv
  import re
  import os

  search_engine = '${search_engine}'
  
  # Parse percolator XML to get PSM info
  perco_psms = {}
  tree = ET.parse('${perco}')
  root = tree.getroot()
  ns = {'p': 'http://per-colator.com/percolator_out/15'}
  
  for psm in root.findall('.//p:psm', ns):
      psm_id = psm.get('{http://per-colator.com/percolator_out/15}psm_id')
      decoy = psm.get('{http://per-colator.com/percolator_out/15}decoy')
      if decoy == 'true':
          continue
      
      q_value = psm.find('p:q_value', ns).text
      pep = psm.find('p:pep', ns).text
      svm_score = psm.find('p:svm_score', ns).text
      peptide_elem = psm.find('p:peptide_seq', ns)
      peptide = peptide_elem.get('seq') if peptide_elem is not None else ''
      
      # Handle both MSGF+ and FragPipe PSM ID formats
      if search_engine == 'fragpipe':
          # FragPipe PIN format: spectrum_charge_rank or similar
          # Try to extract peptide and charge from PSM ID
          parts = psm_id.split('_')
          if len(parts) >= 2:
              charge = parts[-2] if parts[-2].isdigit() else parts[-1]
              key = (peptide, charge)
              perco_psms[key] = {
                  'q_value': q_value,
                  'pep': pep,
                  'svm_score': svm_score,
                  'psm_id': psm_id
              }
      else:
          # MSGF+ format: filename_SII_specidx_rank_scannum_charge_1
          parts = psm_id.rsplit('_SII_', 1)
          if len(parts) == 2:
              filename = parts[0]
              scan_parts = parts[1].split('_')
              if len(scan_parts) >= 4:
                  spec_idx = scan_parts[0]
                  scan_num = scan_parts[2]
                  charge = scan_parts[3]
                  key = (filename, peptide, charge)
                  perco_psms[key] = {
                      'q_value': q_value,
                      'pep': pep,
                      'svm_score': svm_score,
                      'psm_id': psm_id
                  }

  # Read all TSV files and match
  tsv_files = sorted(glob.glob('tsv*'))
  
  # For MSGF+ mode, also parse mzid files
  mzid_mapping = {}
  if search_engine != 'fragpipe':
      mzid_files = sorted(glob.glob('mzid*'))
      for mzid_file in mzid_files:
          try:
              tree = ET.parse(mzid_file)
              root = tree.getroot()
              mzid_ns = {'mzid': 'http://psidev.info/psi/pi/mzIdentML/1.1'}
              
              for spec_data in root.findall('.//mzid:SpectraData', mzid_ns):
                  location = spec_data.get('location', '')
                  mzid_filename = location.split('/')[-1].replace('.mzML', '').replace('.d', '').replace('_filtered', '')
              
              peptide_sequences = {}
              for pep_elem in root.findall('.//mzid:Peptide', mzid_ns):
                  pep_id = pep_elem.get('id')
                  pep_seq = pep_elem.find('mzid:PeptideSequence', mzid_ns)
                  if pep_seq is not None and pep_id:
                      peptide_sequences[pep_id] = pep_seq.text
              
              for spec_result in root.findall('.//mzid:SpectrumIdentificationResult', mzid_ns):
                  spec_id = spec_result.get('spectrumID')
                  for spec_item in spec_result.findall('mzid:SpectrumIdentificationItem', mzid_ns):
                      charge = spec_item.get('chargeState')
                      peptide_ref = spec_item.get('peptide_ref')
                      seq = peptide_sequences.get(peptide_ref)
                      if seq:
                          key = (mzid_filename, seq, charge)
                          mzid_mapping[key] = spec_id
          except Exception as e:
              print(f"Warning: Could not parse mzid file {mzid_file}: {e}")

  # Process TSV files
  all_rows = []
  header = None
  
  for tsv_file in tsv_files:
      with open(tsv_file, 'r') as f:
          reader = csv.DictReader(f, delimiter='\\t')
          if header is None:
              header = list(reader.fieldnames) + ['percolator svm-score', 'PSM q-value', 'PSM PEP', 'peptide q-value', 'peptide PEP', 'TD']
          
          for row in reader:
              spec_file = row.get('#SpecFile', row.get('SpecFile', ''))
              filename = spec_file.split('/')[-1].replace('.mzML', '').replace('.d', '').replace('_filtered', '')
              peptide = row.get('Peptide', '')
              charge = row.get('Charge', '')
              
              # Match based on search engine
              if search_engine == 'fragpipe':
                  key = (peptide, charge)
              else:
                  key = (filename, peptide, charge)
              
              if key in perco_psms:
                  perco_data = perco_psms[key]
                  row['percolator svm-score'] = perco_data['svm_score']
                  row['PSM q-value'] = perco_data['q_value']
                  row['PSM PEP'] = perco_data['pep']
                  row['peptide q-value'] = perco_data['q_value']
                  row['peptide PEP'] = perco_data['pep']
                  row['TD'] = 'target'
                  all_rows.append(row)

  # Write output
  print(f"Matched {len(all_rows)} PSMs with percolator results")
  if all_rows:
      with open('target.tsv', 'w', newline='') as f:
          writer = csv.DictWriter(f, fieldnames=header, delimiter='\\t', extrasaction='ignore')
          writer.writeheader()
          writer.writerows(all_rows)
  else:
      # Write empty file with header
      with open('target.tsv', 'w') as f:
          f.write('\\t'.join(header) + '\\n')
  
  print(f"Matched {len(all_rows)} PSMs from percolator to TSV")
  """
}

// Merge the two channels
mzidtsv_perco_standard
  .mix(mzidtsv_perco_timstof)
  .set { mzidtsv_perco }

mzidtsv_perco
  .combine(spec_lookup)
  .set { prepsm }

process createPSMTable {

  input:
  set val(setname), val(peptype), file('psms'), file('lookup') from prepsm
  val(mzmlcount) from mzmlcount_psm

  output:
  set val(setname), val(peptype), file("${setname}_${peptype}_psmtable.txt") into psmtable

  script:
  // For FragPipe/timsTOF: skip spectra lookup (--dbfile), just filter and create basic table
  if (search_engine == 'fragpipe')
  """
  msstitch conffilt -i psms -o filtpsm --confidence-better lower --confidence-lvl 0.01 --confcolpattern 'PSM q-value'
  msstitch conffilt -i filtpsm -o filtpep --confidence-better lower --confidence-lvl 0.01 --confcolpattern 'peptide q-value'
  # Create PSM table without spectra lookup for timsTOF
  cp filtpep ${setname}_${peptype}_psmtable.txt
  sed 's/\\#SpecFile/SpectraFile/' -i ${setname}_${peptype}_psmtable.txt
  """
  else
  """
  msstitch conffilt -i psms -o filtpsm --confidence-better lower --confidence-lvl 0.01 --confcolpattern 'PSM q-value'
  msstitch conffilt -i filtpsm -o filtpep --confidence-better lower --confidence-lvl 0.01 --confcolpattern 'peptide q-value'
  cp lookup psmlookup
  msstitch psmtable -i filtpep --dbfile psmlookup --addbioset -o ${setname}_${peptype}_psmtable.txt ${params.isobaric ? '--isobaric': ''}
  sed 's/\\#SpecFile/SpectraFile/' -i ${setname}_${peptype}_psmtable.txt
  """
}


variantpsms = Channel.create()
novelpsms = Channel.create()
psmtable
  .tap { setmergepsmtables; peppsms }
  .choice( variantpsms, novelpsms ) { it -> it[1] == 'variant' ? 0 : 1 }

 
setmergepsmtables
  .groupTuple(by: 1)
  .set { psmmerge_in }


process mergeSetPSMtable {
  publishDir "${params.outdir}", mode: 'copy', overwrite: true

  input:
  set val(setnames), val(peptype), file(psmtables) from psmmerge_in

  output:
  file "${peptype}_psmtable.txt" into produced_psmtables

  """
  head -n1 ${psmtables[0]} > ${peptype}_psmtable.txt
  for fn in ${psmtables.join(' ')}; do tail -n+2 \$fn >> ${peptype}_psmtable.txt; done
  """
}


process prePeptideTable {

  input:
  set val(setname), val(peptype), file('psms.txt') from peppsms 

  output:
  set val(setname), val(peptype), file('peptidetable.txt') into peptable

  script:
  """
  msstitch peptides -i psms.txt -o peptidetable.txt --scorecolpattern svm --spectracol 1 \
    ${setisobaric && setisobaric[setname] ? "--isobquantcolpattern plex --minint 0.1 --logisoquant --denompatterns ${setdenoms[setname].join(' ')}" : ''}
  """
}

novelpsms
  .into{novelpsmsFastaBedGFF; novelpsms_specai}

novelprepep = Channel.create()
presai_peptable = Channel.create()
peptable
  .choice( presai_peptable, novelprepep ) { it -> it[1] == 'variant' ? 0 : 1 }
novelprepep
  .join(novelpsmsFastaBedGFF)
  .set { novelFaBdGfPep }

process createFastaBedGFF {
 
  publishDir "${params.outdir}", mode: 'copy', overwrite: true, saveAs: { it == "${setname}_novel_peptides.gff3" ? "${setname}_novel_peptides.gff3" : null}
 
  input:
  set val(setname), val(peptype), file(peptides) , val(psmtype), file(psms) from novelFaBdGfPep
  file gtffile
  file tdb from tdb
 
  output:
  set val(setname), file('novel_peptides.fa') into novelfasta
  set val(setname), file('novel_peptides.bed') into novelbed
  set val(setname), file("${setname}_novel_peptides.gff3") into novelGFF3
  set val(setname), file('novel_peptides.tab.txt') into novelpep
  set val(setname), file('novpep_perco_quant.txt') into novelpep_percoquant
 
  """
  map_novelpeptide2genome.py --input $psms --gtf $gtffile --fastadb $tdb --tab_out novel_peptides.tab.txt --fasta_out novel_peptides.fa --gff3_out ${setname}_novel_peptides.gff3 --bed_out novel_peptides.bed
  sort -k 1b,1 <(tail -n+2 $peptides) |cut -f 1,14-500 > peptable_sorted
  sort -k 2b,2 <(tail -n+2 novel_peptides.tab.txt) > novpep_sorted
  paste <(cut -f 2 novpep_sorted) <(cut -f1,3-500 novpep_sorted) > novpep_pepcols
  join novpep_pepcols peptable_sorted -j 1 -a1 -o auto -e 'NA' -t \$'\\t' > novpep_pqjoin
  # Cut only bare peptide col and q-values/isoquant
  paste <(cut -f 2 novpep_pqjoin) <(cut -f8-500 novpep_pqjoin) > novpep_joined_pepcols
  paste <(head -n1 novel_peptides.tab.txt | cut -f1)  <(cut -f 14-500 $peptides |head -n1) > header
  cat header novpep_joined_pepcols > novpep_perco_quant.txt
  """
}

novelpep
  .into {blastnovelpep; blatnovelpep; annonovelpep; snpnovelpep}
novelfasta
  .into {blastnovelfasta; blatnovelfasta}

process BlastPNovel {

  input:
  set val(setname), file(novelfasta) from blastnovelfasta
  file blastdb

  output:
  set val(setname), file('blastp_out.txt') into novelblast
  
  """
  makeblastdb -in $blastdb -dbtype prot
  blastp -db $blastdb -query $novelfasta -outfmt '6 qseqid sseqid pident qlen slen qstart qend sstart send mismatch positive gapopen gaps qseq sseq evalue bitscore' -num_threads 4 -max_target_seqs 1 -evalue 1000 -out blastp_out.txt
  """
}

novelpsms_specai
  .map { it -> [it[0], it[2]] }
  .join(blastnovelpep)
  .join(novelblast)
  .set { novelblastout }

process ParseBlastpOut {
 
 input:
 set val(setname), file(psms), file(novelpep), file(novelblast) from novelblastout
 file blastdb

 output:
 set val(setname), file('peptable_blastp.txt') into peptable_blastp
 set val(setname), file('single_mismatch_novpeps.txt') into novpeps_singlemis

 """
 parse_BLASTP_out.py --input $novelpep --blastp_result $novelblast --fasta $blastdb --output peptable_blastp.txt
 extract_1mismatch_novpsm.py peptable_blastp.txt $psms single_mismatch_novpeps.txt
 """

}

groupset_mzmls
  .map { it -> [it[0], it[1], it[2]] } // strip fraction, strip
  .groupTuple()
  .tap { var_specaimzmls }
  .join(novpeps_singlemis)
  .set { grouped_saavnov_mzml_psms }
  

process ValidateSingleMismatchNovpeps {
  
  publishDir "${params.outdir}", mode: 'copy', overwrite: true, saveAs: { it == "precursorError.histogram.plot.pdf" ? "${setname}_novel_precursorError_plot.pdf" : it }
  
  input:
  set val(setname), val(samples), file(mzmls), file(psms) from grouped_saavnov_mzml_psms

  output:
  set val(setname), file("${setname}_novel_saav_specai.txt") into singlemis_specai
  file 'precursorError.histogram.plot.pdf' optional true into novel_specai_plot

  """
  mkdir mzmls
  for fn in $mzmls; do ln -s `pwd`/\$fn mzmls/; done
  Rscript /SpectrumAI/SpectrumAI.R mzmls $psms ${setname}_novel_saav_specai.txt || cp $psms singlemis_specai.txt
  """
}

singlemis_specai
  .join(peptable_blastp)
  .set { nov_specaiparse }

process novpepSpecAIOutParse {

  input:
  set val(setname), file(x), file('peptide_table.txt') from nov_specaiparse
  
  output:
  set val(setname), file('novpep_specai.txt') into novpep_singlemisspecai

  """
  parse_spectrumAI_out.py --spectrumAI_out $x --input peptide_table.txt --output novpep_sa
  cut -f 1,8-19 novpep_sa > novpep_specai.txt
  """
}


process BLATNovel {

  input:
  set val(setname), file(novelfasta) from blatnovelfasta
  file genomefa

  output:
  set val(setname), file('blat_out.pslx') into novelblat

  """
  blat $genomefa $novelfasta -t=dnax -q=prot -tileSize=5 -minIdentity=99 -out=pslx blat_out.pslx 
  """
}

novelblat
  .join(blatnovelpep)
  .set { novblatparse }

process parseBLATout {

 input:
 set val(setname), file(novelblat), file(novelpep) from novblatparse

 output:
 set val(setname), file('peptable_blat.txt') into peptable_blat

 """
 parse_BLAT_out.py $novelblat $novelpep peptable_blat.txt

 """
}

process labelnsSNP {
  
  input:
  set val(setname), file(peptable) from snpnovelpep
  file snpfa
  file(snpdb) from snpseqdb

  output:
  set val(setname), file('nssnp.txt') into ns_snp_out

  """
  label_nsSNP_pep.py --input $peptable --nsSNPdb $snpfa --dbfile "$snpdb" --output nssnp.txt --minlen $params.minlen
  """
}

bwfile = Channel.fromPath(params.bigwigs)
novelGFF3
  .combine(bwfile)
  .into { novelGFF3_phast; novelGFF3_phylo; novelGFF3_bams }

process phastcons {
  
  input:
  set val(setname), file(novelgff), file('bigwigs') from novelGFF3_phast

  output:
  set val(setname), file ('phastcons.txt') into phastcons_out

  """
  calculate_phastcons.py $novelgff bigwigs/hg19.100way.phastCons.bw phastcons.txt
  """
}

process phyloCSF {
  
  input:
  set val(setname), file(novelgff), file('bigwigs') from novelGFF3_phylo

  output:
  set val(setname), file('phylocsf.txt') into phylocsf_out

  """
  calculate_phylocsf.py $novelgff bigwigs phylocsf.txt
  """

}


bamFiles = params.bamfiles ? Channel.fromPath(params.bamfiles).map { fn -> [ fn, fn + '.bai' ] } : Channel.empty()

process scanBams {

  when: params.bamfiles

  input:
  set val(setname), file(gff) from novelGFF3_bams
  file bams from bamFiles.collect()
  
  output:
  set val(setname), file('scannedbams.txt') into scannedbams

  """
  ls *.bam > bamfiles.txt
  scan_bams.py  --gff_input $gff --bam_files bamfiles.txt --output scannedbams.txt
  """
}

annoperl = Channel.fromPath("$params.annovar_dir/annotate_variation.pl")
annohumdb = Channel.fromPath("$params.annovar_dir/humandb/")

novelbed
  .combine(annoperl)
  .combine(annohumdb)
  .set { anno_in }

process annovar {
  
  input:
  set val(setname), file(novelbed), file(perlscript), file(humdb) from anno_in

  output:
  set val(setname), file('novpep_annovar.variant_function') into annovar_out

  """
  ./annotate_variation.pl -out novpep_annovar -build hg19 $novelbed humandb/
  """

}

annovar_out
  .join(annonovelpep)
  .set { parseanno }

process parseAnnovarOut {
  
  input:
  set val(setname), file(anno), file(novelpep) from parseanno

  output:
  set val(setname), file('parsed_annovar.txt') into annovar_parsed

  """
  parse_annovar_out.py --input $novelpep --output parsed_annovar.txt --annovar_out $anno 
  """
}


ns_snp_out
  .join(novpep_singlemisspecai)
  .join(peptable_blat)
  .join(annovar_parsed)
  .join(phastcons_out)
  .join(phylocsf_out)
  .join(novelpep_percoquant)
  .set { combined_novelprebam }

combined_novel = (params.bamfiles ? combined_novelprebam.join(scannedbams) : combined_novelprebam)


process combineResults{
  
  publishDir "${params.outdir}", mode: 'copy', overwrite: true
  
  input:
  set val(setname), file(a), file(b), file(c), file(d), file(e), file(f), file(g), file(h) from combined_novel
  
  output:
  set val('nov'), val(setname), file("${setname}_novel_peptides.txt") into novpeps_finished 
  
  script:
  if (!params.bamfiles)
  """
  for fn in $a $b $c $d $e $f $g; do sort -k 1b,1 \$fn > tmpfn; mv tmpfn \$fn; done
  join $a $b -a1 -a2 -o auto -e 'NA' -t \$'\\t' > joined1
  join joined1 $c -a1 -a2 -o auto -e 'NA' -t \$'\\t' > joined2
  join joined2 $d -a1 -a2 -o auto -e 'NA' -t \$'\\t' > joined3
  join joined3 $e -a1 -a2 -o auto -e 'NA' -t \$'\\t' > joined4
  join joined4 $f -a1 -a2 -o auto -e 'NA' -t \$'\\t' > joined5
  join joined5 $g -a1 -a2 -o auto -e 'NA' -t \$'\\t' > joined6
  grep '^Bare peptide' joined6 > ${setname}_novel_peptides.txt
  grep -v '^Bare peptide' joined6 >> ${setname}_novel_peptides.txt
  """

  else
  """
  for fn in $a $b $c $d $e $f $g $h; do sort -k 1b,1 \$fn > tmpfn; mv tmpfn \$fn; done
  join $a $b -a1 -a2 -o auto -e 'NA' -t \$'\\t' > joined1
  join joined1 $c -a1 -a2 -o auto -e 'NA' -t \$'\\t' > joined2
  join joined2 $d -a1 -a2 -o auto -e 'NA' -t \$'\\t' > joined3
  join joined3 $e -a1 -a2 -o auto -e 'NA' -t \$'\\t' > joined4
  join joined4 $f -a1 -a2 -o auto -e 'NA' -t \$'\\t' > joined5
  join joined5 $g -a1 -a2 -o auto -e 'NA' -t \$'\\t' > joined6
  join joined6 $h -a1 -a2 -o auto -e 'NA' -t \$'\\t' > joined7
  grep '^Bare peptide' joined7 > ${setname}_novel_peptides.txt
  grep -v '^Bare peptide' joined7 >> ${setname}_novel_peptides.txt
  """
}

process prepSpectrumAI {

  input:
  set val(setname), val(peptype), file(psms) from variantpsms
  
  output:
  set val(setname), file('specai_in.txt') into var_specai_input
  
  script:
  if (params.saavheader)
  """
  cat <(head -n1 $psms) <(grep $params.saavheader $psms) > saavpsms
  label_sub_pos.py --input_psm saavpsms --output specai_in.txt ${params.splitchar ? "--splitchar ${params.splitchar}" : ''}
  """
  else
  """
  label_sub_pos.py --input_psm $psms --output specai_in.txt
  """
}


var_specaimzmls
  .join(var_specai_input)
  .set { var_specai_inmzml }

process SpectrumAI {

  publishDir "${params.outdir}", mode: 'copy', overwrite: true, saveAs: { it == "precursorError.histogram.plot.pdf" ? "${setname}_variant_precursorError_plot.pdf" : it }

  input:
  set val(setname), val(samples), file(mzmls), file(specai_in) from var_specai_inmzml

  output:
  set val(setname), file("${setname}_variant_specairesult.txt") into specai
  file "precursorError.histogram.plot.pdf" into specai_plot

  """
  mkdir mzmls
  for fn in $mzmls; do ln -s `pwd`/\$fn mzmls/; done
  ls mzmls
  Rscript /SpectrumAI/SpectrumAI.R mzmls $specai_in ${setname}_variant_specairesult.txt
  """
}

specai
  .join(presai_peptable)
  .set { specai_peptable }

process mapVariantPeptidesToGenome {

  publishDir "${params.outdir}", mode: 'copy', overwrite: true

  input:
  set val(setname), file(x), val(peptype), file(peptides) from specai_peptable
  file cosmic
  file dbsnp
  
  output:
  set val('var'), val(setname), file("${setname}_variant_peptides.txt") into varpeps_finished
  file "${setname}_variant_peptides.saav.pep.hg19cor.vcf" into saavvcfs_finished

  """
  ${params.saavheader ? "cat <(head -n1 ${peptides}) <(grep ${params.saavheader} ${peptides}) > saavpeps" : "mv ${peptides} saavpeps" }
  parse_spectrumAI_out.py --spectrumAI_out $x --input saavpeps --output setsaavs
  ${params.saavheader ? "cat setsaavs <(grep -v ${params.saavheader} ${peptides} | sed \$'s/\$/\tNA/') > ${setname}_variant_peptides.txt" : "mv setsaavs ${setname}_variant_peptides.txt"}
  map_cosmic_snp_tohg19.py --input ${setname}_variant_peptides.txt --output ${setname}_variant_peptides.saav.pep.hg19cor.vcf --cosmic_input $cosmic --dbsnp_input $dbsnp
  # Remove PSM-table specific stuff (RT, precursor, etc etc) from variant PEPTIDE table
  cut -f 1,2,14-5000 ${setname}_variant_peptides.txt > pepsfix
  mv pepsfix ${setname}_variant_peptides.txt
  """
}

novpeps_finished
  .concat(varpeps_finished) 
  .groupTuple()
  .set { setmerge_peps }


accession_keymap = ['var': 'Peptide sequence', 'nov': 'Peptide']
acc_removemap = ['nov': 'Bare peptide', 'var': 'Mod.peptide']


process mergeSetPeptidetable {

  publishDir "${params.outdir}", mode: 'copy', overwrite: true

  input:
  set val(peptype), val(setnames), file('peps?') from setmerge_peps

  output:
  file "${peptype}_peptidetable.txt" into produced_peptables

  """
  # build non-changing fields (seq based fields) table:
  fixfields=`head -n1 peps1 |tr -s '\\t' '\\n' | egrep -vn '(Setname|Spectrum|Files|Charge|q-val|plex|${acc_removemap[peptype]})' | cut -f 1 -d ':'`
  fixfields=`echo \$fixfields | sed 's/ /,/g'`
  head -n1 peps1 | cut -f `echo \$fixfields` > fixheader
  count=1; for setn in ${setnames.join(' ')} ; do
  cut -f  `echo \$fixfields` peps\$count | tail -n+2 >> fixpeps
  ((count++))
  done
  if [ ${peptype} == 'nov' ]
  then
     cat fixheader <(sort -u -k1b,1 fixpeps) > temp
     group_novpepToLoci.py  --input temp --output temp.loci --distance 10kb
     head -n1 temp.loci > fixheader
     tail -n+2 temp.loci > fixpeps
  fi
  sort -u -k1b,1 fixpeps > temp
  mv temp fixpeps

  ## Build changing fields table
  touch peptable
  count=1; for setn in ${setnames.join(' ')}; do
    varfields=`head -n1 peps\$count |tr -s '\\t' '\\n' | egrep -n '(${accession_keymap[peptype]}|Spectrum|q-val|plex)' | cut -f 1 -d ':'`
    varfields=`echo \$varfields| sed 's/ /,/g'`
    # first add to header, cut from f2 to remove join-key pep seq field
    head -n1 peps\$count | cut -f `echo \$varfields` | cut -f 2-5000| sed "s/^\\(\\w\\)/\${setn}_\\1/;s/\\(\\s\\)/\\1\${setn}_/g" > varhead
    paste fixheader varhead > newheader && mv newheader fixheader
    # then join the values
    tail -n+2 peps\$count | cut -f `echo \$varfields` | sort -k1b,1 > sortpep; join peptable sortpep -a1 -a2 -o auto -e 'NA' -t \$'\\t' > joined
    mv joined peptable
    ((count++))
  done
  join fixpeps peptable -a1 -a2 -o auto -e 'NA' -t \$'\\t' > fixvarpeps
  cat fixheader fixvarpeps > ${peptype}_peptidetable.txt
  """
}


produced_psmtables
  .concat(produced_peptables)
  .subscribe { println "Pipeline output ready: ${it}" }