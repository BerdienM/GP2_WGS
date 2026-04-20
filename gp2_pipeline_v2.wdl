## GP2 Purine Metabolism Pipeline — WDL Workflow
##
## Runs steps 02-06 for multiple ancestries in parallel on Terra/Cromwell.
## Assumes per-ancestry VCFs and covariate files are already in the workspace bucket
## (produced by notebooks 00_extract and 01_covariate_file).
##
## Steps per ancestry (run in sequence within each scatter):
##   02_annotation      — VEP annotation → Tier 1/2 TSVs
##   00b_group_files    — Tier TSVs → SAIGE group file
##   00c_genome_wide    — per-chr PLINK2 → genome-wide pruned bed
##   00d_regenie_annot  — Tier TSVs → REGENIE annotation files
##   03_rare_variant    — SAIGE-GENE+ null model + SKAT-O
##   04_common_variant  — REGENIE step 1 + step 2
##   06_geneset         — MAGMA gene-level + gene-set analysis
##
## Usage:
##   Upload this WDL and the inputs JSON to Terra.
##   Run via Terra workflow submission — Cromwell handles parallelism and retries.

version 1.0

# ─────────────────────────────────────────────────────────────────────────────
# Main workflow
# ─────────────────────────────────────────────────────────────────────────────

workflow GP2PipelineMultiAncestry {

  input {
    # Ancestries to process
    Array[String] ancestries = ["EUR", "AFR", "EAS", "CAS", "MDE"]

    # Gene set
    String gene_set = "purine_metabolism"

    # Workspace bucket
    String workspace_bucket

    # GP2 data bucket (requester-pays)
    String gp2_bucket = "gs://gp2tier2/release10"
    String billing_project

    # Gene set CSV in workspace bucket
    File gene_set_csv

    # VEP cache and plugin paths (in workspace bucket or persistent disk)
    String vep_cache_dir
    String cadd_snv
    String cadd_indel
    String loftee_dir

    # MAGMA reference
    File magma_gene_loc

    # Per-ancestry covariate files (Map: ancestry -> covariate CSV)
    Map[String, File] covariate_files

    # Runtime parameters
    Int vep_cpus        = 8
    Int vep_mem_gb      = 16
    Int saige_cpus      = 4
    Int saige_mem_gb    = 16
    Int regenie_cpus    = 4
    Int regenie_mem_gb  = 16
    Int magma_cpus      = 2
    Int magma_mem_gb    = 8
    String disk_type    = "SSD"
  }

  # Scatter over ancestries — all run in parallel
  scatter (ancestry in ancestries) {

    # ── Step 02: VEP annotation ───────────────────────────────────────────────
    call Annotate {
      input:
        ancestry         = ancestry,
        gene_set         = gene_set,
        workspace_bucket = workspace_bucket,
        vep_cache_dir    = vep_cache_dir,
        cadd_snv         = cadd_snv,
        cadd_indel       = cadd_indel,
        loftee_dir       = loftee_dir,
        cpus             = vep_cpus,
        mem_gb           = vep_mem_gb,
        disk_type        = disk_type
    }

    # ── Step 00b: Build SAIGE group file ──────────────────────────────────────
    call BuildGroupFile {
      input:
        ancestry         = ancestry,
        gene_set         = gene_set,
        workspace_bucket = workspace_bucket,
        tier1_tsv        = Annotate.tier1_tsv,
        tier2_tsv        = Annotate.tier2_tsv,
        gene_set_csv     = gene_set_csv
    }

    # ── Step 00c: Build genome-wide bfile ─────────────────────────────────────
    call BuildGenomeWide {
      input:
        ancestry         = ancestry,
        workspace_bucket = workspace_bucket,
        gp2_bucket       = gp2_bucket,
        billing_project  = billing_project,
        disk_type        = disk_type
    }

    # ── Step 00d: Build REGENIE annotation files ──────────────────────────────
    call BuildRegenieAnnotations {
      input:
        ancestry         = ancestry,
        gene_set         = gene_set,
        workspace_bucket = workspace_bucket,
        tier1_tsv        = Annotate.tier1_tsv,
        tier2_tsv        = Annotate.tier2_tsv,
        gene_set_csv     = gene_set_csv
    }

    # ── Step 03: SAIGE rare variant analysis ──────────────────────────────────
    call RareVariant {
      input:
        ancestry         = ancestry,
        gene_set         = gene_set,
        workspace_bucket = workspace_bucket,
        group_file       = BuildGroupFile.group_file,
        genome_wide_bed  = BuildGenomeWide.genome_wide_bed,
        genome_wide_bim  = BuildGenomeWide.genome_wide_bim,
        genome_wide_fam  = BuildGenomeWide.genome_wide_fam,
        gene_set_vcf     = BuildGroupFile.gene_set_vcf_gp2id,
        covariate_file   = covariate_files[ancestry],
        cpus             = saige_cpus,
        mem_gb           = saige_mem_gb,
        disk_type        = disk_type
    }

    # ── Step 04: REGENIE common variant analysis ──────────────────────────────
    call CommonVariant {
      input:
        ancestry         = ancestry,
        gene_set         = gene_set,
        workspace_bucket = workspace_bucket,
        anno_file        = BuildRegenieAnnotations.anno_file,
        set_list         = BuildRegenieAnnotations.set_list,
        mask_file        = BuildRegenieAnnotations.mask_file,
        genome_wide_bed  = BuildGenomeWide.genome_wide_bed,
        genome_wide_bim  = BuildGenomeWide.genome_wide_bim,
        genome_wide_fam  = BuildGenomeWide.genome_wide_fam,
        gene_set_vcf     = BuildGroupFile.gene_set_vcf_gp2id,
        genome_wide_psam = BuildGenomeWide.genome_wide_psam,
        covariate_file   = covariate_files[ancestry],
        cpus             = regenie_cpus,
        mem_gb           = regenie_mem_gb,
        disk_type        = disk_type
    }

    # ── Step 06: MAGMA gene-set analysis ──────────────────────────────────────
    call GeneSetAnalysis {
      input:
        ancestry             = ancestry,
        gene_set             = gene_set,
        workspace_bucket     = workspace_bucket,
        regenie_sv_out       = CommonVariant.sv_results,
        genome_wide_bed      = BuildGenomeWide.genome_wide_bed,
        genome_wide_bim      = BuildGenomeWide.genome_wide_bim,
        genome_wide_fam      = BuildGenomeWide.genome_wide_fam,
        gene_set_csv         = gene_set_csv,
        magma_gene_loc       = magma_gene_loc,
        n_samples            = CommonVariant.n_samples,
        cpus                 = magma_cpus,
        mem_gb               = magma_mem_gb,
        disk_type            = disk_type
    }
  }

  output {
    # Per-ancestry outputs
    Array[File] tier1_tsvs        = Annotate.tier1_tsv
    Array[File] tier2_tsvs        = Annotate.tier2_tsv
    Array[File] saige_results     = RareVariant.skat_results
    Array[File] regenie_sv        = CommonVariant.sv_results
    Array[File] regenie_masks     = CommonVariant.mask_results
    Array[File] magma_gene_out    = GeneSetAnalysis.gene_results
    Array[File] magma_geneset_out = GeneSetAnalysis.geneset_results
  }
}

# ─────────────────────────────────────────────────────────────────────────────
# Task: 02 — VEP Annotation
# ─────────────────────────────────────────────────────────────────────────────

task Annotate {
  input {
    String ancestry
    String gene_set
    String workspace_bucket
    String vep_cache_dir
    String cadd_snv
    String cadd_indel
    String loftee_dir
    Int    cpus
    Int    mem_gb
    String disk_type
  }

  String input_vcf = "~{workspace_bucket}/data/extracted/~{gene_set}/~{ancestry}/~{gene_set}_~{ancestry}.vcf.gz"
  String out_dir   = "~{workspace_bucket}/results/~{gene_set}/~{ancestry}/02_annotation"

  command <<<
    set -euo pipefail

    ANCESTRY="~{ancestry}"
    GENE_SET="~{gene_set}"
    INPUT_VCF="~{input_vcf}"
    OUT_DIR="/tmp/annotation"
    mkdir -p $OUT_DIR

    # Copy VCF locally
    gsutil cp ${INPUT_VCF} $OUT_DIR/${GENE_SET}_${ANCESTRY}.vcf.gz
    gsutil cp ${INPUT_VCF}.tbi $OUT_DIR/${GENE_SET}_${ANCESTRY}.vcf.gz.tbi 2>/dev/null || \
    gsutil cp ${INPUT_VCF}.csi $OUT_DIR/${GENE_SET}_${ANCESTRY}.vcf.gz.csi 2>/dev/null || true

    VEP_OUT_VCF=$OUT_DIR/${GENE_SET}_${ANCESTRY}_vep.vcf
    VEP_OUT_TSV=$OUT_DIR/${GENE_SET}_${ANCESTRY}_vep.tsv

    # Run VEP
    vep \
      --input_file $OUT_DIR/${GENE_SET}_${ANCESTRY}.vcf.gz \
      --output_file $VEP_OUT_VCF \
      --tab \
      --stats_file $OUT_DIR/${GENE_SET}_${ANCESTRY}_vep_stats.html \
      --warning_file $OUT_DIR/${GENE_SET}_${ANCESTRY}_vep_warnings.txt \
      --species homo_sapiens \
      --assembly GRCh38 \
      --offline \
      --cache \
      --dir_cache ~{vep_cache_dir} \
      --everything \
      --canonical \
      --vcf \
      --fork ~{cpus} \
      --plugin CADD,~{cadd_snv},~{cadd_indel} \
      --plugin LoF,loftee_path:~{loftee_dir},human_ancestor_fa:false,filter_position:0.05,min_intron_size:15 \
      --fasta ~{vep_cache_dir}/homo_sapiens/105_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz

    # Convert VCF to TSV
    grep -v "^##" $VEP_OUT_VCF | sed 's/^#//' > $VEP_OUT_TSV

    # Tier filtering
    python3 << 'PYEOF'
import csv, sys

tier1_csq = {
    "transcript_ablation","splice_acceptor_variant","splice_donor_variant",
    "stop_gained","frameshift_variant","stop_lost","start_lost"
}
tier2_missense = {"missense_variant","protein_altering_variant","inframe_insertion","inframe_deletion"}
tier2_splice   = {"splice_region_variant","splice_donor_5th_base_variant","splice_polypyrimidine_tract_variant"}
MAF_THRESH = 0.01
CADD_THRESH = 20.0

in_path  = "$VEP_OUT_TSV"
t1_path  = "$OUT_DIR/~{gene_set}_~{ancestry}_tier1_lof.tsv"
t2_path  = "$OUT_DIR/~{gene_set}_~{ancestry}_tier2_missense_splice.tsv"

def get_float(row, *cols):
    for col in cols:
        v = row.get(col,"").strip()
        if v not in ("",".","-"):
            try: return float(v)
            except: pass
    return None

def passes_maf(row):
    af = get_float(row,"MAX_AF","AF")
    return af is None or af < MAF_THRESH

def is_canonical(row):
    return row.get("CANONICAL","") in ("YES","")

t1_fh = t2_fh = None
t1_n = t2_n = 0

with open(in_path) as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    for row in reader:
        if not is_canonical(row) or not passes_maf(row):
            continue
        csq   = set(row.get("Consequence","").split("&"))
        cadd  = get_float(row,"CADD_PHRED","CADD_phred")
        lof_hc = row.get("LoF","").upper() == "HC"
        tier1 = bool(csq & tier1_csq)
        if tier1 and row.get("LoF","") != "":
            tier1 = lof_hc
        missense_ok = bool(csq & tier2_missense) and (cadd is None or cadd > CADD_THRESH)
        splice_ok   = bool(csq & tier2_splice)
        tier2 = (missense_ok or splice_ok) and not tier1
        for flag, path, fh_ref, n_ref in [(tier1,t1_path,"t1","t1_n"),(tier2,t2_path,"t2","t2_n")]:
            if flag:
                if eval(fh_ref) is None:
                    exec(f"{fh_ref} = open(path,'w',newline='')")
                    exec(f"csv.DictWriter(eval(fh_ref),fieldnames=reader.fieldnames,delimiter='\\t').writeheader()")
                exec(f"csv.DictWriter(eval(fh_ref),fieldnames=reader.fieldnames,delimiter='\\t').writerow(row)")
                exec(f"{n_ref} += 1")

for h in [t1_fh, t2_fh]:
    if h: h.close()

print(f"Tier 1: {t1_n}, Tier 2: {t2_n}")
PYEOF

    # Copy outputs to bucket
    gsutil -m cp $OUT_DIR/*.tsv $OUT_DIR/*.html $OUT_DIR/*.txt ~{out_dir}/
    echo "Annotation complete for ~{ancestry}"
  >>>

  output {
    File tier1_tsv = "~{out_dir}/~{gene_set}_~{ancestry}_tier1_lof.tsv"
    File tier2_tsv = "~{out_dir}/~{gene_set}_~{ancestry}_tier2_missense_splice.tsv"
  }

  runtime {
    docker: "ensemblorg/ensembl-vep:release_105"
    cpu: cpus
    memory: "~{mem_gb} GB"
    disks: "local-disk 200 ~{disk_type}"
    preemptible: 1
  }
}

# ─────────────────────────────────────────────────────────────────────────────
# Task: 00b — Build SAIGE group file
# ─────────────────────────────────────────────────────────────────────────────

task BuildGroupFile {
  input {
    String ancestry
    String gene_set
    String workspace_bucket
    File   tier1_tsv
    File   tier2_tsv
    File   gene_set_csv
  }

  String out_dir = "~{workspace_bucket}/results/~{gene_set}/~{ancestry}/03_rare_variant"

  command <<<
    set -euo pipefail

    python3 << 'PYEOF'
import pandas as pd
from pathlib import Path

MISSENSE_CSQ = {"missense_variant","protein_altering_variant","inframe_insertion","inframe_deletion"}
SPLICE_CSQ   = {"splice_region_variant","splice_donor_5th_base_variant","splice_polypyrimidine_tract_variant"}

def make_var_id(row):
    return f"{row['CHROM']}:{row['POS']}:{row['REF']}:{row['ALT']}"

def get_anno(row):
    csq = set(str(row.get('Consequence','')).split('&'))
    return 'missense' if csq & MISSENSE_CSQ else 'splice' if csq & SPLICE_CSQ else 'lof'

t1 = pd.read_csv("~{tier1_tsv}", sep='\t', low_memory=False)
t2 = pd.read_csv("~{tier2_tsv}", sep='\t', low_memory=False)

records = []
for _, row in t1.iterrows():
    gene = str(row.get('SYMBOL','')).strip()
    if gene and gene != 'nan':
        records.append({'gene': gene, 'var_id': make_var_id(row), 'anno': 'lof'})
for _, row in t2.iterrows():
    gene = str(row.get('SYMBOL','')).strip()
    if gene and gene != 'nan':
        records.append({'gene': gene, 'var_id': make_var_id(row), 'anno': get_anno(row)})

df = pd.DataFrame(records).drop_duplicates(subset=['gene','var_id'])

with open('/tmp/group_file.txt','w') as f:
    for gene, grp in df.groupby('gene'):
        f.write(f"{gene}\tvar\t" + '\t'.join(grp['var_id']) + '\n')
        f.write(f"{gene}\tanno\t" + '\t'.join(grp['anno']) + '\n')

print(f"Group file: {df['gene'].nunique()} genes, {len(df)} variants")
PYEOF

    gsutil cp /tmp/group_file.txt ~{out_dir}/~{gene_set}_~{ancestry}_saige_group.txt

    # Also create GP2ID-renamed VCF for SAIGE step 2 and REGENIE
    INPUT_VCF="~{workspace_bucket}/data/extracted/~{gene_set}/~{ancestry}/~{gene_set}_~{ancestry}.vcf.gz"
    gsutil cp $INPUT_VCF /tmp/input.vcf.gz
    gsutil cp $INPUT_VCF.tbi /tmp/input.vcf.gz.tbi 2>/dev/null || true

    bcftools query -l /tmp/input.vcf.gz | awk -F'_' '{print $0, $(NF-1)"_"$NF}' > /tmp/rename.txt
    bcftools reheader --samples /tmp/rename.txt /tmp/input.vcf.gz \
        -o /tmp/~{gene_set}_~{ancestry}_gp2id.vcf.gz
    bcftools index /tmp/~{gene_set}_~{ancestry}_gp2id.vcf.gz

    gsutil cp /tmp/~{gene_set}_~{ancestry}_gp2id.vcf.gz \
        ~{workspace_bucket}/data/extracted/~{gene_set}/~{ancestry}/
    gsutil cp /tmp/~{gene_set}_~{ancestry}_gp2id.vcf.gz.csi \
        ~{workspace_bucket}/data/extracted/~{gene_set}/~{ancestry}/
  >>>

  output {
    File group_file        = "~{out_dir}/~{gene_set}_~{ancestry}_saige_group.txt"
    File gene_set_vcf_gp2id = "~{workspace_bucket}/data/extracted/~{gene_set}/~{ancestry}/~{gene_set}_~{ancestry}_gp2id.vcf.gz"
  }

  runtime {
    docker: "us.gcr.io/broad-gatk/gatk:4.3.0.0"
    cpu: 2
    memory: "8 GB"
    disks: "local-disk 100 SSD"
    preemptible: 2
  }
}

# ─────────────────────────────────────────────────────────────────────────────
# Task: 00c — Build genome-wide bfile
# ─────────────────────────────────────────────────────────────────────────────

task BuildGenomeWide {
  input {
    String ancestry
    String workspace_bucket
    String gp2_bucket
    String billing_project
    String disk_type
  }

  String out_dir = "~{workspace_bucket}/data/genome_wide/~{ancestry}"

  command <<<
    set -euo pipefail

    ANCESTRY="~{ancestry}"
    WORKDIR="/tmp/genome_wide"
    SCRATCH="/tmp/genome_wide/tmp_chr"
    mkdir -p $WORKDIR $SCRATCH

    # Download and QC per chromosome
    for CHR in $(seq 1 22); do
        BUCKET_PREFIX="~{gp2_bucket}/wgs/deepvariant_joint_calling/plink/${ANCESTRY}/chr${CHR}_${ANCESTRY}_release10"
        LOCAL_PREFIX="$SCRATCH/chr${CHR}_${ANCESTRY}_release10"
        QC_PREFIX="$SCRATCH/chr${CHR}_${ANCESTRY}_qc"

        for EXT in pgen pvar psam; do
            gsutil -u ~{billing_project} cp ${BUCKET_PREFIX}.${EXT} ${LOCAL_PREFIX}.${EXT}
        done

        plink2 \
            --pfile $LOCAL_PREFIX \
            --maf 0.01 --geno 0.05 --hwe 1e-4 \
            --make-pgen --out $QC_PREFIX --silent

        rm -f ${LOCAL_PREFIX}.pgen ${LOCAL_PREFIX}.pvar ${LOCAL_PREFIX}.psam
        echo "chr${CHR}: QC done"
    done

    # Merge
    ls $SCRATCH/*_qc.pgen | sed 's/.pgen//' > $WORKDIR/merge_list.txt
    plink2 --pmerge-list $WORKDIR/merge_list.txt pfile \
        --make-pgen --out $WORKDIR/${ANCESTRY}_merged --silent

    # LD prune
    plink2 --pfile $WORKDIR/${ANCESTRY}_merged \
        --indep-pairwise 1000 10 0.02 \
        --out $WORKDIR/${ANCESTRY} --silent
    plink2 --pfile $WORKDIR/${ANCESTRY}_merged \
        --extract $WORKDIR/${ANCESTRY}.prune.in \
        --make-pgen --out $WORKDIR/${ANCESTRY}_genome_wide_pruned --silent

    # Sort and convert to PLINK1 (autosomes only for SAIGE/REGENIE)
    plink2 --pfile $WORKDIR/${ANCESTRY}_genome_wide_pruned \
        --sort-vars --make-pgen \
        --out $WORKDIR/${ANCESTRY}_genome_wide_pruned_sorted --silent
    plink2 --pfile $WORKDIR/${ANCESTRY}_genome_wide_pruned_sorted \
        --chr 1-22 --make-bed \
        --out $WORKDIR/${ANCESTRY}_genome_wide_pruned_bed --silent

    # Copy to bucket
    for EXT in pgen pvar psam; do
        gsutil cp $WORKDIR/${ANCESTRY}_genome_wide_pruned.${EXT} ~{out_dir}/
    done
    for EXT in bed bim fam; do
        gsutil cp $WORKDIR/${ANCESTRY}_genome_wide_pruned_bed.${EXT} ~{out_dir}/
    done
    echo "Genome-wide bfile complete for ${ANCESTRY}"
  >>>

  output {
    File genome_wide_bed  = "~{out_dir}/~{ancestry}_genome_wide_pruned_bed.bed"
    File genome_wide_bim  = "~{out_dir}/~{ancestry}_genome_wide_pruned_bed.bim"
    File genome_wide_fam  = "~{out_dir}/~{ancestry}_genome_wide_pruned_bed.fam"
    File genome_wide_psam = "~{out_dir}/~{ancestry}_genome_wide_pruned.psam"
  }

  runtime {
    docker: "us.gcr.io/broad-gatk/gatk:4.3.0.0"
    cpu: 4
    memory: "16 GB"
    disks: "local-disk 200 ~{disk_type}"
    preemptible: 1
  }
}

# ─────────────────────────────────────────────────────────────────────────────
# Task: 00d — Build REGENIE annotation files
# ─────────────────────────────────────────────────────────────────────────────

task BuildRegenieAnnotations {
  input {
    String ancestry
    String gene_set
    String workspace_bucket
    File   tier1_tsv
    File   tier2_tsv
    File   gene_set_csv
  }

  String out_dir = "~{workspace_bucket}/results/~{gene_set}/~{ancestry}/04_common_variant"

  command <<<
    set -euo pipefail

    python3 << 'PYEOF'
import pandas as pd

MISSENSE_CSQ = {"missense_variant","protein_altering_variant","inframe_insertion","inframe_deletion"}
SPLICE_CSQ   = {"splice_region_variant","splice_donor_5th_base_variant","splice_polypyrimidine_tract_variant"}

def make_var_id(row):
    chrom = str(row['CHROM'])
    if not chrom.startswith('chr'): chrom = f"chr{chrom}"
    return f"{chrom}:{row['POS']}:{row['REF']}:{row['ALT']}"

def get_anno(row):
    csq = set(str(row.get('Consequence','')).split('&'))
    return 'missense' if csq & MISSENSE_CSQ else 'splice' if csq & SPLICE_CSQ else 'lof'

t1 = pd.read_csv("~{tier1_tsv}", sep='\t', low_memory=False)
t2 = pd.read_csv("~{tier2_tsv}", sep='\t', low_memory=False)
gs = pd.read_csv("~{gene_set_csv}")

records = []
for _, row in t1.iterrows():
    gene = str(row.get('SYMBOL','')).strip()
    if gene and gene != 'nan':
        records.append({'gene': gene, 'var_id': make_var_id(row), 'anno': 'lof'})
for _, row in t2.iterrows():
    gene = str(row.get('SYMBOL','')).strip()
    if gene and gene != 'nan':
        records.append({'gene': gene, 'var_id': make_var_id(row), 'anno': get_anno(row)})

df = pd.DataFrame(records).drop_duplicates(subset=['gene','var_id'])

# Anno file
with open('/tmp/anno.txt','w') as f:
    for _, row in df.iterrows():
        f.write(f"{row['var_id']} {row['gene']} {row['anno']}\n")

# Set list
all_vars = pd.concat([
    t1[['CHROM','POS','SYMBOL']].rename(columns={'SYMBOL':'gene'}),
    t2[['CHROM','POS','SYMBOL']].rename(columns={'SYMBOL':'gene'})
])
all_vars['POS'] = pd.to_numeric(all_vars['POS'], errors='coerce')
all_vars = all_vars.dropna(subset=['POS'])

with open('/tmp/setlist.txt','w') as f:
    for gene, grp in df.groupby('gene'):
        chrom = str(grp.iloc[0]['var_id'].split(':')[0]).replace('chr','')
        start = min(int(v.split(':')[1]) for v in grp['var_id'])
        var_list = ','.join(grp['var_id'].tolist())
        f.write(f"{gene}\t{chrom}\t{start}\t{var_list}\n")

# Mask file
masks = [
    ("M1_lof",          "lof"),
    ("M2_missense",     "missense"),
    ("M3_splice",       "splice"),
    ("M4_lof_missense", "lof,missense"),
    ("M5_all",          "lof,missense,splice"),
]
with open('/tmp/masks.txt','w') as f:
    for name, annos in masks:
        f.write(f"{name} {annos}\n")

print(f"Anno: {len(df)} variants, SetList: {df['gene'].nunique()} genes, Masks: {len(masks)}")
PYEOF

    gsutil cp /tmp/anno.txt    ~{out_dir}/~{gene_set}_~{ancestry}_regenie_anno.txt
    gsutil cp /tmp/setlist.txt ~{out_dir}/~{gene_set}_~{ancestry}_regenie_setlist.txt
    gsutil cp /tmp/masks.txt   ~{out_dir}/~{gene_set}_~{ancestry}_regenie_masks.txt
  >>>

  output {
    File anno_file = "~{out_dir}/~{gene_set}_~{ancestry}_regenie_anno.txt"
    File set_list  = "~{out_dir}/~{gene_set}_~{ancestry}_regenie_setlist.txt"
    File mask_file = "~{out_dir}/~{gene_set}_~{ancestry}_regenie_masks.txt"
  }

  runtime {
    docker: "us.gcr.io/broad-gatk/gatk:4.3.0.0"
    cpu: 2
    memory: "8 GB"
    disks: "local-disk 50 SSD"
    preemptible: 2
  }
}

# ─────────────────────────────────────────────────────────────────────────────
# Task: 03 — SAIGE rare variant analysis
# ─────────────────────────────────────────────────────────────────────────────

task RareVariant {
  input {
    String ancestry
    String gene_set
    String workspace_bucket
    File   group_file
    File   genome_wide_bed
    File   genome_wide_bim
    File   genome_wide_fam
    File   gene_set_vcf
    File   covariate_file
    Int    cpus
    Int    mem_gb
    String disk_type
  }

  String out_dir = "~{workspace_bucket}/results/~{gene_set}/~{ancestry}/03_rare_variant"

  command <<<
    set -euo pipefail

    ANCESTRY="~{ancestry}"
    GENE_SET="~{gene_set}"
    WORKDIR="/tmp/saige"
    mkdir -p $WORKDIR

    # Prepare covariate file — GP2ID format, tab-separated
    python3 << 'PYEOF'
import pandas as pd, subprocess

covar = pd.read_csv("~{covariate_file}")
if 'GP2ID' in covar.columns and 'IID' not in covar.columns:
    covar = covar.rename(columns={'GP2ID': 'IID'})
if 'FID' not in covar.columns:
    covar.insert(0, 'FID', covar['IID'])
covar.to_csv('/tmp/saige/covariate_step1.tsv', sep='\t', index=False)

# Remap IDs for step 2 VCF
vcf_ids = subprocess.run(
    "bcftools query -l ~{gene_set_vcf}",
    shell=True, capture_output=True, text=True
).stdout.strip().split('\n')
def gp2id(v): parts=v.split('_'); return '_'.join(parts[-2:])
vcf_map = {gp2id(v): v for v in vcf_ids}
covar2 = covar.copy()
covar2['FID'] = covar2['FID'].map(vcf_map)
covar2['IID'] = covar2['IID'].map(vcf_map)
covar2 = covar2.dropna(subset=['FID','IID'])
covar2.to_csv('/tmp/saige/covariate_step2.tsv', sep='\t', index=False)
print(f"Step1: {len(covar)} samples, Step2: {len(covar2)} samples")
PYEOF

    # Copy VCF locally and create CSI index
    # (WDL localizes files but index must be co-located with VCF)
    cp ~{gene_set_vcf} $WORKDIR/input.vcf.gz
    bcftools index --csi $WORKDIR/input.vcf.gz

    # PLINK1 bed prefix — strip .bed suffix from localized path
    BED_PREFIX=$(echo "~{genome_wide_bed}" | sed 's/\.bed$//')

    # SAIGE step 1 — null model
    Rscript /usr/local/bin/step1_fitNULLGLMM.R \
        --plinkFile=$BED_PREFIX \
        --phenoFile=/tmp/saige/covariate_step1.tsv \
        --phenoCol=PHENO \
        --sampleIDColinphenoFile=FID \
        --covarColList=SEX,age,cohort,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
        --qCovarColList=cohort \
        --traitType=binary \
        --outputPrefix=$WORKDIR/${GENE_SET}_${ANCESTRY} \
        --nThreads=~{cpus} \
        --LOCO=FALSE

    # SAIGE step 2 — SKAT-O gene tests
    Rscript /usr/local/bin/step2_SPAtests.R \
        --vcfFile=$WORKDIR/input.vcf.gz \
        --vcfFileIndex=$WORKDIR/input.vcf.gz.csi \
        --vcfField=GT \
        --GMMATmodelFile=$WORKDIR/${GENE_SET}_${ANCESTRY}.rda \
        --varianceRatioFile=$WORKDIR/${GENE_SET}_${ANCESTRY}.varianceRatio.txt \
        --groupFile=~{group_file} \
        --LOCO=FALSE \
        --is_output_moreDetails=TRUE \
        --MACCutoff_to_CollapseUltraRare=5 \
        --maxMAF_in_groupTest=0.01 \
        --annotation_in_groupTest=lof,missense,splice,lof:missense,lof:missense:splice \
        --SAIGEOutputFile=$WORKDIR/${GENE_SET}_${ANCESTRY}_skat_results

    # Copy results to bucket
    gsutil -m cp $WORKDIR/*.txt $WORKDIR/*skat_results* ~{out_dir}/
    echo "SAIGE complete for ${ANCESTRY}"
  >>>

  output {
    File skat_results = "~{out_dir}/~{gene_set}_~{ancestry}_skat_results"
    File var_ratio    = "~{out_dir}/~{gene_set}_~{ancestry}.varianceRatio.txt"
  }

  runtime {
    docker: "wzhou88/saige:1.3.1"
    cpu: cpus
    memory: "~{mem_gb} GB"
    disks: "local-disk 100 ~{disk_type}"
    preemptible: 1
  }
}

# ─────────────────────────────────────────────────────────────────────────────
# Task: 04 — REGENIE common variant analysis
# ─────────────────────────────────────────────────────────────────────────────

task CommonVariant {
  input {
    String ancestry
    String gene_set
    String workspace_bucket
    File   anno_file
    File   set_list
    File   mask_file
    File   genome_wide_bed
    File   genome_wide_bim
    File   genome_wide_fam
    File   gene_set_vcf
    File   genome_wide_psam
    File   covariate_file
    Int    cpus
    Int    mem_gb
    String disk_type
  }

  String out_dir = "~{workspace_bucket}/results/~{gene_set}/~{ancestry}/04_common_variant"

  command <<<
    set -euo pipefail

    ANCESTRY="~{ancestry}"
    GENE_SET="~{gene_set}"
    WORKDIR="/tmp/regenie"
    mkdir -p $WORKDIR

    BED_PREFIX=$(echo "~{genome_wide_bed}" | sed 's/\.bed$//')
    N_SAMPLES=$(tail -n +2 ~{covariate_file} | wc -l)

    # Prepare covariate file
    python3 << 'PYEOF'
import pandas as pd
covar = pd.read_csv("~{covariate_file}")
if 'GP2ID' in covar.columns and 'IID' not in covar.columns:
    covar = covar.rename(columns={'GP2ID': 'IID'})
if 'FID' not in covar.columns:
    covar.insert(0, 'FID', covar['IID'])
covar.to_csv('/tmp/regenie/covariate.tsv', sep='\t', index=False)
print(f"Samples: {len(covar)}")
PYEOF

    # Filter bfile to analysis samples
    awk 'NR>1 {print $1, $2}' /tmp/regenie/covariate.tsv > /tmp/regenie/sample_ids.txt
    plink2 \
        --bfile $BED_PREFIX \
        --keep /tmp/regenie/sample_ids.txt \
        --maf 0.01 --hwe 1e-15 --geno 0.1 \
        --make-bed --out $WORKDIR/genome_wide_filtered --silent

    # Convert gene-set VCF to PLINK2
    plink2 \
        --vcf ~{gene_set_vcf} \
        --psam ~{genome_wide_psam} \
        --split-par hg38 \
        --set-all-var-ids chr@:#:\$r:\$a \
        --new-id-max-allele-len 1000 \
        --make-pgen \
        --out $WORKDIR/geneset_plink2 --silent

    # REGENIE step 1
    regenie \
        --step 1 \
        --bed $WORKDIR/genome_wide_filtered \
        --phenoFile /tmp/regenie/covariate.tsv \
        --phenoCol PHENO \
        --covarFile /tmp/regenie/covariate.tsv \
        --covarColList SEX,age,cohort,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
        --catCovarList cohort \
        --bsize 1000 \
        --bt --lowmem \
        --lowmem-prefix $WORKDIR/regenie_tmp \
        --out $WORKDIR/${GENE_SET}_${ANCESTRY}_step1 \
        --threads ~{cpus}

    # REGENIE step 2 — gene-level with masks
    regenie \
        --step 2 \
        --pgen $WORKDIR/geneset_plink2 \
        --phenoFile /tmp/regenie/covariate.tsv \
        --phenoCol PHENO \
        --covarFile /tmp/regenie/covariate.tsv \
        --covarColList SEX,age,cohort,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
        --catCovarList cohort \
        --bsize 400 \
        --bt --firth --approx \
        --pred $WORKDIR/${GENE_SET}_${ANCESTRY}_step1_pred.list \
        --anno-file ~{anno_file} \
        --set-list ~{set_list} \
        --mask-def ~{mask_file} \
        --aaf-bins 0.01 \
        --build-mask max \
        --write-mask-snplist \
        --out $WORKDIR/${GENE_SET}_${ANCESTRY}_step2 \
        --threads ~{cpus}

    # REGENIE step 2 — single variants for MAGMA
    regenie \
        --step 2 \
        --pgen $WORKDIR/geneset_plink2 \
        --phenoFile /tmp/regenie/covariate.tsv \
        --phenoCol PHENO \
        --covarFile /tmp/regenie/covariate.tsv \
        --covarColList SEX,age,cohort,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
        --catCovarList cohort \
        --bsize 400 \
        --bt --firth --approx \
        --pred $WORKDIR/${GENE_SET}_${ANCESTRY}_step1_pred.list \
        --out $WORKDIR/${GENE_SET}_${ANCESTRY}_step2_singlevariants \
        --threads ~{cpus}

    echo $N_SAMPLES > $WORKDIR/n_samples.txt

    # Copy to bucket
    gsutil -m cp $WORKDIR/*.regenie $WORKDIR/*.log $WORKDIR/n_samples.txt ~{out_dir}/
    echo "REGENIE complete for ${ANCESTRY}"
  >>>

  output {
    File mask_results = "~{out_dir}/~{gene_set}_~{ancestry}_step2_PHENO.regenie"
    File sv_results   = "~{out_dir}/~{gene_set}_~{ancestry}_step2_singlevariants_PHENO.regenie"
    Int  n_samples    = read_int("~{out_dir}/n_samples.txt")
  }

  runtime {
    docker: "ghcr.io/rgcgithub/regenie/regenie:v3.4.1.gz"
    cpu: cpus
    memory: "~{mem_gb} GB"
    disks: "local-disk 200 ~{disk_type}"
    preemptible: 1
  }
}

# ─────────────────────────────────────────────────────────────────────────────
# Task: 06 — MAGMA gene-set analysis
# ─────────────────────────────────────────────────────────────────────────────

task GeneSetAnalysis {
  input {
    String ancestry
    String gene_set
    String workspace_bucket
    File   regenie_sv_out
    File   genome_wide_bed
    File   genome_wide_bim
    File   genome_wide_fam
    File   gene_set_csv
    File   magma_gene_loc
    Int    n_samples
    Int    cpus
    Int    mem_gb
    String disk_type
  }

  String out_dir = "~{workspace_bucket}/results/~{gene_set}/~{ancestry}/06_geneset"

  command <<<
    set -euo pipefail

    ANCESTRY="~{ancestry}"
    GENE_SET="~{gene_set}"
    WORKDIR="/tmp/magma"
    mkdir -p $WORKDIR

    BED_PREFIX=$(echo "~{genome_wide_bed}" | sed 's/.bed//')

    # Build SNP file from REGENIE single variant output
    python3 << 'PYEOF'
import pandas as pd

reg = pd.read_csv("~{regenie_sv_out}", sep=' ', comment='#')
# Filter to single variants only
if 'TEST' in reg.columns:
    reg = reg[reg['TEST'] == 'ADD']
else:
    reg = reg[reg['ID'].str.contains(':', na=False)]

reg['P'] = 10 ** (-reg['LOG10P'].astype(float))
reg['CHR'] = reg['CHROM'].astype(str).str.replace('chr','',regex=False)
reg = reg[reg['P'].between(0, 1, inclusive='neither')]

snp_df = pd.DataFrame({
    'SNP': reg['ID'],
    'CHR': reg['CHR'],
    'BP' : reg['GENPOS'].astype(int),
    'P'  : reg['P'],
    'N'  : ~{n_samples}
})
snp_df.to_csv('/tmp/magma/snps.txt', sep=' ', index=False)
print(f"SNP file: {len(snp_df)} variants")

# Build gene-set file
gs = pd.read_csv("~{gene_set_csv}")
gs_clean = gs.dropna(subset=['entrez_id']).copy()
gs_clean['entrez_id'] = gs_clean['entrez_id'].astype(int).astype(str)

with open('/tmp/magma/genesets.txt','w') as f:
    all_ids = gs_clean['entrez_id'].unique().tolist()
    f.write(f"purine_metabolism_full\t" + '\t'.join(all_ids) + '\n')
    for subset, grp in gs_clean.groupby('gene_set'):
        ids = grp['entrez_id'].unique().tolist()
        if len(ids) >= 5:
            f.write(f"{subset}\t" + '\t'.join(ids) + '\n')
print("Gene-set file written")
PYEOF

    # MAGMA step 1 — annotation using genome-wide bim
    # WDL localizes bed/bim/fam separately — reconstruct prefix from bed path
    BED_PREFIX=$(echo "~{genome_wide_bed}" | sed 's/\.bed$//')

    magma \
        --annotate window=10,10 \
        --snp-loc ~{genome_wide_bim} \
        --gene-loc ~{magma_gene_loc} \
        --out $WORKDIR/magma_annot

    # MAGMA step 2 — gene analysis
    magma \
        --bfile $BED_PREFIX \
        --gene-annot $WORKDIR/magma_annot.genes.annot \
        --pval /tmp/magma/snps.txt ncol=N \
        --gene-model snp-wise=mean \
        --out $WORKDIR/magma_genes

    # MAGMA step 3 — gene-set analysis
    magma \
        --gene-results $WORKDIR/magma_genes.genes.raw \
        --set-annot /tmp/magma/genesets.txt \
        --out $WORKDIR/magma_geneset

    # Copy to bucket
    gsutil -m cp $WORKDIR/*.out $WORKDIR/*.raw $WORKDIR/*.log ~{out_dir}/
    echo "MAGMA complete for ${ANCESTRY}"
  >>>

  output {
    File gene_results    = "~{out_dir}/magma_genes.genes.out"
    File geneset_results = "~{out_dir}/magma_geneset.gsa.out"
  }

  runtime {
    docker: "gcr.io/broad-cga-aarong-gtex/magma:latest"
    cpu: cpus
    memory: "~{mem_gb} GB"
    disks: "local-disk 50 ~{disk_type}"
    preemptible: 2
  }
}
