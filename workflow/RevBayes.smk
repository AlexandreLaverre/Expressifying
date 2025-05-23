import os

FOLDER = os.path.abspath('../../..')

dataset = config["dataset"]
trait = config["trait"]
clade = config["clade"]

exec_path = "/opt/homebrew/Caskroom/miniforge/base/envs/osx-64/bin/rb"
if not os.path.exists(exec_path):
    # Find executable in the path using whereis. If not found, raise an error.
    split = os.popen(f'whereis rb').read().split()
    if len(split) > 1:
        exec_path = split[1].strip()
    else:
        raise FileNotFoundError(f'rb not found. Please install RevBayes and add it to your path.')
print(f"Found rb at {exec_path}")

GENE_LIST = list(set([i.split(".")[0] for i in os.listdir(f"{FOLDER}/data_processed/filtered_{clade}/{dataset}_{trait}") if i.endswith(".tree") and not i.startswith(".")]))
print(f"Found {len(GENE_LIST)} genes in {clade} for {dataset} ({trait})")

rule all:
    input:
        f"{FOLDER}/data_processed/Switchnodes/{clade}_{dataset}_{trait}.tsv"

rule reconstructed_chronogram:
    input:
        script=f"{FOLDER}/scripts/calib_pl.R",
        tree=f"{FOLDER}/data_processed/filtered_{clade}/{dataset}_{trait}/{{gene}}.tree"
    output:
        tree=f"{FOLDER}/data_processed/filtered_{clade}/{dataset}_{trait}/{{gene}}.timetree.newick"
    shell:
        'Rscript {input.script} {input.tree} {output.tree}'


rule convert_to_RevBayes:
    input:
        script=f"{FOLDER}/scripts/convert_to_nexus.py",
        nuc_tree=f"{FOLDER}/data_processed/filtered_{clade}/{dataset}_{trait}/{{gene}}.tree",
        time_tree=rules.reconstructed_chronogram.output.tree,
        traits=f"{FOLDER}/data_processed/filtered_{clade}/{dataset}_{trait}/{{gene}}.traits.tsv"
    output:
        nuc_tree=f"{FOLDER}/data_processed/filtered_{clade}/{dataset}_{trait}/{{gene}}_nuctree.nex",
        time_tree=f"{FOLDER}/data_processed/filtered_{clade}/{dataset}_{trait}/{{gene}}_timetree.nex",
        nexus_traits=f"{FOLDER}/data_processed/filtered_{clade}/{dataset}_{trait}/{{gene}}.traits.nex"
    shell:
        'python3 {input.script} --input_nuctree {input.nuc_tree} --input_timetree {input.time_tree} --input_traits {input.traits} --output_traits {output.nexus_traits} --output_nuctree {output.nuc_tree} --output_timetree {output.time_tree}'


rule run_RevBayes:
    input:
        exec=exec_path,
        rev_file=f"{FOLDER}/scripts/mcmc_simple_BM_Switchnodes.Rev",
        nuctree=rules.convert_to_RevBayes.output.nuc_tree,
        timetree=rules.convert_to_RevBayes.output.time_tree,
        traits=rules.convert_to_RevBayes.output.nexus_traits
    output:
        log=f"{FOLDER}/data_processed/RevBayes_{clade}/{dataset}_{trait}/{{gene}}/simple_BM_Switchnodes.log.gz"
    params:
        folder=lambda wildcards: f"{FOLDER}/data_processed/RevBayes_{clade}/{dataset}_{trait}/{wildcards.gene}"
    shell:
        'mkdir -p {params.folder};'
        'cp {input.timetree} {params.folder}/tree_time.nex;'
        'cp {input.nuctree} {params.folder}/tree_nuc.nex;'
        'cp {input.traits} {params.folder}/traits.nex;'
        'cd {params.folder} && {input.exec} {input.rev_file};'
        'for f in ./*.log; do if [ -f $f ]; then gzip -f $f; fi; done;'


rule gather_RevBayes_log:
    input:
        genes=expand(f"{FOLDER}/data_processed/RevBayes_{clade}/{dataset}_{trait}/{{gene}}/simple_BM_Switchnodes.log.gz", gene=GENE_LIST),
        script=f"{FOLDER}/scripts/gather_RevBayes_log.py"
    output:
        tsv=f"{FOLDER}/data_processed/Switchnodes/{clade}_{dataset}_{trait}.tsv"
    params:
        folder=f"{FOLDER}/data_processed/RevBayes_{clade}/{dataset}_{trait}"
    shell:
        'python3 {input.script} --folder {params.folder} --output_tsv {output.tsv}'