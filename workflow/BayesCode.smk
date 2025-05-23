import os

FOLDER = os.path.abspath('../../..')

dataset = config["dataset"]
trait = config["trait"]
clade = config["clade"]
BURN_IN = config['BURN_IN']
UNTIL = config['UNTIL']

GENE_LIST = list(set([i.split(".")[0] for i in os.listdir(f"{FOLDER}/data_processed/filtered_{clade}/{dataset}_{trait}") if i.endswith(".tree") and not i.startswith(".")]))
print(f"Found {len(GENE_LIST)} genes in {clade} for {dataset} ({trait})")

rule all:
    input:
        f"{FOLDER}/data_processed/BayesCode_{clade}/{dataset}_{trait}.tsv"

rule run_inference:
    input:
        exec=f"{FOLDER}/utils/BayesCode/bin/nodetraits",
        traits=f"{FOLDER}/data_processed/filtered_{clade}/{dataset}_{trait}/{{gene}}.traits.tsv",
        tree=f"{FOLDER}/data_processed/filtered_{clade}/{dataset}_{trait}/{{gene}}.tree"
    output:
        run=f"{FOLDER}/data_processed/inference_{clade}/{dataset}_{trait}/{{gene}}.run"
    params:
        chain=f"{FOLDER}/data_processed/inference_{clade}/{dataset}_{trait}/{{gene}}",
        until=f"--until {UNTIL}"
    log:
        stdout=f"{FOLDER}/data_processed/inference_{clade}/{dataset}_{trait}/{{gene}}.log"
    shell:
        '{input.exec} {params.until} --uniq_kappa --df 1 --tree {input.tree} --traitsfile {input.traits} {params.chain} &> {log.stdout}'


rule read_trace:
    input:
        exec=f"{FOLDER}/utils/BayesCode/bin/readnodetraits",
        inference=f"{FOLDER}/data_processed/inference_{clade}/{dataset}_{trait}/{{gene}}.run",
        var_within=f"{FOLDER}/data_processed/filtered_{clade}/{dataset}_{trait}/{{gene}}.var_within.tsv"
    output:
        tsv=f"{FOLDER}/data_processed/BayesCode_{clade}/{dataset}_{trait}/{{gene}}.tsv"
    params:
        points=f"--burnin {BURN_IN} --until {UNTIL}",
        chain=lambda wildcards: f"{FOLDER}/data_processed/inference_{clade}/{dataset}_{trait}/{wildcards.gene}"
    log:
        stdout=f"{FOLDER}/data_processed/BayesCode_{clade}/{dataset}_{trait}/{{gene}}.trace.log"
    shell:
        '{input.exec} {params.points} --var_within {input.var_within} --output {output.tsv} {params.chain} &> {log.stdout} && gzip -f {params.chain}.chain && gzip -f {params.chain}.trace'


rule merge_traces:
    input:
        genes=expand(f"{FOLDER}/data_processed/BayesCode_{clade}/{dataset}_{trait}/{{gene}}.tsv", gene=GENE_LIST),
        script=f"{FOLDER}/scripts/merge_traces.py"
    output:
        tsv=f"{FOLDER}/data_processed/BayesCode_{clade}/{dataset}_{trait}.tsv"
    params:
        traces=f"{FOLDER}/data_processed/BayesCode_{clade}/{dataset}_{trait}"
    shell:
        'python3 {input.script} --input {params.traces} --output {output.tsv}'