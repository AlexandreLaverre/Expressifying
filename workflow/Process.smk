import os

configfile: 'config/config.yaml'
BURN_IN, UNTIL = config['BURN_IN'], config['UNTIL']

FOLDER = os.path.abspath('.')
bgee_ortho = f"{FOLDER}/results/gene_expression_orthogroups"
datasets = [i.replace(".csv","") for i in os.listdir(f"{bgee_ortho}/log2TPM") if
            not i.startswith(".") and os.path.isfile(f"{bgee_ortho}/log2TPM/{i}")]
print(datasets)
traits = config["TRAITS"]
oma_data = config["OMA_DATA"]
methods = config["METHODS"]
clades = config["CLADES"]
orthologs_path = f"{FOLDER}/data/gene_orthologies/one2one_orthogroups_{oma_data}_mammals.csv"
heritability_path = f"{FOLDER}/data/heritability-gene-expression/41431_2019_511_MOESM2_ESM.tsv"
mammals_polymorphism_path = f"{FOLDER}/data/science.abn5856.Ne/science.abn5856_table_s1.csv"
mammal_tree_path = f"{FOLDER}/data/science.abl8189.Timescale/Table_S2_RootedTrees/Concatenation_HRA_neutral_241_10miss_rooted.tree"
primate_tree_path = f"{FOLDER}/data/science.abn7829.primates/science.abn7829_data_s3.nw.tree"

exec_path = "/opt/homebrew/Caskroom/miniforge/base/envs/osx-64/bin/rb"
if not os.path.exists(exec_path):
    # Find executable in the path using whereis. If not found, raise an error.
    split = os.popen(f'whereis rb').read().split()
    if len(split) > 1:
        exec_path = split[1].strip()
    else:
        raise FileNotFoundError(f'rb not found. Please install RevBayes and add it to your path.')
print(f"Found rb at {exec_path}")

ruleorder: read_trace > merge_traces
ruleorder: run_RevBayes > gather_RevBayes_log

wildcard_constraints:
    # constrain to only alphanumeric characters and underscore
    trait=r"[a-zA-Z0-9]+",method=r"[a-zA-Z]+",gene=r"[a-zA-Z0-9]+",clade=r"[a-zA-Z]+"


rule all:
    input:
        expand(f"{FOLDER}/results/{{trait}}_{{clade}}/switch_RevBayes.tsv",trait=traits,clade=clades),
        expand(f"{FOLDER}/results/{{trait}}_{{clade}}/merge_{{method}}.tsv",method=methods,trait=traits,clade=clades),


rule mammal_to_primate_tree:
    input:
        script=f"{FOLDER}/scripts/mammal_to_primate_tree.py",
        primate_tree=primate_tree_path,
        mammal_neutral_tree=mammal_tree_path
    output:
        tree=f"{FOLDER}/data_processed/primates.tree"
    log:
        stdout=f"{FOLDER}/data_processed/mammal_to_primate_tree.log"
    shell:
        'python3 {input.script} --mammal_tree {input.mammal_neutral_tree} --primate_tree {input.primate_tree} --output {output.tree} > {log.stdout}'

rule mammal_to_nonprimate_tree:
    input:
        script=f"{FOLDER}/scripts/mammal_to_nonprimate_tree.py",
        primate_tree=primate_tree_path,
        mammal_neutral_tree=mammal_tree_path
    output:
        tree=f"{FOLDER}/data_processed/nonMammalPrimates.tree"
    log:
        stdout=f"{FOLDER}/data_processed/nonMammalPrimates.log"
    shell:
        'python3 {input.script} --mammal_tree {input.mammal_neutral_tree} --primate_tree {input.primate_tree} --output {output.tree} > {log.stdout}'


rule pre_processed_traits:
    input:
        script=f"{FOLDER}/scripts/pre_processed_traits.py",
        input_dS=lambda wildcards: mammal_tree_path if wildcards.clade == "mammals" else f"{FOLDER}/data_processed/{wildcards.clade}.tree",
        input_pS=mammals_polymorphism_path,
        input_traits=f"{bgee_ortho}/{{trait}}/{{dataset}}.csv",
        input_orthogroups=orthologs_path,
        input_heritability=heritability_path
    output:
        tree=f"{FOLDER}/data_processed/{{clade}}/{{dataset}}_{{trait}}.tree",
        traits=f"{FOLDER}/data_processed/{{clade}}/{{dataset}}_{{trait}}.traits.tsv",
        var_within=f"{FOLDER}/data_processed/{{clade}}/{{dataset}}_{{trait}}.var_within.tsv"
    log:
        stdout=f"{FOLDER}/data_processed/{{clade}}/{{dataset}}_{{trait}}.preprocessed.log"
    params:
        heritability=lambda wildcards, input, output: f"--input_heritability {input.input_heritability} " if wildcards.trait in [
            "logTPM", "log2TPM"] else ""
    shell:
        'python3 {input.script} --input_dS {input.input_dS} --input_pS {input.input_pS} '
        '--input_traits {input.input_traits} --input_orthogroups {input.input_orthogroups} {params.heritability}'
        '--output_tree {output.tree} --output_traits {output.traits} --output_var_within {output.var_within} > {log.stdout}'

rule neutrality_index:
    input:
        tree=rules.pre_processed_traits.output.tree,
        traits=rules.pre_processed_traits.output.traits,
        var_within=rules.pre_processed_traits.output.var_within
    output:
        tsv=f"{FOLDER}/data_processed/LnL_{{clade}}/{{dataset}}_{{trait}}.tsv"
    params:
        script=f"{FOLDER}/scripts/neutrality_index.py",
    log:
        stdout=f"{FOLDER}/data_processed/LnL_{{clade}}/{{dataset}}_{{trait}}.processed.log"
    shell:
        'python3 {params.script} --tree {input.tree} --traits {input.traits} --var_within {input.var_within}'
        ' --output {output.tsv} > {log.stdout}'


checkpoint filter_traits:
    input:
        script=f"{FOLDER}/scripts/filter_traits.py",
        tree=rules.pre_processed_traits.output.tree,
        traits=rules.pre_processed_traits.output.traits,
        var_within=rules.pre_processed_traits.output.var_within,
        neutrality_index=rules.neutrality_index.output.tsv
    output:
        dir=directory(f"{FOLDER}/data_processed/filtered_{{clade}}/{{dataset}}_{{trait}}")
    log:
        stdout=f"{FOLDER}/data_processed/filtered_{{clade}}/{{dataset}}_{{trait}}.preprocessed.log"
    shell:
        'python3 {input.script} --input_tree {input.tree} --input_traits {input.traits} --input_var_within {input.var_within} --neutrality_index {input.neutrality_index}'
        ' --output_dir {output.dir} > {log.stdout}'


rule reconstructed_chronogram:
    input:
        script=f"{FOLDER}/scripts/calib_pl.R",
        tree=f"{FOLDER}/data_processed/filtered_{{clade}}/{{dataset}}_{{trait}}/{{gene}}.tree"
    output:
        tree=f"{FOLDER}/data_processed/filtered_{{clade}}/{{dataset}}_{{trait}}/{{gene}}.timetree.newick"
    shell:
        'Rscript {input.script} {input.tree} {output.tree}'


rule convert_to_RevBayes:
    input:
        script=f"{FOLDER}/scripts/convert_to_nexus.py",
        nuc_tree=f"{FOLDER}/data_processed/filtered_{{clade}}/{{dataset}}_{{trait}}/{{gene}}.tree",
        time_tree=rules.reconstructed_chronogram.output.tree,
        traits=f"{FOLDER}/data_processed/filtered_{{clade}}/{{dataset}}_{{trait}}/{{gene}}.traits.tsv"
    output:
        nuc_tree=f"{FOLDER}/data_processed/filtered_{{clade}}/{{dataset}}_{{trait}}/{{gene}}_nuctree.nex",
        time_tree=f"{FOLDER}/data_processed/filtered_{{clade}}/{{dataset}}_{{trait}}/{{gene}}_timetree.nex",
        nexus_traits=f"{FOLDER}/data_processed/filtered_{{clade}}/{{dataset}}_{{trait}}/{{gene}}.traits.nex"
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
        log=f"{FOLDER}/data_processed/RevBayes_{{clade}}/{{dataset}}_{{trait}}/{{gene}}/simple_BM_Switchnodes.log.gz"
    params:
        folder=lambda wildcards: f"{FOLDER}/data_processed/RevBayes_{wildcards.clade}/{wildcards.dataset}_{wildcards.trait}/{wildcards.gene}"
    shell:
        'mkdir -p {params.folder};'
        'cp {input.timetree} {params.folder}/tree_time.nex;'
        'cp {input.nuctree} {params.folder}/tree_nuc.nex;'
        'cp {input.traits} {params.folder}/traits.nex;'
        'cd {params.folder} && {input.exec} {input.rev_file};'
        'for f in ./*.log; do if [ -f $f ]; then gzip -f $f; fi; done;'


def aggregate_genes(wildcards):
    checkpoint_output = checkpoints.filter_traits.get(**wildcards).output[0]
    return expand(f"{FOLDER}/data_processed/RevBayes_{{clade}}/{{dataset}}_{{trait}}/{{gene}}/simple_BM_Switchnodes.log.gz",clade=wildcards.clade,dataset=wildcards.dataset,trait=wildcards.trait,gene=glob_wildcards(os.path.join(checkpoint_output,"{gene}.tree")).gene)


rule gather_RevBayes_log:
    input: aggregate_genes
    output:
        tsv=f"{FOLDER}/data_processed/Switchnodes/{{clade}}_{{dataset}}_{{trait}}.tsv"
    params:
        folder=f"{FOLDER}/data_processed/RevBayes_{{clade}}/{{dataset}}_{{trait}}"
    shell:
        'python3 {FOLDER}/scripts/gather_RevBayes_log.py --folder {params.folder} --output_tsv {output.tsv}'


rule merge_RevBayes_results:
    input:
        script=f"{FOLDER}/scripts/merge_RevBayes_results.py",
        tsv_results=expand(f"{FOLDER}/data_processed/Switchnodes/{{{{clade}}}}_{{dataset}}_{{{{trait}}}}.tsv",dataset=datasets)
    output:
        tsv=f"{FOLDER}/results/{{trait}}_{{clade}}/switch_RevBayes.tsv"
    shell:
        'python3 {input.script} --tsv_results {input.tsv_results} --output {output.tsv}'


rule run_inference:
    input:
        exec=f"{FOLDER}/utils/BayesCode/bin/nodetraits",
        traits=f"{FOLDER}/data_processed/filtered_{{clade}}/{{dataset}}_{{trait}}/{{gene}}.traits.tsv",
        tree=f"{FOLDER}/data_processed/filtered_{{clade}}/{{dataset}}_{{trait}}/{{gene}}.tree"
    output:
        run=f"{FOLDER}/data_processed/inference_{{clade}}/{{dataset}}_{{trait}}/{{gene}}.run"
    params:
        chain=f"{FOLDER}/data_processed/inference_{{clade}}/{{dataset}}_{{trait}}/{{gene}}",
        until=f"--until {UNTIL}"
    log:
        stdout=f"{FOLDER}/data_processed/inference_{{clade}}/{{dataset}}_{{trait}}/{{gene}}.log"
    shell:
        '{input.exec} {params.until} --uniq_kappa --df 1 --tree {input.tree} --traitsfile {input.traits} {params.chain} &> {log.stdout}'


rule read_trace:
    input:
        exec=f"{FOLDER}/utils/BayesCode/bin/readnodetraits",
        inference=f"{FOLDER}/data_processed/inference_{{clade}}/{{dataset}}_{{trait}}/{{gene}}.run",
        var_within=f"{FOLDER}/data_processed/filtered_{{clade}}/{{dataset}}_{{trait}}/{{gene}}.var_within.tsv"
    output:
        tsv=f"{FOLDER}/data_processed/Bayes_{{clade}}/{{dataset}}_{{trait}}/{{gene}}.tsv"
    params:
        points=f"--burnin {BURN_IN} --until {UNTIL}",
        chain=lambda wildcards: f"{FOLDER}/data_processed/inference_{wildcards.clade}/{wildcards.dataset}_{wildcards.trait}/{wildcards.gene}"
    log:
        stdout=f"{FOLDER}/data_processed/Bayes_{{clade}}/{{dataset}}_{{trait}}/{{gene}}.trace.log"
    shell:
        '{input.exec} {params.points} --var_within {input.var_within} --output {output.tsv} {params.chain} &> {log.stdout} && gzip -f {params.chain}.chain && gzip -f {params.chain}.trace'


def aggregate_input(wildcards):
    checkpoint_output = checkpoints.filter_traits.get(**wildcards).output[0]
    return expand(f"{FOLDER}/data_processed/Bayes_{{clade}}/{{dataset}}_{{trait}}/{{gene}}.tsv",clade=wildcards.clade,dataset=wildcards.dataset,trait=wildcards.trait,gene=glob_wildcards(os.path.join(checkpoint_output,"{gene}.tree")).gene)


rule merge_traces:
    input: aggregate_input
    output:
        tsv=f"{FOLDER}/data_processed/Bayes_{{clade}}/{{dataset}}_{{trait}}.tsv"
    params:
        traces=f"{FOLDER}/data_processed/Bayes_{{clade}}/{{dataset}}_{{trait}}"
    shell:
        'python3 {FOLDER}/scripts/merge_traces.py --input {params.traces} --output {output.tsv}'


rule merge_results:
    input:
        script=f"{FOLDER}/scripts/merge_results.py",
        tsv_results=expand(f"{FOLDER}/data_processed/{{{{method}}}}_{{{{clade}}}}/{{dataset}}_{{{{trait}}}}.tsv",dataset=datasets)
    output:
        tsv=f"{FOLDER}/results/{{trait}}_{{clade}}/merge_{{method}}.tsv"
    shell:
        'python3 {input.script} --tsv_results {input.tsv_results} --output {output.tsv}'
