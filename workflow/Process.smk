import os

configfile: 'config/config.yaml'
BURN_IN, UNTIL = config['BURN_IN'], config['UNTIL']

FOLDER = os.path.abspath('.')
SNAKERUN = config["SNAKERUN"] if "SNAKERUN" in config else ""
bgee_ortho = f"{FOLDER}/data/gene_expression_orthogroups_15_2"
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

wildcard_constraints:
    # constrain to only alphanumeric characters and underscore
    trait=r"[a-zA-Z0-9]+",method=r"[a-zA-Z]+",clade=r"[a-zA-Z]+"

rule all:
    input:
        expand(f"{FOLDER}/data_merged/{{trait}}_{{clade}}/switch_RevBayes.tsv",trait=traits,clade=clades),
        expand(f"{FOLDER}/data_merged/{{trait}}_{{clade}}/merge_{{method}}.tsv",method=methods,trait=traits,clade=clades),

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
        tree=f"{FOLDER}/data_processed/NoPrimates.tree"
    log:
        stdout=f"{FOLDER}/data_processed/NoPrimates.log"
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

rule filter_traits:
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

rule run_RevBayes:
    input:
        smk=f"{FOLDER}/workflow/RevBayes.smk",
        dir=rules.filter_traits.output.dir,
    output:
        tsv=f"{FOLDER}/data_processed/Switchnodes/{{clade}}_{{dataset}}_{{trait}}.tsv"
    threads: workflow.cores
    shell:
        "cd {input.dir}; snakemake {SNAKERUN} --rerun-incomplete -k -s {input.smk} -j {threads} --config clade={wildcards.clade} dataset={wildcards.dataset} trait={wildcards.trait}"

rule merge_RevBayes_results:
    input:
        script=f"{FOLDER}/scripts/merge_RevBayes_results.py",
        tsv_results=expand(f"{FOLDER}/data_processed/Switchnodes/{{{{clade}}}}_{{dataset}}_{{{{trait}}}}.tsv",dataset=datasets)
    output:
        tsv=f"{FOLDER}/data_merged/{{trait}}_{{clade}}/switch_RevBayes.tsv"
    shell:
        'python3 {input.script} --tsv_results {input.tsv_results} --output {output.tsv}'

rule run_BayesCode:
    input:
        smk=f"{FOLDER}/workflow/BayesCode.smk",
        dir=rules.filter_traits.output.dir,
    output:
        tsv=f"{FOLDER}/data_processed/BayesCode_{{clade}}/{{dataset}}_{{trait}}.tsv"
    threads: workflow.cores
    shell:
        "cd {input.dir}; snakemake {SNAKERUN} --rerun-incomplete  -k -s {input.smk} -j {threads} --config clade={wildcards.clade} dataset={wildcards.dataset} trait={wildcards.trait} UNTIL={UNTIL} BURN_IN={BURN_IN}"

rule merge_results:
    input:
        script=f"{FOLDER}/scripts/merge_results.py",
        tsv_results=expand(f"{FOLDER}/data_processed/{{{{method}}}}_{{{{clade}}}}/{{dataset}}_{{{{trait}}}}.tsv",dataset=datasets)
    output:
        tsv=f"{FOLDER}/data_merged/{{trait}}_{{clade}}/merge_{{method}}.tsv"
    shell:
        'python3 {input.script} --tsv_results {input.tsv_results} --output {output.tsv}'
