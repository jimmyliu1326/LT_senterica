import pandas as pd

# override working directory
workdir: config['output_dir']

# load samples metadata
samples_table = config["samples"]
samples_meta = pd.read_csv(samples_table).set_index("Sample", drop=False)

def pipeline_output():
    if config['tree'] == 0:
        return(expand("neighbours/{samples}.list", samples=samples_meta.Sample))
    else:
        return("tree.nwk")

def input_assemblies(wildcards):
    samples = wildcards["samples"]
    if samples_meta.loc[samples,'Type'] == "raw":
        return "assemblies/{samples}.fasta"
    else:
        return samples_meta.R1[wildcards.samples]

rule all:
    input:
        pipeline_output()

rule trim:
    input:
        R1_in=lambda wildcards: samples.R1[wildcards.samples],
        R2_in=lambda wildcards: samples.R2[wildcards.samples]
    output:
        R1_out="fastp_res/{samples}_R1_trimmed.fastq",
        R2_out="fastp_res/{samples}_R2_trimmed.fastq",
        json="fastp_res/{samples}.json",
        html="fastp_res/{samples}.html"
    threads: 16
    log:
    "logs/fastp/{samples}.log"
    shell:
        """
        fastp -i {input.R1_in} -I {input.R2_in} -o {output.R1_out} -O {output.R2_out} \
            -w {threads} -j {output.json} -h {output.html}
        """

rule assembly:
    input:
        R1_in="fastp_res/{samples}_R1_trimmed.fastq",
        R2_in="fastp_res/{samples}_R2_trimmed.fastq"
    params:
        outdir="shovill_res/{samples}"
    output:
        assembly="assemblies/{samples}.fasta"
    threads: 30
    log:
        "logs/shovill/{samples}.log"
    shell:
        """
        shovill --outdir {params.outdir} --R1 {input.R1_in} --R2 {input.R2_in} \
            --gsize 4.5M --cpus {threads} --force
        mkdir -p assemblies
        ln -s $(realpath {params.outdir})/contigs.fa {output.assembly}
        """

rule sketch:
    input:
        assemblies=input_assemblies
    output:
        sketches="sketches/{samples}.{k}.msh"
    params:
        k="{k}"
    threads: 1
    log:
        "logs/sketches/{samples}.{k}.log"
    shell:
        "mash sketch -k {params.k} -s 10000 -o {output.sketches} {input.assemblies}"

rule mash_dist:
    input:
        sketches="sketches/{samples}.{k}.msh",
        reference=config['pipeline_dir']+"/database/Senterica_ref.{k}.msh"
    output:
        distances="distances/{samples}.{k}.tsv"
    threads: 8
    log:
        "log/mash_distance/{samples}.{k}.log"
    shell:
        "mash dist -p {threads} {input.reference} {input.sketches} > {output.distances}"

rule amrfinder:
    input:
        assemblies=input_assemblies
    params:
        sample_name="{samples}"
    output:
        amr_res="amrfinder/{samples}.tsv"
    threads: 4
    log:
        "log/amrfinder/{samples}.log"
    shell:
        """
        amrfinder --threads {threads} -O Salmonella -n {input.assemblies} \
            -o {output.amr_res} --name {params.sample_name}
        """

rule vfdb:
    input:
        assemblies=input_assemblies
    output:
        vfdb_res="vfdb/{samples}.tsv"
    threads: 4
    log:
        "log/vfdb/{samples}.log"
    shell:
        """
        abricate --db vfdb --nopath --threads {threads} \
            {input.assemblies} > {output.vfdb_res}
        """

rule plasmidfinder:
    input:
        assemblies=input_assemblies
    output:
        plasmidfinder_res="plasmidfinder/{samples}.tsv"
    threads: 4
    log:
        "log/plasmidfinder/{samples}.log"
    shell:
        """
        abricate --db plasmidfinder --nopath --threads {threads} \
            {input.assemblies} > {output.plasmidfinder_res}
        """

rule molecular_linkage:
    input:
        distances=expand("distances/{sample}.{k}.tsv", k=[13,17,21,25,29], sample=samples_meta.Sample),
        plasmidfinder_res="plasmidfinder/{samples}.tsv",
        amrfinder_res="amrfinder/{samples}.tsv",
        vfdb_res="vfdb/{samples}.tsv"
    params:
        R_script=config["pipeline_dir"]+"/src/MLFinder.R",
        sample_name="{samples}",
        distances_directory="distances",
        output_directory=config['output_dir'],
        annotations=config['annotations']
    output:
        neighbours="neighbours/{samples}.list",
        references="references/{samples}.list",
        cluster_plot="plots/{samples}.png"
    threads: 16
    log:
        "log/molecular_linkage/{samples}.log"
    shell:
        """
        {params.R_script} {params.distances_directory} {params.sample_name} {params.output_directory} \
            {threads} {input.amrfinder_res} {input.vfdb_res} {input.plasmidfinder_res} {params.annotations}
        """

rule consolidate:
    input:
       neighbours=expand("neighbours/{samples}.list", samples=samples_meta.Sample),
       references=expand("references/{samples}.list", samples=samples_meta.Sample)
    output:
        all_neighbours="all_neighbours.list",
        all_references="all_references.list",
        all_genomes="all_genomes.list"
    threads: 4
    shell:
        """
        cat {input.neighbours} {input.references} | sort | uniq > {output.all_genomes}
        cat {input.neighbours} | sort | uniq > {output.all_neighbours}
        cat {input.references} | sort | uniq > {output.all_references}
        """

rule download_genomes:
    input:
        all_genomes="all_genomes.list"
    params:
        asm_acc="{genomes}"
    output:
        genomes="public_genomes/{genomes}.fa"
    threads: 4
    shell:
        """
        filename=$(cat {input.all_genomes} | grep "^{params.asm_acc}" | cut -d',' -f1)
        source=$(echo $filename | cut -d'_' -f1)
        url=$(cat {input.all_genomes} | grep "^{params.asm_acc}" | cut -d',' -f2)

        if [[ $source == "GCA" ]]; then
            var=$(echo $url | sed 's/.*GCA_//')
            wget -q -nc $url"/GCA_"${var}"_genomic.fna.gz" -O {output.genomes}.gz
            #echo $url"/GCA_"${var}"_genomic.fna.gz"
            # check filesize
            size=$(stat -c %s {output.genomes}.gz)

            # retrieve again if failed
            while [[ $size -eq 0 ]]; do
                rm {output.genomes}gz
                wget -q -nc $url"/GCA_"$var"_genomic.fna.gz" -O {output.genomes}.gz
                size=$(stat -c %s {output.genomes})
                sleep 2s
            done
        
            # unzip file
            gunzip {output.genomes}
        else
            wget -q -nc -c $url -O $OUT_DIR/LT_Senterica_tmp/genomes/$filename.fa

            # check filesize
            size=$(stat -c %s $OUT_DIR/LT_Senterica_tmp/genomes/$filename.fa)

            # retrieve again if failed
            while [[ $size -eq 0 ]]; do
                rm $OUT_DIR/LT_Senterica_tmp/genomes/$filename.fa
                wget -q -nc -c $url -O $OUT_DIR/LT_Senterica_tmp/genomes/$filename.fa
                size=$(stat -c %s $OUT_DIR/LT_Senterica_tmp/genomes/$filename.fa)
                sleep 2s
            done
        fi
        """

rule create_phame_query_input:
    input:
        query_assemblies=input_assemblies
    output:
        phame_query_input="phame/work_dir/{samples}.fa"
    threads: 1
    shell:
        "ln -s $(realpath {input.query_assemblies}) {output.phame_query_input}"

rule create_phame_public_input:
    input:
        public_assemblies="public_genomes/{genomes}.fa"
    output:
        phame_public_input="phame/work_dir/{genomes}.fa"
    threads: 1
    shell:
        "ln -s $(realpath {input.public_assemblies}) {output.phame_public_input}"

rule create_phame_reference:
    input:
        all_references="all_references.list"
    output:
        reference="phame/ref_dir/reference.fa"
    threads: 1
    shell:
        """
        reference=$(cat {input.all_references} | head -n 1 | cut -f1 -d',')
        mv public_genomes/$reference.fa {output.reference}
        """

#def input_genomes_for_phame():
#    df=pd.read_csv(checkpoints.consolidate.get().output[0], header = None)
#    return df[0].tolist()

rule phame:
    input:
        query_assemblies=expand("phame/work_dir/{samples}.fa", samples=samples_meta.Sample),
        #public_assemblies=expand("phame/work_dir/{genomes}.fa", genomes=input_genomes_for_phame()),
        reference="phame/ref_dir/reference.fa"
    params:
        ctl_file="phame.ctl",
        create_ctl_script=config['pipeline_dir']+"/src/generatePhameCtl.sh",
        phame_input_dir="phame/work_dir",
        phame_ref_dir="phame/ref_dir",
        phame_ref_file="reference.fa"
    output:
        phame_res=directory("phame/work_dir/results"),
        tree="tree.nwk"
    threads:config['threads']
    shell:
        """
        {params.create_ctl_script} {params.phame_input_dir} {params.phame_ref_dir} {params.phame_ref_file} {threads} .
        phame {params.ctl_file}
        cp {output.phame_res}/trees/*.fasttree {output.tree}
        """