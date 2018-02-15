import re
def get_sample_input(wildcards):
    units = pd.read_table(config["units"], index_col=["sample"], dtype=str).sort_index(level=0)
    bam_list = []
    if units.loc[[wildcards.sample]]['unit'].notnull().values.any():
        # has units
        bam_list.append("bam/merge/{sample}.bam".format(sample=wildcards.sample))
    else:
        # no units
        bam_list.append("bam/{sample}.bam".format(sample=wildcards.sample))
    return bam_list

def get_chip_input(wildcards):
    units = pd.read_table(config["units"], index_col=["sample"], dtype=str).sort_index(level=0)
    bam_list = []
    if units.loc[[wildcards.sample]]['unit'].notnull().values.any():
        # has units
        num = re.search(r"([\d+])", units.loc[[wildcards.sample]]['type'].values[0])
        bam_list.append("bam/merge/{sample}.bam".format(sample=units[units[units['unit'].isnull()]['type'].str.contains("input"+num.group(0))].reset_index()['sample'].values[0]))
    else:
        # no units
        num = re.search(r"([\d+])", units.loc[[wildcards.sample]]['type'].values[0])
        bam_list.append("bam/{sample}.bam".format(sample=units[units[units['unit'].isnull()]['type'].str.contains("input"+num.group(0))].reset_index()['sample'].values[0]))
    return bam_list

def get_chip_psuerep_input(wildcards):
    units = pd.read_table(config["units"], index_col=["sample"], dtype=str).sort_index(level=0)
    bam_list = []
    num = re.search(r"([\d+])", units.loc[[wildcards.sample]]['type'].values[0])
    bam_list.append("dedup/split/{sample}-{pseurep}.bed".format(sample=units[units[units['unit'].isnull()]['type'].str.contains("input"+num.group(0))].reset_index()['sample'].values[0], pseurep=wildcards.pseurep))
    return bam_list


def get_bdg_input(wildcards):
    units = pd.read_table(config["units"], index_col=["sample"], dtype=str).sort_index(level=0)
    bam_list = []
    if units.loc[[wildcards.sample]]['unit'].notnull().values.any():
        # has units
        bam_list.append("units/{sample}-{unit}_treat_pileup.bdg".format(sample=wildcards.sample))
    else:
        # no units
        bam_list.append("dedup/{sample}.bam".format(sample=wildcards.sample))
    return bam_list

rule call_peaks_sample:
    input:
        treat = get_sample_input,
        control = get_chip_input
    output:
        peak = "macs/{sample}_peaks.narrowPeak",
        bdg = "macs/{sample}_treat_pileup.bdg",
        ctrl = "macs/{sample}_control_lambda.bdg"
    log:
        "logs/macs2/call_peaks/{sample}.log"
    conda:
        "../envs/macs2.yaml"
    params:
        name = "{sample}",
        extra = "{}".format(config["params"]["macs2"])
    shell:
        """
        macs2 callpeak -t {input.treat} -c {input.control} -f BAMPE -g {params.extra} --outdir macs -n {params.name} --nomodel --keep-dup all -B --call-summits -q 0.01 &> {log}
        """

rule call_peaks_units:
    input:
        treat = "bam/units/{sample}-{unit}.bam",
        control = get_chip_input
    output:
        peak = "macs/units/{sample}-{unit}_peaks.narrowPeak",
        bdg = "macs/units/{sample}-{unit}_treat_pileup.bdg",
        ctrl = "macs/units/{sample}-{unit}_control_lambda.bdg"
    log: "logs/macs2/call_peaks/{sample}-{unit}.log"
    conda: "../envs/macs2.yaml"
    params:
        name = "{sample}-{unit}",
        extra = "{}".format(config["params"]["macs2"])
    shell:
        """
        macs2 callpeak -t {input.treat} -c {input.control} -f BAMPE -g {params.extra} --outdir macs/units -n {params.name} --nomodel --keep-dup all -B --call-summits -p 0.01 &> {log}
        """

rule call_peaks_pseudoreps:
    input:
        treat = "dedup/split/{sample}-{pseurep}.bed",
        control = get_chip_psuerep_input
    output:
        peak = "macs/pseudoreps/{sample}-{pseurep}_peaks.narrowPeak",
        bdg = "macs/pseudoreps/{sample}-{pseurep}_treat_pileup.bdg",
        ctrl = "macs/pseudoreps/{sample}-{pseurep}_control_lambda.bdg",
        summits = "macs/pseudoreps/{sample}-{pseurep}_summits.bed",
        xls = "macs/pseudoreps/{sample}-{pseurep}_peaks.xls"
    log: "logs/macs2/call_peaks/{sample}-{pseurep}.log"
    conda: "../envs/macs2.yaml"
    params:
        name = "{sample}-{pseurep}",
        extra = "{}".format(config["params"]["macs2"])
    shell:
        """
        macs2 callpeak -t {input.treat} -c {input.control} -f BEDPE -g {params.extra} --outdir macs/pseudoreps -n {params.name} --nomodel --keep-dup all -B --call-summits -q 0.01 &> {log}
        """

rule macs_bigwig_sample:
    input: "macs/{sample}_treat_pileup.bdg"
    output: "bw/macs2/{sample}_treat_pileup.bw"
    log:
        "logs/macs2/bw/{sample}.log"
    conda:
        "../envs/macs2.yaml"
    params:
        genome = "{}".format(config["params"]["genome"])
    log:
        "logs/macs2/bw/{sample}.log"
    shell:
        """
        LC_COLLATE=C sort -k 1,1 -k 2,2n {input} > {input}.sorted && bedGraphToBigWig {input}.sorted {params.genome} {output} && rm {input}.sorted &> {log}
        """

rule macs_bigwig_unit:
    input: "macs/units/{sample}-{unit}_treat_pileup.bdg"
    output: "bw/macs2/units/{sample}-{unit}_treat_pileup.bw"
    log: "logs/macs2/bw/{sample}-{unit}.log"
    conda: "../envs/macs2.yaml"
    params:
        genome = "{}".format(config["params"]["genome"])
    log:
        "logs/macs2/bw/{sample}-{unit}.log"
    shell:
        """
        LC_COLLATE=C sort -k1,1 -k2,2n {input} > {input}.sorted && bedGraphToBigWig {input}.sorted {params.genome} {output} && rm {input}.sorted &> {log}
        """

