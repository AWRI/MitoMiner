# preflight
include: "rules/00-common.smk"


# rules
include: "rules/00-intermediateFiles.smk"
include: "rules/01-map.smk"
include: "rules/02-assemble.smk"
include: "rules/03-circularise.smk"

rule all:
    input:
        finalOutput


rule clean:
    shell:
        'rm -rf data/processing/errorLogs/out/summaries/data/ref/S288C.NC_001133.9.fna.*'
