# WGS-pipeline

[![Gitpod ready-to-code](https://img.shields.io/badge/Gitpod-ready--to--code-908a85?logo=gitpod)](https://gitpod.io/#https://github.com/itaysol/Gene-pipeline)

---

## pipeline structure

from bohra (https://github.com/MDU-PHL/bohra): ![pipeline](Bohra.workflow.png)

| step | name | nodes on figure | tools involved |
| --- | --- | --- | --- |
| 1 | QC | [Paired end reads] -> [Quality check reads] | fastp [https://github.com/opengene/fastp], mutliQC [https://multiqc.info/] |
| 2 | reads species ID | [Quality check reads] -> [Species identification (kraken2)] | kraken2 [https://github.com/DerrickWood/kraken2], bracken [https://github.com/jenniferlu717/Bracken], KrakenTools [https://github.com/jenniferlu717/KrakenTools], multiQC [https://multiqc.info/] |
| 3 | assembly | [Quality check reads] -> [assembly with Shovill (skesa, spades)] | shovill [https://github.com/tseemann/shovill] |
| 4 | assembly species ID | [assembly with Shovill (skesa, spades)] -> [Quality check assemblies] -> [mlst] | kraken2, bracken, KrakenTools, seqkit [https://github.com/shenwei356/seqkit], mlst [https://github.com/tseemann/mlst], rMLST [https://github.com/Kincekara/rMLST] |
| 5 | resistome | [Quality check assemblies] -> [AMR genes] | abricate [https://github.com/tseemann/abricate], amrfinderplus [https://github.com/ncbi/amr], hamronization [https://github.com/pha4ge/hAMRonization] |
| 6 | virulome | [Quality check assemblies] -> [Virulence genes] | abricate, hamronization |
