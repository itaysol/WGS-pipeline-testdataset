# Gene-pipeline

[![Gitpod ready-to-code](https://img.shields.io/badge/Gitpod-ready--to--code-908a85?logo=gitpod)](https://gitpod.io/#https://github.com/motroy/Gene-pipeline)

---

## pipeline structure

from bohra (https://github.com/MDU-PHL/bohra): ![pipeline](https://github.com/MDU-PHL/bohra/blob/master/workflow.png?raw=true)

| step | name | nodes on figure | tools involved |
| --- | --- | --- | --- |
| 1 | QC | [Paired end reads] -> [Quality check reads] | fastp [https://github.com/opengene/fastp], mutliQC [https://multiqc.info/], seqkit [https://github.com/shenwei356/seqkit] |
| 2 | species ID | [Quality check reads] -> [Species identification (kraken2)] | kraken2 [https://github.com/DerrickWood/kraken2], bracken [https://github.com/jenniferlu717/Bracken], KrakenTools [https://github.com/jenniferlu717/KrakenTools], multiQC [https://multiqc.info/] |
| 3 | assembly | [Quality check reads] -> [assembly with Shovill (skesa, spades)] | shovill [] |
