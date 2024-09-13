## Clinseq barcodes

Each sample+preparation+capture item should have a corresponding barcode with the format `PROJECT-SDID-TYPE-SAMPLEID-PREPID-CAPTUREID` where:

* `PROJECT` is a two-letter short project designator. One of `AL` (alascca), `LB` (liquid biopspy) and `OT` (other)
* `SDID` is an identifier for a single individual. It must match the pattern `P-[a-zA-Z0-9]+` (*NOTE:* This necessitates an additional "-" within this field).
* `TYPE` is the sample type, one of `T` (tumor), `N` (normal) and `CFDNA` (ctDNA)
* `SAMPLEID` identifies a single biological sample, for example piece of a tumor or a single tube of plasma. It must match the pattern `[a-zA-Z0-9]+`.
* `PREPID` specifies the library preparation kit used. It must be a two-letter shortname followed by a string matching `[0-9]+`, which can be used to indicate the date on which the prep was performed. The date string should *preferably* be in the format `YYYYMMDDHHMM`. For example, `201701241540` would indicate year 2017, January 24th, at 15:40.
* `CAPTUREID` specifies the capture that was performed on the library (if any). It must match either `WGS` (indicating that no capture was performed), or else a two-letter shortname indicating the capture kit used, followed by a string matching `[0-9]+`, which can be used to indicate the date on which the capture was performed. The date should *preferably* be in the format `YYYYMMDDHHMM`.

**NOTE:** The combination `SDID-TYPE-SAMPLEID` must uniquely identify a single sample.

**NOTE:** A clinseq barcode is not garuanteed to uniquely specify a single sample+library+capture item, but in practice it should be unique if precise preparation and capture times are included within the `PREPID` and `CAPTUREID` fields.

It is always ideal to have the full sample string with prep and capture dates. example - PB-P-00462065-CFDNA-04055058-KH20221214-C420221214

### Allowed Project IDs

* `AL` = `ALASCCA` 
* `LB` = `LIQBIO` 
* `OT` = `OTHERS` #use this for all extra projects 
* `PB` = `PROBIO` 
* `PSFF` = `PSFF` 
* `UL` = `ULLEN` 
* `iPCM` = `IPCM` 
* `CRCR` = `CRC Reflex`
* `SARC` = `Sarcoma WGS`
* `CPC` = `Clinical prostate cancer`
* `BM` = `Benchmarking`

### Allowed Prep IDs

Autoseq knows about the following preparation methods: 

* `BN` = `BIOO_NEXTFLEX`
* `KH` = `KAPA_HYPERPREP`
* `KP` = `KAPA_HYPERPLUS`
* `TD` = `THRUPLEX_DNASEQ`
* `TP` = `THRUPLEX_PLASMASEQ`
* `TF` = `THRUPLEX_FD`
* `TS` = `TRUSEQ_RNA`
* `NN` = `NEBNEXT_RNA`
* `VI` = `VILO_RNA`            

### Allowed Capture IDs

Autoseq knows about the following capture kits:

* `CS` = `clinseq_v3_targets`
* `CZ` = `clinseq_v4`
* `EX` = `EXOMEV3`
* `EO` = `EXOMEV1`
* `RF` = `fusion_v1`
* `CC` = `core_design`
* `CD` = `discovery_coho`
* `CB` = `big_design`
* `AL` = `alascca_targets`
* `TT` = `test-regions`
* `CP` = `progression`
* `CM` = `monitor`
* `PC` = `probio_comprehensive`
* `PB` = `probio_biomarker_signature`
* `PA` = `pancancer`
* `C2` = `probio_comprehensive2`
* `C3` = `probio_comprehensive3`
* `C4` = `probio_comprehensive4`
* `PN` = `pancancer2`
* `PE` = `pancancer2_enzymatic`
* `P2` = `probio_biomarkersignature2`
* `S2` = `probio_snvindel2`
* `B2` = `probio_baseline2`
* `PS` = `probio_snvindel`
* `P3` = `probio_biomarkersignature2`
* `N3` = `pancancer3`
