class config:

    ####################################
    # DEVELOPER CONFIGURABLE VARIABLES #
    ####################################

    # TESTING MODE
    TESTING = False

    # DB file name
    DB = 'Indexing.db'

    # Acceptable file formats. Files not following this format will be ignored
    UCSC_ACCEPTABLE_FILE_FORMATS = ["bw", "bb", "bed", "bigPsl"]
    LOCAL_ACCEPTABLE_FILE_FORMATS = ["bed", "bed.gz"]

    # UCSC Genomes to download
    UCSC_BIG_DATA_LINK = "http://hgdownload.soe.ucsc.edu"
    UCSC_SQL_DB_HOST = "genome-mysql.soe.ucsc.edu"
    UCSC_SQL_DB_USER = "genome"
    UCSC_API = "https://api.genome.ucsc.edu"

    ###############################
    # USER CONFIGURABLE VARIABLES #
    ###############################

    # Max Number of Intervals per Index
    MAX_INTERVALS_PER_INDEX = 1000000

    # Number of proccesses available to server
    AVAILABLE_PROCCESSES = 12

    # Max file size for setup (intervals)
    max_setup_file_size = 100000*155 # conversion: ~155 intervals per KB, 1000000 KB per GB
    # Timeout on a file download in seconds
    timeout_file_download = 60*5
    timeout_file_processing = 60*10 # 10 minutes

    # List genomes from UCSC to download and index
    '''
    example:
        UCSC_GENOMES = ["hg19", "rn6", "petMar2"]
    '''
    UCSC_GENOMES = ["petMar2"]

    # List local genomes to download and index 
    '''
    example:
        LOCAL_GENOMES = [
                         {
                         "project_name": "test",
                         "reference_genome": "hg19",
                         "metadata_path": "",
                         "data_path": "roadmap_sort/"
                         },
                        ...]
    '''
    LOCAL_GENOMES = [
                      {
                        "project_name": "fantom5",
                        "reference_genome": "hg19",
                        "metadata_path": "",
                        "data_path": "fantom_sort/"
                        }
                    ]


    # Metadata file for local directories
    # Format: .csv format with columns = ["file_name", "short_name", "long_name", "short_info", "long_info"]
    # Example  "metadata.csv"

    # Uncomment UCSC public hubs to index, see https://api.genome.ucsc.edu/list/publicHubs for further info on each hub 
    '''
    example:
        UCSC_HUBS = [
            # "ALFA Hub",
            "Bird assemblies",
            # "Blueprint Hub",
            ...
            ""
            ]
    '''

    UCSC_HUBS = [
                # "ALFA Hub",
                # "Bird assemblies",
                # "Blueprint Hub",
                # "BrainEpigenomeHub",
                # "Breast Cancer lncRNA",
                # "Broad Improved Canine Annotation v1",
                # "C_elegans_isolates",
                # "CADD",
                # "Cancer Genomics Tracks",
                # "CEMT (CEEHRC)",
                # "CESAR Gene Mappings",
                # "ChIP-seq data track HUBs from MSC cells from GSE79815",
                "Coloc segments",
                # "Cotney Lab Human Craniofacial Epigenomics",
                # "Cotney Lab Human Embryonic Heart Hub",
                # "COVID-19 Gene Annotation",
                # "CPTAC Hub v1",
                # "Croc and Bird Hub",
                # "DANIO-CODE Track Hub",
                # "DASHR small ncRNA Hub",
                # "DASHR v2.0 Hub",
                # "dbRIP_Hub",
                # "dbVar Hub",
                # "Digital genomic footprinting (Vierstra et al., 2020)",
                # "DNA Methylation",
                # "ENCODE Analysis Hub",
                # "ENCODE DNA Trackhub",
                # "ENCODE integrative Trackhub",
                # "ENCODE RNA Trackhub",
                # "Ensembl Regulatory Build",
                # "EPD Viewer Hub",
                # "Exp/Meth VNTR hub",
                # "FaceBase Hub",
                # "FANTOM5",
                # "FANTOM5 CAGE RECLU DATA",
                # "Fish assemblies",
                # "GENCODE Annotation Updates",
                # "GeneHancer",
                # "GRC Genome Issues under Review",
                # "GTEx Analysis Hub",
                # "GTEx RNA-seq Signal Hub",
                # "Human cellular microRNAome",
                # "Human cellular microRNAome barCharts",
                # "Human Global Reference Genomes",
                # "Human Islet lncRNAs",
                # "human p53 Binding And Expression Resource (BAER) hub",
                # "IDEAS roadmap 20states",
                # "JASPAR TFBS",
                # "LIBD Human DLPFC Development",
                # "LNCipedia 5.2",
                # "Mammal assemblies",
                # "McGill EMC (CEEHRC)",
                # "MGI Alleles and Phenotypes",
                # "miRcode microRNA sites",
                # "mm9.SMC1.ChIAPET",
                # "Mouse strain assemblies",
                # "ORF Predictions",
                # "Peptide evidences CNIO",
                # "Peterhof_yeasts",
                # "PhyloCSF",
                # "Plant assemblies",
                # "Polytract Repeats",
                # "Primate assemblies",
                # "Primate x4 NeuroDiff and Human CRISPRa",
                # "Principal Splice Isoforms APPRIS",
                # "Promoterome CAGE and nucleosome positioning",
                # "Purple sea star Assembly Hub",
                # "RefSeqFE Hub",
                # "refTSS",
                # "ReMap 2020 Regulatory Atlas",
                # "ReMap2018 Regulatory Atlas",
                # "rfam12_ncRNA",
                # "Roadmap Epigenomics Data Complete Collection at Wash U VizHub",
                # "Roadmap Epigenomics Integrative Analysis Hub",
                # "Seq-data on mm9 NS5 cells",
                # "Splice Site Usage Hub",
                # "Synonymous Constraint",
                # "T_cell_ATAC_ChIP_Pipkin",
                # "TOBIAS footprint prediction",
                # "Translation Initiation Sites (TIS)",
                # "UCD Methylation",
                # "UCSC Repeat Browser 2020 Update",
                # "Ultraconserved Elements",
                # "Umap and Bismap mappability",
                # "UMassMed ZHub",
                # "UniBind 2.0 permissive Hub",
                # "UniBind 2018 Hub",
                # "UniBind 2021 Robust Hub",
                # "UniProt Features",
                # "Vertebrate assemblies",
                # "VGP primary",
                # "Vista Enhancers",
                # "Wasp spider hub",
                # "WormBase",
                # "Xmac",
                # "ZebrafishGenomics",
                ""   # leave this
                ]
