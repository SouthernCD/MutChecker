def get_versions():
    return versions[0]["number"]


versions = [
    {
        "number": "0.0.5",
        "features": [
            "1. zip the output files if use clean flag",
        ],
    },    
    {
        "number": "0.0.4",
        "features": [
            "1. fix bug in miniprot frame shift",
        ],
    },    
    {
        "number": "0.0.3",
        "features": [
            "1. add bamstat subcommand",
        ],
    },    
    {
        "number": "0.0.2",
        "features": [
            "1. let reseq can parse whole genome",
        ],
    },
    {
        "number": "0.0.1",
        "features": [
            "1. init",
            "2. build mutchecker for resequencing data",
        ],
    },
]
