process GET_ACCESSIONS {
    input:
    file sra_data
    
    output:
    path "acc_list.csv"
    
    script:
    """
    #!/usr/bin/env python3

    import pandas as pd
    
    df = pd.read_csv('${sra_data}',index_col=0)
    pd.Series(df.index).to_csv('acc_list.csv', index=False, header=False)
    """ 
}

process GET_AMPLICON_SCHEME {

    input:
    val sra_accession
    file sra_data    

    output:
    val sra_accession
    path('primer_scheme.txt')


    script:
    """
    #!/usr/bin/env python3

    import pandas as pd
    
    df = pd.read_csv('${sra_data}',index_col=0)
    scheme = df.loc['${sra_accession}','amplicon_PCR_primer_scheme']
    
    scheme = str(scheme)

    if scheme == 'nan':
        primer_scheme = 'unknown'
    elif 'QIAseq' in scheme or 'v3' in scheme:
        primer_scheme = 'ARTICv3'
    elif 'V5.3' in scheme:
        primer_scheme = 'ARTICv5.3.2'
    elif 'V4.1' or 'v4.1' in scheme:
        primer_scheme = 'ARTICv4.1'
    elif 'SNAP' in scheme:
        primer_scheme = 'snap_primers'
    else:
        primer_scheme = 'unknown'

    with open('primer_scheme.txt', 'w') as f:
            f.write(primer_scheme)
    """
}


process FASTERQ_DUMP {
    container "ncbi/sra-tools"
    //errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }
    //maxRetries 3
    
    disk '1TB'
    input:
    val accession


    output:
    tuple val(accession), path("*.fastq.gz")

    script:
    """
    #!/bin/sh
    prefetch ${accession}
    fasterq-dump --split-files ${accession} --disk-limit-tmp '64G'
    gzip *.fastq
    """
}