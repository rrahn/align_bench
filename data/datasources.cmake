cmake_minimum_required (VERSION 3.14)

include (pairalign_datasources)

declare_datasource (FILE sim_reads_n1K_rl150.fq
                    URL ${CMAKE_SOURCE_DIR}/data/sim_reads_n1K_rl150.fq
                    URL_HASH SHA256=26a9b8eb173b8d2c5cf2810dd2612a91e193924b17a733ca66a5e5309e2334b4)

declare_datasource (FILE sim_reads_sanger_n1k_min400_max800.fq
                    URL ${CMAKE_SOURCE_DIR}/data/sim_reads_sanger_n1k_min400_max800.fq
                    URL_HASH SHA256=f0d58f3d3028935219621b1ef5fa090abbcf5dfb10a688fcab0c2e5224699a5e)

# Data sources for benchmarking
declare_datasource (FILE sim_as500_n1k.fasta
                    URL ${CMAKE_SOURCE_DIR}/data/sim_as500_n1k.fasta
                    URL_HASH SHA256=0d2d3ef5cc1e675eb2159e1a92d882438e70309b0388afd495379319e7a1ae73)

declare_datasource (FILE uniprot_sprot.fasta.gz
                    URL https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
                    URL_HASH SHA256=2562cf025013d6e70ebbe612a5c41ddd0cc381dba8df5f7c6ad36f0c0592f579)
