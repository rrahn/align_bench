cmake_minimum_required (VERSION 3.14)

include (pairalign_datasources)

declare_datasource (FILE sim_reads_n1K_rl150.fq
                    URL ${CMAKE_SOURCE_DIR}/data/sim_reads_n1K_rl150.fq
                    URL_HASH SHA256=26a9b8eb173b8d2c5cf2810dd2612a91e193924b17a733ca66a5e5309e2334b4)

# Data sources for benchmarking

# declare_datasource (FILE Ash1_v2.2.fa.gz
#                     URL ftp://ftp.ccb.jhu.edu/pub/data/Homo_sapiens/Ash1/v2.2/Assembly/Ash1_v2.2.fa.gz)

# declare_datasource (FILE ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.jst
#                     URL ftp://ftp.imp.fu-berlin.de/pub/rmaerker/jstmap/current/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.jst
#                     URL_HASH SHA256=9f787f1f0b14375a9da695d42457da13a7cf8a4ee977ff5b0d1867717b2acd17)

# declare_datasource (FILE ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.jst.ibf
#                     URL ftp://ftp.imp.fu-berlin.de//pub/rmaerker/jstmap/current/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.jst.ibf
#                     URL_HASH SHA256=c1cc172b447bb0b4d95b62a1b948e5c23e8bb82b743c38aba0fac071b6399a33)
