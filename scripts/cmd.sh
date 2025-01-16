1542  ls
 1543  tar -xvzf RNAgunjan.tar.gz 
 1544  cd /home/rgitbt/RNAgunjan/try/pipeline
 1545  code .
 1546  wget -P hg38/ https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
 1547  gunzip hg38/hg38.fa.gz 
 1548  samtools faidx hg38/hg38.fa 
 1549  conda env list
 1550  conda activate samtools
 1551  fastp
 1552  samtools faidx hg38/hg38.fa 
 1553  wget -P hg38/ https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
 1554  wget -P hg38/ https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx
 1555  ls -lh hg38/
 1556  conda install bioconda::gatk4
 1557  conda install bioconda::bwa-mem2
 1558  ls
 1559  bwa mem -t 4 -R "@RG\tID:HG008-N-D\tPL:ILLUMINA\tSM:HG008-N-D" hg38/hg38.fa trim/PA220KH-lib09-P19-Tumor_S2_L001_R1.fastq.gz raw/PA220KH-lib09-P19-Tumor_S2_L001_R2_001.fastq.gz > Tumor_S2_L001.sam
 1560  bwa-mem2 -t 4 -R "@RG\tID:HG008-N-D\tPL:ILLUMINA\tSM:HG008-N-D" hg38/hg38.fa trim/PA220KH-lib09-P19-Tumor_S2_L001_R1.fastq.gz raw/PA220KH-lib09-P19-Tumor_S2_L001_R2_001.fastq.gz > Tumor_S2_L001.sam
 1561  bwa-mem2
 1562  bwa-mem2 -R "@RG\tID:HG008-N-D\tPL:ILLUMINA\tSM:HG008-N-D" hg38/hg38.fa trim/PA220KH-lib09-P19-Tumor_S2_L001_R1.fastq.gz raw/PA220KH-lib09-P19-Tumor_S2_L001_R2_001.fastq.gz > Tumor_S2_L001.sam
 1563  bwa-mem2 hg38/hg38.fa trim/PA220KH-lib09-P19-Tumor_S2_L001_R1.fastq.gz raw/PA220KH-lib09-P19-Tumor_S2_L001_R2_001.fastq.gz > Tumor_S2_L001.sam
 1564  bwa-mem2 index hg38/hg38.fa
 1565  df -h
 1566  conda install bioconda::bwa
 1567  bwa index hg38/hg38.fa
 1568  bwa mem -t 4 -R '@RG\tID:sample\tSM:sample\tPL:Illumina'  trim/PA220KH-lib09-P19-Tumor_S2_L001_R1.fastq.gz trim/PA220KH-lib09-P19-Tumor_S2_L001_R2.fastq.gz > PA220KH-lib09-P19-Tumor_S2_L001.sam
 1569  bwa mem -t 4 -R '@RG\tID:sample\tSM:sample\tPL:Illumina' hg38/hg38.fa trim/PA220KH-lib09-P19-Tumor_S2_L001_R1.fastq.gz trim/PA220KH-lib09-P19-Tumor_S2_L001_R2.fastq.gz > PA220KH-lib09-P19-Tumor_S2_L001.sam
 1570  ls
 1571  ls -lh
 1572  bwa mem -t 4 -R '@RG\tID:sample\tSM:sample\tPL:Illumina' hg38/hg38.fa trim/PA221MH-lib09-P19-Norm_S1_L001_R1.fastq.gz trim/PA221MH-lib09-P19-Norm_S1_L001_R2.fastq.gz > PA221MH-lib09-P19-Norm_S1_L001.sam
 1573  gatk MarkDuplicatesSpark -I PA220KH-lib09-P19-Tumor_S2_L001.sam  -O PA220KH_sorted_dedup_reads.bam
 1574  ls
 1575  gatk MarkDuplicatesSpark -I PA220KH-lib09-P19-Tumor_S2_L001.sam -O PA220KH_sorted_dedup_reads.bam
 1576  ls
 1577  ls hg38/
 1578  ls
 1579  mkdir sambam
 1580  rm fastp.*
 1581  rm PA22*
 1582  ls
 1583  rm Tumor_S2_L001.sam 
 1584  rm PA22*
 1585  ls
 1586  ls trim/
 1587  mkdir trim/subset
 1588  conda install bioconda::seqtk
 1589  gatk
 1590  gatk --list
 1591  seqtk 
 1592  seqtk sample -s100 trim/PA220KH-lib09-P19-Tumor_S2_L001_R1.fastq.gz 1000000 > PA220KH-Tumor_S2_L001_R1.fastq.gz
 1593  ls
 1594  ls -lh raw/
 1595  ls -lh trim/
 1596  seqtk sample -s100 trim/PA220KH-lib09-P19-Tumor_S2_L001_R2.fastq.gz 1000000 > PA220KH-Tumor_S2_L001_R2.fastq.gz
 1597  seqtk sample -s100 trim/PA221MH-lib09-P19-Norm_S1_L001_R1.fastq.gz 1000000 > PA221MH-Norm_S1_L001_R1.fastq.gz
 1598  seqtk sample -s100 trim/PA221MH-lib09-P19-Norm_S1_L001_R2.fastq.gz 1000000 > PA221MH-Norm_S1_L001_R2.fastq.gz
 1599  mv PA22* trim/subset/
 1600  ls -lh raw/
 1601  ls -lh trim/
 1602  ls -lh trim/subset/
 1603  bwa mem -t 4 -R '@RG\tID:sample\tSM:sample\tPL:Illumina' hg38/hg38.fa trim/subset/PA220KH-Tumor_S2_L001_R1.fastq.gz trim/subset/PA220KH-Tumor_S2_L001_R2.fastq.gz > sambam/PA220KH-Tumor_S2.sam
 1604  ls -lh sambam/
 1605  head -n 2 sambam/PA220KH-Tumor_S2.sam 
 1606  head -n1 2 sambam/PA220KH-Tumor_S2.sam 
 1607  head -n12 sambam/PA220KH-Tumor_S2.sam 
 1608  samtools view -Sb sambam/PA220KH-Tumor_S2.sam | samtools sort -o tumor_sorted.bam
 1609  ls -lh
 1610  gatk MarkDuplicates -I tumor_sorted.bam -O tumor_dedup.bam -M tumor.metrics.txt
 1611  bwa mem -t 4 -R '@RG\tID:sample\tSM:sample\tPL:Illumina' hg38/hg38.fa trim/subset/PA221MH-Norm_S1_L001_R1.fastq.gz trim/subset/PA221MH-Norm_S1_L001_R2.fastq.gz > sambam/PA221MH-Norm_S1.sam
 1612  samtools view -Sb sambam/PA221MH-Norm_S1.sam | samtools sort -o normal_sorted.bam
 1613  ls -lh
 1614  gatk MarkDuplicates -I normal_sorted.bam -O normal_dedup.bam -M normal.metrics.txt
 1615  ls -lh
 1616  gatk BaseRecalibrator -I normal_dedup.bam -R hg38/hg38.fa --known-sites hg38/Homo_sapiens_assembly38.dbsnp138.vcf -O normal_data.table
 1617  ls hg38/
 1618  gatk BaseRecalibrator -I normal_dedup.bam -R hg38/hg38.fa  --known-sites hg38/Homo_sapiens_assembly38.dbsnp138.vcf -O normal_data.table
 1619  ls hg38/
 1620  ls
 1621  samtools faidx hg38/hg38.fa 
 1622  ls hg38/
 1623  java -jar picard.jar CreateSequenceDictionary R=hg38/hg38.fa O=hg38/hg38.dict
 1624  gatk picard CreateSequenceDictionary R=hg38/hg38.fa O=hg38/hg38.dict
 1625  gatk CreateSequenceDictionary R=hg38/hg38.fa O=hg38/hg38.dict
 1626  ls hg38/
 1627  gatk BaseRecalibrator -I normal_dedup.bam -R hg38/hg38.fa  --known-sites hg38/Homo_sapiens_assembly38.dbsnp138.vcf -O normal_data.table
 1628  ls
 1629  head -n 12 normal_data.table 
 1630  gatk BaseRecalibrator -I tumor_dedup.bam -R hg38/hg38.fa  --known-sites hg38/Homo_sapiens_assembly38.dbsnp138.vcf -O tumor_data.table
 1631  cat tumor_data.table 
 1632  gatk ApplyBQSR -I tumor_dedup.bam -R hg38/hg38.fa --bqsr-recal-file tumor_data.table -O tumor_sorted_dedup_bqsr_reads.bam 
 1633  gatk ApplyBQSR -I normal_dedup.bam -R hg38/hg38.fa --bqsr-recal-file normal_data.table -O normal_sorted_dedup_bqsr_reads.bam 
 1634  gatk CollectAlignmentSummaryMetrics R= hg38/hg38.fa I=normal_sorted_dedup_bqsr_reads.bam  O=normal_alignment_metrics.txt
 1635  head -n 12 normal_alignment_metrics.txt 
 1636  gatk CollectInsertSizeMetrics INPUT=normal_sorted_dedup_bqsr_reads.bam OUTPUT=normal_insert_size_metrics.txt HISTOGRAM_FILE=normal_insert_size_histogram.pdf
 1637  gatk CollectAlignmentSummaryMetrics R= hg38/hg38.fa I=tumor_sorted_dedup_bqsr_reads.bam  O=tumor_alignment_metrics.txt
 1638  gatk CollectInsertSizeMetrics INPUT=normal_sorted_dedup_bqsr_reads.bam OUTPUT=tumor_insert_size_metrics.txt HISTOGRAM_FILE=tumor_insert_size_histogram.pdf
 1639  mkdir mutect
 1640  wget https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/somatic-hg38/af-only-gnomad.hg38.vcf.gz -P mutect/
 1641  wget https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/somatic-hg38/af-only-gnomad.hg38.vcf.gz.tbi  -P mutect/
 1642  wget https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/exome_calling_regions.v1.1.interval_list -P mutect/
 1643  gatk CreateSomaticPanelOfNormals   -vcfs normal.vcf.gz   -O norm_pon.vcf.gz
 1644  gatk CreateSomaticPanelOfNormals   -vcfs normal_for_pon.vcf.gz   -O mutect/norm_pon.vcf.gz
 1645  Mutect2 -R hg38/hg38.fa      -I sambam/tumor_sorted_dedup_bqsr_reads.bam      -I sambam/normal_sorted_dedup_bqsr_reads.bam      -tumor PA220KH-lib09-P19-Tumor_S2      -normal PA221MH-lib09-P19-Norm_S1      --germline-resource mutect/af-only-gnomad.hg38.vcf.gz      --panel-of-normals mutect/normal_for_pon.vcf.gz      -O somatic_variants_mutect2.vcf.gz      --f1r2-tar-gz somatic.tar.gz 
 1646  gatk Mutect2 -R hg38/hg38.fa      -I sambam/tumor_sorted_dedup_bqsr_reads.bam      -I sambam/normal_sorted_dedup_bqsr_reads.bam      -tumor PA220KH-lib09-P19-Tumor_S2      -normal PA221MH-lib09-P19-Norm_S1      --germline-resource mutect/af-only-gnomad.hg38.vcf.gz      --panel-of-normals mutect/normal_for_pon.vcf.gz      -O somatic_variants_mutect2.vcf.gz      --f1r2-tar-gz somatic.tar.gz 
 1647  head -n 2 sambam/normal_sorted_dedup_bqsr_reads.bam 
 1648  gatk Mutect2 -R hg38/hg38.fa      -I sambam/tumor_sorted_dedup_bqsr_reads.bam      -I sambam/normal_sorted_dedup_bqsr_reads.bam      -tumor PA220KH-Tumor_S2      -normal PA221MH-Norm_S1      --germline-resource mutect/af-only-gnomad.hg38.vcf.gz      --panel-of-normals mutect/normal_for_pon.vcf.gz      -O somatic_variants_mutect2.vcf.gz      --f1r2-tar-gz somatic.tar.gz 
 1649  gatk Mutect2 -R hg38/hg38.fa      -I sambam/tumor_sorted_dedup_bqsr_reads.bam      -I sambam/normal_sorted_dedup_bqsr_reads.bam      -tumor tumor      -normal normal      --germline-resource mutect/af-only-gnomad.hg38.vcf.gz      --panel-of-normals mutect/normal_for_pon.vcf.gz      -O somatic_variants_mutect2.vcf.gz      --f1r2-tar-gz somatic.tar.gz 
 1650  gatk Mutect2 -R hg38/hg38.fa      -I sambam/tumor_sorted_dedup_bqsr_reads.bam      -I sambam/normal_sorted_dedup_bqsr_reads.bam    --germline-resource mutect/af-only-gnomad.hg38.vcf.gz      --panel-of-normals mutect/normal_for_pon.vcf.gz      -O somatic_variants_mutect2.vcf.gz      --f1r2-tar-gz somatic.tar.gz 
 1651  head somatic_variants_mutect2.vcf.gz  -n 12
 1652  head somatic_variants_mutect2.vcf.gz.stats -n 12
 1653  head somatic_variants_mutect2.vcf.gz.stats 
 1654  ls
 1655  gatk GetPileupSummaries     --java-options '-Xmx50G' --tmp-dir ~/RNAgunjan/project/tmp/     -I sambam/tumor_sorted_dedup_bqsr_reads.bam     -V mutect/af-only-gnomad.hg38.vcf.gz     -L mutect/exome_calling_regions.v1.1.interval_list     -O T_getpileupsummaries.table
 1656  mkdir tmp
 1657  gatk GetPileupSummaries     --java-options '-Xmx50G' --tmp-dir tmp/     -I sambam/tumor_sorted_dedup_bqsr_reads.bam     -V mutect/af-only-gnomad.hg38.vcf.gz     -L mutect/exome_calling_regions.v1.1.interval_list     -O T_getpileupsummaries.table
 1658  gatk GetPileupSummaries     --java-options '-Xmx50G' --tmp-dir tmp/     -I sambam/normal_sorted_dedup_bqsr_reads.bam     -V mutect/af-only-gnomad.hg38.vcf.gz     -L mutect/exome_calling_regions.v1.1.interval_list     -O N_getpileupsummaries.table
 1659  gatk CalculateContamination     -I T_getpileupsummaries.table     -matched N_getpileupsummaries.table     -O somatic_pair_calculatecontamination.table 
 1660  gatk LearnReadOrientationModel     -I somatic_f1r2.tar.gz     -O read-orientation-model.tar.gz
 1661  gatk LearnReadOrientationModel     -I somatic.tar.gz     -O read-orientation-model.tar.gz
 1662  FilterMutectCalls         -V somatic_variants_mutect2.vcf.gz         -R hg38/hg38.fa         --contamination-table pair_calculatecontamination.table         --ob-priors read-orientation-model.tar.gz         -O somatic_variants_filtered_mutect2.vcf
 1663  gatk FilterMutectCalls         -V somatic_variants_mutect2.vcf.gz         -R hg38/hg38.fa         --contamination-table pair_calculatecontamination.table         --ob-priors read-orientation-model.tar.gz         -O somatic_variants_filtered_mutect2.vcf
 1664  ls
 1665  gatk FilterMutectCalls         -V somatic_variants_mutect2.vcf.gz         -R hg38/hg38.fa         --contamination-table somatic_pair_calculatecontamination.table         --ob-priors read-orientation-model.tar.gz         -O somatic_variants_filtered_mutect2.vcf
 1666  gatk FuncotatorDataSourceDownloader     --data-source-version v1.8     --reference-version hg38     --output ~/RNAgunjan/project/
 1667  gatk FuncotatorDataSourceDownloader     --reference-version hg38     --output ~/RNAgunjan/project/
 1668  ahtop
 1669  htop
 1670  find -h
 1671  find --help
 1672  find -n functotator* 
 1673  find -name functotator* 
 1674  find -name functotator* ~
 1675  find ~ -name "functotator"
 1676  find ~ -name "functotator*"
 1677  find ~ -name "functotator_prepackaged_sources"
 1678  find ~ -name "HG38"
 1679  find ~ -name "hg38"
 1680  sudo apt-get remove anydesk
 1681  sudo apt-get purge anydesk
 1682  clear
 1683  zip -r Sourav-proj/
 1684  zip -r Sourav-proj/ Sourav-proj.zip
 1685  zip -r Sourav-proj Sourav-proj
 1686  du -h
 1687  ls
 1688  df -h
 1689  ls -lh
 1690  ls -lh Sourav-proj
 1691  free -m
 1692  history
