#!/usr/bin/env python3
"""
Simple Data Download Script for FerroScore-Radio
Manual download instructions included
"""

import os

print("=" * 70)
print("FerroScore-Radio: Data Download Instructions")
print("=" * 70)

print("""
由于网络限制，自动下载经常失败。请使用以下方法之一：

═══════════════════════════════════════════════════════════════════════
方法 1: 浏览器直接下载（推荐）
═══════════════════════════════════════════════════════════════════════

请打开浏览器，访问以下链接直接下载：

1. TCGA 基因表达数据（~2 GB）
   https://toil-xenagoatskins.s3.amazonaws.com/tcga_RSEM_gene_tpm.gz

2. TCGA 临床信息（~5 MB）
   https://toil-xenagoatskins.s3.amazonaws.com/TCGA_phenotype_denseDataOnlyDownload.tsv.gz

3. TCGA 生存数据（~2 MB）
   https://xenagoatskins.s3.amazonaws.com/survivalTCGA_pancan.txt

下载后放到这个文件夹：
   D:\Research\FerroScore-Radio\data\raw\

═══════════════════════════════════════════════════════════════════════
方法 2: 使用 IDM/迅雷等下载工具
═══════════════════════════════════════════════════════════════════════

把上面的链接复制到下载工具中，可以断点续传。

═══════════════════════════════════════════════════════════════════════
方法 3: 从 GDC Data Portal 下载
═══════════════════════════════════════════════════════════════════════

1. 访问 https://portal.gdc.cancer.gov/
2. 点击 "Repository"
3. 选择：
   - Cases: TCGA
   - Files: RNA-seq, Gene Expression Quantification
4. 下载数据

═══════════════════════════════════════════════════════════════════════
方法 4: 使用 R 代码（需要翻墙）
═══════════════════════════════════════════════════════════════════════

在 R 中运行：

library(TCGAbiolinks)

# 下载表达数据
query <- GDCquery(
  project = "TCGA-PAAD",  # 或 TCGA-LUAD, TCGA-BRCA 等
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

GDCdownload(query)
data <- GDCprepare(query)

═══════════════════════════════════════════════════════════════════════
下载完成后
═══════════════════════════════════════════════════════════════════════

文件应该放在：
   D:\Research\FerroScore-Radio\data\raw\

需要的文件：
   - tcga_RSEM_gene_tpm.gz (或解压后的文件)
   - TCGA_phenotype_denseDataOnlyDownload.tsv.gz
   - survivalTCGA_pancan.txt

然后运行：
   python 02_data_preprocessing.py

═══════════════════════════════════════════════════════════════════════
""")

# 创建目录
os.makedirs('../data/raw', exist_ok=True)
print(f"\n✓ 已创建目录: {os.path.abspath('../data/raw')}")
print("\n请把下载的文件放到上面的文件夹中。")
