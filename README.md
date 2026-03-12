# FerroScore-Immuno

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

> **首个专门针对免疫治疗反应的铁死亡评分系统**
> 
> FerroScore-Immuno: A machine learning-derived ferroptosis signature for predicting immunotherapy response across cancers

## 🎯 核心创新

1. **Ferro-Immuno联合基因集**：铁死亡基因 + 免疫相关基因 + TME基因
2. **免疫治疗特异性预测**：专门预测免疫治疗反应，而非一般预后
3. **泛癌种统一模型**：一个模型适用于多种癌症
4. **联合用药指导**：预测免疫治疗+铁死亡诱导剂的协同效应

## 📊 工作流程

```
TCGA/GTEx 数据 → 基因表达矩阵 → FerroScore算法 → 机器学习模型 → 免疫治疗反应预测
     ↓                ↓               ↓              ↓              ↓
  数据下载        预处理          评分计算        模型训练        验证分析
```

## 🚀 快速开始

### 安装依赖

```bash
pip install -r requirements.txt
```

### 运行分析

```bash
# Step 1: 下载数据
python code/01_data_download.py

# Step 2: 数据预处理
python code/02_data_preprocessing.py

# Step 3: 计算FerroScore
python code/03_ferroscore_algorithm.py

# Step 4: 机器学习模型
python code/04_model_training.py

# Step 5: 验证分析
python code/05_validation.py

# Step 6: 可视化
python code/06_visualization.py
```

## 📁 项目结构

```
FerroScore-Radio/
├── data/                   # 数据文件夹
│   ├── raw/               # 原始数据 (TCGA, GTEx)
│   ├── processed/         # 处理后数据
│   └── external/          # 外部验证队列
├── gene_sets/             # 基因集
│   └── ferro_radio_genes.txt
├── code/                  # 分析代码
├── results/               # 结果输出
│   ├── tables/           # 数据表格
│   └── figures/          # 图表
├── shiny_app/             # 在线工具
└── manuscript/            # 手稿
```

## 🧬 Ferro-Immuno 基因集

| 类别 | 基因数 | 功能 |
|------|--------|------|
| 铁死亡 Driver | 9 | 促进铁死亡 |
| 铁死亡 Suppressor | 15 | 抑制铁死亡 |
| 抗原呈递 | 12 | MHC分子、抗原加工 |
| T细胞浸润 | 14 | 趋化因子、T细胞标志物 |
| 免疫检查点 | 8 | PD-1/PD-L1、CTLA-4等 |
| 肿瘤微环境 | 10 | TME相关基因 |
| 炎症因子 | 6 | 细胞因子、炎症信号 |

**总计：约80个核心基因**

## 📈 预期结果

- **FerroScore**: 铁死亡活性评分
- **Immune Score**: 免疫浸润评分  
- **FerroImmuno Score**: 综合免疫治疗敏感性评分
- **机器学习模型**: 预测免疫治疗反应 (AUC目标 > 0.70)

## 📚 参考文献

1. FerrDb: http://www.zhounan.org/ferrdb/
2. UCSC Xena: https://xena.ucsc.edu/
3. TCGA: https://portal.gdc.cancer.gov/

## 👥 作者

- 霍悉尼 (医学生)
- 对虾 (AI助手) 🦐

## 📄 许可

MIT License

## 🙏 致谢

感谢 SenScoreR 文章的方法学参考
