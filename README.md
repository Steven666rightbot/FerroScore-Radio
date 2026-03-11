# FerroScore-Radio

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

> **首个专门针对放疗反应的铁死亡评分系统**
> 
> FerroScore-Radio: A machine learning-derived ferroptosis signature for predicting radiotherapy response across cancers

## 🎯 核心创新

1. **Ferro-Radio联合基因集**：铁死亡基因 + DNA修复基因 + ROS基因
2. **放疗特异性预测**：专门预测放疗反应，而非一般预后
3. **泛癌种统一模型**：一个模型适用于多种癌症
4. **联合用药指导**：预测放疗+铁死亡诱导剂的协同效应

## 📊 工作流程

```
TCGA/GTEx 数据 → 基因表达矩阵 → FerroScore算法 → 机器学习模型 → 放疗反应预测
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

## 🧬 Ferro-Radio 基因集

| 类别 | 基因数 | 功能 |
|------|--------|------|
| 铁死亡 Driver | 9 | 促进铁死亡 |
| 铁死亡 Suppressor | 15 | 抑制铁死亡 |
| DNA修复 | 17 | 放疗响应 |
| ROS相关 | 13 | 放疗诱导ROS |
| 细胞周期/凋亡 | 16 | 放疗细胞效应 |
| 缺氧相关 | 4 | 放疗抵抗 |

**总计：约80个核心基因**

## 📈 预期结果

- **FerroScore**: 铁死亡活性评分
- **DDR Score**: DNA损伤修复评分  
- **FerroRadio Score**: 综合放疗敏感性评分
- **机器学习模型**: 预测放疗反应 (AUC目标 > 0.70)

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
