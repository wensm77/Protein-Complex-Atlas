# Multimer-Surface
该项目的脚本用于预测蛋白质的表面残基，使用 **SaProt** 进行特征提取，并通过接入 **Transformer** 模型和预测头实现残基的预测。

## 安装

在开始使用前，请确保运行环境与  [SaProt](https://github.com/westlake-repl/SaProt)  项目保持一致。我们将提供一个 `environment.sh` 脚本来帮助配置环境。

### 环境配置

请下载并运行我们提供的 `environment.sh` 脚本来配置所需的运行环境。

```bash
conda create -n SaProt python=3.10
conda activate SaProt

bash environment.sh
```

## 使用

我们选择使用[SaProt_650M_AF2](https://huggingface.co/westlake-repl/SaProt_650M_AF2)，SaProt的GitHub首页也提供了其他参数供下载


### 特征提取
我们选择使用SaProt将特征保存到本地供后续使用。
使用脚本`python prep_data.py`并修改其中的`data_csv`值。

### 推理
`python inf.py` 
