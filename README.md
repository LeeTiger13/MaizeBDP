# MaizeBDP

MaizeBDP is a maize breeding decision platform prototype composed of two major modules:

- `MIKA/`: a hybrid maize breeding assistant that routes questions to either a RAG pipeline or a maize pedigree knowledge graph.
- `MaizeBDP_Function_Module/`: analytical scripts for genotype pretreatment, genetic background analysis, breeding guidance, and genomic selection.

This public repository shares the platform architecture, script organization, and runnable examples. It is intended primarily as a reusable framework reference. Private production resources, including internal databases, vector stores, pretrained local model assets, and API credentials, are intentionally not distributed.

## Repository Structure

```text
MaizeBDP/
|-- MIKA/
|   |-- main.py
|   |-- mika_langgraph.py
|   |-- utils.py
|   |-- language_adapter.py
|   |-- config.yaml
|   `-- Dockerfile
`-- MaizeBDP_Function_Module/
    |-- Data_Pretreatment/
    |-- Genetic_Background/
    |-- Breeding_Guidance/
    |-- Genomic_Selection/
    `-- tools/
```

## Platform Overview

`MIKA` uses a dual-path question answering architecture:

1. A cloud LLM rewrites multi-turn questions and decides whether the query should go to RAG or KG.
2. KG questions are answered from a Neo4j maize pedigree graph.
3. RAG questions are answered by combining a vector retrieval backend, reranking, and a local or remote answer model.

`MaizeBDP_Function_Module` contains task-oriented scripts that can be integrated into a web platform or run independently from the command line. These scripts mainly process VCF-based genotype inputs and generate tables, figures, or prediction results for downstream breeding analysis.

## What Is Not Included

The following assets are project-specific and are not included in the public release:

- production API keys and endpoint credentials
- internal Neo4j data and graph construction pipelines
- vector database files, indexed literature corpora, and retrieval services
- local model assets, commercial LLM API settings, and reranker checkpoints
- proprietary or unpublished breeding datasets
- platform-side reference panels required by some function modules

You can still reuse the architecture by replacing these resources with your own equivalents.

## Requirements

### Core Environment

- Python 3.11 recommended for `MIKA`
- Python 3.9+ recommended for `MaizeBDP_Function_Module`
- Git
- Docker optional, for containerized chatbot deployment

### Python Packages

For `MIKA`, install packages equivalent to:

- `fastapi`
- `uvicorn[standard]`
- `langgraph`
- `langchain-core`
- `pydantic`
- `pyyaml`
- `requests`
- `httpx`
- `openai`
- `pymilvus`
- `neo4j`
- `sentence-transformers`
- `torch`
- `typing_extensions`

Install `ollama` only if you choose `ollama` as the provider for your answer model or embedding model. If you use only OpenAI-compatible APIs, the `openai` package is sufficient for those parts.

For `MaizeBDP_Function_Module`, common packages include:

- `pandas`
- `numpy`
- `matplotlib`
- `scikit-learn`
- `lightgbm`

### External Bioinformatics Tools

Several scripts depend on command-line tools available from `PATH`:

- `plink`
- `bcftools`
- `bgzip`
- `tabix`
- `Rscript`
- `admixture`

Some submodules also require matching reference panel files that are not included in this repository.

## Quick Start

### 1. Clone the repository

```cmd
git clone <your-repo-url>
cd MaizeBDP_git
```

### 2. Prepare a Python environment

Example with `conda`:

```cmd
conda create -n maizebdp python=3.11 -y
conda activate maizebdp
```

### 3. Install Python dependencies

For the chatbot:

```cmd
pip install langgraph langchain-core fastapi "uvicorn[standard]" pydantic pyyaml httpx requests pymilvus neo4j ollama openai sentence-transformers typing_extensions torch
```

For the function module:

```cmd
pip install pandas numpy matplotlib scikit-learn lightgbm
```

### 4. Configure the chatbot

The committed file `MIKA/config.yaml` is a public template only.

Create a private local config:

```cmd
copy MIKA\config.yaml MIKA\config.local.yaml
```

Then edit `MIKA/config.local.yaml` with your own:

- `routing_model`: routing LLM API key, endpoint, and model
- `answer_model`: provider, endpoint, credentials, and model for final answer generation
- `embedding_model`: provider, endpoint, credentials, and model for embedding generation
- `retrieval_backend`: retrieval API endpoint, timeout, and collection name
- `knowledge_graph`: Neo4j connection information
- `assets`: reranker path
- `search`: retrieval and generation parameters

These settings should be adjusted according to the user's own server environment. In practice, service addresses, model names, authentication information, storage paths, and retrieval backend settings will usually need to be modified before deployment.

The application now loads configuration in this order:

1. `MIKA_CONFIG` environment variable, if provided
2. `MIKA/config.local.yaml`
3. `MIKA/config.yaml`

You can also override selected values through environment variables:

- `MAIZEBDP_CLOUD_API_KEY`
- `MAIZEBDP_CLOUD_BASE_URL`
- `MAIZEBDP_CLOUD_MODEL`
- `MAIZEBDP_ANSWER_API_BASE_URL`
- `MAIZEBDP_ANSWER_API_KEY`
- `MAIZEBDP_ANSWER_MODEL`
- `MAIZEBDP_EMBEDDING_API_BASE_URL`
- `MAIZEBDP_EMBEDDING_API_KEY`
- `MAIZEBDP_EMBEDDING_MODEL`
- `MAIZEBDP_RETRIEVAL_API_BASE_URL`
- `MAIZEBDP_NEO4J_URI`
- `MAIZEBDP_NEO4J_USER`
- `MAIZEBDP_NEO4J_PASSWORD`

### 5. Start the chatbot API

```cmd
cd MIKA
python main.py
```

Default service endpoint:

- `POST http://127.0.0.1:8082/chat`

Example request body:

```json
{
  "question": "What is the information of variety Zheng58?",
  "session_id": "demo-session-001"
}
```

## Reproducing the MIKA Architecture with Your Own Resources

If you want to implement a chatbot similar to `MIKA` without using our private backend resources, the recommended replacement strategy is:

1. Replace the cloud routing model with any OpenAI-compatible chat model that can handle query rewriting and JSON routing.
2. Build your own Neo4j maize pedigree graph with nodes such as `Variety` and relations such as `is_maternal` and `is_paternal`.
3. Build your own vector corpus and expose it through a retrieval backend that can return relevant text chunks for RAG.
4. Replace the answer model with any local or remote model that supports the chat completion and embedding workflow you want to use.
5. Replace the reranker path with your own cross-encoder checkpoint.
6. Keep the overall orchestration logic in `mika_langgraph.py` if you want the same dual-mode design.

In other words, this repository is best understood as an architecture and workflow reference rather than a one-click fully reproducible service release. You can swap Milvus for another vector store, or provide your own compatible retrieval service, as long as it supports the same role in the pipeline. Likewise, you do not have to use Ollama for local deployment: you can also connect any suitable commercial or self-hosted LLM API, as long as the surrounding code and configuration are adjusted to match your own environment.

## Function Module Usage

The scripts under `MaizeBDP_Function_Module` can be run independently. A few representative examples are shown below.

### Genotype Pretreatment

```cmd
cd MaizeBDP_Function_Module\Data_Pretreatment\code
python data_pretreat_report_multiple.py -i ..\test_data\input.maize_20k_v5.30sample.vcf -s 0.05 -m 0.05 -n 0.2 -a 0.05
```

### Genomic Selection

```cmd
cd MaizeBDP_Function_Module\Genomic_Selection
python DH_gs_main.py --genotype_train test_data\test_training_geno.vcf --phenotype_train test_data\test_training_pheno.txt --genotype_test test_data\test_testing_geno.vcf --models rrBLUP
```

### Marker Extraction

```cmd
cd MaizeBDP_Function_Module\tools
python extract_marker_multiple.py -i ..\Breeding_Guidance\variety_pattern_recommend\test_data\input.maize_20k_v5.query.vcf -r chr1:1-1000000
```

## Notes on Portability

- Public configuration files have been sanitized. Do not commit `config.local.yaml`.
- `Rscript` and `admixture` are now resolved from system `PATH` by default.
- If needed, you can explicitly set executable paths through environment variables:
  - `set RSCRIPT_BIN=C:\path\to\Rscript.exe`
  - `set ADMIXTURE_BIN=C:\path\to\admixture.exe`
- Several scripts still assume the presence of project-specific reference panel files. Those inputs must be prepared separately for full reproduction.
- Some bundled `test_data` VCF files may retain historical source metadata in VCF headers; these headers are not used by the scripts themselves.

## Docker

`MIKA/Dockerfile` is retained as an optional deployment method for packaging the chatbot API service. It is not required to use this repository.

If you do not use Docker, you can run the framework directly in a normal Python environment.

If you do use Docker, note that the image only packages the API service itself. You still need to adapt the project to your own server environment and provide your own:

- model endpoints or external LLM APIs
- vector retrieval backend
- Neo4j service
- reranker files
- private datasets and credentials

So in this project, Docker should be understood as a convenience option for deployment, not as a complete reproduction of the full private backend.

## Citation and Reuse

If you reuse this codebase or adapt the architecture in academic work, please cite the corresponding MaizeBDP project or paper when available.

## Authors
- Hu Li lihu0628@qq.com
- Chunchao Wang wangchunchao@caas.cn
- Junjie Fu fujunjie@caas.cn

Copyright (c) Institute of Crop Sciences, Chinese Academy of Agricultural Sciences. 2024-2026. All rights reserved.
