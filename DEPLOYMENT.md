# MaizeBDP Deployment Guide

This document explains what is required to run the public MaizeBDP server code. It is separate from `README.md`: the README describes the project and its modules, while this file describes deployment requirements, external services, and installation steps.

## 1. What is included

This repository contains the server-side implementation of MaizeBDP:

- genomic analysis scripts for genotype pretreatment, genetic-background analysis, breeding guidance, and genomic selection;
- the MIKA API and its LangGraph-based orchestration code;
- representative test data and example commands.

MaizeBDP is a data-dependent decision-support platform. It is not a standalone language model and it is not a self-contained package that becomes fully functional after installing Python dependencies. The complete MIKA service also requires a knowledge graph, a literature retrieval backend, model services, and the corresponding domain data.

## 2. Choose a deployment mode

| Mode | What runs locally | GPU requirement | Use case |
|---|---|---|---|
| Genomic analysis only | PLINK/bcftools/ADMIXTURE, R and Python scripts | Not required | Run genotype analysis independently |
| MIKA with remote services | MIKA API; Ollama, Neo4j and Milvus are remote | Not required locally | Lightweight server deployment |
| MIKA with a model API | MIKA API plus user-provided Neo4j and retrieval services | Not required locally | Use a cloud or institutional LLM API |
| Fully local MIKA | MIKA, Ollama, Qwen3:14B, Neo4j and retrieval services | L20-class GPU recommended | Self-contained institutional deployment |

The current reference deployment runs `qwen3:14b` through Ollama on an NVIDIA L20. This is a reference configuration for stable local inference, not an absolute minimum. A smaller quantized model may run on a lower-memory GPU, but response speed, context length, and concurrency may be reduced. Users who do not want to maintain a local GPU can replace the local answer and embedding models with compatible remote APIs.

## 3. Resource requirements

The following values apply to the current 20K-marker examples and the present MIKA implementation.

| Component | Practical minimum | Recommended |
|---|---|---|
| Genomic analysis | 4 CPU cores, 8 GB RAM, 20 GB free storage | 8 CPU cores, 16 GB RAM, 50 GB free storage |
| MIKA API with remote model/KG/retrieval services | 4 CPU cores, 8-16 GB RAM, 30 GB free storage | 8 CPU cores, 16-32 GB RAM, 50-100 GB free storage |
| Fully local MIKA with Qwen3:14B | 8 CPU cores, 32 GB RAM, 100 GB free storage | 16 CPU cores, 64 GB RAM, 200 GB free storage; L20-class GPU recommended |

Genomic selection requires more memory as the number of training samples and markers increases. For large training populations, 32 GB RAM or more is recommended. The required storage also increases with the size of VCF files, temporary PLINK files, literature documents, vector indexes, model caches, and analysis outputs.

## 4. Software requirements

### Operating system

Linux x86_64 is the recommended environment. Windows users should use WSL2 or Docker for the command-line bioinformatics tools.

### MIKA

- Python 3.11;
- Docker, if the API is deployed using `MIKA/Dockerfile`;
- FastAPI and Uvicorn;
- LangGraph and LangChain Core;
- PyTorch and sentence-transformers;
- PyMilvus and the Neo4j Python driver;
- Ollama client or an OpenAI-compatible client;
- a reranker checkpoint, such as `ms-marco-MiniLM-L-6-v2`.

### Genomic analysis

- Python 3.9 or later;
- PLINK;
- bcftools, bgzip and tabix;
- ADMIXTURE for population-structure projection;
- R and `Rscript` for genomic-selection workflows;
- R packages `rrBLUP`, `BGLR`, `cvTools`, and `optparse` where the corresponding genomic-selection model is selected.

The required executables must be available through `PATH`, or their paths must be supplied through the environment variables documented in `README.md`.

## 5. External services and data

The MIKA API connects to the following services through `MIKA/config.local.yaml`:

| Service | Role |
|---|---|
| Routing model | Determines whether a question uses RAG or KG |
| Answer model | Generates the final response |
| Embedding model | Creates query/document embeddings |
| Neo4j | Stores the maize variety and pedigree knowledge graph |
| Milvus or a compatible retrieval backend | Stores and searches literature vectors |
| Reranker model | Re-ranks retrieved literature passages |

The public repository does not include the production Neo4j graph, indexed literature corpus, Milvus database, private reference panels, model weights, or service credentials. Users must provide equivalent resources or connect to institutional services. The genomic-analysis examples can be run with the included test data, but complete MIKA functionality cannot be reproduced from this repository alone.

`MIKA/variety_alias.yaml` is included because it is part of the query-normalization layer. It maps English or abbreviated variety names to the names used in the Chinese knowledge graph. The Chinese values in this file are intentional and are the only non-English data file retained in the public repository.

## 6. Installation

### 6.1 Clone the repository

```text
git clone <your-repository-url>
cd MaizeBDP
```

### 6.2 Create the MIKA environment

```text
conda create -n maizebdp python=3.11 -y
conda activate maizebdp
python -m pip install langgraph langchain-core fastapi "uvicorn[standard]" pydantic pyyaml httpx requests pymilvus neo4j ollama openai sentence-transformers typing_extensions torch
```

### 6.3 Configure services

```text
cp MIKA/config.yaml MIKA/config.local.yaml
```

Edit `MIKA/config.local.yaml` and provide the routing model, answer model, embedding model, retrieval backend, Neo4j connection, and reranker path. Do not commit this local configuration file or any credentials.

On Windows, use `copy MIKA\config.yaml MIKA\config.local.yaml` instead. The genomic command-line tools are still best run in Linux, WSL2, or Docker.

### 6.4 Start the MIKA API

```text
cd MIKA
python main.py
```

The default endpoint is:

```text
POST http://127.0.0.1:8082/chat
```

### 6.5 Run a genomic-analysis example

Each function module contains its own `test_data` directory and example command. The command-line tools listed in Section 4 must be installed before running these examples. The example workflows create temporary files such as filtered VCFs, PLINK files, distance matrices, and ADMIXTURE outputs in the working directory.

## 7. Basic environment checks

Before running a module, check that the required executables are visible:

```text
python --version
plink --version
bcftools --version
Rscript --version
admixture
```

For MIKA, first check that the configured Neo4j, retrieval backend, and model endpoints are reachable. A successful Python installation alone does not confirm that the MIKA service is ready.

## 8. Scaling and reproducibility

The resource values above are practical values for the current release and example datasets. Actual requirements depend on marker count, sample count, the selected genomic-selection model, literature-corpus size, model quantization, and the number of concurrent MIKA requests. For a production deployment, record the peak memory, storage use, inference latency, and concurrent-request capacity on the target server and report those measurements together with the software versions.
