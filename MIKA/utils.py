import os
import yaml
import logging
import torch
import requests
from pathlib import Path
from typing import List, Dict, Any, Optional
from collections import deque
from neo4j import GraphDatabase
from sentence_transformers import CrossEncoder
from openai import OpenAI
from language_adapter import LanguageAdapter

try:
    import ollama
except ImportError:
    ollama = None

# ========== Milvus Proxy ==========
class RemoteMilvusClient:
    """
    Using docker container to connect to Milvus-Lite
    """
    def __init__(self, proxy_url: str, collection_name: str, timeout: int = 30):
        self.proxy_url = proxy_url.rstrip('/')
        self.collection_name = collection_name
        self.timeout = timeout
        logger.info(f"Linking Milvus Proxy: {proxy_url} | Collection: {collection_name}")
    
    def query(self, collection_name: str = None, filter_expr: str = None, 
              output_fields: List[str] = None, limit: int = 10, **kwargs):
        try:
            resp = requests.post(
                f"{self.proxy_url}/query",
                json={
                    "collection": collection_name or self.collection_name,
                    "filter": filter_expr,  
                    "output_fields": output_fields or ["*"],
                    "limit": limit
                },
                timeout=self.timeout
            )
            resp.raise_for_status()
            result = resp.json()
            if result.get("status") == "ok":
                return result.get("data", [])
            else:
                logger.warning(f"Query failed: {result}")
                return []
        except requests.RequestException as e:
            logger.error(f"Linking Milvus Proxy fail: {e}")
            return []
    
    def search(self, collection_name: str = None, data: List[List[float]] = None,
           limit: int = 10, output_fields: List[str] = None, **kwargs):
        try:
            vector = data[0] if data and isinstance(data[0], list) else (data[0] if data else [])
            
            resp = requests.post(
                f"{self.proxy_url}/search",
                json={
                    "collection": collection_name or self.collection_name,
                    "vector": vector,
                    "limit": limit,
                    "output_fields": output_fields or ["text"]
                },
                timeout=self.timeout
            )
            resp.raise_for_status()
            result = resp.json()
            
            if result.get("status") == "ok":
                return [result.get("data", [])]
            else:
                logger.warning(f"❌ Search failed: {result}")
                return [[]]
        except requests.RequestException as e:
            logger.error(f"❌ Connect Milvus Proxy failed: {e}")
            return [[]]
    
    def has_collection(self, collection_name: str = None) -> bool:
        try:
            resp = requests.get(f"{self.proxy_url}/health", timeout=5)
            return resp.status_code == 200
        except:
            return False
# ========== Milvus ==========

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("MaizeBDP MIKA")


def _normalize_config(cfg: dict) -> dict:
    """Normalize legacy config keys into the public framework-oriented schema."""
    cfg.setdefault("routing_model", cfg.pop("cloud_api", {}))

    legacy_local = cfg.pop("local_settings", {})
    cfg.setdefault("answer_model", {})
    cfg.setdefault("embedding_model", {})
    cfg.setdefault("retrieval_backend", {})
    cfg.setdefault("knowledge_graph", {})
    cfg.setdefault("assets", {})
    cfg.setdefault("search", {})

    if legacy_local:
        cfg["answer_model"].setdefault("provider", "ollama")
        cfg["answer_model"].setdefault("api_base_url", legacy_local.get("ollama_host"))
        cfg["answer_model"].setdefault("model", legacy_local.get("llm_model"))

        cfg["embedding_model"].setdefault("provider", "ollama")
        cfg["embedding_model"].setdefault("api_base_url", legacy_local.get("ollama_host"))
        cfg["embedding_model"].setdefault("model", legacy_local.get("embedding_model"))

        cfg["retrieval_backend"].setdefault("api_base_url", legacy_local.get("milvus_proxy_url"))
        cfg["retrieval_backend"].setdefault("timeout", legacy_local.get("milvus_timeout"))
        cfg["retrieval_backend"].setdefault("collection_name", legacy_local.get("collection_name"))

        cfg["knowledge_graph"].setdefault("uri", legacy_local.get("neo4j_uri"))
        cfg["knowledge_graph"].setdefault("auth", legacy_local.get("neo4j_auth"))

        cfg["assets"].setdefault("reranker_path", legacy_local.get("rerank_path"))

        cfg["search"].setdefault("query_topk", legacy_local.get("query_topk"))
        cfg["search"].setdefault("final_topk", legacy_local.get("final_topk"))
        cfg["search"].setdefault("token_limit", legacy_local.get("token_limit"))

    cfg["routing_model"].setdefault("api_key", "")
    cfg["routing_model"].setdefault("base_url", "")
    cfg["routing_model"].setdefault("model", "")

    cfg["answer_model"].setdefault("provider", "ollama")
    cfg["answer_model"].setdefault("api_key", "")
    cfg["answer_model"].setdefault("api_base_url", "http://127.0.0.1:11434")
    cfg["answer_model"].setdefault("model", "")

    cfg["embedding_model"].setdefault("provider", cfg["answer_model"]["provider"])
    cfg["embedding_model"].setdefault("api_key", cfg["answer_model"].get("api_key", ""))
    cfg["embedding_model"].setdefault("api_base_url", cfg["answer_model"].get("api_base_url", "http://127.0.0.1:11434"))
    cfg["embedding_model"].setdefault("model", "")

    cfg["retrieval_backend"].setdefault("api_base_url", "http://127.0.0.1:8081")
    cfg["retrieval_backend"].setdefault("timeout", 30)
    cfg["retrieval_backend"].setdefault("collection_name", "")

    cfg["knowledge_graph"].setdefault("uri", "bolt://127.0.0.1:7687")
    cfg["knowledge_graph"].setdefault("auth", ["neo4j", "password"])

    cfg["assets"].setdefault("reranker_path", "")

    cfg["search"].setdefault("query_topk", 10)
    cfg["search"].setdefault("final_topk", 5)
    cfg["search"].setdefault("token_limit", 2048)
    return cfg


def load_config() -> dict:
    """Load public template config and allow local/private overrides."""
    base_dir = Path(__file__).resolve().parent
    candidates = [
        os.environ.get("MIKA_CONFIG"),
        str(base_dir / "config.local.yaml"),
        str(base_dir / "config.yaml"),
    ]

    config_path = next((Path(path) for path in candidates if path and Path(path).exists()), None)
    if config_path is None:
        raise FileNotFoundError(
            "No configuration file found. Create MIKA/config.local.yaml or provide MIKA_CONFIG."
        )

    with config_path.open("r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f)

    cfg = _normalize_config(cfg)

    cfg["routing_model"]["api_key"] = os.environ.get(
        "MAIZEBDP_CLOUD_API_KEY",
        cfg["routing_model"].get("api_key", ""),
    )
    cfg["routing_model"]["base_url"] = os.environ.get(
        "MAIZEBDP_CLOUD_BASE_URL",
        cfg["routing_model"].get("base_url", ""),
    )
    cfg["routing_model"]["model"] = os.environ.get(
        "MAIZEBDP_CLOUD_MODEL",
        cfg["routing_model"].get("model", ""),
    )

    cfg["answer_model"]["api_base_url"] = os.environ.get(
        "MAIZEBDP_ANSWER_API_BASE_URL",
        cfg["answer_model"].get("api_base_url", "http://127.0.0.1:11434"),
    )
    cfg["answer_model"]["model"] = os.environ.get(
        "MAIZEBDP_ANSWER_MODEL",
        cfg["answer_model"].get("model", ""),
    )
    cfg["answer_model"]["api_key"] = os.environ.get(
        "MAIZEBDP_ANSWER_API_KEY",
        cfg["answer_model"].get("api_key", ""),
    )

    cfg["embedding_model"]["api_base_url"] = os.environ.get(
        "MAIZEBDP_EMBEDDING_API_BASE_URL",
        cfg["embedding_model"].get("api_base_url", cfg["answer_model"]["api_base_url"]),
    )
    cfg["embedding_model"]["model"] = os.environ.get(
        "MAIZEBDP_EMBEDDING_MODEL",
        cfg["embedding_model"].get("model", ""),
    )
    cfg["embedding_model"]["api_key"] = os.environ.get(
        "MAIZEBDP_EMBEDDING_API_KEY",
        cfg["embedding_model"].get("api_key", cfg["answer_model"].get("api_key", "")),
    )

    cfg["retrieval_backend"]["api_base_url"] = os.environ.get(
        "MAIZEBDP_RETRIEVAL_API_BASE_URL",
        cfg["retrieval_backend"].get("api_base_url", "http://127.0.0.1:8081"),
    )
    cfg["knowledge_graph"]["uri"] = os.environ.get(
        "MAIZEBDP_NEO4J_URI",
        cfg["knowledge_graph"].get("uri", "bolt://127.0.0.1:7687"),
    )

    neo4j_user = os.environ.get("MAIZEBDP_NEO4J_USER")
    neo4j_password = os.environ.get("MAIZEBDP_NEO4J_PASSWORD")
    if neo4j_user and neo4j_password:
        cfg["knowledge_graph"]["auth"] = [neo4j_user, neo4j_password]

    alias_path = cfg.get("variety_alias_path")
    if alias_path and not Path(alias_path).is_absolute():
        cfg["variety_alias_path"] = str((config_path.parent / alias_path).resolve())

    logger.info("Loaded configuration from %s", config_path)
    return cfg


# parameters
cfg = load_config()

# Client initialization
routing_client = OpenAI(
    api_key=cfg['routing_model']['api_key'],
    base_url=cfg['routing_model']['base_url']
)

answer_provider = cfg["answer_model"]["provider"]
embedding_provider = cfg["embedding_model"]["provider"]


def _make_ollama_client(base_url: str):
    if ollama is None:
        raise ImportError("The `ollama` package is required when provider is set to `ollama`.")
    return ollama.Client(host=base_url)


answer_client = None
if answer_provider == "ollama":
    answer_client = _make_ollama_client(cfg["answer_model"]["api_base_url"])
elif answer_provider == "openai_compatible":
    answer_client = OpenAI(
        api_key=cfg["answer_model"]["api_key"],
        base_url=cfg["answer_model"]["api_base_url"],
    )

embedding_client = None
if embedding_provider == "ollama":
    embedding_client = _make_ollama_client(cfg["embedding_model"]["api_base_url"])
elif embedding_provider == "openai_compatible":
    embedding_client = OpenAI(
        api_key=cfg["embedding_model"]["api_key"],
        base_url=cfg["embedding_model"]["api_base_url"],
    )

proxy_url = cfg['retrieval_backend'].get('api_base_url', 'http://localhost:8081')
collection_name = cfg['retrieval_backend']['collection_name']
milvus_client = RemoteMilvusClient(
    proxy_url=proxy_url,
    collection_name=collection_name,
    timeout=cfg['retrieval_backend'].get('timeout', 30)
)
neo4j_driver = GraphDatabase.driver(
    cfg['knowledge_graph']['uri'],
    auth=tuple(cfg['knowledge_graph']['auth'])
)
reranker = CrossEncoder(cfg['assets']['reranker_path'],
                        device="cuda" if torch.cuda.is_available() else "cpu")

# Language Adapter
language_adapter = LanguageAdapter(cfg)


def generate_embedding(text: str) -> List[float]:
    """Generate an embedding using the configured embedding provider."""
    model = cfg["embedding_model"]["model"]
    if embedding_provider == "ollama":
        return embedding_client.embeddings(model=model, prompt=text)["embedding"]
    if embedding_provider == "openai_compatible":
        resp = embedding_client.embeddings.create(model=model, input=text)
        return resp.data[0].embedding
    raise ValueError(f"Unsupported embedding provider: {embedding_provider}")


def generate_answer(messages: List[Dict[str, str]], temperature: float = 0.3) -> str:
    """Generate a chat response using the configured answer-model provider."""
    model = cfg["answer_model"]["model"]
    if answer_provider == "ollama":
        resp = answer_client.chat(
            model=model,
            messages=messages,
            options={"temperature": temperature},
        )
        return resp["message"]["content"]
    if answer_provider == "openai_compatible":
        resp = answer_client.chat.completions.create(
            model=model,
            messages=messages,
            temperature=temperature,
            max_tokens=cfg["search"]["token_limit"],
        )
        return resp.choices[0].message.content.strip()
    raise ValueError(f"Unsupported answer-model provider: {answer_provider}")


# Hybrid Retrieval for RAG
def hybrid_retrieval(query: str) -> List[Dict[str, Any]]:
    """
    Hybrid Retrieval for RAG
    """
    output_fields = ["text", "source_id", "chunk_id"]
    # 1. keywords
    keywords = [word for word in query.split() if len(word) > 4]
    filter_expr = " OR ".join([f'text LIKE "%{kw}%"' for kw in keywords])
    
    kw_hits = milvus_client.query(
        collection_name=cfg['retrieval_backend']['collection_name'],
        filter_expr=filter_expr,
        output_fields=output_fields,
        limit=cfg['search']['query_topk']
    )

    # 2. semantic
    emb = generate_embedding(query)

    sem_hits = milvus_client.search(
        collection_name=cfg['retrieval_backend']['collection_name'],
        data=[emb],
        limit=cfg['search']['query_topk'],
        output_fields=output_fields
    )[0]

    if not sem_hits:
        return []

    # merge + dedup
    seen = {}
    for h in kw_hits + sem_hits:
        key = (h["source_id"], h["chunk_id"])
        seen[key] = h

    unique_hits = list(seen.values())

    texts = [h["text"] for h in unique_hits]
    scores = reranker.predict([[query, t] for t in texts])

    for h, s in zip(unique_hits, scores):
        h["rerank_score"] = float(s)

    return sorted(
        unique_hits,
        key=lambda x: x["rerank_score"],
        reverse=True
    )[:cfg['search']['final_topk']]
