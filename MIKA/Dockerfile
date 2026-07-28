FROM python:3.11-slim

WORKDIR /app

RUN apt-get update && apt-get install -y --no-install-recommends git && \
    rm -rf /var/lib/apt/lists/*

# 配置 pip 使用国内镜像 + 超时设置
RUN pip config set global.index-url https://pypi.tuna.tsinghua.edu.cn/simple && \
    pip config set global.timeout 300 && \
    pip config set global.retries 10

# 先装 torch CPU-only（最大的包）
RUN pip install --no-cache-dir torch --index-url https://download.pytorch.org/whl/cpu \
    -i https://pypi.tuna.tsinghua.edu.cn/simple

# 再装其他包
RUN pip install --no-cache-dir \
    langgraph langchain-core \
    fastapi "uvicorn[standard]" \
    pydantic pyyaml httpx requests \
    pymilvus neo4j ollama openai \
    sentence-transformers typing_extensions \
    -i https://pypi.tuna.tsinghua.edu.cn/simple

EXPOSE 8082

CMD ["python", "main.py"]
